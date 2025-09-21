// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

#include <random>
#include <limits>
#include <cstdint>

/*
 * our_lm_solve
 * ------------
 * Fit linear models y ~ X for multiple response columns (Y),
 * compute a contrast on coefficients, and obtain permutation p-values.
 *
 * PARAMETERS (key)
 *  X  (n x p)           : design matrix
 *  Y  (n x m)           : response matrix (columns modeled independently)
 *  contrast (length p)  : linear contrast on coefficients
 *  perm_groups (list)   : NULL => full shuffle; else list of 1-based integer
 *                         vectors with within-group permutations
 *  n_randomizations     : number of permutations
 *  alternative          : "two-sided" | "greater" | "less"
 *  return_residuals     : include n x m residual matrix (NA where y was NA)
 *  return_sampled_fits  : include per-permutation β (B x p) per Y column
 *  use_huber            : use Huber M-estimation via IRLS (robust regression)
 *  huber_k, huber_maxit, huber_tol : Huber parameters
 *  na_mode              : "drop" | "impute_weak" (impute to mean/zero, small weight)
 *  na_weight            : weight for imputed rows in "impute_weak"
 *  na_center            : "mean" | "zero"       (imputation center)
 *  illcond_rcond        : rcond(X'X) threshold; below -> treat as ill-conditioned
 *  pinv_tol             : tolerance for pinv
 *  n_cores              : per-column OpenMP parallelism (>1 enables parallel if available)
 *
 * RETURNS
 *  coef (p x m), stat (m), p_value (m), residuals (opt), sampled_fits (opt list)
 */

/*** ======================= Utilities & helpers ======================= ***/

// ---- SERIAL (R RNG) shuffles (fast when n_cores == 1) ----
static inline void shuffle_vec_in_place_R(arma::vec& v) {
  for (arma::uword i = v.n_elem; i > 1; --i) {
    arma::uword j = (arma::uword) std::floor(R::unif_rand() * static_cast<double>(i));
    if (j >= i) j = i - 1;
    std::swap(v[i - 1], v[j]);
  }
}
static inline void shuffle_selected_in_place_R(arma::vec& y, const arma::uvec& idx) {
  arma::uword k = idx.n_elem;
  for (arma::uword i = k; i > 1; --i) {
    arma::uword j = (arma::uword) std::floor(R::unif_rand() * static_cast<double>(i));
    if (j >= i) j = i - 1;
    double tmp = y[idx[i - 1]];
    y[idx[i - 1]] = y[idx[j]];
    y[idx[j]]     = tmp;
  }
}
static inline void permute_by_groups_R(arma::vec& y, const std::vector<std::vector<arma::uword>>& groups) {
  for (const auto& g : groups) {
    if (g.size() < 2) continue;
    for (std::size_t i = g.size(); i > 1; --i) {
      std::size_t j = (std::size_t) std::floor(R::unif_rand() * (double)i);
      if (j >= i) j = i - 1;
      std::swap(y[g[i - 1]], y[g[j]]);
    }
  }
}

// ---- PARALLEL (C++ RNG) shuffles (thread-safe when n_cores > 1) ----
template <class URNG>
static inline void shuffle_vec_in_place_rng(arma::vec& v, URNG& rng) {
  for (arma::uword i = v.n_elem; i > 1; --i) {
    std::uniform_int_distribution<arma::uword> dist(0, i - 1);
    arma::uword j = dist(rng);
    std::swap(v[i - 1], v[j]);
  }
}
template <class URNG>
static inline void shuffle_selected_in_place_rng(arma::vec& y, const arma::uvec& idx, URNG& rng) {
  arma::uword k = idx.n_elem;
  for (arma::uword i = k; i > 1; --i) {
    std::uniform_int_distribution<arma::uword> dist(0, i - 1);
    arma::uword j = dist(rng);
    double tmp = y[idx[i - 1]];
    y[idx[i - 1]] = y[idx[j]];
    y[idx[j]]     = tmp;
  }
}
template <class URNG>
static inline void permute_by_groups_rng(arma::vec& y, const std::vector<std::vector<arma::uword>>& groups, URNG& rng) {
  for (const auto& g : groups) {
    if (g.size() < 2) continue;
    for (std::size_t i = g.size(); i > 1; --i) {
      std::uniform_int_distribution<std::size_t> dist(0, i - 1);
      std::size_t j = dist(rng);
      std::swap(y[g[i - 1]], y[g[j]]);
    }
  }
}

// Build 0-based groups once (from 1-based perm_groups in full data)
static std::vector<std::vector<arma::uword>>
  build_groups0_full(const List& perm_groups_full, arma::uword n) {
    std::vector<std::vector<arma::uword>> out;
    if (perm_groups_full.size() == 0) return out;
    out.reserve(perm_groups_full.size());
    for (int g = 0; g < perm_groups_full.size(); ++g) {
      IntegerVector grp = perm_groups_full[g];
      std::vector<arma::uword> mapped; mapped.reserve(grp.size());
      for (int a = 0; a < grp.size(); ++a) {
        int i1 = grp[a];
        if (i1 >= 1 && i1 <= (int)n) mapped.push_back((arma::uword)(i1 - 1));
      }
      out.push_back(std::move(mapped));
    }
    return out;
  }

// Map full-data (0-based) groups to subset (0-based) groups using idx_fin
static std::vector<std::vector<arma::uword>>
  map_groups_full_to_subset(const std::vector<std::vector<arma::uword>>& groups_full0,
                            const arma::uvec& idx_fin,
                            arma::uword n_full) {
    std::vector<int> pos(n_full, -1);
    for (arma::uword t = 0; t < idx_fin.n_elem; ++t) pos[(std::size_t)idx_fin[t]] = (int)t;
    
    std::vector<std::vector<arma::uword>> out;
    out.reserve(groups_full0.size());
    for (const auto& g : groups_full0) {
      std::vector<arma::uword> h; h.reserve(g.size());
      for (arma::uword full_idx : g) {
        int p = pos[(std::size_t)full_idx];
        if (p >= 0) h.push_back((arma::uword)p);
      }
      out.push_back(std::move(h));
    }
    return out;
  }

// Symmetrize & check finiteness
static inline bool symmetrize_and_check(arma::mat& A) {
  A = arma::symmatu(A);
  return A.is_finite();
}

// pinv helper (avoid ternary overload mismatch)
static inline arma::mat pinv_safe(const arma::mat& A, double tol) {
  arma::mat P;
  if (tol > 0.0) P = arma::pinv(A, tol);
  else           P = arma::pinv(A);
  return P;
}

// inv(X'X) via SPD, else pinv; also return Xt
static inline bool inv_xtx(const arma::mat& X, arma::mat& invXtX, arma::mat& Xt, double pinv_tol=0.0) {
  Xt = X.t();
  arma::mat XtX = Xt * X;
  if (!symmetrize_and_check(XtX)) return false;
  if (inv_sympd(invXtX, XtX)) return true;
  invXtX = pinv_safe(XtX, pinv_tol);
  return invXtX.is_finite();
}

// OLS from precomputed inv(X'X), X'
static inline arma::vec coef_from_inv(const arma::mat& invXtX, const arma::mat& Xt, const arma::vec& y) {
  return invXtX * (Xt * y);
}

// Contrast statistic
static inline double contrast_stat(const arma::vec& beta, const arma::vec& contrast) {
  return dot(beta, contrast);
}

// Tally permutation
static inline void tally_perm(double obs, double perm, int alt, int& ge, int& le, int& ge_abs) {
  if (alt == 0) { if (std::fabs(perm) >= std::fabs(obs)) ge_abs++; }
  else if (alt == 1) { if (perm >= obs) ge++; }
  else { if (perm <= obs) le++; }
}

// Safe MAD scale (NA/Inf-aware, positive fallback)
static inline double robust_scale_mad_safe(const arma::vec& r) {
  arma::uvec idx = arma::find_finite(r);
  if (idx.n_elem == 0) return 1e-8;
  arma::vec rf = r.elem(idx);
  double med = arma::median(rf);
  arma::vec af = arma::abs(rf - med);
  double mad = (af.n_elem > 0) ? arma::median(af) : 0.0;
  double s = 1.4826 * mad;
  if (!(s > 0.0)) {
    s = std::sqrt(arma::mean(arma::square(rf))) + 1e-12;
    if (!(s > 0.0)) s = 1e-8;
  }
  return s;
}

// Weighted normal equations with small adaptive ridge
static arma::vec wls_solve_ridge(const arma::mat& X,
                                 const arma::vec& y,
                                 const arma::vec& w,
                                 double ridge_eps = 1e-10,
                                 double pinv_tol = 0.0) {
  arma::mat XtWX = X.t() * (X.each_col() % w);
  arma::vec XtWy = X.t() * (y % w);
  
  if (!symmetrize_and_check(XtWX)) {
    arma::vec out(X.n_cols); out.fill(arma::datum::nan); return out;
  }
  
  arma::mat Minv;
  if (inv_sympd(Minv, XtWX)) return Minv * XtWy;
  
  double tr = arma::trace(XtWX);
  double p  = static_cast<double>(X.n_cols);
  double lambda = ridge_eps * ((tr > 0.0) ? tr / p : 1.0);
  XtWX.diag() += lambda;
  
  if (inv_sympd(Minv, XtWX)) return Minv * XtWy;
  
  arma::mat P = pinv_safe(XtWX, pinv_tol);
  if (!P.is_finite()) { arma::vec out(X.n_cols); out.fill(arma::datum::nan); return out; }
  return P * XtWy;
}

// Huber IRLS (optionally base weights); init at LS; safe MAD; ridge in inner WLS
static arma::vec huber_irls(const arma::mat& X, const arma::vec& y,
                            double k, int maxit, double tol,
                            const arma::mat* invXtX_opt = nullptr,
                            const arma::mat* Xt_opt     = nullptr,
                            const arma::vec* base_w     = nullptr,
                            double pinv_tol = 0.0) {
  const arma::uword p = X.n_cols;
  arma::vec beta(p, arma::fill::zeros);
  
  // Initialize at (weighted) LS
  if (base_w) {
    beta = wls_solve_ridge(X, y, *base_w, 1e-10, pinv_tol);
  } else if (invXtX_opt && Xt_opt) {
    arma::vec b0 = (*invXtX_opt) * ((*Xt_opt) * y);
    if (!b0.is_finite()) { beta.fill(arma::datum::nan); return beta; }
    beta = b0;
  } else {
    arma::mat Xt, invXtX;
    if (!inv_xtx(X, invXtX, Xt, pinv_tol)) { beta.fill(arma::datum::nan); return beta; }
    beta = invXtX * (Xt * y);
  }
  
  for (int it = 0; it < maxit; ++it) {
    if (!beta.is_finite()) return beta;
    arma::vec r = y - X * beta;
    if (!r.is_finite()) { beta.fill(arma::datum::nan); return beta; }
    
    double s = robust_scale_mad_safe(r);
    if (!(s > 0.0)) break;
    
    const double ks = k * s;
    arma::vec w_hub(r.n_elem, arma::fill::ones);
    for (arma::uword i = 0; i < r.n_elem; ++i) {
      double ar = std::fabs(r[i]);
      if (ar > ks) w_hub[i] = ks / ar;
    }
    arma::vec w_eff = base_w ? (*base_w) % w_hub : w_hub;
    
    arma::vec beta_new = wls_solve_ridge(X, y, w_eff, 1e-10, pinv_tol);
    if (!beta_new.is_finite()) { beta.fill(arma::datum::nan); return beta; }
    
    double denom = arma::norm(beta, 2) + 1e-12;
    double rel_change = arma::norm(beta_new - beta, 2) / denom;
    beta = beta_new;
    if (rel_change < tol) break;
  }
  return beta;
}

static inline double mean_finite(const arma::vec& y) {
  arma::uvec idx = arma::find_finite(y);
  if (idx.n_elem == 0) return 0.0;
  return arma::mean(y.elem(idx));
}

// Ill-conditioning via rcond(X'X) or k < p
static inline bool is_ill_conditioned(const arma::mat& X, double rcond_thresh) {
  const arma::uword k = X.n_rows, p = X.n_cols;
  if (k < p) return true;
  arma::mat XtX = X.t() * X;
  if (!symmetrize_and_check(XtX)) return true;
  double rc = arma::rcond(XtX);
  if (!(rc > 0.0)) return true;
  return (rc < rcond_thresh);
}

// Deterministic per-column RNG from base seed
static inline std::mt19937_64 make_rng_for_column(std::uint64_t base, arma::uword j) {
  std::uint64_t x = base ^ (0x9e3779b97f4a7c15ULL + j + (j<<6) + (j>>2));
  return std::mt19937_64(x);
}

/*** ================= Main exported solver (parallel over Y columns) ================= ***/

// [[Rcpp::export]]
Rcpp::List our_lm_solve(const arma::mat& X,
                        const arma::mat& Y,
                        const arma::vec& contrast,
                        Rcpp::Nullable<Rcpp::List> perm_groups = R_NilValue,
                        int n_randomizations = 100,
                        std::string alternative = "two-sided",
                        bool return_residuals = true,
                        bool return_sampled_fits = false,
                        bool use_huber = false,
                        double huber_k = 1.345,
                        int huber_maxit = 8,
                        double huber_tol = 1e-6,
                        std::string na_mode = "drop",
                        double na_weight = 1e-4,
                        std::string na_center = "mean",
                        double illcond_rcond = 1e-12,
                        double pinv_tol = 0.0,
                        int n_cores = 1) {
  RNGScope scope;
  
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const arma::uword m = Y.n_cols;
  
  if (Y.n_rows != n) stop("X and Y must have the same number of rows.");
  if (contrast.n_elem != p) stop("contrast length must equal ncol(X).");
  if (n_randomizations < 0) stop("n_randomizations must be >= 0.");
  
  int alt = 0;
  if (alternative == "two-sided") alt = 0;
  else if (alternative == "greater") alt = 1;
  else if (alternative == "less") alt = 2;
  else stop("alternative must be 'two-sided','greater', or 'less'.");
  
  const bool impute_mode = (na_mode == "impute_weak");
  const bool center_mean = (na_center == "mean");
  
  // Precompute full-design OLS pieces
  arma::mat invXtX_full, Xt_full;
  if (!inv_xtx(X, invXtX_full, Xt_full, pinv_tol)) {
    stop("Design X is ill-conditioned or non-finite in full data.");
  }
  arma::mat B_full = invXtX_full * Xt_full;  // p x n
  
  // Permutation groups
  List perm_full = perm_groups.isNotNull() ? List(perm_groups) : List();
  const bool full_permute_requested = (perm_full.size() == 0);
  std::vector<std::vector<arma::uword>> groups_full0;
  if (!full_permute_requested) groups_full0 = build_groups0_full(perm_full, n);
  
  // Outputs
  arma::mat coef_obs(p, m); coef_obs.fill(arma::datum::nan);
  arma::vec stat_obs(m);    stat_obs.fill(arma::datum::nan);
  arma::vec pvals(m);       pvals.fill(arma::datum::nan);
  arma::mat resid_out;
  if (return_residuals) { resid_out.set_size(n, m); resid_out.fill(arma::datum::nan); }
  std::vector<arma::mat> sampled_list; 
  if (return_sampled_fits) sampled_list.resize(m);
  
  // Parallel?
  bool parallel_mode = false;
#ifdef _OPENMP
  if (n_cores > 1) { omp_set_num_threads(n_cores); parallel_mode = true; }
#endif
  std::uint64_t base_seed = 0xD1B54A32D192ED03ULL
  ^ (std::uint64_t) std::floor(R::unif_rand() * std::numeric_limits<uint32_t>::max())
    ^ (((std::uint64_t) std::floor(R::unif_rand() * std::numeric_limits<uint32_t>::max())) << 32);
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(parallel_mode)
#endif
    for (int jj = 0; jj < static_cast<int>(m); ++jj) {
      const arma::uword j = static_cast<arma::uword>(jj);
      
      // single-column uvec for residual writes
      arma::uvec cj(1); cj[0] = j;
      
      std::mt19937_64 rng;
      if (parallel_mode) rng = make_rng_for_column(base_seed, j);
      
      const bool need_beta_perm = (use_huber || return_sampled_fits);
      const bool stat_only      = !need_beta_perm;
      
      arma::vec yj = Y.col(j);
      arma::uvec idx_fin = arma::find_finite(yj);
      bool all_finite = (idx_fin.n_elem == n);
      
      arma::vec yc;
      arma::uword k = n;
      arma::mat Xc;                         // persist subset design for lifetime of loop
      const arma::mat* Xobs = &X;           // pointer to current design (X or Xc)
      arma::mat invXtX, Xt;
      arma::mat B;                          // mapping y -> beta (p x n or p x k)
      arma::vec alpha;                      // mapping y -> stat (alpha^T y), only when stat_only
      
      bool use_full_shuffle = full_permute_requested;
      std::vector<std::vector<arma::uword>> groups_sub;
      arma::uvec eligible_idx;
      
      // ----------- NA handling: all-finite, drop, impute_weak -----------
      if (all_finite) {
        yc = yj;
        invXtX = invXtX_full;  Xt = Xt_full;  B = B_full;
        Xobs = &X;
        if (!full_permute_requested) { use_full_shuffle = false; groups_sub = groups_full0; }
        if (stat_only) alpha = B.t() * contrast;  // n-vector
        
      } else if (!impute_mode || na_weight <= 0.0) {  // DROP
        if (idx_fin.n_elem == 0) { pvals[j] = NA_REAL; continue; }
        Xc = X.rows(idx_fin);              // keep Xc alive
        yc = yj.elem(idx_fin);
        k  = Xc.n_rows;
        Xobs = &Xc;
        
        if (is_ill_conditioned(Xc, illcond_rcond)) {
          pvals[j]=NA_REAL; stat_obs[j]=NA_REAL; coef_obs.col(j).fill(arma::datum::nan);
          continue;
        }
        if (!inv_xtx(Xc, invXtX, Xt, pinv_tol)) {
          pvals[j]=NA_REAL; stat_obs[j]=NA_REAL; coef_obs.col(j).fill(arma::datum::nan);
          continue;
        }
        B = invXtX * Xt;  // p x k
        if (!full_permute_requested) {
          groups_sub = map_groups_full_to_subset(groups_full0, idx_fin, n);
          use_full_shuffle = false;
        } else {
          eligible_idx = arma::regspace<uvec>(0, k - 1);
        }
        if (stat_only) alpha = B.t() * contrast; // k-vector
        
      } else {  // IMPUTE_WEAK
        yc = yj;
        double mu = center_mean ? mean_finite(yj) : 0.0;
        
        arma::vec w_base(n); w_base.fill(na_weight);
        for (arma::uword t = 0; t < idx_fin.n_elem; ++t) w_base[idx_fin[t]] = 1.0;
        arma::vec sqrtw = arma::sqrt(w_base);
        
        if (idx_fin.n_elem < n) {
          yc.fill(mu);
          for (arma::uword t = 0; t < idx_fin.n_elem; ++t) yc[idx_fin[t]] = yj[idx_fin[t]];
        }
        
        if (use_huber) {
          // Effective sample check (observed or near-observed rows)
          arma::uword neff = arma::accu(w_base > 1e-8);
          if (neff < p + 1) {  // not enough information
            coef_obs.col(j).fill(arma::datum::nan); stat_obs[j] = NA_REAL; pvals[j] = NA_REAL; 
            continue;
          }
          
          arma::vec beta = huber_irls(X, yc, huber_k, huber_maxit, huber_tol, nullptr, nullptr, &w_base, pinv_tol);
          if (!beta.is_finite()) { pvals[j]=NA_REAL; stat_obs[j]=NA_REAL; coef_obs.col(j).fill(arma::datum::nan); continue; }
          
          double sobs = dot(beta, contrast);
          coef_obs.col(j) = beta; stat_obs[j] = sobs;
          
          if (return_residuals) {
            arma::vec rc = yc - X * beta;                // length n
            resid_out.submat(idx_fin, cj) = rc.elem(idx_fin);
          }
          
          if (n_randomizations == 0) { pvals[j] = NA_REAL; continue; }
          
          // Prepare groups: restrict permutations to observed rows
          bool can_permute = false;
          if (full_permute_requested) {
            eligible_idx = idx_fin; 
            can_permute = (eligible_idx.n_elem >= 2);
          } else {
            std::vector<std::vector<arma::uword>> groups_filt;
            groups_filt.reserve(groups_full0.size());
            std::vector<char> elig(n, 0);
            for (arma::uword t = 0; t < idx_fin.n_elem; ++t) elig[(std::size_t)idx_fin[t]] = 1;
            for (const auto& g : groups_full0) {
              std::vector<arma::uword> h; h.reserve(g.size());
              for (arma::uword full_idx : g) if (elig[(std::size_t)full_idx]) h.push_back(full_idx);
              if (h.size() >= 2) can_permute = true;
              groups_filt.push_back(std::move(h));
            }
            groups_sub = std::move(groups_filt);
            use_full_shuffle = false;
          }
          if (!can_permute) { pvals[j] = 1.0; continue; }
          
          int ge=0, le=0, ge_abs=0; arma::vec yperm = yc;
          arma::mat Mcoef; if (return_sampled_fits) Mcoef.set_size(n_randomizations, p);
          bool failed=false;
          
          for (int r=0; r<n_randomizations; ++r) {
            yperm = yc;
            if (use_full_shuffle) {
              if (parallel_mode) shuffle_selected_in_place_rng(yperm, eligible_idx, rng);
              else               shuffle_selected_in_place_R (yperm, eligible_idx);
            } else {
              if (parallel_mode) permute_by_groups_rng(yperm, groups_sub, rng);
              else               permute_by_groups_R  (yperm, groups_sub);
            }
            
            arma::vec bperm = huber_irls(X, yperm, huber_k, huber_maxit, huber_tol, nullptr, nullptr, &w_base, pinv_tol);
            if (!bperm.is_finite()) { failed = true; break; }
            
            double sperm = dot(bperm, contrast);
            tally_perm(sobs, sperm, alt, ge, le, ge_abs);
            if (return_sampled_fits) Mcoef.row(r) = bperm.t();
          }
          
          if (failed) pvals[j] = NA_REAL;
          else {
            double pval = (alt==0) ? (ge_abs+1.0)/(n_randomizations+1.0)
              : (alt==1) ? (ge+1.0)/(n_randomizations+1.0)
              : (le+1.0)/(n_randomizations+1.0);
            pvals[j] = pval;
            if (return_sampled_fits) sampled_list[j] = std::move(Mcoef);
          }
          continue;
        }
        
        // Weighted OLS fast map (full X)
        arma::mat Xw = X; Xw.each_col() %= sqrtw;
        arma::mat Xt_w = Xw.t();
        arma::mat invXtXw;
        if (!inv_xtx(Xw, invXtXw, Xt_w, pinv_tol)) {
          pvals[j]=NA_REAL; stat_obs[j]=NA_REAL; coef_obs.col(j).fill(arma::datum::nan); continue;
        }
        arma::mat Bw = invXtXw * Xt_w;     // p x n
        arma::vec ycw = sqrtw % yc;
        arma::vec beta = Bw * ycw;
        
        if (!beta.is_finite()) { pvals[j]=NA_REAL; stat_obs[j]=NA_REAL; coef_obs.col(j).fill(arma::datum::nan); continue; }
        
        double sobs = dot(beta, contrast);
        coef_obs.col(j) = beta; stat_obs[j] = sobs;
        
        if (return_residuals) {
          arma::vec rc = yc - X * beta;    // length n
          resid_out.submat(idx_fin, cj) = rc.elem(idx_fin);
        }
        
        if (n_randomizations == 0) { pvals[j] = NA_REAL; continue; }
        
        // Prepare groups (restrict permutations to observed rows)
        bool can_permute = false;
        if (full_permute_requested) { can_permute = (idx_fin.n_elem >= 2); }
        else {
          for (const auto& g : groups_full0) {
            std::size_t cnt=0; for (auto u: g) if (u < n && arma::is_finite(yc[u])) ++cnt;
            if (cnt >= 2) { can_permute = true; break; }
          }
        }
        if (!can_permute) { pvals[j] = 1.0; continue; }
        
        // stat-only alpha for weighted case: alpha = (Bw^T c) ∘ sqrt(w)
        arma::vec alpha_w = Bw.t() * contrast;           // n
        arma::vec alpha_w_unweighted = alpha_w % sqrtw;
        
        int ge=0, le=0, ge_abs=0; arma::vec yperm = yc;
        arma::mat Mcoef; if (return_sampled_fits) Mcoef.set_size(n_randomizations, p);
        bool failed=false;
        
        for (int r=0; r<n_randomizations; ++r) {
          yperm = yc;
          if (full_permute_requested) {
            if (parallel_mode) shuffle_selected_in_place_rng(yperm, idx_fin, rng);
            else               shuffle_selected_in_place_R (yperm, idx_fin);
          } else {
            std::vector<std::vector<arma::uword>> groups_filt;
            groups_filt.reserve(groups_full0.size());
            std::vector<char> elig(n, 0);
            for (arma::uword t = 0; t < idx_fin.n_elem; ++t) elig[(std::size_t)idx_fin[t]] = 1;
            for (const auto& g : groups_full0) {
              std::vector<arma::uword> h; h.reserve(g.size());
              for (arma::uword full_idx : g) if (elig[(std::size_t)full_idx]) h.push_back(full_idx);
              groups_filt.push_back(std::move(h));
            }
            if (parallel_mode) permute_by_groups_rng(yperm, groups_filt, rng);
            else               permute_by_groups_R  (yperm, groups_filt);
          }
          
          if (return_sampled_fits) {
            arma::vec bw = Bw * (sqrtw % yperm);
            if (!bw.is_finite()) { failed = true; break; }
            double sperm = dot(bw, contrast);
            tally_perm(sobs, sperm, alt, ge, le, ge_abs);
            Mcoef.row(r) = bw.t();
          } else {
            double sperm = dot(alpha_w_unweighted, yperm);
            tally_perm(sobs, sperm, alt, ge, le, ge_abs);
          }
        }
        
        if (failed) pvals[j] = NA_REAL;
        else {
          double pval = (alt==0) ? (ge_abs+1.0)/(n_randomizations+1.0)
            : (alt==1) ? (ge+1.0)/(n_randomizations+1.0)
            : (le+1.0)/(n_randomizations+1.0);
          pvals[j] = pval;
          if (return_sampled_fits) sampled_list[j] = std::move(Mcoef);
        }
        continue;  // end IMPUTE_WEAK
      }
      
      // ----------- Observed fit (unweighted: NA-free or DROP) -----------
      arma::vec beta;
      if (use_huber) {
        beta = huber_irls(*Xobs, yc, huber_k, huber_maxit, huber_tol, &invXtX, &Xt, nullptr, pinv_tol);
      } else {
        beta = coef_from_inv(invXtX, Xt, yc);
      }
      if (!beta.is_finite()) { pvals[j]=NA_REAL; stat_obs[j]=NA_REAL; coef_obs.col(j).fill(arma::datum::nan); continue; }
      
      double sobs = dot(beta, contrast);
      coef_obs.col(j) = beta; stat_obs[j] = sobs;
      
      if (return_residuals) {
        if (all_finite) { resid_out.col(j) = yc - (*Xobs) * beta; }
        else {
          arma::vec rc = yc - (*Xobs) * beta;     // length k
          resid_out.submat(idx_fin, cj) = rc;     // scatter into full n
        }
      }
      
      // ----------- Permutations (unweighted) -----------
      if (n_randomizations == 0) { pvals[j] = NA_REAL; continue; }
      
      bool can_permute = false;
      if (all_finite) {
        if (full_permute_requested)  can_permute = (n >= 2);
        else for (const auto& g : groups_sub) if (g.size() >= 2) { can_permute = true; break; }
      } else {
        if (full_permute_requested)  can_permute = (k >= 2);
        else for (const auto& g : groups_sub) if (g.size() >= 2) { can_permute = true; break; }
      }
      if (!can_permute) { pvals[j] = 1.0; continue; }
      
      if (!need_beta_perm) {
        if (alpha.n_elem == 0) alpha = B.t() * contrast;  // n- or k-vector
      }
      
      int ge=0, le=0, ge_abs=0; arma::vec yperm = yc;
      arma::mat Mcoef; if (return_sampled_fits) Mcoef.set_size(n_randomizations, p);
      bool failed=false;
      
      for (int r=0; r<n_randomizations; ++r) {
        yperm = yc;
        if (all_finite) {
          if (full_permute_requested) {
            if (parallel_mode) shuffle_vec_in_place_rng(yperm, rng);
            else               shuffle_vec_in_place_R (yperm);
          } else {
            if (parallel_mode) permute_by_groups_rng(yperm, groups_sub, rng);
            else               permute_by_groups_R  (yperm, groups_sub);
          }
        } else {
          if (full_permute_requested) {
            if (eligible_idx.n_elem == 0) eligible_idx = arma::regspace<uvec>(0, k - 1);
            if (parallel_mode) shuffle_selected_in_place_rng(yperm, eligible_idx, rng);
            else               shuffle_selected_in_place_R (yperm, eligible_idx);
          } else {
            if (parallel_mode) permute_by_groups_rng(yperm, groups_sub, rng);
            else               permute_by_groups_R  (yperm, groups_sub);
          }
        }
        
        if (need_beta_perm) {
          arma::vec bperm;
          if (use_huber) {
            bperm = huber_irls(*Xobs, yperm, huber_k, huber_maxit, huber_tol, &invXtX, &Xt, nullptr, pinv_tol);
          } else {
            bperm = B * yperm;
          }
          if (!bperm.is_finite()) { failed = true; break; }
          double sperm = dot(bperm, contrast);
          tally_perm(sobs, sperm, alt, ge, le, ge_abs);
          if (return_sampled_fits) Mcoef.row(r) = bperm.t();
        } else {
          double sperm = dot(alpha, yperm);
          tally_perm(sobs, sperm, alt, ge, le, ge_abs);
        }
      }
      
      if (failed) pvals[j] = NA_REAL;
      else {
        double pval = (alt==0) ? (ge_abs+1.0)/(n_randomizations+1.0)
          : (alt==1) ? (ge+1.0)/(n_randomizations+1.0)
          : (le+1.0)/(n_randomizations+1.0);
        pvals[j] = pval;
        if (return_sampled_fits) sampled_list[j] = std::move(Mcoef);
      }
    } // end column loop
    
    // Return
    Rcpp::List out = Rcpp::List::create(
      _["coef"]    = coef_obs,
      _["stat"]    = stat_obs,
      _["p_value"] = pvals
    );
    if (return_residuals) out["residuals"] = resid_out;
    if (return_sampled_fits) {
      Rcpp::List L(m);
      for (arma::uword j = 0; j < m; ++j) L[j] = Rcpp::wrap(sampled_list[j]);
      out["sampled_fits"] = L;
    }
    return out;
}

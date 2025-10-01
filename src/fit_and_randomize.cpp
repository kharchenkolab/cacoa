// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]


/*
 fit_and_randomize: fast OLS / robust (winsor / Huber IRLS) with permutation z-scores
 
 WHAT IT DOES
 ------------
 Fits each column of Y on a design matrix X, computes a linear contrast of coefficients,
 generates a permutation null (full or blocked), and returns:
 • observed coefficients (coef) and observed contrast (stat)
 • z-score (z_score) derived from permutation p-value with a direction reflecting which tail
 of the permutation distribution the observed statistic falls in (above/below median)
 • optional residuals, per-permutation coefficient matrices, and per-permutation statistics
 
 ROBUSTNESS
 ----------
 robust = "none"   : OLS (fastest). Permutations use a "stat-only" path (no refit) unless
 return_sampled_fits=TRUE (then refit to collect coefficients).
 robust = "winsor" : One extra pass: clip residuals at ±k·MAD and refit OLS on adjusted y.
 Much faster than IRLS; often close to Huber.
 robust = "huber"  : IRLS with WLS inner solves; slower but fully robust.
 
 NA HANDLING
 -----------
 na_mode = "drop"        : Drop NA rows per Y column (fastest, exact).
 na_mode = "impute_weak" : Impute missing y to mean/zero and give them tiny weight (na_weight);
 only observed rows are permuted in this mode.
 
 PERMUTATIONS
 ------------
 perm_groups = NULL => full permutations.
 Otherwise a list of 1-based integer vectors; permutations are performed **within** each group.
 For na_mode = "impute_weak", only observed entries of y are permuted (within their groups).
 
 PARALLELISM
 -----------
 OpenMP across Y columns; set n_cores > 1 to enable. n_cores = 1 uses R RNG for permutations
 and avoids OpenMP for maximum compatibility.
 
 RETURNS (see the function signature doc below):
 - coef (p x m), stat (m), z_score (m)
 - residuals (n x m) if return_residuals
 - sampled_fits: list of (n_randomizations x p) matrices if return_sampled_fits
 - sampled_stats: (n_randomizations x m) if return_sampled_stats
 
 NOTE
 ----
 • We use "add-one" permutation p-value counts, so p is never 0 or 1.
 • Z-score is derived from p via qnorm in the direction specified above and alternative.
 */

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

#include <random>
#include <limits>
#include <cstdint>
#include <unordered_map>
#include <sstream>
#include <string>

/*** ============================= RNG helpers ============================= ***/

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

/*** =========================== Group mappers ============================ ***/

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

/*** ====================== Linear algebra small utils ===================== ***/

static inline bool symmetrize_and_check(arma::mat& A) {
  A = arma::symmatu(A);
  return A.is_finite();
}
static inline arma::mat pinv_safe(const arma::mat& A, double tol) {
  arma::mat P;
  if (tol > 0.0) P = arma::pinv(A, tol);
  else           P = arma::pinv(A);
  return P;
}
static inline bool inv_xtx(const arma::mat& X, arma::mat& invXtX, arma::mat& Xt, double pinv_tol=0.0) {
  Xt = X.t();
  arma::mat XtX = Xt * X;
  if (!symmetrize_and_check(XtX)) return false;
  if (inv_sympd(invXtX, XtX)) return true;
  invXtX = pinv_safe(XtX, pinv_tol);
  return invXtX.is_finite();
}
static inline arma::vec coef_from_inv(const arma::mat& invXtX, const arma::mat& Xt, const arma::vec& y) {
  return invXtX * (Xt * y);
}
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
static arma::vec wls_solve_ridge(const arma::mat& X,
                                 const arma::vec& y,
                                 const arma::vec& w,
                                 double ridge_eps = 1e-10,
                                 double pinv_tol  = 0.0) {
  arma::mat XtWX = X.t() * (X.each_col() % w);
  arma::vec XtWy = X.t() * (y % w);
  if (!symmetrize_and_check(XtWX)) { arma::vec out(X.n_cols); out.fill(arma::datum::nan); return out; }
  
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
static arma::vec huber_irls(const arma::mat& X, const arma::vec& y,
                            double k, int maxit, double tol,
                            const arma::mat* invXtX_opt = nullptr,
                            const arma::mat* Xt_opt     = nullptr,
                            const arma::vec* base_w     = nullptr,
                            double pinv_tol = 0.0) {
  const arma::uword p = X.n_cols;
  arma::vec beta(p, arma::fill::zeros);
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
static inline bool is_ill_conditioned(const arma::mat& X, double rcond_thresh) {
  const arma::uword k = X.n_rows, p = X.n_cols;
  if (k < p) return true;
  arma::mat XtX = X.t() * X;
  if (!symmetrize_and_check(XtX)) return true;
  double rc = arma::rcond(XtX);
  if (!(rc > 0.0)) return true;
  return (rc < rcond_thresh);
}
static inline std::mt19937_64 make_rng_for_column(std::uint64_t base, arma::uword j) {
  std::uint64_t x = base ^ (0x9e3779b97f4a7c15ULL + j + (j<<6) + (j>>2));
  return std::mt19937_64(x);
}

/*** ============== Residual winsorization (fast robust) ================== ***/

static inline arma::vec winsor_fit_unweighted(const arma::mat& B,
                                              const arma::mat& X,
                                              const arma::vec& y,
                                              double k) {
  arma::vec beta0 = B * y;
  arma::vec r     = y - X * beta0;
  double s = robust_scale_mad_safe(r);
  if (!(s > 0.0)) return beta0;
  const double ks = k * s;
  for (arma::uword i = 0; i < r.n_elem; ++i) {
    double v = r[i];
    if (v >  ks) r[i] =  ks;
    else if (v < -ks) r[i] = -ks;
  }
  arma::vec ytil = X * beta0 + r;
  return B * ytil;
}
static inline arma::vec winsor_fit_weighted(const arma::mat& Bw,
                                            const arma::mat& X,
                                            const arma::vec& y,
                                            const arma::uvec& idx_obs,
                                            const arma::vec& sqrtw,
                                            double k) {
  arma::vec beta0 = Bw * (sqrtw % y);
  arma::vec r     = y - X * beta0;
  
  arma::vec r_obs = r.elem(idx_obs);
  double s = robust_scale_mad_safe(r_obs);
  if (!(s > 0.0)) return beta0;
  
  const double ks = k * s;
  for (arma::uword t = 0; t < idx_obs.n_elem; ++t) {
    arma::uword i = idx_obs[t];
    double v = r[i];
    if (v >  ks) r[i] =  ks;
    else if (v < -ks) r[i] = -ks;
  }
  arma::vec ytil = X * beta0 + r;
  return Bw * (sqrtw % ytil);
}

/*** ======================== Perm tally + z-score ======================== ***/

static inline void tally_perm(double obs, double perm, int alt, int& ge, int& le, int& ge_abs) {
  if (alt == 0) { if (std::fabs(perm) >= std::fabs(obs)) ge_abs++; }
  else if (alt == 1) { if (perm >= obs) ge++; }
  else { if (perm <= obs) le++; }
}
// Replace the old helper with this version
static inline double z_from_p(double p, int alt, double obs, double med_perm) {
  // Handle boundary / invalid p explicitly
  if (!std::isfinite(p)) return NA_REAL;
  
  // If add-one p-value hits the upper boundary (no evidence in the tested tail),
  // return 0 instead of NA/Inf.
  if (p >= 1.0) return 0.0;
  
  // With add-one, p should never be 0, but guard anyway:
  if (p <= 0.0) {
    // Map to a very large finite z with the correct direction.
    // (Users almost never hit this with add-one.)
    const double p_eff = 1e-16;
    if (alt == 0) {
      double zmag = R::qnorm(1.0 - p_eff/2.0, 0.0, 1.0, 1, 0);
      double sgn  = (obs >= med_perm) ? 1.0 : -1.0;
      return sgn * zmag;
    } else if (alt == 1) {
      return R::qnorm(1.0 - p_eff, 0.0, 1.0, 1, 0);
    } else {
      return -R::qnorm(1.0 - p_eff, 0.0, 1.0, 1, 0);
    }
  }
  
  // Regular case: 0 < p < 1
  if (alt == 0) { // two-sided: sign by location vs permutation median
    double zmag = R::qnorm(1.0 - p/2.0, 0.0, 1.0, 1, 0);
    double sgn  = (obs >= med_perm) ? 1.0 : -1.0;
    return sgn * zmag;
  } else if (alt == 1) { // greater: upper-tail
    return R::qnorm(1.0 - p, 0.0, 1.0, 1, 0);
  } else {               // less: lower-tail
    return -R::qnorm(1.0 - p, 0.0, 1.0, 1, 0);
  }
}


/*** =========================== Main entry =============================== ***/

/*
 * fit_and_randomize(
 *   X, Y, contrast,
 *   perm_groups = NULL, n_randomizations = 100, alternative = "two-sided",
 *   return_residuals = TRUE, return_sampled_fits = FALSE, return_sampled_stats = FALSE,
 *   robust = "none", huber_k = 1.345, huber_maxit = 8, huber_tol = 1e-6,
 *   na_mode = "drop", na_weight = 1e-4, na_center = "mean",
 *   illcond_rcond = 1e-12, pinv_tol = 0.0, n_cores = 1
 * )
 *
 * INPUTS
 *  - X          : (n x p) double; design matrix (include intercept if desired).
 *  - Y          : (n x m) double; responses (columns). NAs allowed.
 *  - contrast   : (p) double; linear contrast vector over coefficients.
 *  - perm_groups: NULL => full permutations; else list of integer vectors (1-based row sets).
 *  - n_randomizations: number of permutations for each column (>= 0).
 *  - alternative: "two-sided" | "greater" | "less" (defines counting and z-score direction).
 *  - return_residuals    : if TRUE, return n x m residual matrix (NaN at dropped/NA rows).
 *  - return_sampled_fits : if TRUE, return list(m) of (n_randomizations x p) permuted betas.
 *  - return_sampled_stats: if TRUE, return (n_randomizations x m) matrix of permutation stats
 *                          used to form the z-scores (permutation analog of the contrast).
 *  - robust     : "none" | "winsor" | "huber".
 *  - huber_*    : IRLS parameters (used by "huber"; k is also used as winsor clamp).
 *  - na_mode    : "drop" (drop NA rows) | "impute_weak" (tiny weights on imputed values).
 *  - na_weight  : weight for imputed points (e.g., 1e-4).
 *  - na_center  : "mean" or "zero" imputation center (only used for impute_weak).
 *  - illcond_rcond: rcond threshold to treat X'X as ill-conditioned in the drop-NA branch.
 *  - pinv_tol   : tolerance for pinv fallback (0 => Armadillo default).
 *  - n_cores    : OpenMP threads across Y columns (1 => no OpenMP, serial & R RNG).
 *
 * RETURNS
 *  List with components:
 *   - coef    : (p x m) matrix of observed coefficients (NaN if column invalid).
 *   - stat    : (m) vector of observed contrast (contrast' * coef).
 *   - z_score : (m) vector; z-score derived from permutation p-value. For "two-sided",
 *               the sign indicates whether obs is above/below the permutation median.
 *               For "greater", positive z => obs in upper tail; for "less", negative z
 *               => obs in lower tail.
 *  - p_value : (m) vector; permutation p-value (add-one, never 0 or 1).
 *   - residuals     : (n x m) matrix if return_residuals=TRUE (NaN at NA/dropped rows).
 *   - sampled_fits  : list(m) of (n_randomizations x p) matrices if return_sampled_fits=TRUE.
 *   - sampled_stats : (n_randomizations x m) matrix if return_sampled_stats=TRUE.
 */
// [[Rcpp::export]]
Rcpp::List fit_and_randomize(const arma::mat& X,
                             const arma::mat& Y,
                             const arma::vec& contrast,
                             Rcpp::Nullable<Rcpp::List> perm_groups = R_NilValue,
                             int n_randomizations = 100,
                             std::string alternative = "two-sided",
                             bool return_residuals = true,
                             bool return_sampled_fits = false,
                             bool return_sampled_stats = false,
                             std::string robust = "none",
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
  if      (alternative == "two-sided") alt = 0;
  else if (alternative == "greater")   alt = 1;
  else if (alternative == "less")      alt = 2;
  else stop("alternative must be 'two-sided','greater','less'.");
  
  const bool rob_huber   = (robust == "huber");
  const bool rob_winsor  = (robust == "winsor");
  const bool impute_mode = (na_mode == "impute_weak");
  const bool center_mean = (na_center == "mean");
  
  // Precompute OLS factors for full X (used by OLS, winsor, and Huber init when no NA)
  arma::mat invXtX_full, Xt_full;
  if (!inv_xtx(X, invXtX_full, Xt_full, pinv_tol)) {
    stop("Design X is ill-conditioned or non-finite in full data.");
  }
  arma::mat B_full = invXtX_full * Xt_full;  // p x n
  
  // Permutation blocks
  List perm_full = perm_groups.isNotNull() ? List(perm_groups) : List();
  const bool full_permute_requested = (perm_full.size() == 0);
  std::vector<std::vector<arma::uword>> groups_full0;
  if (!full_permute_requested) groups_full0 = build_groups0_full(perm_full, n);
  
  // Allocate outputs
  arma::mat coef_obs(p, m); coef_obs.fill(arma::datum::nan);
  arma::vec stat_obs(m);    stat_obs.fill(arma::datum::nan);
  arma::vec zscore(m);      zscore.fill(arma::datum::nan);
  arma::vec pvalue(m);      pvalue.fill(arma::datum::nan);
  
  arma::mat resid_out;
  if (return_residuals) { resid_out.set_size(n, m); resid_out.fill(arma::datum::nan); }
  
  std::vector<arma::mat> sampled_list;
  if (return_sampled_fits) sampled_list.resize(m);
  
  arma::mat sampled_stats_out;
  if (return_sampled_stats && n_randomizations > 0) {
    sampled_stats_out.set_size(n_randomizations, m);
    sampled_stats_out.fill(arma::datum::nan);
  }
  
  // Parallel seeds
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
      
      std::mt19937_64 rng;
      if (parallel_mode) rng = make_rng_for_column(base_seed, j);
      
      const bool need_beta_perm = (rob_huber || rob_winsor || return_sampled_fits);
      const bool stat_only      = !need_beta_perm;
      
      arma::vec yj = Y.col(j);
      arma::uvec idx_fin = arma::find_finite(yj);
      bool all_finite = (idx_fin.n_elem == n);
      
      arma::vec yc; arma::uword k = n; arma::mat Xc; const arma::mat* Xobs = &X;
      arma::mat invXtX, Xt, B; arma::vec alpha;
      std::vector<std::vector<arma::uword>> groups_sub;
      arma::uvec eligible_idx;
      
      // ------- NA handling branches -------
      if (all_finite) {
        yc = yj; invXtX = invXtX_full; Xt = Xt_full; B = B_full; Xobs = &X;
        if (!full_permute_requested) groups_sub = groups_full0;
        if (stat_only) alpha = B.t() * contrast;
        
      } else if (!impute_mode || na_weight <= 0.0) {
        if (idx_fin.n_elem == 0) continue;
        Xc = X.rows(idx_fin); yc = yj.elem(idx_fin); k = Xc.n_rows; Xobs = &Xc;
        
        if (is_ill_conditioned(Xc, illcond_rcond)) continue;
        if (!inv_xtx(Xc, invXtX, Xt, pinv_tol))    continue;
        B = invXtX * Xt;
        
        if (!full_permute_requested) groups_sub = map_groups_full_to_subset(groups_full0, idx_fin, n);
        else                         eligible_idx = arma::regspace<uvec>(0, k - 1);
        
        if (stat_only) alpha = B.t() * contrast;
        
      } else {
        // impute_weak
        yc = yj;
        double mu = center_mean ? mean_finite(yj) : 0.0;
        arma::vec w_base(n); w_base.fill(na_weight);
        for (arma::uword t = 0; t < idx_fin.n_elem; ++t) w_base[idx_fin[t]] = 1.0;
        arma::vec sqrtw = arma::sqrt(w_base);
        
        if (idx_fin.n_elem < n) { yc.fill(mu); for (arma::uword t = 0; t < idx_fin.n_elem; ++t) yc[idx_fin[t]] = yj[idx_fin[t]]; }
        
        arma::mat Xw = X; Xw.each_col() %= sqrtw;
        arma::mat Xt_w = Xw.t(), invXtXw;
        if (!inv_xtx(Xw, invXtXw, Xt_w, pinv_tol)) continue;
        arma::mat Bw = invXtXw * Xt_w;  // p x n
        
        // observed fit
        if (rob_huber) {
          arma::uword neff = arma::accu(w_base > 1e-8);
          if (neff < p + 1) continue;
          arma::vec beta = huber_irls(X, yc, huber_k, huber_maxit, huber_tol, nullptr, nullptr, &w_base, pinv_tol);
          if (!beta.is_finite()) continue;
          
          double sobs = dot(beta, contrast);
          coef_obs.col(j) = beta; stat_obs[j] = sobs;
          
          if (return_residuals) {
            arma::vec rc = yc - X * beta;
            arma::vec colj_full(n); colj_full.fill(arma::datum::nan);
            if (idx_fin.n_elem > 0) colj_full.elem(idx_fin) = rc.elem(idx_fin);
            resid_out.col(j) = colj_full;
          }
          
          if (n_randomizations == 0) continue;
          
          // build groups restricted to observed rows
          bool can_permute = false;
          arma::uvec obs_idx = idx_fin;
          std::vector<std::vector<arma::uword>> groups_filt;
          if (full_permute_requested) {
            can_permute = (obs_idx.n_elem >= 2);
          } else {
            groups_filt.reserve(groups_full0.size());
            std::vector<char> elig(n, 0);
            for (arma::uword t = 0; t < obs_idx.n_elem; ++t) elig[(std::size_t)obs_idx[t]] = 1;
            for (const auto& g : groups_full0) {
              std::vector<arma::uword> h; h.reserve(g.size());
              for (arma::uword full_idx : g) if (elig[(std::size_t)full_idx]) h.push_back(full_idx);
              if (h.size() >= 2) can_permute = true;
              groups_filt.push_back(std::move(h));
            }
          }
          if (!can_permute) { zscore[j] = NA_REAL; pvalue[j] = NA_REAL; continue; }
          
          int ge=0, le=0, ge_abs=0;
          arma::mat Mcoef; if (return_sampled_fits) Mcoef.set_size(n_randomizations, p);
          arma::vec stat_perm(n_randomizations); stat_perm.fill(arma::datum::nan);
          arma::vec yperm = yc;
          bool failed=false;
          
          for (int r=0; r<n_randomizations; ++r) {
            yperm = yc;
            if (full_permute_requested) {
              if (parallel_mode) shuffle_selected_in_place_rng(yperm, obs_idx, rng);
              else               shuffle_selected_in_place_R (yperm, obs_idx);
            } else {
              if (parallel_mode) permute_by_groups_rng(yperm, groups_filt, rng);
              else               permute_by_groups_R  (yperm, groups_filt);
            }
            arma::vec bperm = huber_irls(X, yperm, huber_k, huber_maxit, huber_tol, nullptr, nullptr, &w_base, pinv_tol);
            if (!bperm.is_finite()) { failed = true; break; }
            double sperm = dot(bperm, contrast);
            stat_perm[r] = sperm;
            tally_perm(sobs, sperm, alt, ge, le, ge_abs);
            if (return_sampled_fits) Mcoef.row(r) = bperm.t();
          }
          if (failed) { zscore[j] = NA_REAL; pvalue[j] = NA_REAL; }
          else {
            double pval = (alt==0) ? (ge_abs+1.0)/(n_randomizations+1.0)
              : (alt==1) ? (ge+1.0)/(n_randomizations+1.0)
              : (le+1.0)/(n_randomizations+1.0);
            double med = arma::median(stat_perm.elem(find_finite(stat_perm)));
            zscore[j] = z_from_p(pval, alt, sobs, med);
            pvalue[j] = pval;
            if (return_sampled_fits) sampled_list[j] = std::move(Mcoef);
            if (return_sampled_stats) sampled_stats_out.col(j) = stat_perm;
          }
          continue;
        }
        
        // winsor or OLS in impute_weak
        arma::vec ycw = sqrtw % yc;
        arma::vec beta = rob_winsor ? winsor_fit_weighted(Bw, X, yc, idx_fin, sqrtw, huber_k)
          : Bw * ycw;
        if (!beta.is_finite()) continue;
        
        double sobs = dot(beta, contrast);
        coef_obs.col(j) = beta; stat_obs[j] = sobs;
        
        if (return_residuals) {
          arma::vec rc = yc - X * beta;
          arma::vec colj_full(n); colj_full.fill(arma::datum::nan);
          if (idx_fin.n_elem > 0) colj_full.elem(idx_fin) = rc.elem(idx_fin);
          resid_out.col(j) = colj_full;
        }
        
        if (n_randomizations == 0) continue;
        
        if (idx_fin.n_elem < 2) { zscore[j] = NA_REAL; pvalue[j] = NA_REAL; continue; }
        
        int ge=0, le=0, ge_abs=0;
        arma::vec yperm = yc;
        arma::mat Mcoef; if (return_sampled_fits) Mcoef.set_size(n_randomizations, p);
        arma::vec stat_perm(n_randomizations); stat_perm.fill(arma::datum::nan);
        
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
          arma::vec bperm = rob_winsor ? winsor_fit_weighted(Bw, X, yperm, idx_fin, sqrtw, huber_k)
            : Bw * (sqrtw % yperm);
          double sperm = dot(bperm, contrast);
          stat_perm[r] = sperm;
          tally_perm(sobs, sperm, alt, ge, le, ge_abs);
          if (return_sampled_fits) Mcoef.row(r) = bperm.t();
        }
        {
          double pval = (alt==0) ? (ge_abs+1.0)/(n_randomizations+1.0)
            : (alt==1) ? (ge+1.0)/(n_randomizations+1.0)
            : (le+1.0)/(n_randomizations+1.0);
          double med = arma::median(stat_perm.elem(find_finite(stat_perm)));
          pvalue[j] = pval;
          zscore[j] = z_from_p(pval, alt, sobs, med);
          if (return_sampled_fits) sampled_list[j] = std::move(Mcoef);
          if (return_sampled_stats) sampled_stats_out.col(j) = stat_perm;
        }
        continue;
      } // end impute_weak
      
      // ---------- Observed fit (no impute) ----------
      arma::vec beta;
      if      (rob_huber)  beta = huber_irls(*Xobs, yc, huber_k, huber_maxit, huber_tol, &invXtX, &Xt, nullptr, pinv_tol);
      else if (rob_winsor) beta = winsor_fit_unweighted(B, *Xobs, yc, huber_k);
      else                 beta = coef_from_inv(invXtX, Xt, yc);
      
      if (!beta.is_finite()) continue;
      
      double sobs = dot(beta, contrast);
      coef_obs.col(j) = beta; stat_obs[j] = sobs;
      
      if (return_residuals) {
        if (all_finite) {
          resid_out.col(j) = yc - (*Xobs) * beta;
        } else {
          arma::vec rc = yc - (*Xobs) * beta; // length k
          arma::vec colj_full(n); colj_full.fill(arma::datum::nan);
          if (idx_fin.n_elem == rc.n_elem) {
            colj_full.elem(idx_fin) = rc;
          } else {
            for (arma::uword t = 0; t < idx_fin.n_elem && t < rc.n_elem; ++t)
              colj_full[idx_fin[t]] = rc[t];
          }
          resid_out.col(j) = colj_full;
        }
      }
      
      if (n_randomizations == 0) continue;
      
      bool can_permute = false;
      if (all_finite) {
        if (full_permute_requested)  can_permute = (n >= 2);
        else for (const auto& g : groups_sub) if (g.size() >= 2) { can_permute = true; break; }
      } else {
        if (full_permute_requested)  can_permute = (k >= 2);
        else for (const auto& g : groups_sub) if (g.size() >= 2) { can_permute = true; break; }
      }
      if (!can_permute) { zscore[j] = NA_REAL; pvalue[j] = NA_REAL; continue; }
      
      int ge=0, le=0, ge_abs=0;
      arma::vec yperm = yc;
      arma::mat Mcoef; if (return_sampled_fits) Mcoef.set_size(n_randomizations, p);
      arma::vec stat_perm(n_randomizations); stat_perm.fill(arma::datum::nan);
      
      if (!rob_huber && !rob_winsor && !return_sampled_fits) {
        alpha = B.t() * contrast; // sperm = alpha' yperm
      }

      bool failed = false;
      
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
        
        if (rob_huber) {
          arma::vec bperm = huber_irls(*Xobs, yperm, huber_k, huber_maxit, huber_tol, &invXtX, &Xt, nullptr, pinv_tol);
          if (!bperm.is_finite()) { failed = true; break; }
          double sperm = dot(bperm, contrast);
          stat_perm[r] = sperm;
          tally_perm(sobs, sperm, alt, ge, le, ge_abs);
          if (return_sampled_fits) Mcoef.row(r) = bperm.t();
        } else if (rob_winsor) {
          arma::vec bperm = winsor_fit_unweighted(B, *Xobs, yperm, huber_k);
          double sperm = dot(bperm, contrast);
          stat_perm[r] = sperm;
          tally_perm(sobs, sperm, alt, ge, le, ge_abs);
          if (return_sampled_fits) Mcoef.row(r) = bperm.t();
        } else {
          if (return_sampled_fits) {
            arma::vec bperm = B * yperm;
            double sperm = dot(bperm, contrast);
            stat_perm[r] = sperm;
            tally_perm(sobs, sperm, alt, ge, le, ge_abs);
            Mcoef.row(r) = bperm.t();
          } else {
            double sperm = dot(alpha, yperm);
            stat_perm[r] = sperm;
            tally_perm(sobs, sperm, alt, ge, le, ge_abs);
          }
        }
      }

      if (rob_huber && failed) { 
        zscore[j] = NA_REAL; 
        pvalue[j] = NA_REAL;   
        continue;
      }
      
      double pval = (alt==0) ? (ge_abs+1.0)/(n_randomizations+1.0)
        : (alt==1) ? (ge+1.0)/(n_randomizations+1.0)
        : (le+1.0)/(n_randomizations+1.0);
      double med = arma::median(stat_perm.elem(find_finite(stat_perm)));
      pvalue[j] = pval;
      zscore[j] = z_from_p(pval, alt, sobs, med);
      if (return_sampled_fits)  sampled_list[j] = std::move(Mcoef);
      if (return_sampled_stats) sampled_stats_out.col(j) = stat_perm;
    } // end Y columns
    
    // Assemble return
    Rcpp::List out = Rcpp::List::create(
      _["coef"]    = coef_obs,
      _["stat"]    = stat_obs,
      _["z_score"] = zscore,
      _["p_value"] = pvalue
      
    );
    if (return_residuals)     out["residuals"]     = resid_out;
    if (return_sampled_fits) {
      Rcpp::List L(m);
      for (arma::uword j = 0; j < m; ++j) L[j] = Rcpp::wrap(sampled_list[j]);
      out["sampled_fits"] = L;
    }
    if (return_sampled_stats && n_randomizations > 0) out["sampled_stats"] = sampled_stats_out;
    
    return out;
}

// --- FL helpers ---

// Residualize columns of M on Z using an SVD projector (rank aware)
// Returns R = (I - P_Z) M, where P_Z projects onto col(Z)
// Robust when Z is rank-deficient or ill-conditioned
static inline arma::mat residualize_unweighted(const arma::mat& M, const arma::mat& Z, double rtol = 1e-12) {
  if (Z.n_elem == 0) return M;
  arma::mat U, V; arma::vec s;
  if (!arma::svd_econ(U, s, V, Z, 'l')) return M;
  double tol = rtol * (s.n_elem ? s.max() : 1.0);
  arma::uvec keep = arma::find(s > tol);
  if (keep.n_elem == 0) return M;
  arma::mat Ur = U.cols(keep);
  return M - Ur * (Ur.t() * M);
}

// Weighted residuals via SVD on Zw = diag(sw) Z
// y_imp has NAs weak-imputed already; sw = sqrt(weights); weight=1 obs, =sqrt(na_weight) on NAs
static inline arma::vec residualize_weighted_svd(const arma::mat& Z,
                                                 const arma::vec& y_imp,
                                                 const arma::vec& sw,
                                                 double rtol = 1e-12)
{
  if (Z.n_elem == 0) return y_imp;
  arma::mat Zw = Z; Zw.each_col() %= sw;          // Zw = diag(sw) Z
  arma::vec yw = sw % y_imp;                      // yw = diag(sw) y
  arma::mat U, V; arma::vec s;
  if (!arma::svd_econ(U, s, V, Zw, 'l')) return y_imp;

  const double tol = rtol * (s.n_elem ? s.max() : 1.0);
  arma::uvec keep = arma::find(s > tol);
  if (keep.n_elem == 0) return y_imp;

  arma::mat Ur = U.cols(keep);
  arma::mat Vr = V.cols(keep);
  arma::vec  si = 1.0 / s.elem(keep);

  // Minimum-norm WLS solution: beta = Vr * diag(si) * Ur^T * yw
  arma::vec beta = Vr * (si % (Ur.t() * yw));

  // residuals on original scale
  return y_imp - Z * beta;
}

// Residualize a single y on Z using weak imputation: closer to R's qr.resid()
//   • NAs in y are replaced with mu (mean or zero),
//   • weighted LS with weight=1 for observed rows, na_weight for NA rows,
//   • then put back NA markers in the residual vector so downstream code
//     knows which rows were originally missing.
static inline arma::mat residualize_weighted_column_fullZ(const arma::mat& Z,
                                                          const arma::vec& y,
                                                          double na_weight,
                                                          bool center_mean,
                                                          double rtol = 1e-12)
{
  arma::vec y_imp = y;
  arma::uvec fin = arma::find_finite(y);
  const double mu = (center_mean && fin.n_elem>0) ? arma::mean(y.elem(fin)) : 0.0;
  y_imp.elem(arma::find_nonfinite(y)).fill(mu);

  if (Z.n_elem == 0) {
    arma::vec r = y_imp;
    r.elem(arma::find_nonfinite(y)).fill(arma::datum::nan);
    return arma::mat(r);
  }
  // weights: 1 for observed, tiny for NA rows
  arma::vec w(y.n_elem, arma::fill::ones);
  w.elem(arma::find_nonfinite(y)).fill(na_weight);
  arma::vec sw = arma::sqrt(w);

  arma::vec r = residualize_weighted_svd(Z, y_imp, sw, rtol);
  r.elem(arma::find_nonfinite(y)).fill(arma::datum::nan); // restore NA markers
  return arma::mat(r);
}

// key for a row mask; group Y columns by identical NA pattern
static inline std::string mask_key(const arma::uvec& keep_idx, arma::uword n) {
  std::ostringstream oss; oss << n << ':';
  for (arma::uword t=0; t<keep_idx.n_elem; ++t) { if (t) oss << ','; oss << keep_idx[t]; }
  return oss.str();
}

// Remap core-groups to a reduced core-space after rows are dropped.
// groups_core: list of 1-based indices in ORIGINAL core space
// keep_core_pos: 0-based positions (in original core space) that survive for this subset
static inline Rcpp::List subset_core_groups_1based(const Rcpp::List& groups_core,
                                                   const arma::uvec& keep_core_pos,
                                                   arma::uword n_core_orig) {
  std::vector<int> map_old_to_new(n_core_orig, -1);
  for (arma::uword i=0; i<keep_core_pos.n_elem; ++i) map_old_to_new[ keep_core_pos[i] ] = (int)i;

  Rcpp::List out(groups_core.size());
  for (int g=0; g<groups_core.size(); ++g) {
    Rcpp::IntegerVector grp = groups_core[g]; // 1-based indices
    std::vector<int> kept; kept.reserve(grp.size());
    for (int v : grp) {
      int old0 = v - 1;
      if (old0 >= 0 && old0 < (int)n_core_orig) {
        int nw = map_old_to_new[old0]; // new position or -1
        if (nw >= 0) kept.push_back(nw + 1); // back to 1-based
      }
    }
    out[g] = Rcpp::IntegerVector(kept.begin(), kept.end());
  }
  return out;
}

// ---------- Main entry: preprocess FL, then call fit_and_randomize ----------

//' Freedman–Lane driver with NA-aware residualization and core-space batching
//'
//' @description
//' `cpp_fl()` prepares inputs for permutation inference under the
//' Freedman–Lane scheme and then calls `fit_and_randomize()`.
//' It handles missingness per response column, reuses factorizations,
//' and respects core-only fitting and permutation blocks.
//'
//' @details
//' Workflow:
//' 1) Residualizes `X` against `Z` on the full data once via QR.
//' 2) Depending on `na_mode`:
//'    - `"drop"`: splits `Y` columns into those with and without `NA`.
//'       * No-NA columns: residualize against full `Z` and fit in one batch (using core rows).
//'       * NA columns: group columns by identical row mask; for each group,
//'         re-residualize `X` and the group's `Y` columns on the reduced set of rows,
//'         intersect with `core_rows`, remap `perm_groups_core` into that reduced
//'         core space, and fit.
//'    - `"impute_weak"`: for each `Y` column, compute `Z`-residuals using weighted LS
//'         (weight 1 for observed rows and `na_weight` for NA rows), then restore the
//'         original `NA` flags in the residuals. After slicing to `core_rows`, delegate
//'         to `fit_and_randomize(…, na_mode = "impute_weak")`, which may apply its own
//'         tiny-weight logic and permute only observed rows.
//'
//' `perm_groups_core` must be a list of integer vectors containing 1-based positions
//' in the *core* row space (as produced by `permutationGroups(...)[["core"]]`).
//' If `perm_groups_core` is an empty list, a full permutation of eligible rows is used.
//'
//' @param X Numeric matrix (n × p). Core regressors including the contrasted term(s).
//' @param Z Numeric matrix (n × k) of nuisance regressors; use a 0×0 matrix if none.
//' @param Y Numeric matrix (n × m) of responses (each column a response).
//' @param contrast Numeric vector (length p) specifying the contrast over columns of `X`.
//' @param core_rows_opt Logical vector (length n) selecting rows to keep for fitting
//'   (the “core” rows). If `NULL`, all rows are used.
//' @param perm_groups_core_opt List of integer vectors (1-based indices) defining
//'   permutation blocks in the core row space.
//' @param n_randomizations Integer; number of random permutations.
//' @param alternative Character: `"two-sided"`, `"greater"`, or `"less"`.
//' @param return_residuals Logical; include residuals in the result.
//' @param return_sampled_stats Logical; include permutation statistics matrix.
//' @param robust Character: `"none"`, `"huber"`, or `"winsor"`; passed through.
//' @param huber_k,huber_maxit,huber_tol Robust-fitting parameters; passed through.
//' @param na_mode Character: `"drop"` or `"impute_weak"`. See Details.
//' @param na_weight Double; tiny weight for NA rows when `na_mode = "impute_weak"`.
//' @param na_center Character: centering for constructing imputed values in the weighted
//'   residualization step; `"mean"` or `"zero"`.
//' @param illcond_rcond,pinv_tol,n_cores Numeric tuning and threading; passed through.
//'
//' @returns A list with:
//' * `stat`       — numeric vector (length m): observed contrast statistics per Y column.
//' * `p_value`    — numeric vector (length m).
//' * `z_score`    — numeric vector (length m).
//' * `residuals`  — numeric matrix (|core_rows| × m) if requested; residuals are `NA` for rows
//'                   not used under `"drop"`.
//' * `sampled_stats` — numeric matrix (`n_randomizations` × m) if requested.
//'
//' @section Invariants & errors:
//' * Requires `nrow(X) == nrow(Y)` and `length(contrast) == ncol(X)`.
//' * For `"drop"`, columns with `< 2` eligible rows in the (reduced) core space are skipped;
//'   their outputs are `NA`.
//' * Column space of `X` is not trimmed automatically; any rank/conditioning checks
//'   are handled inside `fit_and_randomize()`.
//'
//' @section Performance notes:
//' * Reuses one full-data residualization for `X`.
//' * Batches all no-NA columns of `Y` into a single call.
//' * Groups NA columns by identical row mask to reduce Z/X residualization and fitting.
//'
//' @examples
//' \dontrun{
//' out <- cpp_fl(X = x$X,
//'               Z = if (is.null(x$Z)) matrix(0,0,0) else x$Z,
//'               Y = Y,
//'               contrast = x$contrast.X,
//'               core_rows_opt = x$core.rows,
//'               perm_groups_core_opt = x$perm.groups.core,
//'               n_randomizations = 1000,
//'               alternative = "two-sided",
//'               na_mode = "drop",
//'               robust = "none")
//' }
// [[Rcpp::export]]
Rcpp::List cpp_fl_old(const arma::mat& X,
                  const arma::mat& Z,
                  const arma::mat& Y,
                  const arma::vec& contrast,
                  Rcpp::Nullable<Rcpp::LogicalVector> core_rows_opt = R_NilValue,
                  Rcpp::Nullable<Rcpp::List> perm_groups_core_opt   = R_NilValue, // ignored (full shuffles)
                  int n_randomizations = 1000,
                  std::string alternative = "two-sided",
                  bool return_residuals = true,
                  bool return_sampled_stats = true,
                  std::string robust = "none",
                  double huber_k = 1.345,
                  int huber_maxit = 8,
                  double huber_tol = 1e-6,
                  std::string na_mode = "drop",
                  double na_weight = 1e-4,
                  std::string na_center = "mean",
                  double illcond_rcond = 1e-12,
                  double pinv_tol = 0.0,
                  int n_cores = 1)
{

  // ---- input checks ----
  if (X.n_rows != Y.n_rows) Rcpp::stop("X and Y must have the same number of rows.");
  if (contrast.n_elem != X.n_cols) Rcpp::stop("contrast length must equal ncol(X).");

  const arma::uword n = X.n_rows, m = Y.n_cols;

  arma::uvec core_idx;
  if (core_rows_opt.isNotNull()) {
    Rcpp::LogicalVector cr(core_rows_opt);
    if ((arma::uword)cr.size() != n) Rcpp::stop("core.rows length mismatch.");
    std::vector<arma::uword> tmp; tmp.reserve(cr.size());
    for (int i = 0; i < cr.size(); ++i) if (cr[i]) tmp.push_back((arma::uword)i);
    core_idx = arma::conv_to<arma::uvec>::from(tmp);
  } else {
    core_idx = arma::regspace<arma::uvec>(0, n-1); // all rows are core
  }
  const arma::uword n_core = core_idx.n_elem;

  // map full row index -> core position (or -1 if not a core row)
  std::vector<int> pos_in_core(n, -1);
  for (arma::uword i=0; i<n_core; ++i) pos_in_core[ core_idx[i] ] = (int)i;


  // residualize X on FULL Z once 
  arma::mat Xr_full = (Z.n_elem == 0) ? X : residualize_unweighted(X, Z);

  // outputs
  arma::vec stat(m); stat.fill(arma::datum::nan);
  arma::vec pval(m); pval.fill(arma::datum::nan);
  arma::vec z(m);    z.fill(arma::datum::nan);
  arma::mat resid;
  if (return_residuals) { resid.set_size(n_core, m); resid.fill(arma::datum::nan); }
  arma::mat sampled_stats;
  if (return_sampled_stats && n_randomizations > 0) {
    sampled_stats.set_size(n_randomizations, m);
    sampled_stats.fill(arma::datum::nan);
  }

  // helper to key uvecs
  auto key_from_uvec = [](const arma::uvec& idx)->std::string {
    std::string s; s.reserve(idx.n_elem * 6);
    for (arma::uword i=0; i<idx.n_elem; ++i) { if (i) s.push_back(','); s += std::to_string(idx[i]); }
    return s;
  };

  // Empty permutation groups => full shuffles
  Rcpp::List empty_groups;

  // ========================== branch: impute_weak ==========================
  if (na_mode == "impute_weak") {
    // Residualize EVERY y-column on FULL Z using weak WLS; restore NA markers
    arma::mat Yr_full(n, m);
    const bool center_mean_flag = (na_center == "mean");
    for (arma::uword j=0; j<m; ++j)
      Yr_full.col(j) = residualize_weighted_column_fullZ(Z, Y.col(j), na_weight, center_mean_flag);

    arma::mat Xr_use = Xr_full.rows(core_idx);
    arma::mat Yr_use = Yr_full.rows(core_idx);

    // Delegate to fitter for all columns at once, with FULL SHUFFLES
    Rcpp::List fit = fit_and_randomize(
      Xr_use, Yr_use, contrast,
      empty_groups, n_randomizations, alternative,
      return_residuals, false, return_sampled_stats,
      robust, huber_k, huber_maxit, huber_tol,
      "impute_weak", na_weight, na_center,
      illcond_rcond, pinv_tol, n_cores
    );

    stat = Rcpp::as<arma::vec>(fit["stat"]);
    pval = Rcpp::as<arma::vec>(fit["p_value"]);
    z    = Rcpp::as<arma::vec>(fit["z_score"]);
    if (return_residuals)     resid  = Rcpp::as<arma::mat>(fit["residuals"]);
    if (return_sampled_stats && n_randomizations > 0)
                              sampled_stats = Rcpp::as<arma::mat>(fit["sampled_stats"]);

    Rcpp::List out = Rcpp::List::create(
      Rcpp::_["stat"]    = stat,
      Rcpp::_["p_value"] = pval,
      Rcpp::_["z_score"] = z
    );
    if (return_residuals) out["residuals"] = resid; else out["residuals"] = R_NilValue;
    if (return_sampled_stats && n_randomizations > 0) out["sampled_stats"] = sampled_stats; else out["sampled_stats"] = R_NilValue;
    return out;
  }

  // ============================ branch: drop ===============================
  // Split Y columns by NA presence (full rows)
  std::vector<arma::uword> cols_no_na, cols_has_na;
  cols_no_na.reserve(m); cols_has_na.reserve(m);
  for (arma::uword j=0; j<m; ++j) {
    arma::uvec fin = arma::find_finite(Y.col(j));
    if (fin.n_elem == n) cols_no_na.push_back(j); else cols_has_na.push_back(j);
  }

  // precompute row indices for submat writes

  arma::uvec all_resid_rows;
  if (return_residuals) all_resid_rows = arma::regspace<arma::uvec>(0, n_core-1);
  arma::uvec all_stat_rows;
  if (return_sampled_stats && n_randomizations > 0) all_stat_rows = arma::regspace<arma::uvec>(0, n_randomizations-1);

  // 1) batch fit for columns with NO NA
  if (!cols_no_na.empty()) {
    arma::uvec idx = arma::conv_to<arma::uvec>::from(cols_no_na);
    arma::mat Ysub = Y.cols(idx);
    arma::mat Yr_full = (Z.n_elem == 0) ? Ysub : residualize_unweighted(Ysub, Z);

    arma::mat Xr_use = Xr_full.rows(core_idx);
    arma::mat Yr_use = Yr_full.rows(core_idx);

    Rcpp::List fit = fit_and_randomize(
      Xr_use, Yr_use, contrast,
      empty_groups, n_randomizations, alternative,        // FULL SHUFFLES
      return_residuals, false, return_sampled_stats,
      robust, huber_k, huber_maxit, huber_tol,
      "drop", na_weight, na_center,
      illcond_rcond, pinv_tol, n_cores
    );

    arma::vec s  = Rcpp::as<arma::vec>(fit["stat"]);
    arma::vec pv = Rcpp::as<arma::vec>(fit["p_value"]);
    arma::vec zz = Rcpp::as<arma::vec>(fit["z_score"]);
    stat.elem(idx) = s; pval.elem(idx) = pv; z.elem(idx) = zz;

    if (return_residuals) {
      arma::mat R = Rcpp::as<arma::mat>(fit["residuals"]);
      resid.submat(all_resid_rows, idx) = R;
    }
    if (return_sampled_stats && n_randomizations > 0) {
      arma::mat S = Rcpp::as<arma::mat>(fit["sampled_stats"]);
      sampled_stats.submat(all_stat_rows, idx) = S;
    }
  }

  // 2) columns WITH NA: group by identical observed-row mask, then fit per group
  if (!cols_has_na.empty()) {
    std::unordered_map<std::string, std::vector<arma::uword>> groups;
    groups.reserve(cols_has_na.size());
    for (arma::uword j : cols_has_na) {
      arma::uvec keep_full = arma::find_finite(Y.col(j));
      groups[ key_from_uvec(keep_full) ].push_back(j);
    }
    for (auto &kv : groups) {
      const std::vector<arma::uword>& cols = kv.second;
      arma::uvec keep_full = arma::find_finite(Y.col(cols[0]));
      if (keep_full.n_elem == 0) continue;

      arma::mat Z_sub = (Z.n_elem == 0) ? arma::mat(0,0) : Z.rows(keep_full);
      arma::mat X_sub = X.rows(keep_full);
      arma::mat Xr_sub = (Z_sub.n_elem == 0) ? X_sub : residualize_unweighted(X_sub, Z_sub);

      arma::uvec cols_idx = arma::conv_to<arma::uvec>::from(cols);
      arma::mat Y_sub = Y.submat(keep_full, cols_idx);
      arma::mat Yr_sub = (Z_sub.n_elem == 0) ? Y_sub : residualize_unweighted(Y_sub, Z_sub);

      // intersect with core rows
      std::vector<arma::uword> keep_core_pos_vec, local_rows_vec;
      keep_core_pos_vec.reserve(keep_full.n_elem); local_rows_vec.reserve(keep_full.n_elem);
      for (arma::uword i=0; i<keep_full.n_elem; ++i) {
        int core_pos = pos_in_core[ keep_full[i] ];
        if (core_pos >= 0) { keep_core_pos_vec.push_back((arma::uword)core_pos); local_rows_vec.push_back(i); }
      }
      arma::uvec keep_core_pos = arma::conv_to<arma::uvec>::from(keep_core_pos_vec);
      arma::uvec local_rows    = arma::conv_to<arma::uvec>::from(local_rows_vec);
      if (local_rows.n_elem == 0) continue;

      arma::mat Xr_use = Xr_sub.rows(local_rows);
      arma::mat Yr_use = Yr_sub.rows(local_rows);

      // FULL SHUFFLES (ignore blocks)
      Rcpp::List fit = fit_and_randomize(
        Xr_use, Yr_use, contrast,
        empty_groups, n_randomizations, alternative,
        return_residuals, false, return_sampled_stats,
        robust, huber_k, huber_maxit, huber_tol,
        "drop", na_weight, na_center,
        illcond_rcond, pinv_tol, n_cores
      );

      arma::vec s  = Rcpp::as<arma::vec>(fit["stat"]);
      arma::vec pv = Rcpp::as<arma::vec>(fit["p_value"]);
      arma::vec zz = Rcpp::as<arma::vec>(fit["z_score"]);
      for (arma::uword k=0; k<cols_idx.n_elem; ++k) {
        stat[cols_idx[k]] = s[k];
        pval[cols_idx[k]] = pv[k];
        z[cols_idx[k]]    = zz[k];
      }

      if (return_residuals) {
        arma::mat R = Rcpp::as<arma::mat>(fit["residuals"]); // |local_rows| × |cols|
        for (arma::uword r=0; r<local_rows.n_elem; ++r) {
          arma::uword core_pos = keep_core_pos[r];
          arma::uvec rr(1); rr[0] = core_pos;
          resid.submat(rr, cols_idx) = R.row(r);
        }
      }
      if (return_sampled_stats && n_randomizations > 0) {
        arma::mat S = Rcpp::as<arma::mat>(fit["sampled_stats"]);
        sampled_stats.submat(all_stat_rows, cols_idx) = S;
      }
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["stat"]    = stat,
    Rcpp::_["p_value"] = pval,
    Rcpp::_["z_score"] = z
  );
  if (return_residuals) out["residuals"] = resid; else out["residuals"] = R_NilValue;
  if (return_sampled_stats && n_randomizations > 0) out["sampled_stats"] = sampled_stats; else out["sampled_stats"] = R_NilValue;
  return out;
}



/// FWL-optimized FL implementation 

// ---- helpers specific to fl_fwl_cpp ----------------------------------

// FNV-1a hash of a 0/1 mask (for NA pattern grouping)
static inline std::uint64_t fnv1a64_mask(const std::vector<unsigned char>& mask) {
  const std::uint64_t FNV_OFFSET = 1469598103934665603ULL;
  const std::uint64_t FNV_PRIME  = 1099511628211ULL;
  std::uint64_t h = FNV_OFFSET;
  for (unsigned char b : mask) { h ^= (std::uint64_t)b; h *= FNV_PRIME; }
  return h;
}

// Orthonormal basis of col(Z) via econ QR, truncated to numerical rank
static inline bool qr_basis(const arma::mat& Z, arma::mat& Q, double rank_tol = 1e-12) {
  if (Z.n_cols == 0) { Q.set_size(Z.n_rows, 0); return true; }
  arma::mat R;
  bool ok = arma::qr_econ(Q, R, Z);  // Q: n x k
  if (!ok || !Q.is_finite() || !R.is_finite()) { Q.reset(); return false; }
  arma::uword r = arma::rank(R, rank_tol);
  if (r == 0) { Q.set_size(Z.n_rows, 0); return true; }
  if (r < Q.n_cols) Q = Q.cols(0, r - 1);
  return true;
}

// Apply residual maker M = I - QQ' to a matrix
static inline arma::mat project_out_Q(const arma::mat& Q, const arma::mat& B) {
  if (Q.n_cols == 0) return B;
  return B - Q * (Q.t() * B);
}

// Conditioning check for core design after residualization
static inline bool is_ill_conditioned_mat(const arma::mat& Xr, double rcond_thresh) {
  if (Xr.n_rows < Xr.n_cols) return true;
  arma::mat XtX = Xr.t() * Xr;
  if (!XtX.is_finite()) return true;
  double rc = arma::rcond(arma::symmatu(XtX));
  if (!(rc > 0.0)) return true;
  return (rc < rcond_thresh);
}

// Parse core_rows (NULL | logical | integer 1-based) → arma::uvec (0-based)
static inline arma::uvec parse_core_rows(SEXP core_rows, arma::uword n) {
  if (Rf_isNull(core_rows)) return arma::regspace<arma::uvec>(0, n - 1);
  if (Rf_isLogical(core_rows)) {
    Rcpp::LogicalVector L(core_rows);
    if ((arma::uword)L.size() != n) Rcpp::stop("core_rows (logical) must have length nrow(X).");
    std::vector<arma::uword> v; v.reserve(n);
    for (int i = 0; i < L.size(); ++i) if (L[i] == TRUE) v.push_back((arma::uword)i);
    return arma::uvec(v);
  }
  Rcpp::IntegerVector I(core_rows);
  std::vector<arma::uword> v; v.reserve(I.size());
  for (int a = 0; a < I.size(); ++a) {
    int ii = I[a];
    if (ii < 1 || ii > (int)n) Rcpp::stop("core_rows indices must be in [1,n].");
    v.push_back((arma::uword)(ii - 1));
  }
  return arma::uvec(v);
}

/*** ===================================================================== ***/
/***                           fl_fwl_cpp (FWL)                           ***/
/*** ===================================================================== ***/

// [[Rcpp::export]]
Rcpp::List fl_fwl_cpp(const arma::mat& X,
                      const arma::mat& Z,
                      const arma::mat& Y,
                      const arma::vec& contrast,
                      SEXP core_rows = R_NilValue,
                      int n_randomizations = 100,
                      std::string alternative = "two-sided",
                      std::string robust = "none",
                      double huber_k = 1.345,
                      int huber_maxit = 8,
                      double huber_tol = 1e-6,
                      std::string na_mode = "drop",       // "drop" | "impute_weak"
                      double na_weight = 1e-4,
                      std::string na_center = "mean",     // used if impute_weak
                      double illcond_rcond = 1e-12,
                      double pinv_tol = 0.0,
                      int n_cores = 1,
                      bool return_residuals = true,
                      bool return_sampled_fits = false,
                      bool return_sampled_stats = false) {
  RNGScope scope;
  
  const arma::uword n = X.n_rows, p = X.n_cols;
  if (Y.n_rows != n) stop("X and Y must have the same number of rows.");
  if (contrast.n_elem != p) stop("contrast length must equal ncol(X).");
  
  const bool use_drop   = (na_mode == "drop");
  const bool use_impute = (na_mode == "impute_weak");
  if (!use_drop && !use_impute) stop("na_mode must be 'drop' or 'impute_weak'.");
  
  // core rows (0-based indices in [0,n))
  arma::uvec idx_core = parse_core_rows(core_rows, n);
  const arma::uword nc = idx_core.n_elem;
  if (nc < p + 1) stop("Not enough core rows (|core_rows| < p+1).");
  
  // map global row -> position within core (or -1 if not in core)
  std::vector<int> pos_in_core((size_t)n, -1);
  for (arma::uword t = 0; t < nc; ++t) pos_in_core[(size_t)idx_core[t]] = (int)t;
  
  const arma::uword m = Y.n_cols;
  
  // Group Y columns by identical NA mask over **all n rows**
  struct GInfo {
    std::vector<arma::uword> cols;
    std::vector<unsigned char> mask_all; // length n; 1=finite, 0=NA
  };
  std::unordered_map<std::uint64_t, GInfo> groups;
  groups.reserve((size_t)std::max<arma::uword>(8u, m / 8u));
  
  for (arma::uword j = 0; j < m; ++j) {
    std::vector<unsigned char> mask; mask.reserve(n);
    arma::vec colj = Y.col(j);
    for (arma::uword i = 0; i < n; ++i) mask.push_back( arma::is_finite(colj[i]) ? 1u : 0u );
    std::uint64_t key = fnv1a64_mask(mask);
    auto it = groups.find(key);
    if (it == groups.end()) {
      GInfo g; g.cols.push_back(j); g.mask_all = std::move(mask);
      groups.emplace(key, std::move(g));
    } else {
      it->second.cols.push_back(j);
    }
  }
  
  // Allocate outputs (core-row space for residuals)
  arma::mat Coef(p, m);   Coef.fill(arma::datum::nan);
  arma::vec Stat(m);      Stat.fill(arma::datum::nan);
  arma::vec Zscore(m);    Zscore.fill(arma::datum::nan);
  arma::vec Pval(m);      Pval.fill(arma::datum::nan);
  arma::mat Resid;        if (return_residuals) { Resid.set_size(nc, m); Resid.fill(arma::datum::nan); }
  std::vector<arma::mat> SampledFits; if (return_sampled_fits) SampledFits.resize(m);
  arma::mat SampledStats; if (return_sampled_stats && n_randomizations > 0) { SampledStats.set_size(n_randomizations, m); SampledStats.fill(arma::datum::nan); }
  
  // Process NA-pattern groups
  for (auto & kv : groups) {
    const GInfo& g = kv.second;
    const std::vector<unsigned char>& mask = g.mask_all;
    const std::vector<arma::uword>& J = g.cols;
    
    // column index vector for Armadillo submat ops
    arma::uvec Jv(J.size());
    for (size_t a = 0; a < J.size(); ++a) Jv[a] = J[a];
    
    if (use_drop) {
      // ---- 1) finite rows (over ALL rows), residualize there ----
      std::vector<arma::uword> pos_fin; pos_fin.reserve(n);
      for (arma::uword i = 0; i < n; ++i) if (mask[i]) pos_fin.push_back(i);
      arma::uword n_fin = (arma::uword)pos_fin.size();
      if (n_fin < p + 1) continue; // nowhere to fit even before core filtering
      
      arma::uvec sel_fin = arma::uvec(pos_fin);             // indices in [0,n)
      arma::mat Xf = X.rows(sel_fin);
      arma::mat Zf = (Z.n_cols > 0) ? Z.rows(sel_fin) : arma::mat(n_fin, 0);
      arma::mat Yf = Y.submat(sel_fin, Jv);                 // n_fin x |J|
      
      arma::mat Xr_fin, Yr_fin;
      if (Zf.n_cols > 0) {
        arma::mat Qz;
        if (!qr_basis(Zf, Qz, pinv_tol)) continue;          // numerical failure → skip
        Xr_fin = project_out_Q(Qz, Xf);
        Yr_fin = project_out_Q(Qz, Yf);
      } else {
        Xr_fin = Xf; Yr_fin = Yf;
      }
      
      // ---- 2) now apply core_rows: intersect(sel_fin, idx_core) ----
      std::vector<arma::uword> pos_corefin; pos_corefin.reserve(n_fin);
      std::vector<arma::uword> pos_corefin_corepos; pos_corefin_corepos.reserve(n_fin);
      // map global->position in sel_fin
      std::vector<int> pos_in_fin((size_t)n, -1);
      for (arma::uword r = 0; r < n_fin; ++r) pos_in_fin[(size_t)sel_fin[r]] = (int)r;
      
      for (arma::uword r = 0; r < n_fin; ++r) {
        arma::uword i_glob = sel_fin[r];
        int pc = pos_in_core[(size_t)i_glob];
        if (pc >= 0) {
          pos_corefin.push_back(r);                  // row index in Xr_fin / Yr_fin
          pos_corefin_corepos.push_back((arma::uword)pc); // row index in Resid (0..nc-1)
        }
      }
      if (pos_corefin.size() < (size_t)(p + 1)) continue;   // not enough rows to fit on core
      
      arma::uvec sel_corefin = arma::uvec(pos_corefin);                  // rows in Xr_fin/Yr_fin
      arma::uvec sel_corepos = arma::uvec(pos_corefin_corepos);          // rows in Resid
      
      arma::mat Xr = Xr_fin.rows(sel_corefin);
      arma::mat Yr = Yr_fin.rows(sel_corefin);
      
      if (is_ill_conditioned_mat(Xr, illcond_rcond)) continue;
      
      // Engine on residualized subset; permutations are full (perm_groups=NULL)
      Rcpp::List ans = fit_and_randomize(
        Xr, Yr, contrast,
        R_NilValue, n_randomizations, alternative,
        return_residuals, return_sampled_fits, return_sampled_stats,
        robust, huber_k, huber_maxit, huber_tol,
        "drop", na_weight, na_center,
        illcond_rcond, pinv_tol, n_cores
      );
      
      arma::mat B = ans["coef"];        // p x |J|
      arma::vec s = ans["stat"];        // |J|
      arma::vec z = ans["z_score"];     // |J|
      arma::vec pv = ans["p_value"];    // |J|
      for (size_t a = 0; a < J.size(); ++a) {
        arma::uword j = J[a];
        Coef.col(j) = B.col(a);
        Stat[j]     = s[a];
        Zscore[j]   = z[a];
        Pval[j]     = pv[a];
      }
      if (return_residuals) {
        arma::mat Rg = ans["residuals"];           // |core∩finite| x |J|
        Resid.submat(sel_corepos, Jv) = Rg;        // write back at core positions
      }
      if (return_sampled_fits) {
        Rcpp::List L = ans["sampled_fits"];
        for (size_t a = 0; a < J.size(); ++a) SampledFits[J[a]] = Rcpp::as<arma::mat>(L[a]);
      }
      if (return_sampled_stats && n_randomizations > 0) {
        arma::mat SS = ans["sampled_stats"];       // n_perm x |J|
        SampledStats.cols(Jv) = SS;
      }
      
    } else { // -------- impute_weak: residualize on ALL rows, then apply core_rows --------
      arma::vec w(n); w.fill(na_weight);
      for (arma::uword i = 0; i < n; ++i) if (mask[i]) w[i] = 1.0;
      arma::vec sqrtw = arma::sqrt(w);
      
      arma::mat Xw = X.each_col() % sqrtw;
      arma::mat Zw = (Z.n_cols > 0) ? Z.each_col() % sqrtw : arma::mat(n, 0);
      
      arma::mat Yw(n, J.size());
      for (size_t a = 0; a < J.size(); ++a) {
        arma::uword j = J[a];
        arma::vec y = Y.col(j);
        if (na_center == "mean") {
          double mu = 0.0; arma::uword cnt = 0;
          for (arma::uword i = 0; i < n; ++i) if (mask[i]) { mu += y[i]; ++cnt; }
          mu = (cnt > 0) ? mu / (double)cnt : 0.0;
          for (arma::uword i = 0; i < n; ++i) if (!mask[i]) y[i] = mu;
        } else {
          for (arma::uword i = 0; i < n; ++i) if (!mask[i]) y[i] = 0.0;
        }
        Yw.col(a) = y % sqrtw;
      }
      
      arma::mat Xr_full, Yr_full;
      if (Zw.n_cols > 0) {
        arma::mat Qw;
        if (!qr_basis(Zw, Qw, pinv_tol)) continue;
        Xr_full = project_out_Q(Qw, Xw);
        Yr_full = project_out_Q(Qw, Yw);
      } else {
        Xr_full = Xw; Yr_full = Yw;
      }
      
      // Now apply core_rows only
      arma::mat Xr = Xr_full.rows(idx_core);                // nc x p
      arma::mat Yr = Yr_full.rows(idx_core);                // nc x |J|
      
      if (is_ill_conditioned_mat(Xr, illcond_rcond)) continue;
      
      Rcpp::List ans = fit_and_randomize(
        Xr, Yr, contrast,
        R_NilValue, n_randomizations, alternative,
        return_residuals, return_sampled_fits, return_sampled_stats,
        robust, huber_k, huber_maxit, huber_tol,
        "drop", na_weight, na_center,
        illcond_rcond, pinv_tol, n_cores
      );
      
      arma::mat B = ans["coef"];        // p x |J|
      arma::vec s = ans["stat"];        // |J|
      arma::vec z = ans["z_score"];     // |J|
      arma::vec pv = ans["p_value"];    // |J|
      for (size_t a = 0; a < J.size(); ++a) {
        arma::uword j = J[a];
        Coef.col(j) = B.col(a);
        Stat[j]     = s[a];
        Zscore[j]   = z[a];
        Pval[j]     = pv[a];
      }
      if (return_residuals) {
        arma::mat Rg = ans["residuals"];            // nc x |J|
        Resid.cols(Jv) = Rg;                        // directly into core space
      }
      if (return_sampled_fits) {
        Rcpp::List L = ans["sampled_fits"];
        for (size_t a = 0; a < J.size(); ++a) SampledFits[J[a]] = Rcpp::as<arma::mat>(L[a]);
      }
      if (return_sampled_stats && n_randomizations > 0) {
        arma::mat SS = ans["sampled_stats"];        // n_perm x |J|
        SampledStats.cols(Jv) = SS;
      }
    }
  } // groups
  
  // Return
  Rcpp::List out = Rcpp::List::create(
    _["coef"]    = Coef,
    _["stat"]    = Stat,
    _["z_score"] = Zscore,
    _["p_value"] = Pval
  );
  if (return_residuals)     out["residuals"]     = Resid;
  if (return_sampled_fits) {
    Rcpp::List L(m); for (arma::uword j = 0; j < m; ++j) L[j] = Rcpp::wrap(SampledFits[j]);
    out["sampled_fits"] = L;
  }
  if (return_sampled_stats && n_randomizations > 0) out["sampled_stats"] = SampledStats;
  return out;
}



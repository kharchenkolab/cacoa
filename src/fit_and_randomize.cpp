// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * UNIFIED OPTIMIZED STATISTICS MODULE
 * ===================================
 * * This module provides high-performance tools for linear modeling with:
 * 1. Robust estimation (Huber, Winsorization).
 * 2. Complex NA handling (Drop-NA or Weighted Imputation).
 * 3. Permutation-based inference (Z-scores/P-values).
 * 4. Partial regression (Frisch-Waugh-Lovell) with nuisance covariates.
 *
 * The implementation is provided by two top-level functions:
 *  - fit_and_randomize (core function implementing fits and testing)
 *  - fl_fwl_cpp (Freedman-Lane procedure)
 */

#include <RcppArmadillo.h>
#include <omp.h>
#include <random>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>

// Disable Armadillo's internal OpenMP to avoid thread oversubscription.
// We manage threads explicitly at the column level.
#define ARMA_DONT_USE_OPENMP 

using namespace Rcpp;
using namespace arma;

/*** ===================================================================== ***/
/*** HELPER FUNCTIONS                                                      ***/
/*** ===================================================================== ***/

// --- RNG & Permutation Helpers ---

static inline std::mt19937_64 make_rng_for_column(std::uint64_t base, arma::uword j) {
  // Hash the base seed with the column index to get a deterministic seed per column
  std::uint64_t x = base ^ (0x9e3779b97f4a7c15ULL + j + (j<<6) + (j>>2));
  return std::mt19937_64(x);
}

// Standard Fisher-Yates shuffle
template <class URNG>
static inline void shuffle_vec_in_place(arma::vec& v, URNG& rng) {
  for (arma::uword i = v.n_elem; i > 1; --i) {
    std::uniform_int_distribution<arma::uword> dist(0, i - 1);
    std::swap(v[i - 1], v[dist(rng)]);
  }
}

// Shuffle only specific indices (preserves NAs in place)
template <class URNG>
static inline void shuffle_selected_in_place(arma::vec& y, const arma::uvec& idx, URNG& rng) {
  for (arma::uword i = idx.n_elem; i > 1; --i) {
    std::uniform_int_distribution<arma::uword> dist(0, i - 1);
    arma::uword j = dist(rng);
    std::swap(y[idx[i - 1]], y[idx[j]]);
  }
}

// Permute within blocks/groups
template <class URNG>
static inline void permute_by_groups(arma::vec& y, const std::vector<std::vector<arma::uword>>& groups, URNG& rng) {
  for (const auto& g : groups) {
    if (g.size() < 2) continue;
    for (std::size_t i = g.size(); i > 1; --i) {
      std::uniform_int_distribution<std::size_t> dist(0, i - 1);
      std::swap(y[g[i - 1]], y[g[dist(rng)]]);
    }
  }
}

// --- Linear Algebra Helpers ---

static inline bool inv_xtx_safe(const arma::mat& X, arma::mat& invXtX, arma::mat& Xt, double pinv_tol) {
  Xt = X.t();
  arma::mat XtX = Xt * X;
  XtX = arma::symmatu(XtX); // Ensure exact symmetry
  if (!XtX.is_finite()) return false;
  
  // Try Cholesky/fast inverse first
  if (inv_sympd(invXtX, XtX)) return true;
  
  // Fallback to Moore-Penrose pseudo-inverse
  if (pinv_tol > 0.0) invXtX = arma::pinv(XtX, pinv_tol);
  else                invXtX = arma::pinv(XtX);
  
  return invXtX.is_finite();
}

static inline double robust_scale_mad(const arma::vec& r) {
  arma::uvec idx = arma::find_finite(r);
  if (idx.n_elem == 0) return 1e-8;
  arma::vec rf = r.elem(idx);
  double med = arma::median(rf);
  arma::vec af = arma::abs(rf - med);
  double mad = (af.n_elem > 0) ? arma::median(af) : 0.0;
  double s = 1.4826 * mad;
  return (s > 1e-12) ? s : 1e-8;
}

static inline bool is_ill_conditioned(const arma::mat& X, double rcond_thresh) {
  if (X.n_rows < X.n_cols) return true;
  arma::mat XtX = X.t() * X;
  if (!XtX.is_finite()) return true;
  double rc = arma::rcond(arma::symmatu(XtX));
  return (rc < rcond_thresh);
}

// FNV-1a Hash for NA masks
static inline std::uint64_t hash_na_mask(const arma::vec& y) {
  const std::uint64_t FNV_OFFSET = 1469598103934665603ULL;
  const std::uint64_t FNV_PRIME  = 1099511628211ULL;
  std::uint64_t h = FNV_OFFSET;
  for (const double& val : y) {
    unsigned char b = arma::is_finite(val) ? 1u : 0u;
    h ^= (std::uint64_t)b;
    h *= FNV_PRIME;
  }
  return h;
}

static inline double z_from_p(double p, int alt, double obs, double med_perm) {
  if (!std::isfinite(p)) return arma::datum::nan;
  if (p >= 1.0) return 0.0;
  if (p <= 0.0) p = 1e-16; // Guard
  
  if (alt == 0) { // two-sided
    double z = R::qnorm(1.0 - p/2.0, 0.0, 1.0, 1, 0);
    return (obs >= med_perm) ? z : -z;
  } else if (alt == 1) { // greater
    return R::qnorm(1.0 - p, 0.0, 1.0, 1, 0);
  } else { // less
    return -R::qnorm(1.0 - p, 0.0, 1.0, 1, 0);
  }
}

// --- FWL Helpers (QR Decomposition) ---

static inline bool qr_basis(const arma::mat& Z, arma::mat& Q, double rank_tol = 1e-12) {
  if (Z.n_cols == 0) { Q.set_size(Z.n_rows, 0); return true; }
  arma::mat R;
  bool ok = arma::qr_econ(Q, R, Z);
  if (!ok || !Q.is_finite() || !R.is_finite()) { Q.reset(); return false; }
  arma::uword r = arma::rank(R, rank_tol);
  if (r == 0) { Q.set_size(Z.n_rows, 0); return true; }
  if (r < Q.n_cols) Q = Q.cols(0, r - 1);
  return true;
}

static inline arma::mat project_out_Q(const arma::mat& Q, const arma::mat& B) {
  if (Q.n_cols == 0) return B;
  return B - Q * (Q.t() * B);
}

/*** ===================================================================== ***/
/*** STRUCTS & LOGIC                                                       ***/
/*** ===================================================================== ***/

struct Config {
  std::string robust, na_mode;
  double huber_k, huber_tol, na_weight, pinv_tol, illcond_rcond;
  int huber_maxit, n_randomizations, alt_code; // 0=two, 1=gr, 2=less
  bool center_mean, ret_res, ret_fits, ret_stats;
};

// Represents a group of columns that share the same NA pattern 
struct DesignGroup {
  bool valid_design = false;
  bool can_permute = false;
  
  // Data Indices
  arma::uvec obs_indices; 
  std::vector<std::vector<arma::uword>> perm_groups_sub;

  // Linear Algebra Cache (computed once per group)
  arma::mat X_sub;    // Design matrix (weighted if impute)
  arma::mat B;        // Projection (invXtX * Xt)
  arma::mat invXtX;   // Cached for Huber
  arma::mat Xt;       // Cached for Huber
  arma::vec weights;  // For impute_weak
};

// Represents a specific task (Column J belongs to DesignGroup G)
struct Job {
  arma::uword col_idx;
  int group_idx; 
};

// --- Fitting Routines ---

static arma::vec huber_irls(const arma::mat& X, const arma::vec& y, 
                            const DesignGroup& g, const Config& cfg) {
  arma::vec beta = g.invXtX * (g.Xt * y); // Init OLS
  if (!beta.is_finite()) return beta;
  
  for (int it = 0; it < cfg.huber_maxit; ++it) {
    arma::vec r = y - X * beta;
    double s = robust_scale_mad(r);
    if (s <= 1e-12) break;
    
    const double ks = cfg.huber_k * s;
    arma::vec w_hub(r.n_elem);
    for(arma::uword i=0; i<r.n_elem; ++i) {
      double ar = std::abs(r[i]);
      w_hub[i] = (ar > ks) ? (ks / ar) : 1.0;
    }
    
    // Combine with imputation weights if present
    if (cfg.na_mode == "impute_weak") w_hub %= g.weights;
    
    // Weighted Ridge Solve
    arma::mat XtWX = X.t() * (X.each_col() % w_hub);
    arma::vec XtWy = X.t() * (y % w_hub);
    
    // Add mild ridge for stability in IRLS
    double tr = arma::trace(XtWX);
    double lambda = 1e-10 * ((tr > 0.0) ? tr / X.n_cols : 1.0);
    XtWX.diag() += lambda;

    arma::mat Minv;
    if (!inv_sympd(Minv, XtWX)) Minv = arma::pinv(XtWX);
    
    arma::vec beta_new = Minv * XtWy;
    if (!beta_new.is_finite()) return beta; 
    
    double change = arma::norm(beta_new - beta) / (arma::norm(beta) + 1e-12);
    beta = beta_new;
    if (change < cfg.huber_tol) break;
  }
  return beta;
}

static inline arma::vec winsor_fit(const DesignGroup& g, const arma::vec& y, double k) {
  arma::vec beta = g.B * y;
  arma::vec r = y - g.X_sub * beta;
  double s = robust_scale_mad(r);
  if (s <= 1e-12) return beta;
  
  const double ks = k * s;
  r.clamp(-ks, ks); 
  // Refit on y_clean = X*beta + r_clipped
  return g.B * (g.X_sub * beta + r);
}

/*** ===================================================================== ***/
/*** MAIN EXPORT: FIT_AND_RANDOMIZE                                        ***/
/*** ===================================================================== ***/

/*
 fit_and_randomize: Fast OLS / Robust (Winsor / Huber IRLS) with permutation Z-scores
 
 WHAT IT DOES
 ------------
 Fits each column of Y on a design matrix X, computes a linear contrast of coefficients,
 generates a permutation null (full or blocked), and returns observed stats and Z-scores.
 
 OPTIMIZED ARCHITECTURE
 ----------------------
 This implementation uses a "Group & Flatten" strategy to maximize parallelism:
 1. Columns are grouped by their NA pattern (Mask).
 2. Linear algebra (X'X inverse) is computed once per Group in parallel.
 3. A "Job List" is created (mapping Column -> Group).
 4. Execution is "flattened" into a single parallel loop over columns, preventing
    thread starvation if groups vary significantly in size.
 
 ROBUSTNESS
 ----------
 robust = "none"   : OLS. Uses a fast "stat-only" path for permutations if sampled fits aren't requested.
 robust = "winsor" : Refits OLS on residuals clipped at ±k·MAD. Fast and resistant to outliers.
 robust = "huber"  : Iterative Reweighted Least Squares (IRLS). Fully robust but slower.
 
 NA HANDLING
 -----------
 na_mode = "drop"        : Drop NA rows per Y column (exact).
 na_mode = "impute_weak" : Impute missing Y to mean/zero with tiny weight (na_weight).
                           Allows keeping X fixed size, but only observed slots are permuted.
 
 PERMUTATIONS
 ------------
 perm_groups = NULL => Full permutations.
 Otherwise, a list of 1-based integer vectors. Permutations occur **within** these groups.
 
 INPUTS
 ------
 - X: (n x p) Design matrix.
 - Y: (n x m) Response matrix.
 - contrast: (p) Contrast vector.
 - perm_groups: List of groups for restricted permutation (optional).
 - n_randomizations: Number of permutations per column.
 - alternative: "two-sided" | "greater" | "less".
 - robust: "none", "winsor", "huber".
 - na_mode: "drop" or "impute_weak".
 - n_cores: Number of OpenMP threads.
 
 RETURNS
 -------
 List containing:
 - coef (p x m): Observed coefficients.
 - stat (m): Observed contrast statistic.
 - z_score (m): Z-score derived from permutation P-value.
 - p_value (m): Permutation P-value (add-one).
 - residuals (n x m): (Optional) Residual matrix.
 - sampled_fits: (Optional) List of permutation matrices.
 - sampled_stats: (Optional) Matrix of permutation statistics.
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
  
  // 1. Configuration & Validation
  Config cfg;
  cfg.robust = robust; cfg.na_mode = na_mode; cfg.huber_k = huber_k;
  cfg.huber_maxit = huber_maxit; cfg.huber_tol = huber_tol;
  cfg.na_weight = na_weight; cfg.center_mean = (na_center == "mean");
  cfg.pinv_tol = pinv_tol; cfg.illcond_rcond = illcond_rcond;
  cfg.n_randomizations = n_randomizations;
  cfg.ret_res = return_residuals; cfg.ret_fits = return_sampled_fits;
  cfg.ret_stats = return_sampled_stats;
  
  if (alternative == "two-sided") cfg.alt_code = 0;
  else if (alternative == "greater") cfg.alt_code = 1;
  else cfg.alt_code = 2;

  arma::uword n = X.n_rows, p = X.n_cols, m = Y.n_cols;
  if (Y.n_rows != n) stop("X and Y dimension mismatch");

  // Parse Permutation Groups
  std::vector<std::vector<arma::uword>> global_perm_groups;
  if (perm_groups.isNotNull()) {
    Rcpp::List pg(perm_groups);
    for (int i = 0; i < pg.size(); ++i) {
      IntegerVector g = pg[i];
      std::vector<arma::uword> idxs; 
      for(int x : g) if(x >= 1 && x <= (int)n) idxs.push_back(x - 1);
      global_perm_groups.push_back(std::move(idxs));
    }
  }

  // 2. Group Columns by NA Pattern (Serial)
  // We use a hash map to quickly identify columns with identical missingness
  std::unordered_map<std::uint64_t, std::vector<arma::uword>> map_mask;
  for (arma::uword j = 0; j < m; ++j) {
    map_mask[hash_na_mask(Y.col(j))].push_back(j);
  }

  // 3. Prepare "DesignGroups" and "Jobs"
  // Separating the Design (Matrix) from the Job (Column) enables flattened parallelism
  std::vector<DesignGroup> designs;
  designs.reserve(map_mask.size());
  
  std::vector<Job> jobs;
  jobs.reserve(m);
  
  int group_counter = 0;
  struct RawGroup { std::vector<arma::uword> cols; arma::uvec obs; };
  std::vector<RawGroup> raw_groups;
  
  for (auto& kv : map_mask) {
    if (kv.second.empty()) continue;
    arma::uword first_col = kv.second[0];
    arma::uvec obs = arma::find_finite(Y.col(first_col));
    
    raw_groups.push_back({kv.second, obs});
    
    for (arma::uword c : kv.second) {
      jobs.push_back({c, group_counter});
    }
    group_counter++;
  }
  designs.resize(group_counter);

  // 4. Compute Linear Algebra (Parallel over Groups)
  // This computes (X'X)^-1 just once per unique NA pattern
  #ifdef _OPENMP
  if (n_cores > 1) omp_set_num_threads(n_cores);
  #endif

  #pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < group_counter; ++i) {
    DesignGroup& g = designs[i];
    const RawGroup& raw = raw_groups[i];
    g.obs_indices = raw.obs;
    bool is_drop = (cfg.na_mode == "drop");
    
    // Setup Design Matrix X_sub
    if (is_drop) {
      if (raw.obs.n_elem >= p + 1) {
        g.X_sub = X.rows(raw.obs);
        // Map global perm groups to this subset
        if (global_perm_groups.empty()) {
          g.can_permute = (raw.obs.n_elem >= 2);
        } else {
          std::vector<int> glob_to_sub(n, -1);
          for(uword k=0; k<raw.obs.n_elem; ++k) glob_to_sub[raw.obs[k]] = k;
          
          bool any_pair = false;
          for(const auto& g_full : global_perm_groups) {
            std::vector<uword> g_sub;
            for(uword idx : g_full) if(glob_to_sub[idx] != -1) g_sub.push_back(glob_to_sub[idx]);
            if (g_sub.size() >= 2) any_pair = true;
            g.perm_groups_sub.push_back(std::move(g_sub));
          }
          g.can_permute = any_pair;
        }
      }
    } else { 
      // impute_weak: Use Weighted Least Squares
      g.weights.set_size(n); g.weights.fill(cfg.na_weight);
      g.weights.elem(raw.obs).fill(1.0);
      arma::vec sqrtw = arma::sqrt(g.weights);
      g.X_sub = X.each_col() % sqrtw;
      
      // Impute Perm Logic: permute only observed slots amongst themselves
      if (global_perm_groups.empty()) {
        g.can_permute = (raw.obs.n_elem >= 2);
      } else {
        std::vector<char> is_obs(n, 0);
        for(uword k : raw.obs) is_obs[k] = 1;
        bool any_pair = false;
        for(const auto& g_full : global_perm_groups) {
          std::vector<uword> g_filt;
          for(uword idx : g_full) if(is_obs[idx]) g_filt.push_back(idx);
          if (g_filt.size() >= 2) any_pair = true;
          g.perm_groups_sub.push_back(std::move(g_filt));
        }
        g.can_permute = any_pair;
      }
    }

    // Matrix Inversion / Factorization
    if ((is_drop && raw.obs.n_elem >= p + 1) || (!is_drop)) {
       if (!is_ill_conditioned(g.X_sub, cfg.illcond_rcond) &&
            inv_xtx_safe(g.X_sub, g.invXtX, g.Xt, cfg.pinv_tol)) {
         g.B = g.invXtX * g.Xt;
         g.valid_design = true;
       }
    }
  }

  // 5. Run Fits & Permutations (Flattened Parallelism over Columns)
  
  arma::mat Coef(p, m);   Coef.fill(arma::datum::nan);
  arma::vec Stat(m);      Stat.fill(arma::datum::nan);
  arma::vec Zscore(m);    Zscore.fill(arma::datum::nan);
  arma::vec Pval(m);      Pval.fill(arma::datum::nan);
  arma::mat Resid;        if (cfg.ret_res) { Resid.set_size(n, m); Resid.fill(arma::datum::nan); }
  
  // Thread-safe storage for list outputs
  std::vector<arma::mat> SampledFits(m);
  arma::mat SampledStats; 
  if (cfg.ret_stats && n_randomizations > 0) {
     SampledStats.set_size(n_randomizations, m); 
     SampledStats.fill(arma::datum::nan); 
  }

  std::uint64_t base_seed = 0xD1B54A32D192ED03ULL + (std::uint64_t)std::time(0);
  bool is_huber  = (cfg.robust == "huber");
  bool is_winsor = (cfg.robust == "winsor");

  #pragma omp parallel for schedule(dynamic)
  for (size_t k = 0; k < jobs.size(); ++k) {
    const Job& job = jobs[k];
    const DesignGroup& grp = designs[job.group_idx];
    arma::uword j = job.col_idx;
    
    if (!grp.valid_design) continue;

    std::mt19937_64 rng = make_rng_for_column(base_seed, j);

    // --- A. Prep Data ---
    arma::vec y_full = Y.col(j);
    arma::vec y_work;

    if (cfg.na_mode == "impute_weak") {
      y_work = y_full;
      double mu = 0.0;
      if (cfg.center_mean) {
         arma::vec vals = y_full.elem(grp.obs_indices);
         if (vals.n_elem > 0) mu = arma::mean(vals);
      }
      for(uword r=0; r<n; ++r) if(!std::isfinite(y_work[r])) y_work[r] = mu;
      
      // Weight Y for OLS/Winsor (Huber weights are passed separately)
      if (!is_huber) y_work %= arma::sqrt(grp.weights); 
    } else {
      y_work = y_full.elem(grp.obs_indices);
    }

    // --- B. Fit Observed ---
    arma::vec beta;
    if (is_huber) {
      if (cfg.na_mode == "impute_weak") beta = huber_irls(X, y_work, grp, cfg);
      else beta = huber_irls(grp.X_sub, y_work, grp, cfg);
    } else if (is_winsor) {
      beta = winsor_fit(grp, y_work, cfg.huber_k);
    } else {
      beta = grp.B * y_work;
    }

    if (!beta.is_finite()) continue;
    
    double stat_obs = arma::dot(beta, contrast);
    Coef.col(j) = beta;
    Stat[j] = stat_obs;

    // --- C. Residuals ---
    if (cfg.ret_res) {
      arma::vec r_out(n); r_out.fill(arma::datum::nan);
      if (cfg.na_mode == "impute_weak") {
         arma::vec y_raw = Y.col(j);
         double mu = 0.0;
         if (cfg.center_mean && grp.obs_indices.n_elem > 0) mu = arma::mean(y_raw.elem(grp.obs_indices));
         arma::uvec na_idx = find_nonfinite(y_raw);
         y_raw.elem(na_idx).fill(mu);
         r_out = y_raw - X * beta;
         r_out.elem(na_idx).fill(arma::datum::nan);
      } else {
         arma::vec r_sub = y_work - grp.X_sub * beta;
         r_out.elem(grp.obs_indices) = r_sub;
      }
      Resid.col(j) = r_out;
    }

    if (cfg.n_randomizations == 0 || !grp.can_permute) continue;

    // --- D. Permutations ---
    arma::vec stats_perm(cfg.n_randomizations);
    arma::mat fits_perm; 
    if (cfg.ret_fits) fits_perm.set_size(cfg.n_randomizations, p);

    arma::vec y_perm = y_work;
    // For OLS/Stat-only optimization
    arma::vec alpha; 
    if (!is_huber && !is_winsor && !cfg.ret_fits) {
       alpha = grp.B.t() * contrast;
    }

    int ge = 0, le = 0, ge_abs = 0;

    for (int r = 0; r < cfg.n_randomizations; ++r) {
      y_perm = y_work; // reset

      // Permute
      if (grp.perm_groups_sub.empty()) {
        if (cfg.na_mode == "impute_weak") shuffle_selected_in_place(y_perm, grp.obs_indices, rng);
        else shuffle_vec_in_place(y_perm, rng);
      } else {
        permute_by_groups(y_perm, grp.perm_groups_sub, rng);
      }

      // Fit Permuted Data
      double s_perm = 0.0;
      arma::vec b_perm;

      if (is_huber) {
         if (cfg.na_mode == "impute_weak") b_perm = huber_irls(X, y_perm, grp, cfg);
         else b_perm = huber_irls(grp.X_sub, y_perm, grp, cfg);
         s_perm = arma::dot(b_perm, contrast);
      } else if (is_winsor) {
         b_perm = winsor_fit(grp, y_perm, cfg.huber_k);
         s_perm = arma::dot(b_perm, contrast);
      } else {
         if (cfg.ret_fits) {
            b_perm = grp.B * y_perm;
            s_perm = arma::dot(b_perm, contrast);
         } else {
            s_perm = arma::dot(alpha, y_perm);
         }
      }

      stats_perm[r] = s_perm;
      
      if (std::abs(s_perm) >= std::abs(stat_obs)) ge_abs++;
      if (s_perm >= stat_obs) ge++;
      if (s_perm <= stat_obs) le++;
      
      if (cfg.ret_fits) fits_perm.row(r) = b_perm.t();
    }

    // --- E. Stats ---
    double pval = (cfg.alt_code==0)? (ge_abs+1.0) : (cfg.alt_code==1)? (ge+1.0) : (le+1.0);
    pval /= (cfg.n_randomizations + 1.0);

    arma::uvec valid_p = arma::find_finite(stats_perm);
    double med = (valid_p.n_elem > 0) ? arma::median(stats_perm.elem(valid_p)) : 0.0;

    Pval[j] = pval;
    Zscore[j] = z_from_p(pval, cfg.alt_code, stat_obs, med);
    
    if (cfg.ret_stats) SampledStats.col(j) = stats_perm;
    if (cfg.ret_fits)  SampledFits[j] = fits_perm;
  }

  // 6. Wrap Output
  Rcpp::List out = Rcpp::List::create(
    _["coef"] = Coef, _["stat"] = Stat, _["z_score"] = Zscore, _["p_value"] = Pval
  );
  if (cfg.ret_res) out["residuals"] = Resid;
  if (cfg.ret_stats && n_randomizations > 0) out["sampled_stats"] = SampledStats;
  if (cfg.ret_fits) {
    Rcpp::List L(m);
    for(uword j=0; j<m; ++j) L[j] = Rcpp::wrap(SampledFits[j]);
    out["sampled_fits"] = L;
  }
  return out;
}

/*** ===================================================================== ***/
/*** MAIN EXPORT: FL_FWL_CPP                        ***/
/*** ===================================================================== ***/

// FNV-1a hash of a 0/1 mask (for NA pattern grouping in FWL)
static inline std::uint64_t fnv1a64_mask_char(const std::vector<unsigned char>& mask) {
  const std::uint64_t FNV_OFFSET = 1469598103934665603ULL;
  const std::uint64_t FNV_PRIME  = 1099511628211ULL;
  std::uint64_t h = FNV_OFFSET;
  for (unsigned char b : mask) { h ^= (std::uint64_t)b; h *= FNV_PRIME; }
  return h;
}

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

/*
 fl_fwl_cpp: Frisch-Waugh-Lovell Partial Regression with Nuisance Covariates
 
 WHAT IT DOES
 ------------
 Performs a partial regression of Y on X, controlling for Z.
 Model: Y ~ X + Z
 1. Projects Z out of X and Y (creating X_resid, Y_resid).
 2. Subsets the data to `core_rows` (optional validation set).
 3. Calls `fit_and_randomize` on the residualized data to infer X effects.
 
 ALGORITHM
 ---------
 1. Groups Y columns by NA pattern.
 2. For each group:
    a. Perform QR decomposition on Z (subsetted to finite rows).
    b. Calculate Residuals: Xr = (I - Q_z Q_z')X, Yr = (I - Q_z Q_z')Y.
    c. Subset Xr and Yr to `core_rows`.
    d. Pass Xr, Yr to `fit_and_randomize`.
 
 INPUTS
 ------
 - X: (n x p) Design matrix of interest.
 - Z: (n x k) Nuisance covariate matrix.
 - Y: (n x m) Response matrix.
 - contrast: (p) Contrast vector for X.
 - core_rows: Indices [1-based] or Logical vector indicating rows to use for final fit.
 - na_mode: "drop" or "impute_weak".
 - n_cores: Number of OpenMP threads.
 
 RETURNS
 -------
 List containing:
 - coef, stat, z_score, p_value: Inference on X (controlled for Z).
 - partial_core (nc x m): Y residuals (Y - Y_hat_Z) on the core rows.
 - residuals (nc x m): Full Model residuals (Y - Y_hat_X - Y_hat_Z) on core rows.
 */
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
                      std::string na_mode = "drop",      
                      double na_weight = 1e-4,
                      std::string na_center = "mean",    
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
  
  // Parse Core Rows
  arma::uvec idx_core = parse_core_rows(core_rows, n);
  const arma::uword nc = idx_core.n_elem;
  if (nc < p + 1) stop("Not enough core rows (|core_rows| < p+1).");
  
  // Map global row index -> position within core output (or -1)
  std::vector<int> pos_in_core((size_t)n, -1);
  for (arma::uword t = 0; t < nc; ++t) pos_in_core[(size_t)idx_core[t]] = (int)t;
  
  const arma::uword m = Y.n_cols;
  
  // Group Y columns by NA mask
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
    std::uint64_t key = fnv1a64_mask_char(mask);
    auto it = groups.find(key);
    if (it == groups.end()) {
      GInfo g; g.cols.push_back(j); g.mask_all = std::move(mask);
      groups.emplace(key, std::move(g));
    } else {
      it->second.cols.push_back(j);
    }
  }
  
  // Allocate outputs
  arma::mat Coef(p, m);   Coef.fill(arma::datum::nan);
  arma::vec Stat(m);      Stat.fill(arma::datum::nan);
  arma::vec Zscore(m);    Zscore.fill(arma::datum::nan);
  arma::vec Pval(m);      Pval.fill(arma::datum::nan);
  arma::mat Resid;        if (return_residuals) { Resid.set_size(nc, m); Resid.fill(arma::datum::nan); }
  arma::mat PartialCore;  PartialCore.set_size(nc, m); PartialCore.fill(arma::datum::nan);
  
  std::vector<arma::mat> SampledFits; if (return_sampled_fits) SampledFits.resize(m);
  arma::mat SampledStats; if (return_sampled_stats && n_randomizations > 0) { 
    SampledStats.set_size(n_randomizations, m); SampledStats.fill(arma::datum::nan); 
  }
  
  // Iterate Groups
  for (auto & kv : groups) {
    const GInfo& g = kv.second;
    const std::vector<unsigned char>& mask = g.mask_all;
    const std::vector<arma::uword>& J = g.cols;
    
    arma::uvec Jv(J.size());
    for (size_t a = 0; a < J.size(); ++a) Jv[a] = J[a];
    
    // --- Case 1: Drop NA ---
    if (use_drop) {
      std::vector<arma::uword> pos_fin; pos_fin.reserve(n);
      for (arma::uword i = 0; i < n; ++i) if (mask[i]) pos_fin.push_back(i);
      arma::uword n_fin = (arma::uword)pos_fin.size();
      if (n_fin < p + 1) continue;
      
      arma::uvec sel_fin = arma::uvec(pos_fin);
      arma::mat Xf = X.rows(sel_fin);
      arma::mat Zf = (Z.n_cols > 0) ? Z.rows(sel_fin) : arma::mat(n_fin, 0);
      arma::mat Yf = Y.submat(sel_fin, Jv);
      
      // Project Z out
      arma::mat Xr_fin, Yr_fin;
      if (Zf.n_cols > 0) {
        arma::mat Qz;
        if (!qr_basis(Zf, Qz, pinv_tol)) continue;
        Xr_fin = project_out_Q(Qz, Xf);
        Yr_fin = project_out_Q(Qz, Yf);
      } else {
        Xr_fin = Xf; Yr_fin = Yf;
      }
      
      // Subset to Core Rows
      std::vector<arma::uword> pos_corefin;
      std::vector<arma::uword> pos_corefin_corepos;
      
      for (arma::uword r = 0; r < n_fin; ++r) {
        arma::uword i_glob = sel_fin[r];
        int pc = pos_in_core[(size_t)i_glob];
        if (pc >= 0) {
          pos_corefin.push_back(r); 
          pos_corefin_corepos.push_back((arma::uword)pc);
        }
      }
      if (pos_corefin.size() < (size_t)(p + 1)) continue;
      
      arma::uvec sel_corefin = arma::uvec(pos_corefin);
      arma::uvec sel_corepos = arma::uvec(pos_corefin_corepos);
      
      arma::mat Xr = Xr_fin.rows(sel_corefin);
      arma::mat Yr = Yr_fin.rows(sel_corefin);
      
      if (is_ill_conditioned(Xr, illcond_rcond)) continue;
      
      // Call fit_and_randomize on the residualized data
      Rcpp::List ans = fit_and_randomize(
        Xr, Yr, contrast,
        R_NilValue, n_randomizations, alternative,
        return_residuals, return_sampled_fits, return_sampled_stats,
        robust, huber_k, huber_maxit, huber_tol,
        "drop", na_weight, na_center,
        illcond_rcond, pinv_tol, n_cores
      );
      
      // Map results back
      arma::mat B = ans["coef"]; 
      arma::vec s = ans["stat"]; 
      arma::vec z = ans["z_score"]; 
      arma::vec pv = ans["p_value"];
      
      for (size_t a = 0; a < J.size(); ++a) {
        arma::uword j = J[a];
        Coef.col(j) = B.col(a);
        Stat[j]     = s[a];
        Zscore[j]   = z[a];
        Pval[j]     = pv[a];
      }
      
      if (return_residuals) {
        arma::mat Rg = ans["residuals"]; 
        Resid.submat(sel_corepos, Jv) = Rg; 
        PartialCore.submat(sel_corepos, Jv) = Yr;
      } else {
        PartialCore.submat(sel_corepos, Jv) = Yr;
      }
      
      if (return_sampled_fits) {
        Rcpp::List L = ans["sampled_fits"];
        for (size_t a = 0; a < J.size(); ++a) SampledFits[J[a]] = Rcpp::as<arma::mat>(L[a]);
      }
      if (return_sampled_stats && n_randomizations > 0) {
        arma::mat SS = ans["sampled_stats"];
        SampledStats.cols(Jv) = SS;
      }
      
    } else { 
      // --- Case 2: Impute Weak ---
      arma::vec w(n); w.fill(na_weight);
      for (arma::uword i = 0; i < n; ++i) if (mask[i]) w[i] = 1.0;
      arma::vec sqrtw = arma::sqrt(w);
      
      arma::mat Xw = X.each_col() % sqrtw;
      arma::mat Zw = (Z.n_cols > 0) ? Z.each_col() % sqrtw : arma::mat(n, 0);
      
      arma::mat Yw(n, J.size());
      for (size_t a = 0; a < J.size(); ++a) {
        arma::uword j = J[a];
        arma::vec y = Y.col(j);
        // Impute
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
      
      arma::mat Xr = Xr_full.rows(idx_core);
      arma::mat Yr = Yr_full.rows(idx_core);
      
      if (is_ill_conditioned(Xr, illcond_rcond)) continue;
      
      Rcpp::List ans = fit_and_randomize(
        Xr, Yr, contrast,
        R_NilValue, n_randomizations, alternative,
        return_residuals, return_sampled_fits, return_sampled_stats,
        robust, huber_k, huber_maxit, huber_tol,
        "drop", na_weight, na_center,
        illcond_rcond, pinv_tol, n_cores
      );
      
      arma::mat B = ans["coef"];
      arma::vec s = ans["stat"];
      arma::vec z = ans["z_score"];
      arma::vec pv = ans["p_value"];
      
      for (size_t a = 0; a < J.size(); ++a) {
        arma::uword j = J[a];
        Coef.col(j) = B.col(a);
        Stat[j]     = s[a];
        Zscore[j]   = z[a];
        Pval[j]     = pv[a];
      }
      if (return_residuals) {
        arma::mat Rg = ans["residuals"]; 
        Resid.cols(Jv) = Rg;
        PartialCore.cols(Jv) = Yr;
      } else {
        PartialCore.cols(Jv) = Yr;
      }
      
      if (return_sampled_fits) {
        Rcpp::List L2 = ans["sampled_fits"];
        for (size_t a = 0; a < J.size(); ++a) SampledFits[J[a]] = Rcpp::as<arma::mat>(L2[a]);
      }
      if (return_sampled_stats && n_randomizations > 0) {
        arma::mat SS = ans["sampled_stats"];
        SampledStats.cols(Jv) = SS;
      }
    }
  }
  
  Rcpp::List out = Rcpp::List::create(
    _["coef"]        = Coef,
    _["stat"]        = Stat,
    _["z_score"]     = Zscore,
    _["p_value"]     = Pval
  );
  if (return_residuals)     out["residuals"]     = Resid;
  if (return_sampled_fits) {
    Rcpp::List L(m); for (arma::uword j = 0; j < m; ++j) L[j] = Rcpp::wrap(SampledFits[j]);
    out["sampled_fits"] = L;
  }
  if (return_sampled_stats && n_randomizations > 0) out["sampled_stats"] = SampledStats;
  out["partial_core"] = PartialCore;
  
  return out;
}


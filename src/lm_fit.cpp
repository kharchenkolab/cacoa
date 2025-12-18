// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "lm_common.h"

/*
 * UNIFIED OPTIMIZED STATISTICS MODULE
 * ===================================
 * This module implements high-performance linear modeling and permutation-based inference
 * for high-dimensional data (e.g., gene expression, distance matrices).
 * 
 * The implementation is headed by two top-level functions: 
 *  fit_and_randomize implementing core fit and statistics
 *  fl_fwl_cpp implementing Freedman-Lane workflow
 *
 * ARCHITECTURE & OPTIMIZATIONS
 * ----------------------------
 * 1. NA Pattern Grouping:
 * Columns of Y are grouped by their missingness (NA) pattern. Expensive linear algebra 
 * factorizations (e.g., (X'X)^-1) are computed only once per group, drastically 
 * reducing overhead for sparse or incomplete data.
 *
 * 2. Flattened Parallelism:
 * Instead of nesting parallel loops (which causes thread starvation when groups vary in size),
 * the workload is flattened into a single list of (Column, Group) jobs. This ensures 
 * perfect load balancing across OpenMP threads.
 *
 * 3. Unified Randomization Engine:
 * Supports a layered randomization logic to handle complex designs:
 * - Standard Mode: Shuffles rows (observations) directly.
 * - Graph/MRQAP Mode: Shuffles nodes (samples) and maps them to edges (rows) using 
 * a provided topology map ('pair_indices').
 * - Block Constraints: Both modes support stratified permutation (shuffling within 
 * defined blocks, e.g., Batches) via 'perm_groups'.
 *
 * 4. Freedman-Lane Procedure (FWL):
 * Implements partial regression to control for nuisance covariates (Z) by projecting 
 * them out of the data before permutation testing, preserving the statistical 
 * validity of the test for the variable of interest (X).
 */


/*** ===================================================================== ***/
/*** MATH & STATS HELPERS                                                  ***/
/*** ===================================================================== ***/

static inline bool is_ill_conditioned(const arma::mat& X, double rcond_thresh) {
  if (X.n_rows < X.n_cols) return true;
  arma::mat XtX = X.t() * X;
  return (!XtX.is_finite() || arma::rcond(arma::symmatu(XtX)) < rcond_thresh);
}


/*** ===================================================================== ***/
/*** CORE LOGIC & EXPORTS                                                  ***/
/*** ===================================================================== ***/


// Represents a group of columns sharing the same NA pattern (Design reuse)
struct DesignGroup {
  bool valid_design = false;
  bool can_permute = false;
  arma::uvec obs_indices; 
  arma::mat X_sub, B, invXtX, Xt;       
  arma::vec weights;
  
  // Permutation config for this group
  std::vector<arma::uvec> perm_blocks;
  arma::uword n_units_for_perm; 
};

struct Job { arma::uword col_idx; int group_idx; };

// --- Iteratively Reweighted Least Squares (Huber) ---
static arma::vec huber_irls(const arma::mat& X, const arma::vec& y, 
                            const DesignGroup& g, const Config& cfg) {
  arma::vec beta = g.invXtX * (g.Xt * y); 
  if (!beta.is_finite()) return beta;
  
  for (int it = 0; it < cfg.huber_maxit; ++it) {
    arma::vec r = y - X * beta;
    double s = robust_scale_mad(r);
    if (s <= 1e-12) break;
    
    double ks = cfg.huber_k * s;
    arma::vec w = arma::abs(r);
    w.transform([&](double val){ return (val > ks) ? (ks/val) : 1.0; });
    if (cfg.na_mode == "impute_weak") w %= g.weights;
    
    arma::mat XtWX = X.t() * (X.each_col() % w);
    arma::vec XtWy = X.t() * (y % w);
    // Mild ridge stabilization
    double tr = arma::trace(XtWX);
    XtWX.diag() += 1e-10 * ((tr > 0.0) ? tr/X.n_cols : 1.0);
    
    arma::mat Minv;
    if (!inv_sympd(Minv, XtWX)) Minv = arma::pinv(XtWX);
    
    arma::vec beta_new = Minv * XtWy;
    if (!beta_new.is_finite()) return beta; 
    
    if (arma::norm(beta_new - beta)/(arma::norm(beta)+1e-12) < cfg.huber_tol) { beta=beta_new; break; }
    beta = beta_new;
  }
  return beta;
}

// --- Winsorized Fit ---
static inline arma::vec winsor_fit(const DesignGroup& g, const arma::vec& y, double k) {
  arma::vec beta = g.B * y;
  arma::vec r = y - g.X_sub * beta;
  double s = robust_scale_mad(r);
  if (s <= 1e-12) return beta;
  r.clamp(-k*s, k*s);
  return g.B * (g.X_sub * beta + r);
}

// -------------------------------------------------------------------------
// FIT_AND_RANDOMIZE
// -------------------------------------------------------------------------
/*
 * fit_and_randomize
 * =================
 * Fits a linear model (Y ~ X) for each column of Y, estimates a linear contrast
 * of coefficients, and computes permutation-based statistics (P-values and Z-scores)
 * using a highly optimized, parallelized engine.
 *
 * ALGORITHMIC FEATURES
 * --------------------
 * 1. NA Pattern Grouping:
 * Columns of Y are automatically grouped by their missingness pattern. The expensive
 * linear algebra factorization (X'X)^-1 is computed only once per group, providing
 * massive speedups for datasets with sporadic missing values.
 *
 * 2. Flattened Parallelism:
 * Work is distributed across 'n_cores' using a dynamic load-balancing strategy
 * that flattens the nested loop (Groups -> Columns) into a single job list.
 * This prevents thread starvation when groups vary significantly in size.
 *
 * 3. Unified Randomization Engine:
 * Supports two fundamental modes of permutation, switched automatically by 'pair_indices':
 * - Standard Mode (Row Shuffling): Permutes the rows of the design matrix directly.
 * Used for standard independent sampling designs.
 * - Graph Mode (Node Shuffling): Permutes the underlying biological units (Samples/Nodes)
 * and maps them to the observation rows (Edges/Pairs) using the topology map.
 * Used for pairwise designs (e.g., distance matrices) to preserve geometric dependencies (MRQAP).
 *
 * 4. Stratified Permutation (Blocking):
 * Both modes support restricted randomization via 'perm_groups'.
 *
 * PARAMETERS
 * ----------
 * @param X (arma::mat)
 * The design matrix (n_rows x n_preds). Must include an intercept column if desired.
 *
 * @param Y (arma::mat)
 * The response matrix (n_rows x n_features). Each column is fit independently.
 *
 * @param contrast (arma::vec)
 * A linear contrast vector (length = n_preds) defining the statistic of interest.
 * The observed statistic is calculated as: stat = dot(contrast, beta).
 *
 * @param perm_groups (Rcpp::Nullable<Rcpp::List>)
 * A list of integer vectors defining the blocks of exchangeable units for stratified permutation.
 * - In Standard Mode (pair_indices = NULL): These are ROW indices [1-based].
 * Rows are only swapped with other rows in the same block.
 * - In Graph Mode (pair_indices != NULL): These are SAMPLE indices [1-based].
 * Samples are only swapped with other samples in the same block (e.g., Batch).
 * If NULL, an unrestricted global permutation is performed.
 *
 * @param pair_indices (Rcpp::Nullable<arma::umat>)
 * An (N x 2) matrix defining the graph topology (Node -> Edge map).
 * - If provided, the function switches to GRAPH/MRQAP Randomization.
 * - Row 'k' of X/Y corresponds to the pair of samples (pair_indices[k, 0], pair_indices[k, 1]).
 * - Randomization involves shuffling the sample IDs and reconstructing the row order.
 *
 * @param n_randomizations (int)
 * The number of permutations to perform per column.
 *
 * @param alternative (std::string)
 * The alternative hypothesis for P-value calculation:
 * - "two-sided": P = (|perm| >= |obs|). Z-score preserves the sign of (obs - median).
 * - "greater":   P = (perm >= obs).
 * - "less":      P = (perm <= obs).
 *
 * @param robust (std::string)
 * The estimation method:
 * - "none":   Ordinary Least Squares (OLS). Fastest.
 * - "winsor": OLS on residuals winsorized at 'huber_k' * MAD. Fast approximation of robust regression.
 * - "huber":  Iteratively Reweighted Least Squares (IRLS) using Huber weights. Fully robust but slower.
 *
 * @param huber_k (double)
 * The tuning constant for Huber/Winsorization (defaults to 1.345).
 *
 * @param huber_maxit (int)
 * Maximum iterations for Huber IRLS.
 *
 * @param huber_tol (double)
 * Convergence tolerance for Huber IRLS.
 *
 * @param na_mode (std::string)
 * How to handle missing values (NAs) in Y:
 * - "drop": Rows with NAs are excluded from the fit (exact OLS on subset).
 * *NOTE*: In Graph Mode (if pair_indices != NULL), 'drop' fits the observed statistic on the valid subset,
 * but utilizes a mean-imputed "Clean Y" vector for generating the null distribution
 * to prevent permutations from pulling NAs into valid slots.
 * - "impute_weak": NAs are replaced by the mean (or 0) and assigned a negligible weight ('na_weight').
 * Maintains constant vector size, which is numerically stable for Graph/MRQAP.
 *
 * @param na_weight (double)
 * The weight assigned to imputed observations in "impute_weak" mode (typically 1e-4).
 *
 * @param na_center (std::string)
 * "mean" (impute with observed mean) or "zero".
 *
 * @param illcond_rcond (double)
 * Reciprocal condition number threshold to detect singular designs. Groups failing this check are skipped.
 *
 * @param pinv_tol (double)
 * Tolerance for the pseudo-inverse (if Cholesky decomposition fails).
 *
 * @param n_cores (int)
 * Number of OpenMP threads to use.
 *
 * @param return_residuals (bool)
 * If TRUE, returns the (n x m) matrix of residuals from the observed fit.
 *
 * @param return_sampled_fits (bool)
 * If TRUE, returns a list of matrices containing the estimated coefficients for every permutation.
 *
 * @param return_sampled_stats (bool)
 * If TRUE, returns the (n_perm x m) matrix of permutation statistics.
 *
 * OUTPUT (Rcpp::List)
 * -------------------
 * coef : (p x m) matrix of observed coefficients.
 * stat : (m) vector of observed contrast statistics.
 * z_score : (m) vector of Z-scores derived from the permutation P-values.
 * (Quantile function of 1 - P/2 for two-sided).
 * p_value : (m) vector of permutation P-values (calculated with add-one smoothing).
 * residuals : (Optional) (n x m) matrix of residuals. NaN where input was NA (in drop mode).
 * sampled_stats : (Optional) (n_perm x m) matrix of null statistics.
 * sampled_fits : (Optional) List of (n_perm x p) matrices containing null coefficients.
 */
// [[Rcpp::export]]
Rcpp::List fit_and_randomize(const arma::mat& X, const arma::mat& Y, const arma::vec& contrast,
                             Rcpp::Nullable<Rcpp::List> perm_groups = R_NilValue,
                             Rcpp::Nullable<arma::umat> pair_indices = R_NilValue,
                             int n_randomizations = 100,
                             std::string alternative = "two-sided",
                             bool return_residuals = true, bool return_sampled_fits = false, bool return_sampled_stats = false,
                             std::string robust = "none", double huber_k = 1.345, int huber_maxit = 8, double huber_tol = 1e-6,
                             std::string na_mode = "drop", double na_weight = 1e-4, std::string na_center = "mean",
                             double illcond_rcond = 1e-12, double pinv_tol = 0.0, int n_cores = 1) {
  
  Config cfg; 
  cfg.robust = robust; cfg.na_mode = na_mode; cfg.huber_k = huber_k; cfg.huber_maxit = huber_maxit; cfg.huber_tol = huber_tol;
  cfg.na_weight = na_weight; cfg.center_mean = (na_center == "mean"); cfg.pinv_tol = pinv_tol; cfg.illcond_rcond = illcond_rcond;
  cfg.n_randomizations = n_randomizations; cfg.ret_res = return_residuals; cfg.ret_fits = return_sampled_fits; cfg.ret_stats = return_sampled_stats;
  cfg.alt_code = (alternative=="two-sided")?0 : (alternative=="greater")?1 : 2;

  arma::uword n = X.n_rows, p = X.n_cols, m = Y.n_cols;
  if (Y.n_rows != n) stop("X and Y dimension mismatch");

  // 1. Setup Randomization Logic
  bool is_graph = pair_indices.isNotNull();
  PairLookup pair_mapper; arma::umat pairs_mat; arma::uword n_units = n;
  
  if (is_graph) {
      pairs_mat = Rcpp::as<arma::umat>(pair_indices);
      if (pairs_mat.min() > 0) pairs_mat -= 1; // 0-based correction
      n_units = pairs_mat.max() + 1;           // n_units = n_samples
      pair_mapper.init(pairs_mat, n_units);
  }

  // Parse Constraints (Blocks)
  std::vector<arma::uvec> blocks;
  if (perm_groups.isNotNull()) {
    Rcpp::List pg(perm_groups);
    for (int i=0; i<pg.size(); ++i) {
      IntegerVector g = pg[i]; arma::uvec ug = as<arma::uvec>(g);
      if(ug.max() > 0) ug -= 1; blocks.push_back(std::move(ug));
    }
  } else {
    // Default: Global shuffling
    blocks.push_back(arma::regspace<arma::uvec>(0, n_units - 1));
  }

  // 2. Group Columns by NA Pattern (Optimization)
  std::unordered_map<std::uint64_t, std::vector<arma::uword>> map_mask;
  for (arma::uword j=0; j<m; ++j) map_mask[hash_vec_mask(Y.col(j))].push_back(j);
  
  std::vector<DesignGroup> designs; designs.reserve(map_mask.size());
  std::vector<Job> jobs; jobs.reserve(m);
  int grp_cnt = 0;
  
  struct RawGroup { std::vector<arma::uword> cols; arma::uvec obs; };
  std::vector<RawGroup> raw_groups;

  for (auto& kv : map_mask) {
    arma::uword first = kv.second[0];
    raw_groups.push_back({kv.second, arma::find_finite(Y.col(first))});
    for(auto c : kv.second) jobs.push_back({c, grp_cnt});
    grp_cnt++;
  }
  designs.resize(grp_cnt);

  // 3. Pre-calculate Designs (Parallel)
  #ifdef _OPENMP
  if (n_cores > 1) omp_set_num_threads(n_cores);
  #endif

  #pragma omp parallel for schedule(dynamic)
  for (int i=0; i<grp_cnt; ++i) {
    DesignGroup& g = designs[i]; const RawGroup& raw = raw_groups[i];
    g.obs_indices = raw.obs;
    bool is_drop = (cfg.na_mode == "drop");
    
    // Matrix factorization
    if (is_drop) {
      if (raw.obs.n_elem >= p + 1) {
        g.X_sub = X.rows(raw.obs);
        g.can_permute = (raw.obs.n_elem >= 2);
      }
    } else { 
      g.weights.set_size(n); g.weights.fill(cfg.na_weight); g.weights.elem(raw.obs).fill(1.0);
      g.X_sub = X.each_col() % arma::sqrt(g.weights);
      g.can_permute = (raw.obs.n_elem >= 2);
    }

    if (g.X_sub.n_rows > 0 && !is_ill_conditioned(g.X_sub, cfg.illcond_rcond) && 
        inv_xtx_safe(g.X_sub, g.invXtX, g.Xt, cfg.pinv_tol)) {
      g.B = g.invXtX * g.Xt; g.valid_design = true;
    }

    // Configure permutation blocks for this NA pattern
    if (is_graph) {
        // Graph Mode: Shuffle Samples. NAs handled by subsetting resulting edges.
        g.perm_blocks = blocks; g.n_units_for_perm = n_units; 
    } else {
        if (is_drop) {
            // Standard Mode + Drop: Map global rows to subset indices to avoid NaN poisoning
            g.perm_blocks = subset_blocks(blocks, g.obs_indices, n);
            // Fallback if subsetting leaves no valid blocks
            if (g.perm_blocks.empty() && g.obs_indices.n_elem > 0) 
                g.perm_blocks.push_back(arma::regspace<arma::uvec>(0, g.obs_indices.n_elem - 1));
            g.n_units_for_perm = g.obs_indices.n_elem;
        } else {
            // Impute Mode: Shuffle full N rows
            g.perm_blocks = blocks; g.n_units_for_perm = n;
        }
    }
  }

  // 4. Execution (Flattened Parallelism)
  arma::mat Coef(p, m); Coef.fill(datum::nan);
  arma::vec Stat(m); Stat.fill(datum::nan);
  arma::vec Pval(m); Pval.fill(datum::nan);
  arma::mat Resid; if (cfg.ret_res) { Resid.set_size(n, m); Resid.fill(datum::nan); }
  
  std::vector<arma::mat> SampledFits(m);
  arma::mat SampledStats; 
  if (cfg.ret_stats && n_randomizations > 0) { SampledStats.set_size(n_randomizations, m); SampledStats.fill(datum::nan); }

  std::uint64_t seed = 0xD1B54A32D192ED03ULL + (std::uint64_t)std::time(0);
  bool is_huber = (cfg.robust == "huber"), is_winsor = (cfg.robust == "winsor");

  #pragma omp parallel for schedule(dynamic)
  for (size_t k=0; k<jobs.size(); ++k) {
    const Job& job = jobs[k]; const DesignGroup& grp = designs[job.group_idx];
    arma::uword j = job.col_idx;
    if (!grp.valid_design) continue;

    std::mt19937_64 rng = make_rng(seed, j);

    // Prepare Vectors
    arma::vec y_raw = Y.col(j);
    arma::vec y_work; 
    
    // y_clean: Dense vector for permutation (NAs filled) to prevent NaN pull-in
    arma::vec y_clean = y_raw; 
    double mu = 0.0;
    if (grp.obs_indices.n_elem > 0) mu = arma::mean(y_raw.elem(grp.obs_indices));
    y_clean.elem(find_nonfinite(y_clean)).fill(mu); 

    // y_work: Vector for OBSERVED fit
    if (cfg.na_mode == "drop") y_work = y_raw.elem(grp.obs_indices);
    else {
        y_work = y_clean;
        if (!is_huber) y_work %= arma::sqrt(grp.weights); 
    }

    // Observed Fit
    arma::vec beta;
    if (is_huber) beta = (cfg.na_mode=="drop") ? huber_irls(grp.X_sub, y_work, grp, cfg) : huber_irls(X, y_work, grp, cfg);
    else if (is_winsor) beta = winsor_fit(grp, y_work, cfg.huber_k);
    else beta = grp.B * y_work;

    if (!beta.is_finite()) continue;
    double stat_obs = arma::dot(beta, contrast);
    Coef.col(j) = beta; Stat[j] = stat_obs;

    if (cfg.ret_res) {
      arma::vec r_out(n); r_out.fill(datum::nan);
      if (cfg.na_mode == "drop") r_out.elem(grp.obs_indices) = y_work - grp.X_sub * beta;
      else {
          r_out = y_clean - X * beta;
          r_out.elem(find_nonfinite(y_raw)).fill(datum::nan); 
      }
      Resid.col(j) = r_out;
    }

    if (cfg.n_randomizations == 0 || !grp.can_permute) continue;

    // Permutation Loop
    arma::vec stats_perm(cfg.n_randomizations);
    arma::mat fits_perm; if (cfg.ret_fits) fits_perm.set_size(cfg.n_randomizations, p);
    
    // Fast path alpha for OLS
    arma::vec alpha; if (!is_huber && !is_winsor && !cfg.ret_fits) alpha = grp.B.t() * contrast;

    int ge=0, le=0, ge_abs=0;
    for (int r=0; r<cfg.n_randomizations; ++r) {
        // Generate indices
        arma::uvec perm_idx = generate_permutation(rng, grp.perm_blocks, grp.n_units_for_perm, is_graph, pair_mapper, pairs_mat);
        
        // Apply indices to data
        arma::vec y_perm;
        if (cfg.na_mode == "drop") {
            if (is_graph) {
                // Graph mode: Global Shuffle -> Subset to Observed. Use 'y_clean' to be safe.
                // Fix: Breaking chained subsetting for older compilers
                arma::vec y_global_perm = y_clean.elem(perm_idx);
                y_perm = y_global_perm.elem(grp.obs_indices); 
            } else {
                // Standard mode: Shuffle directly on the observed subset
                y_perm = y_work.elem(perm_idx);
            }
        } else {
            // Impute mode: Shuffle global
            y_perm = y_work.elem(perm_idx);
        }

        // Fit
        double s_perm = 0.0; arma::vec b_perm;
        if (is_huber) {
             b_perm = (cfg.na_mode=="drop") ? huber_irls(grp.X_sub, y_perm, grp, cfg) : huber_irls(X, y_perm, grp, cfg);
             s_perm = arma::dot(b_perm, contrast);
        } else if (is_winsor) {
             b_perm = winsor_fit(grp, y_perm, cfg.huber_k);
             s_perm = arma::dot(b_perm, contrast);
        } else {
             if (cfg.ret_fits) { b_perm = grp.B * y_perm; s_perm = arma::dot(b_perm, contrast); }
             else s_perm = arma::dot(alpha, y_perm);
        }

        // Stats
        stats_perm[r] = s_perm;
        if (std::abs(s_perm) >= std::abs(stat_obs)) ge_abs++;
        if (s_perm >= stat_obs) ge++;
        if (s_perm <= stat_obs) le++;
        if (cfg.ret_fits) fits_perm.row(r) = b_perm.t();
    }

    double pval = (cfg.alt_code==0)?(ge_abs+1.0):(cfg.alt_code==1)?(ge+1.0):(le+1.0);
    pval /= (cfg.n_randomizations + 1.0);
    
    arma::uvec valid_p = find_finite(stats_perm);
    double med = (valid_p.n_elem > 0) ? arma::median(stats_perm.elem(valid_p)) : 0.0;
    
    Pval[j] = pval;
    if (cfg.ret_stats) SampledStats.col(j) = stats_perm;
    if (cfg.ret_fits) SampledFits[j] = fits_perm;
  }

  Rcpp::List out = Rcpp::List::create(_["coef"]=Coef, _["stat"]=Stat, _["p_value"]=Pval);
  if (cfg.ret_res) out["residuals"] = Resid;
  if (cfg.ret_stats && n_randomizations > 0) out["sampled_stats"] = SampledStats;
  if (cfg.ret_fits) {
    Rcpp::List L(m); for(uword j=0; j<m; ++j) L[j]=Rcpp::wrap(SampledFits[j]);
    out["sampled_fits"] = L;
  }
  return out;
}

// -------------------------------------------------------------------------
// FL_FWL_CPP
// -------------------------------------------------------------------------
/*
 * Performs the Freedman-Lane procedure for partial regression (Y ~ X | Z).
 * This function controls for nuisance covariates (Z) by projecting them out of 
 * both the design matrix (X) and the response (Y) before performing inference 
 * on the variable of interest.
 *
 * ALGORITHM
 * ---------
 * 1. NA Pattern Grouping:
 * Columns of Y are grouped by missingness pattern to minimize redundant 
 * QR decompositions of Z.
 *
 * 2. Nuisance Projection (Residualization):
 * For each group, the nuisance matrix Z (subsetted to valid rows) is decomposed 
 * (QR). Both X and Y are projected onto the orthogonal complement of Z:
 * X_resid = (I - P_z) * X
 * Y_resid = (I - P_z) * Y
 *
 * 3. Core Subset Selection:
 * If 'core_rows' is provided, the residualized matrices are subsetted to these 
 * specific rows (e.g., for split-sample validation schemes).
 *
 * 4. Inference:
 * The residualized (and potentially subsetted) matrices are passed to 
 * 'fit_and_randomize' for parameter estimation and permutation testing.
 *
 * PARAMETERS
 * ----------
 * @param X (arma::mat)
 * The primary design matrix (n_rows x n_preds) containing variables of interest.
 *
 * @param Z (arma::mat)
 * The nuisance design matrix (n_rows x n_nuis) to be projected out.
 *
 * @param Y (arma::mat)
 * The response matrix (n_rows x n_features).
 *
 * @param contrast (arma::vec)
 * Linear contrast vector for X.
 *
 * @param core_rows (SEXP: IntegerVector, LogicalVector, or NULL)
 * Optional subset of rows to use for the final test statistic. 
 * - If NULL: All rows are used.
 * - If provided: The residualization (Step 2) uses ALL valid data to learn Z effects,
 * but the final fit (Step 4) uses only the specified 'core_rows'.
 *
 * @param perm_groups (Rcpp::Nullable<Rcpp::List>)
 * Passed to fit_and_randomize. Defines randomization blocks (Rows or Samples).
 *
 * @param pair_indices (Rcpp::Nullable<arma::umat>)
 * Passed to fit_and_randomize. Defines graph topology for MRQAP.
 *
 * @param na_mode (std::string)
 * Controls how NAs are handled during the Z projection step:
 * - "drop": Z is projected out using only rows where Y is observed.
 * - "impute_weak": Z is projected out using all rows (requires Z to be complete).
 *
 * @param [Pass-through Parameters]
 * n_randomizations, alternative, robust, huber_k, huber_maxit, huber_tol,
 * na_weight, na_center, illcond_rcond, pinv_tol, n_cores, 
 * return_residuals, return_sampled_fits, return_sampled_stats.
 * (See 'fit_and_randomize' documentation for details).
 *
 * OUTPUT (Rcpp::List)
 * -------------------
 * coef : (p x m) Estimated coefficients for X (after controlling for Z).
 * stat : (m) Contrast statistics.
 * z_score : (m) Permutation Z-scores.
 * p_value : (m) Permutation P-values.
 * partial_core : (n_core x m) The Y residuals (Y - Y_hat_Z) restricted to core_rows.
 * residuals : (Optional) (n_core x m) The full model residuals (Y - Y_hat_X - Y_hat_Z).
 * sampled_stats : (Optional) Null distribution statistics.
 * sampled_fits : (Optional) Null distribution coefficients.
 */
// [[Rcpp::export]]
Rcpp::List fl_fwl_cpp(const arma::mat& X, const arma::mat& Z, const arma::mat& Y, const arma::vec& contrast,
                      SEXP core_rows = R_NilValue, 
                      Rcpp::Nullable<Rcpp::List> core_perm_groups = R_NilValue,
                      Rcpp::Nullable<arma::umat> core_pair_indices = R_NilValue,
                      int n_randomizations = 100,
                      std::string alternative = "two-sided", std::string robust = "none",
                      double huber_k = 1.345, int huber_maxit = 8, double huber_tol = 1e-6,
                      std::string na_mode = "drop", double na_weight = 1e-4, std::string na_center = "mean",
                      double illcond_rcond = 1e-12, double pinv_tol = 0.0, int n_cores = 1,
                      bool return_residuals = true, bool return_sampled_fits = false, bool return_sampled_stats = false) {
  
  arma::uword n = X.n_rows, m = Y.n_cols;
  arma::uvec idx_core = parse_core_rows(core_rows, n);
  
  // 1. FAST PATH: No Nuisance Variables (Z is empty)
  if (Z.n_cols == 0) {
    arma::mat X_sub = X.rows(idx_core);
    arma::mat Y_sub = Y.rows(idx_core);
    Rcpp::List out = fit_and_randomize(X_sub, Y_sub, contrast, core_perm_groups, core_pair_indices,
                                       n_randomizations, alternative, return_residuals, return_sampled_fits, return_sampled_stats,
                                       robust, huber_k, huber_maxit, huber_tol, na_mode, na_weight, na_center,
                                       illcond_rcond, pinv_tol, n_cores);
    out["partial_core"] = Y_sub;
    return out;
  }
  
  // 2. FREEDMAN-LANE PATH (Has Z)
  std::unordered_map<std::uint64_t, std::vector<arma::uword>> groups;
  for (arma::uword j=0; j<m; ++j) {
    groups[hash_vec_mask(Y.col(j))].push_back(j);
  }
  
  arma::mat Coef(X.n_cols, m); Coef.fill(datum::nan);
  arma::vec Stat(m), Pval(m); Stat.fill(datum::nan); Pval.fill(datum::nan);
  arma::mat Resid; if(return_residuals) { Resid.set_size(idx_core.n_elem, m); Resid.fill(datum::nan); }
  arma::mat PartialCore; PartialCore.set_size(idx_core.n_elem, m); PartialCore.fill(datum::nan);
  arma::mat SampledStats; if(return_sampled_stats) { SampledStats.set_size(n_randomizations, m); SampledStats.fill(datum::nan); }
  std::vector<arma::mat> SampledFits(m); 
  
  std::vector<int> global_to_core(n, -1);
  for(uword k=0; k < idx_core.n_elem; ++k) global_to_core[idx_core[k]] = k;
  
  for (auto & kv : groups) {
    std::vector<arma::uword> J = kv.second;
    arma::uvec Jv = arma::conv_to<arma::uvec>::from(J);
    bool use_drop = (na_mode == "drop");
    
    // A. Identify Rows
    // 'valid_rows' are strictly those with FINITE data.
    // In 'impute' mode, we will use ALL rows, but we need 'valid_rows' to calculate the mean.
    arma::vec y_rep = Y.col(J[0]);
    arma::uvec valid_rows = arma::find_finite(y_rep);
    
    // We need 'obs_rows' for the projection matrix.
    // Drop: obs_rows = valid_rows
    // Impute: obs_rows = 0..n-1 (all)
    arma::uvec obs_rows;
    if (use_drop) obs_rows = valid_rows;
    else obs_rows = arma::regspace<arma::uvec>(0, n-1);
    
    if (obs_rows.n_elem < X.n_cols) continue;
    
    // B. Calculate Residuals (With Temp Imputation for Weak Mode)
    arma::mat Z_sub = Z.rows(obs_rows);
    arma::mat Qz;
    if (!qr_basis(Z_sub, Qz, pinv_tol)) continue; 
    
    // Handle Y for projection
    // If impute_weak, we must fill NaNs in Y before projecting, 
    // otherwise Q'*Y propagates NaNs everywhere.
    arma::mat Y_for_proj;
    if (use_drop) {
      Y_for_proj = Y.submat(obs_rows, Jv);
    } else {
      // Extract full columns
      Y_for_proj = Y.cols(Jv); 
      // Impute NaNs with Mean of Valid Data (or 0)
      // Since all cols in J share the same mask, we iterate cols
      for(uword c=0; c < Y_for_proj.n_cols; ++c) {
        arma::vec col = Y_for_proj.col(c);
        double mu = 0.0;
        if (valid_rows.n_elem > 0) {
          // valid_rows indices are global. We need valid indices relative to Y_for_proj (which is full size n)
          // Since Y_for_proj is n rows, valid_rows works directly.
          mu = (na_center == "mean") ? arma::mean(col.elem(valid_rows)) : 0.0;
        }
        // Fill non-finite values
        col.elem(find_nonfinite(col)).fill(mu);
        Y_for_proj.col(c) = col;
      }
    }
    
    arma::mat X_resid_full = project_out_Q(Qz, X.rows(obs_rows));
    arma::mat Y_resid_full = project_out_Q(Qz, Y_for_proj);
    
    // C. Construct Core Subset Matrices
    // If 'impute_weak', we must RE-MASK the values that were originally missing to NaN.
    // fit_and_randomize needs NaNs to trigger the weighting logic.
    
    arma::mat X_fin(idx_core.n_elem, X.n_cols); X_fin.fill(datum::nan);
    arma::mat Y_fin(idx_core.n_elem, J.size()); Y_fin.fill(datum::nan);
    
    bool has_data = false;
    
    // Map from obs_rows index -> global index -> core index
    for(uword k=0; k < obs_rows.n_elem; ++k) {
      uword glob_idx = obs_rows[k];
      int core_idx = global_to_core[glob_idx];
      
      if (core_idx != -1) {
        // Check if this row was originally valid
        bool was_valid = arma::is_finite(y_rep[glob_idx]);
        
        if (use_drop) {
          // In Drop mode, obs_rows ARE valid_rows. So always valid.
          X_fin.row(core_idx) = X_resid_full.row(k);
          Y_fin.row(core_idx) = Y_resid_full.row(k);
          has_data = true;
        } else {
          // In Impute mode, obs_rows includes NAs.
          X_fin.row(core_idx) = X_resid_full.row(k);
          
          if (was_valid) {
            Y_fin.row(core_idx) = Y_resid_full.row(k);
          } else {
            // It was originally NA. We imputed it for projection.
            // Now RE-MASK it to NaN so fit_and_randomize applies epsilon weight.
            Y_fin.row(core_idx).fill(datum::nan);
          }
          has_data = true; // Even if NaN, we pass it (as missing/weighted)
        }
      }
    }
    
    if (!has_data) continue;
    
    // D. Call Fitter
    Rcpp::List res = fit_and_randomize(X_fin, Y_fin, contrast, core_perm_groups, core_pair_indices,
                                       n_randomizations, alternative, return_residuals, return_sampled_fits, return_sampled_stats,
                                       robust, huber_k, huber_maxit, huber_tol, na_mode, na_weight, na_center,
                                       illcond_rcond, pinv_tol, n_cores);
    
    // E. Unpack
    arma::mat B = res["coef"];
    arma::vec S = res["stat"];
    arma::vec Pv = res["p_value"];
    
    for(uword i=0; i<J.size(); ++i) {
      uword col = J[i];
      Coef.col(col) = B.col(i);
      Stat(col) = S(i); Pval(col) = Pv(i);
    }
    
    // Fill PartialCore and Resid
    if (return_residuals || true) {
      arma::mat Rr = res["residuals"]; 
      for(uword c=0; c < J.size(); ++c) {
        uword col = J[c];
        Resid.col(col) = Rr.col(c);
        
        // For PartialCore: In 'impute' mode, we might want the IMPUTED residual 
        // rather than NaN for plotting? 
        // Standard behavior: Return what was used. If it was NaN (masked), return NaN.
        PartialCore.col(col) = Y_fin.col(c);
      }
    }
    
    if(return_sampled_stats) {
      arma::mat SS = res["sampled_stats"];
      SampledStats.cols(Jv) = SS;
    }
    if(return_sampled_fits) {
      Rcpp::List L = res["sampled_fits"];
      for(uword i=0; i<J.size(); ++i) SampledFits[J[i]] = Rcpp::as<arma::mat>(L[i]);
    }
  }
  
  Rcpp::List out = Rcpp::List::create(_["coef"]=Coef, _["stat"]=Stat, _["p_value"]=Pval);
  if(return_residuals) out["residuals"] = Resid;
  out["partial_core"] = PartialCore;
  if(return_sampled_stats && n_randomizations>0) out["sampled_stats"] = SampledStats;
  if(return_sampled_fits) {
    Rcpp::List L(m); for(uword j=0; j<m; ++j) L[j]=Rcpp::wrap(SampledFits[j]);
    out["sampled_fits"] = L;
  }
  return out;
}

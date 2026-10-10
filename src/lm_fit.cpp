// [[Rcpp::depends(RcppArmadillo)]]

#include "lm_common.h"

/*
 * PER-COLUMN LINEAR-MODEL FITTER WITH PERMUTATION INFERENCE
 * =========================================================
 * fit_and_randomize() fits Y[, j] ~ X (| Z) for every column j of Y, estimates a linear contrast of the X
 * coefficients and tests it against relabelings supplied from R (perm_matrix). It is the kernel behind the
 * composition (CoDA), cell-density and cluster-free DE tests.
 *
 *  - NA-pattern grouping: columns of Y sharing a missingness pattern share one factorization of the design.
 *  - Nuisance covariates (Z): when Z is given, X and Y are residualized on Z per NA pattern (Frisch-Waugh-Lovell),
 *    and the relabelings act on the residualized response (Freedman-Lane). Every sample stays in the fit.
 *  - Missing responses: "drop" removes the rows for that column; "impute_weak" keeps them with weight na_weight
 *    and the response filled by the observed mean (or zero). Under a relabeling the weights stay with the response
 *    rows and the design rows are permuted instead.
 *  - Robust fits: Huber IRLS or one-step winsorization of the residuals (both re-estimated per relabeling).
 *  - Test statistic: "t" (contrast estimate divided by its standard error from the fit's own residual variance;
 *    pivotal, the default) or "coef" (the raw contrast estimate). The raw estimate and its standard error are
 *    returned in either case as the effect size.
 *  - Flattened parallelism: (column, pattern) jobs on the sccore thread pool; results do not depend on n_cores.
 */

/*** ===================================================================== ***/
/*** HELPERS                                                               ***/
/*** ===================================================================== ***/

static inline bool is_ill_conditioned(const arma::mat& X, double rcond_thresh) {
  if (X.n_rows < X.n_cols) return true;
  arma::mat XtX = X.t() * X;
  return (!XtX.is_finite() || arma::rcond(arma::symmatu(XtX)) < rcond_thresh);
}

// One NA pattern: the design seen by its columns and the permutation units.
struct DesignGroup {
  bool valid_design = false;
  bool can_permute = false;
  arma::uvec obs_indices;              // rows with an observed response
  arma::mat X_sub, B, invXtX, Xt;      // design of the observed fit (rows = units; residualized on Z; sqrt(w)-scaled under impute_weak)
  arma::mat X_full;                    // impute_weak: unscaled (residualized) design over all rows, permuted row-wise under relabeling
  arma::mat Qz;                        // orthonormal basis of Z over the units (empty without Z)
  arma::vec weights;                   // impute_weak: 1 for observed rows, na_weight for missing ones
  int df = 0;                          // residual degrees of freedom: n_obs - rank(X) - rank(Z)
  int rank = 0;                        // rank of the design over the units (may be below p when a level is absent)
  std::vector<arma::uvec> perm_blocks; // permutation cells in unit index space
  arma::uvec perm_units_global;        // global row ids of the units (for perm_matrix induction)
};

struct Job { arma::uword col_idx; int group_idx; };

// Standard error of the contrast: sqrt(sigma2 * c' (X'WX)^-1 c) with sigma2 = sum(w r^2) / df (w empty: unweighted).
// `invXtX`, when given, is (X'WX)^-1 of the design passed in.
static inline double contrast_se(const arma::mat& X, const arma::vec& y, const arma::vec& beta, const arma::vec& w,
                                 const arma::vec& c, int df, const arma::mat* invXtX = nullptr) {
  if (df < 1) return datum::nan;
  arma::vec r = y - X * beta;
  double rss = w.n_elem ? arma::dot(w, arma::square(r)) : arma::dot(r, r);
  double cvar;
  if (invXtX) cvar = arma::as_scalar(c.t() * (*invXtX) * c);
  else {
    arma::mat M = w.n_elem ? arma::mat(X.t() * (X.each_col() % w)) : arma::mat(X.t() * X);
    arma::vec v;
    if (!arma::solve(v, M, c, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx)) v = arma::pinv(M) * c;
    cvar = arma::dot(c, v);
  }
  double se2 = rss / df * cvar;
  return (se2 >= 0.0 && std::isfinite(se2)) ? std::sqrt(se2) : datum::nan;
}

// Studentized contrast. A zero standard error means a perfectly fitted (e.g. constant) response: the statistic is 0
// when the effect is at rounding level relative to the response scale `ysc`, NaN otherwise.
static inline double studentize(double effect, double se, double ysc) {
  if (se > 0.0 && std::isfinite(se)) return effect / se;
  return (std::abs(effect) <= 1e-10 * std::max(1.0, ysc)) ? 0.0 : datum::nan;
}

// --- Ordinary least squares on an explicit design ---
static inline arma::vec ols_free(const arma::mat& X, const arma::vec& y) {
  arma::mat XtX = X.t() * X; arma::vec Xty = X.t() * y; arma::vec b;
  if (!arma::solve(b, XtX, Xty, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx)) b = arma::pinv(XtX) * Xty;
  return b;
}

// --- Huber IRLS; `base_w` (may be empty) multiplies the robust weights; the final weights are returned in w_out ---
static arma::vec huber_irls(const arma::mat& X, const arma::vec& y, const arma::vec& base_w, const arma::vec& beta0,
                            const Config& cfg, arma::vec& w_out) {
  arma::vec beta = beta0;
  w_out = base_w.n_elem ? base_w : arma::vec(y.n_elem, arma::fill::ones);
  if (!beta.is_finite()) return beta;
  for (int it = 0; it < cfg.huber_maxit; ++it) {
    arma::vec r = y - X * beta;
    double s = robust_scale_mad(r);
    if (s <= 1e-12) break;
    double ks = cfg.huber_k * s;
    arma::vec w = arma::abs(r);
    w.transform([&](double val){ return (val > ks) ? (ks/val) : 1.0; });
    if (base_w.n_elem) w %= base_w;
    arma::mat XtWX = X.t() * (X.each_col() % w);
    arma::vec XtWy = X.t() * (y % w);
    double tr = arma::trace(XtWX);
    XtWX.diag() += 1e-10 * ((tr > 0.0) ? tr/X.n_cols : 1.0);    // mild ridge stabilization
    arma::mat Minv;
    if (!inv_sympd(Minv, XtWX)) Minv = arma::pinv(XtWX);
    arma::vec beta_new = Minv * XtWy;
    if (!beta_new.is_finite()) return beta;
    w_out = w;
    if (arma::norm(beta_new - beta)/(arma::norm(beta)+1e-12) < cfg.huber_tol) { beta = beta_new; break; }
    beta = beta_new;
  }
  return beta;
}

// --- One-step winsorized fit: OLS on the pseudo-response X beta0 + clamp(r); the pseudo-response is returned in ystar ---
static inline arma::vec winsor_fit(const arma::mat& X, const arma::vec& y, const arma::vec& beta0, double k,
                                   const arma::mat* B, arma::vec& ystar) {
  arma::vec r = y - X * beta0;
  double s = robust_scale_mad(r);
  if (s <= 1e-12) { ystar = y; return beta0; }
  r.clamp(-k*s, k*s);
  ystar = X * beta0 + r;
  return B ? arma::vec((*B) * ystar) : ols_free(X, ystar);
}

/*** ===================================================================== ***/
/*** FIT_AND_RANDOMIZE                                                     ***/
/*** ===================================================================== ***/
/*
 * @param X          design of interest (n x p); must carry the intercept unless Z does
 * @param Y          responses (n x m), NA allowed
 * @param contrast   contrast over the columns of X; effect = contrast' beta
 * @param Z          optional nuisance design (n x q, complete). When given, X and Y are residualized on Z per NA
 *                   pattern before fitting and permuting (Freedman-Lane / FWL); rank(Z) is removed from the df.
 * @param perm_groups permutation cells (list of 1-based row index vectors) the relabelings were drawn in
 * @param perm_matrix n x B matrix of 1-based relabelings in design convention (permuted design = X[p, ]); required
 *                   whenever permutations are requested; induced onto the observed rows of each NA pattern
 * @param statistic  "t" (default) or "coef"
 * @param alternative "two-sided" / "greater" / "less" (on the chosen statistic)
 * @param robust     "none", "huber" or "winsor"; huber_k, huber_maxit, huber_tol tune them
 * @param na_mode    "drop" or "impute_weak"; na_weight, na_center ("mean" or anything else = zero) for the latter
 * @param illcond_rcond, pinv_tol  conditioning guard and pseudo-inverse tolerance
 * @param return_residuals / return_sampled_fits / return_sampled_stats  optional outputs
 *
 * Returns a list: coef (p x m), effect (m; contrast' beta), se (m), df (m), stat (m; the chosen statistic),
 * p_value (m; add-one smoothing), residuals (n x m, NaN where missing), sampled_stats (B x m),
 * sampled_effects (B x m; only when statistic = "t"), sampled_fits (list of B x p).
 */
// [[Rcpp::export]]
Rcpp::List fit_and_randomize(const arma::mat& X, const arma::mat& Y, const arma::vec& contrast,
                             Rcpp::Nullable<Rcpp::NumericMatrix> Z = R_NilValue,
                             Rcpp::Nullable<Rcpp::List> perm_groups = R_NilValue,
                             int n_randomizations = 100,
                             std::string alternative = "two-sided",
                             bool return_residuals = true, bool return_sampled_fits = false, bool return_sampled_stats = false,
                             std::string robust = "none", double huber_k = 1.345, int huber_maxit = 8, double huber_tol = 1e-6,
                             std::string na_mode = "drop", double na_weight = 1e-4, std::string na_center = "mean",
                             double illcond_rcond = 1e-12, double pinv_tol = 0.0, int n_cores = 1,
                             Rcpp::Nullable<Rcpp::IntegerMatrix> perm_matrix = R_NilValue,
                             std::string statistic = "t") {

  Config cfg;
  cfg.robust = robust; cfg.na_mode = na_mode; cfg.huber_k = huber_k; cfg.huber_maxit = huber_maxit; cfg.huber_tol = huber_tol;
  cfg.na_weight = na_weight; cfg.center_mean = (na_center == "mean"); cfg.pinv_tol = pinv_tol; cfg.illcond_rcond = illcond_rcond;
  cfg.n_randomizations = n_randomizations; cfg.ret_res = return_residuals; cfg.ret_fits = return_sampled_fits; cfg.ret_stats = return_sampled_stats;
  cfg.alt_code = (alternative=="two-sided")?0 : (alternative=="greater")?1 : 2;
  if (statistic != "t" && statistic != "coef") Rcpp::stop("statistic must be \"t\" or \"coef\"");
  const bool use_t = (statistic == "t");
  if (robust != "none" && robust != "huber" && robust != "winsor") Rcpp::stop("robust must be \"none\", \"huber\" or \"winsor\"");
  if (na_mode != "drop" && na_mode != "impute_weak") Rcpp::stop("na_mode must be \"drop\" or \"impute_weak\"");

  const arma::uword n = X.n_rows, p = X.n_cols, m = Y.n_cols;
  if (Y.n_rows != n) stop("X and Y dimension mismatch");
  if (contrast.n_elem != p) stop("contrast must have one entry per column of X");

  // nuisance design
  arma::mat Zm;
  if (Z.isNotNull()) {
    Rcpp::NumericMatrix zz(Z);
    if (zz.ncol() > 0) {
      Zm = Rcpp::as<arma::mat>(zz);
      if (Zm.n_rows != n) stop("Z must have one row per row of X");
      if (!Zm.is_finite()) stop("Z must be complete");
    }
  }
  const bool has_Z = Zm.n_cols > 0;

  // 1. Permutation cells (1-based row index vectors) in which the R-drawn permutations were drawn
  std::vector<arma::uvec> blocks;
  if (perm_groups.isNotNull()) {
    Rcpp::List pg(perm_groups);
    for (int i=0; i<pg.size(); ++i) {
      IntegerVector g = pg[i]; arma::uvec ug = as<arma::uvec>(g);
      if (ug.n_elem && ug.max() > 0) ug -= 1; blocks.push_back(std::move(ug));
    }
  } else {
    blocks.push_back(arma::regspace<arma::uvec>(0, n - 1));   // one cell: every row exchangeable
  }

  // Permutations supplied from R: n x B, 1-based, design convention. Every column of Y sees the same relabeling b.
  const bool has_P = perm_matrix.isNotNull();
  arma::umat Pm;
  if (has_P) {
    Rcpp::IntegerMatrix pm(perm_matrix);
    if ((arma::uword)pm.nrow() != n) Rcpp::stop("perm_matrix must have %d rows (one per unit), got %d", (int)n, pm.nrow());
    Pm.set_size(pm.nrow(), pm.ncol());
    for (int b = 0; b < pm.ncol(); ++b) for (int i = 0; i < pm.nrow(); ++i) {
      int v = pm(i, b) - 1;
      if (v < 0 || v >= (int)n) Rcpp::stop("perm_matrix has an index outside 1..%d", (int)n);
      Pm(i, b) = (arma::uword)v;
    }
    n_randomizations = (int)Pm.n_cols; cfg.n_randomizations = n_randomizations;
  } else if (n_randomizations > 0) Rcpp::stop("perm_matrix is required for permutations (draw it with drawPermutations())");

  // 2. Group columns by NA pattern
  std::unordered_map<std::uint64_t, std::vector<arma::uword>> map_mask;
  for (arma::uword j=0; j<m; ++j) map_mask[hash_vec_mask(Y.col(j))].push_back(j);

  std::vector<DesignGroup> designs;
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

  const bool is_drop = (cfg.na_mode == "drop");
  const bool is_huber = (cfg.robust == "huber"), is_winsor = (cfg.robust == "winsor");

  // 3. Per-pattern designs (parallel)
  cacoa::parallelFor(0, grp_cnt, [&](int i) {
    DesignGroup& g = designs[i]; const RawGroup& raw = raw_groups[i];
    g.obs_indices = raw.obs;
    arma::uword qz = 0;
    if (raw.obs.n_elem < 2) return;

    if (is_drop) {
      g.X_sub = X.rows(raw.obs);
      if (has_Z) {
        if (!qr_basis(Zm.rows(raw.obs), g.Qz)) return;
        qz = g.Qz.n_cols;
        if (raw.obs.n_elem <= qz + 1) return;
        g.X_sub = project_out_Q(g.Qz, g.X_sub);
      }
      g.perm_blocks = subset_blocks(blocks, g.obs_indices, n);
      if (g.perm_blocks.empty() && g.obs_indices.n_elem > 0) g.perm_blocks.push_back(arma::regspace<arma::uvec>(0, g.obs_indices.n_elem - 1));
      g.perm_units_global = g.obs_indices;
    } else {
      g.weights.set_size(n); g.weights.fill(cfg.na_weight); g.weights.elem(raw.obs).fill(1.0);
      g.X_full = X;
      if (has_Z) {
        if (!qr_basis(Zm, g.Qz)) return;
        qz = g.Qz.n_cols;
        if (raw.obs.n_elem <= qz + 1) return;
        g.X_full = project_out_Q(g.Qz, X);
      }
      g.X_sub = g.X_full.each_col() % arma::sqrt(g.weights);
      g.perm_blocks = blocks; g.perm_units_global = arma::regspace<arma::uvec>(0, n - 1);
    }
    g.can_permute = (raw.obs.n_elem >= 2);

    // columns that vanish on these units (an absent level, or a direction entirely absorbed by Z) are made exactly zero,
    // so that the estimability check below sees them as missing rather than as numerical noise
    { double sc = 1e-300; for (arma::uword k = 0; k < p; ++k) sc = std::max(sc, arma::norm(X.col(k)));
      for (arma::uword k = 0; k < p; ++k) if (arma::norm(g.X_sub.col(k)) <= 1e-8 * sc) g.X_sub.col(k).zeros(); }

    // Factorization. A full-rank design uses (X'X)^-1; a rank-deficient one (a factor level or a covariate pattern
    // absent among this pattern's units) uses the pseudo-inverse, provided the contrast stays estimable (c in the
    // row space of X): then c'beta is unique and c' (X'X)^+ c is its variance.
    if (g.X_sub.n_rows == 0) return;
    if (!is_ill_conditioned(g.X_sub, cfg.illcond_rcond) && inv_xtx_safe(g.X_sub, g.invXtX, g.Xt, cfg.pinv_tol)) {
      g.B = g.invXtX * g.Xt; g.rank = (int)p; g.valid_design = true;
    } else {
      arma::mat Xp;
      if (!arma::pinv(Xp, g.X_sub) || !Xp.is_finite()) return;
      arma::mat Proj = Xp * g.X_sub;                                   // projector onto the row space of X
      double cn = arma::norm(contrast);
      if (cn <= 0.0 || arma::norm(contrast - Proj * contrast) > 1e-8 * cn) return;   // not estimable
      g.B = Xp; g.Xt = g.X_sub.t(); g.invXtX = Xp * Xp.t(); g.rank = (int)arma::rank(g.X_sub); g.valid_design = true;
    }
    g.df = (int)raw.obs.n_elem - g.rank - (int)qz;
    if (g.df < 1) { g.valid_design = false; return; }
  }, n_cores, false);

  // 4. Fits and permutations (flattened parallelism over columns)
  arma::mat Coef(p, m); Coef.fill(datum::nan);
  arma::vec Effect(m), Se(m), Stat(m), Pval(m); Effect.fill(datum::nan); Se.fill(datum::nan); Stat.fill(datum::nan); Pval.fill(datum::nan);
  arma::ivec Df(m); Df.fill(0);
  arma::mat Resid; if (cfg.ret_res) { Resid.set_size(n, m); Resid.fill(datum::nan); }
  std::vector<arma::mat> SampledFits(m);
  arma::mat SampledStats, SampledEffects;
  if (cfg.ret_stats && n_randomizations > 0) {
    SampledStats.set_size(n_randomizations, m); SampledStats.fill(datum::nan);
    if (use_t) { SampledEffects.set_size(n_randomizations, m); SampledEffects.fill(datum::nan); }
  }
  const arma::vec no_w;

  cacoa::parallelFor(0, (int)jobs.size(), [&](int k) {
    const Job& job = jobs[k]; const DesignGroup& grp = designs[job.group_idx];
    const arma::uword j = job.col_idx;
    if (!grp.valid_design) return;
    Df[j] = grp.df;

    // response for this column: observed rows (drop) or all rows with the missing ones filled (impute_weak);
    // residualized on Z when Z is given
    arma::vec y_raw = Y.col(j);
    const arma::uvec miss = arma::find_nonfinite(y_raw);
    auto fill_value = [&](const arma::vec& v) { return (cfg.center_mean && grp.obs_indices.n_elem) ? arma::mean(v.elem(grp.obs_indices)) : 0.0; };
    arma::vec y_work, y_clean;
    if (is_drop) {
      y_work = y_raw.elem(grp.obs_indices);
      if (has_Z) y_work = project_out_Q(grp.Qz, y_work);
    } else {
      y_clean = y_raw; y_clean.elem(miss).fill(fill_value(y_raw));
      if (has_Z) {                                   // project the filled vector, then re-fill the missing rows from the residualized observed ones
        y_clean = project_out_Q(grp.Qz, y_clean);
        arma::vec tmp = y_clean; tmp.elem(miss).fill(datum::nan);
        y_clean.elem(miss).fill(fill_value(tmp));
      }
      y_work = y_clean;
      if (!is_huber) y_work %= arma::sqrt(grp.weights);
    }
    const bool huber_imp = is_huber && !is_drop;
    const arma::mat& Xobs = huber_imp ? grp.X_full : grp.X_sub;     // design matching y_work
    const arma::vec& wobs = huber_imp ? grp.weights : no_w;         // explicit weights only for huber under impute_weak

    // observed fit
    arma::vec beta, w_fin, ystar;
    if (is_huber) {
      arma::vec b0 = huber_imp ? arma::vec(grp.B * arma::vec(y_work % arma::sqrt(grp.weights))) : arma::vec(grp.B * y_work);
      beta = huber_irls(Xobs, y_work, wobs, b0, cfg, w_fin);
    } else if (is_winsor) beta = winsor_fit(Xobs, y_work, arma::vec(grp.B * y_work), cfg.huber_k, &grp.B, ystar);
    else beta = grp.B * y_work;
    if (!beta.is_finite()) return;
    const double effect = arma::dot(beta, contrast);
    const double se = is_huber ? contrast_se(Xobs, y_work, beta, w_fin, contrast, grp.df)
                               : contrast_se(Xobs, is_winsor ? ystar : y_work, beta, no_w, contrast, grp.df, &grp.invXtX);
    const double ysc = y_work.n_elem ? arma::abs(y_work).max() : 0.0;
    const double stat_obs = use_t ? studentize(effect, se, ysc) : effect;
    Coef.col(j) = beta; Effect[j] = effect; Se[j] = se; Stat[j] = stat_obs;

    if (cfg.ret_res) {
      arma::vec r_out(n); r_out.fill(datum::nan);
      if (is_drop) r_out.elem(grp.obs_indices) = y_work - grp.X_sub * beta;
      else { r_out = y_clean - grp.X_full * beta; r_out.elem(miss).fill(datum::nan); }
      Resid.col(j) = r_out;
    }

    if (cfg.n_randomizations == 0 || !grp.can_permute || !std::isfinite(stat_obs)) return;

    // permutation loop
    arma::vec stats_perm(cfg.n_randomizations), effects_perm(cfg.n_randomizations);
    arma::mat fits_perm; if (cfg.ret_fits) fits_perm.set_size(cfg.n_randomizations, p);
    arma::vec alpha; if (!is_huber && !is_winsor && !cfg.ret_fits && !use_t) alpha = grp.B.t() * contrast;   // raw-coefficient fast path
    const arma::vec sw = is_drop ? arma::vec() : arma::vec(arma::sqrt(grp.weights));

    int ge=0, le=0, ge_abs=0;
    for (int r=0; r<cfg.n_randomizations; ++r) {
      arma::uvec q = induced_perm(Pm.col(r), grp.perm_units_global, grp.perm_blocks);   // design convention over this pattern's units
      double eff_p = 0.0, s_p = 0.0; arma::vec b_p;
      if (is_drop) {
        arma::vec y_perm = y_work.elem(inverse_perm(q));                                  // response aligned with X_sub[q, ]
        if (is_huber) {
          arma::vec wf; b_p = huber_irls(grp.X_sub, y_perm, no_w, arma::vec(grp.B * y_perm), cfg, wf);
          eff_p = arma::dot(b_p, contrast);
          s_p = use_t ? studentize(eff_p, contrast_se(grp.X_sub, y_perm, b_p, wf, contrast, grp.df), ysc) : eff_p;
        } else if (is_winsor) {
          arma::vec ys; b_p = winsor_fit(grp.X_sub, y_perm, arma::vec(grp.B * y_perm), cfg.huber_k, &grp.B, ys);
          eff_p = arma::dot(b_p, contrast);
          s_p = use_t ? studentize(eff_p, contrast_se(grp.X_sub, ys, b_p, no_w, contrast, grp.df, &grp.invXtX), ysc) : eff_p;
        } else if (alpha.n_elem) {
          eff_p = s_p = arma::dot(alpha, y_perm);
        } else {
          b_p = grp.B * y_perm; eff_p = arma::dot(b_p, contrast);
          s_p = use_t ? studentize(eff_p, contrast_se(grp.X_sub, y_perm, b_p, no_w, contrast, grp.df, &grp.invXtX), ysc) : eff_p;
        }
      } else {
        // impute_weak: the weights stay with the response rows; the (residualized) design rows are permuted instead
        arma::mat Xq = grp.X_full.rows(q);
        arma::mat Xqw = Xq.each_col() % sw; arma::vec yw = y_clean % sw;
        if (is_huber) {
          arma::vec wf; b_p = huber_irls(Xq, y_clean, grp.weights, ols_free(Xqw, yw), cfg, wf);
          eff_p = arma::dot(b_p, contrast);
          s_p = use_t ? studentize(eff_p, contrast_se(Xq, y_clean, b_p, wf, contrast, grp.df), ysc) : eff_p;
        } else if (is_winsor) {
          arma::vec ys; b_p = winsor_fit(Xqw, yw, ols_free(Xqw, yw), cfg.huber_k, nullptr, ys);
          eff_p = arma::dot(b_p, contrast);
          s_p = use_t ? studentize(eff_p, contrast_se(Xqw, ys, b_p, no_w, contrast, grp.df), ysc) : eff_p;
        } else {
          b_p = ols_free(Xqw, yw); eff_p = arma::dot(b_p, contrast);
          s_p = use_t ? studentize(eff_p, contrast_se(Xqw, yw, b_p, no_w, contrast, grp.df), ysc) : eff_p;
        }
      }
      stats_perm[r] = s_p; effects_perm[r] = eff_p;
      if (std::isfinite(s_p)) {
        if (std::abs(s_p) >= std::abs(stat_obs)) ge_abs++;
        if (s_p >= stat_obs) ge++;
        if (s_p <= stat_obs) le++;
      }
      if (cfg.ret_fits) fits_perm.row(r) = b_p.t();
    }

    double pval = (cfg.alt_code==0)?(ge_abs+1.0):(cfg.alt_code==1)?(ge+1.0):(le+1.0);
    Pval[j] = pval / (cfg.n_randomizations + 1.0);
    if (cfg.ret_stats) { SampledStats.col(j) = stats_perm; if (use_t) SampledEffects.col(j) = effects_perm; }
    if (cfg.ret_fits) SampledFits[j] = fits_perm;
  }, n_cores, false);

  if (cfg.ret_fits) for (uword j = 0; j < m; ++j) if (SampledFits[j].n_elem == 0) { SampledFits[j].set_size(std::max(n_randomizations, 0), p); SampledFits[j].fill(datum::nan); }   // failed / unpermuted columns

  Rcpp::List out = Rcpp::List::create(_["coef"]=Coef, _["effect"]=Effect, _["se"]=Se, _["df"]=Df, _["stat"]=Stat, _["p_value"]=Pval);
  if (cfg.ret_res) out["residuals"] = Resid;
  if (cfg.ret_stats && n_randomizations > 0) { out["sampled_stats"] = SampledStats; if (use_t) out["sampled_effects"] = SampledEffects; }
  if (cfg.ret_fits) {
    Rcpp::List L(m); for(uword j=0; j<m; ++j) L[j]=Rcpp::wrap(SampledFits[j]);
    out["sampled_fits"] = L;
  }
  return out;
}

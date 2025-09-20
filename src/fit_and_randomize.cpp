// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ---------- utilities ----------

// In-place Fisher–Yates shuffle (R RNG)
static inline void shuffle_vec_in_place(arma::vec& v) {
  for (arma::uword i = v.n_elem; i > 1; --i) {
    arma::uword j = (arma::uword) std::floor(R::unif_rand() * static_cast<double>(i)); // [0,i)
    if (j >= i) j = i - 1;
    std::swap(v[i - 1], v[j]);
  }
}

// Permute within groups (0-based). Groups <2 unchanged.
static inline void permute_by_groups(arma::vec& y, const std::vector<std::vector<arma::uword>>& groups) {
  for (const auto& g : groups) {
    if (g.size() < 2) continue;
    std::vector<double> tmp; tmp.reserve(g.size());
    for (arma::uword k : g) tmp.push_back(y[k]);
    for (std::size_t i = tmp.size(); i > 1; --i) {
      std::size_t j = (std::size_t) std::floor(R::unif_rand() * (double)i);
      if (j >= i) j = i - 1;
      std::swap(tmp[i - 1], tmp[j]);
    }
    for (std::size_t i = 0; i < g.size(); ++i) y[g[i]] = tmp[i];
  }
}

// Build 0-based groups once (full data)
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

// Map full 1-based groups to subset positions (0-based) for NA columns
static std::vector<std::vector<arma::uword>>
  map_groups_to_subset(const IntegerVector& pos_full, const List& perm_groups_full) {
    std::vector<std::vector<arma::uword>> out;
    if (perm_groups_full.size() == 0) return out;
    out.reserve(perm_groups_full.size());
    for (int g = 0; g < perm_groups_full.size(); ++g) {
      IntegerVector grp = perm_groups_full[g];
      std::vector<arma::uword> mapped; mapped.reserve(grp.size());
      for (int a = 0; a < grp.size(); ++a) {
        int i1 = grp[a];
        if (i1 <= 0 || i1 > pos_full.size()) continue;
        int pos = pos_full[i1 - 1]; // -1 if dropped
        if (pos >= 0) mapped.push_back((arma::uword)pos);
      }
      out.push_back(std::move(mapped));
    }
    return out;
  }

// inv(X'X) via SPD fast path; else pseudo-inverse; also return Xt
static inline bool inv_xtx(const arma::mat& X, arma::mat& invXtX, arma::mat& Xt) {
  Xt = X.t();
  arma::mat XtX = Xt * X;                   // p x p
  bool ok = inv_sympd(invXtX, XtX);         // fast if SPD
  if (!ok) invXtX = pinv(XtX);              // MP pseudo-inverse
  return ok;
}

// OLS coef from precomputed inv(X'X) and X'
static inline arma::vec coef_from_inv(const arma::mat& invXtX, const arma::mat& Xt, const arma::vec& y) {
  return invXtX * (Xt * y);
}

static inline double contrast_stat(const arma::vec& beta, const arma::vec& contrast) {
  return dot(beta, contrast);
}

static inline void tally_perm(double obs, double perm, int alt,
                              int& ge, int& le, int& ge_abs) {
  if (alt == 0) { // two-sided
    if (std::fabs(perm) >= std::fabs(obs)) ge_abs++;
  } else if (alt == 1) { // greater
    if (perm >= obs) ge++;
  } else { // less
    if (perm <= obs) le++;
  }
}

// Robust scale via MAD (median absolute deviation)
// Median absolute deviation with NA/NaN/Inf safety.
// Returns a strictly positive scale (fallback to RMS + tiny eps if needed).
static inline double robust_scale_mad_safe(const arma::vec& r) {
  arma::uvec idx = arma::find_finite(r);
  if (idx.n_elem == 0) return 1e-8;                // nothing finite → tiny scale
  
  arma::vec rf = r.elem(idx);
  double med = arma::median(rf);                    // finite by construction
  
  arma::vec af = arma::abs(rf - med);
  double mad = (af.n_elem > 0) ? arma::median(af) : 0.0;
  
  double s = 1.4826 * mad;
  if (!(s > 0.0)) {
    // fallback: RMS of finite residuals + tiny epsilon
    s = std::sqrt(arma::mean(arma::square(rf))) + 1e-12;
    if (!(s > 0.0)) s = 1e-8;
  }
  return s;
}


// One Huber IRLS fit: returns beta
// - Uses MAD scale each iteration
// - Starts at OLS via invXtX/Xt if provided (fast), else from scratch
static arma::vec huber_irls(const arma::mat& X,
                            const arma::vec& y,
                            double k, int maxit, double tol,
                            const arma::mat* invXtX_opt = nullptr,
                            const arma::mat* Xt_opt     = nullptr) {
  arma::uword p = X.n_cols;
  arma::vec beta(p, fill::zeros);
  
  if (invXtX_opt && Xt_opt) {
    beta = (*invXtX_opt) * ((*Xt_opt) * y); // OLS start
  } else {
    // OLS from scratch
    arma::mat Xt = X.t();
    arma::mat XtX = Xt * X;
    arma::vec Xty = Xt * y;
    arma::mat M;
    if (inv_sympd(M, XtX)) beta = M * Xty; else beta = pinv(XtX) * Xty;
  }
  
  for (int it = 0; it < maxit; ++it) {
    arma::vec r = y - X * beta;
    double s = robust_scale_mad(r);
    if (!(s > 0.0)) break;
    
    double ks = k * s;
    arma::vec w(r.n_elem, fill::ones);
    for (arma::uword i = 0; i < r.n_elem; ++i) {
      double ar = std::fabs(r[i]);
      if (ar > ks) w[i] = ks / ar; // Huber weight
      // else 1.0
    }
    
    // Weighted normal equations: (X' W X) beta = X' W y
    arma::mat XtWX = X.t() * (X.each_col() % w); // p x p
    arma::vec XtWy = X.t() * (y % w);            // p x 1
    arma::mat Minv;
    arma::vec beta_new;
    if (inv_sympd(Minv, XtWX)) beta_new = Minv * XtWy; else beta_new = pinv(XtWX) * XtWy;
    
    double denom = arma::norm(beta, 2) + 1e-12;
    double rel_change = arma::norm(beta_new - beta, 2) / denom;
    beta = beta_new;
    if (rel_change < tol) break;
  }
  return beta;
}

// ---------- main ----------

//' Fast multi-response OLS / Huber-IRLS + (blocked) permutations + contrast p-values
 //'
 //' Optimized for NA-free columns (reuses full-design precomputations).
 //' Full permutation (no blocks) uses in-place shuffle.
 //'
 //' @param X (n x p) design matrix
 //' @param Y (n x m) response matrix; NAs allowed and dropped per column
 //' @param contrast length-p vector (contrast over coefficients of X)
 //' @param perm_groups optional list of 1-based integer vectors; permutation blocks.
 //'        If NULL, all rows form one big block (full permutation).
 //' @param n_randomizations number of permutations (>=1 recommended)
 //' @param alternative "two-sided", "greater", or "less"
 //' @param return_residuals return (n x m) residual matrix (NA where dropped)
 //' @param return_sampled_fits if TRUE, per-column list of (B x p) coef matrices
 //' @param use_huber if TRUE, use Huber M-estimation (IRLS) instead of OLS
 //' @param huber_k Huber tuning constant (default 1.345 ~ 95% Gaussian efficiency)
 //' @param huber_maxit maximum IRLS iterations (e.g., 5–10 is typical)
 //' @param huber_tol relative coefficient tolerance to stop early
 //'
 //' @return list: coef (p x m), stat (m), p_value (m), residuals?, sampled_fits?
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
                         double huber_tol = 1e-6) {
   RNGScope scope; // R RNG
   
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
   
   // Precompute full-design pieces (NA-free fast path)
   arma::mat invXtX_full, Xt_full;
   inv_xtx(X, invXtX_full, Xt_full);
   
   // Full-data permutation groups
   List perm_full = perm_groups.isNotNull() ? List(perm_groups) : List();
   const bool full_permute_requested = (perm_full.size() == 0);
   std::vector<std::vector<arma::uword>> groups_full0; // 0-based, blocked mode
   bool any_movable_full = (n >= 2);
   if (!full_permute_requested) {
     groups_full0 = build_groups0_full(perm_full, n);
     any_movable_full = false;
     for (const auto& g : groups_full0) if (g.size() >= 2) { any_movable_full = true; break; }
   }
   
   // Outputs
   arma::mat coef_obs(p, m); coef_obs.fill(arma::datum::nan);
   arma::vec stat_obs(m);    stat_obs.fill(arma::datum::nan);
   arma::vec pvals(m);       pvals.fill(arma::datum::nan);
   arma::mat resid_out;
   if (return_residuals) { resid_out.set_size(n, m); resid_out.fill(arma::datum::nan); }
   Rcpp::List sampled_fits; if (return_sampled_fits) sampled_fits = Rcpp::List(m);
   
   // Per-column loop
   for (arma::uword j = 0; j < m; ++j) {
     arma::vec yj = Y.col(j);
     
     // Detect NA-free column (fast path)
     bool all_finite = true;
     for (arma::uword i = 0; i < n; ++i) { if (!arma::is_finite(yj[i])) { all_finite = false; break; } }
     
     // Observed fit setup
     const arma::mat* Xobs = &X;
     arma::vec yc;
     arma::uword k = n;
     arma::mat invXtX, Xt;            // may alias full precomputes
     bool use_full_shuffle = full_permute_requested;
     std::vector<std::vector<arma::uword>> groups_sub;
     
     if (all_finite) {
       yc = yj;
       invXtX = invXtX_full;  Xt = Xt_full;  // reuse full precomputes
       if (!full_permute_requested) {
         use_full_shuffle = false;
         groups_sub = groups_full0;
         if (!any_movable_full) {
           // nothing to permute; we still compute observed fit below
           pvals[j] = 1.0;
           if (return_sampled_fits) sampled_fits[j] = NumericMatrix(0, p);
         }
       }
     } else {
       // subset rows for this column
       std::vector<arma::uword> keep_rows; keep_rows.reserve(n);
       for (arma::uword i = 0; i < n; ++i) if (arma::is_finite(yj[i])) keep_rows.push_back(i);
       if (keep_rows.empty()) {
         pvals[j] = NA_REAL;
         if (return_sampled_fits) sampled_fits[j] = NumericMatrix(0, p);
         continue;
       }
       arma::uvec idx = arma::uvec(keep_rows);
       arma::mat Xc   = X.rows(idx);
       yc             = yj.elem(idx);
       k              = Xc.n_rows;
       Xobs           = &Xc;
       inv_xtx(Xc, invXtX, Xt);
       
       if (!full_permute_requested) {
         IntegerVector pos_full(n); std::fill(pos_full.begin(), pos_full.end(), -1);
         for (arma::uword t = 0; t < k; ++t) pos_full[(int)idx[t]] = (int)t;
         use_full_shuffle = false;
         groups_sub = map_groups_to_subset(pos_full, perm_full);
         bool any_movable = false; for (const auto& g : groups_sub) if (g.size() >= 2) { any_movable = true; break; }
         if (!any_movable) {
           pvals[j] = 1.0;
           if (return_sampled_fits) sampled_fits[j] = NumericMatrix(0, p);
         }
       }
     }
     
     // Observed coefficients & residuals
     arma::vec beta;
     if (!use_huber) {
       beta = coef_from_inv(invXtX, Xt, yc);                         // OLS
     } else {
       beta = huber_irls(*Xobs, yc, huber_k, huber_maxit, huber_tol, // Huber
                         &invXtX, &Xt);  // pass OLS start cheaply
     }
     double sobs = contrast_stat(beta, contrast);
     
     coef_obs.col(j) = beta;
     stat_obs[j]     = sobs;
     
     if (return_residuals) {
       arma::vec rc = yc - (*Xobs) * beta;
       if (all_finite) {
         resid_out.col(j) = rc;
       } else {
         // restore into full residual matrix
         arma::uword placed = 0;
         for (arma::uword i = 0; i < n; ++i) {
           if (arma::is_finite(yj[i])) { resid_out(i, j) = rc[placed++]; }
         }
       }
     }
     
     // Permutations
     if (n_randomizations == 0) {
       pvals[j] = NA_REAL;
       if (return_sampled_fits) sampled_fits[j] = NumericMatrix(0, p);
       continue;
     }
     if (k < 2) {
       pvals[j] = 1.0;
       if (return_sampled_fits) sampled_fits[j] = NumericMatrix(0, p);
       continue;
     }
     
     int ge = 0, le = 0, ge_abs = 0;
     arma::vec yperm = yc;
     NumericMatrix Mcoef; if (return_sampled_fits) Mcoef = NumericMatrix(n_randomizations, p);
     
     for (int r = 0; r < n_randomizations; ++r) {
       yperm = yc;
       if (use_full_shuffle) {
         shuffle_vec_in_place(yperm);
       } else {
         permute_by_groups(yperm, groups_sub);
       }
       
       arma::vec bperm;
       if (!use_huber) {
         bperm = coef_from_inv(invXtX, Xt, yperm);                                // OLS
       } else {
         bperm = huber_irls(*Xobs, yperm, huber_k, huber_maxit, huber_tol,        // Huber
                            &invXtX, &Xt);
       }
       double sperm = contrast_stat(bperm, contrast);
       
       if (alt == 0)      tally_perm(sobs, sperm, 0, ge, le, ge_abs);
       else if (alt == 1) tally_perm(sobs, sperm, 1, ge, le, ge_abs);
       else               tally_perm(sobs, sperm, 2, ge, le, ge_abs);
       
       if (return_sampled_fits) {
         for (arma::uword q = 0; q < p; ++q) Mcoef(r, q) = bperm[q];
       }
     }
     
     double pval;
     if (alt == 0) pval = (ge_abs + 1.0) / (n_randomizations + 1.0);
     else if (alt == 1) pval = (ge + 1.0) / (n_randomizations + 1.0);
     else               pval = (le + 1.0) / (n_randomizations + 1.0);
     
     pvals[j] = pval;
     if (return_sampled_fits) sampled_fits[j] = Mcoef;
   }
   
   Rcpp::List out = Rcpp::List::create(
     _["coef"]    = coef_obs,
     _["stat"]    = stat_obs,
     _["p_value"] = pvals
   );
   if (return_residuals)    out["residuals"]    = resid_out;
   if (return_sampled_fits) out["sampled_fits"] = sampled_fits;
   return out;
 }

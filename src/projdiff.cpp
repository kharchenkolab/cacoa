// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <limits>
using namespace arma;
using namespace Rcpp;



// determine mean difference between groups of columns, distances on that projection
// g1, g2 - vectors of indices for the two groups (0-based)
// [[Rcpp::export]]
arma::rowvec projdiff(const arma::mat & mat, const arma::ivec & g1, const arma::ivec & g2) {
  // TODO: remove with common expression shifts

  // mat: columns - samples; rows - genes;
  // determine consensus difference
  arma::mat dm(mat.n_rows,g1.n_elem * g2.n_elem);

  for(int i=0;i<g1.n_elem;i++) {
    for(int j=0;j<g2.n_elem;j++) {
      dm.col(i*g2.n_elem + j) = mat.col(g1[i])-mat.col(g2[j]);
    }
  }
  arma::vec dmm=mean(dm,1); // would be nice to do trimming here
  dmm/=sqrt(sum(dmm % dmm));
  return dmm.t() * dm;

}

// fit linear model for density differences
struct fit_density_lm_result {
  arma::mat  KP;
  arma::mat  KPge;
  arma::cube KP_perm;
  int        n_randomizations;
};

inline fit_density_lm_result fit_density_lm_impl(const arma::mat& M,
                                                 const arma::mat& P,
                                                 int n_randomizations) {
  arma::mat K  = arma::solve(M.t() * M, M.t());
  arma::mat KP = K * P;

  arma::mat  KPge(KP.n_rows, KP.n_cols, arma::fill::zeros);
  arma::cube KP_perm(KP.n_rows, KP.n_cols, n_randomizations);

  for (int r = 0; r < n_randomizations; ++r) {
    arma::mat P_prime = P.rows(arma::randperm(P.n_rows));
    arma::mat KP_prime = K * P_prime;
    KP_perm.slice(r) = KP_prime;
    KPge += arma::conv_to<arma::mat>::from(KP_prime >= KP);
  }
  return {KP, KPge, KP_perm, n_randomizations};
}

// [[Rcpp::export]]
Rcpp::List fit_density_lm(const arma::mat& M,
                          const arma::mat& P,
                          int n_randomizations) {
  auto res = fit_density_lm_impl(M, P, n_randomizations);
  return Rcpp::List::create(
    Rcpp::Named("KP") = res.KP,
    Rcpp::Named("KPge") = res.KPge,
    Rcpp::Named("KP_perm") = res.KP_perm,
    Rcpp::Named("n_randomizations") = res.n_randomizations
  );
}



// permute row indices within blocks (0..L-1); length 0 => one big block
static inline uvec permute_within_blocks(const uvec& blocks, uword n) {
  if (blocks.n_elem == 0) {
    uvec idx = regspace<uvec>(0, n - 1);
    return shuffle(idx);
  }
  uvec out(n);
  uvec levs = unique(blocks);
  for (uword k = 0; k < levs.n_elem; ++k) {
    uword lev = levs[k];
    uvec where = find(blocks == lev);
    uvec shuffled = shuffle(where);
    out.elem(where) = shuffled;
  }
  return out;
}

// ---------------- FULL model: Y is n x p, F is n x q ----------------
// [[Rcpp::export]]
List perm_full_contrast_mat(const arma::mat& F,                // n x q
                            const arma::mat& Y,                // n x p
                            const arma::vec& contrastF,        // length q (aligned to cols(F))
                            const arma::uvec& blocks,          // length n (0-based) or length 0
                            const int B) {
  const uword n = F.n_rows, q = F.n_cols, p = Y.n_cols;
  if (Y.n_rows != n) stop("Y and F must have same nrow.");

  // precompute normal-equation pieces once
  mat XtX = F.t() * F;            // q x q
  mat XtX_inv = inv_sympd(XtX);   // robust SPD inverse; falls back if not SPD
  // If XtX might be singular/non-SPD, use pinv(XtX) instead:
  // mat XtX_inv = pinv(XtX);

  mat XtY = F.t() * Y;            // q x p
  mat Beta = XtX_inv * XtY;       // q x p

  // observed stats: c' * Beta  (1 x p) -> vector length p
  rowvec cF = contrastF.t();      // 1 x q
  rowvec s_obs_row = cF * Beta;   // 1 x p
  vec stat_obs = s_obs_row.t();   // p

  // permutations
  mat stats_perm(B, p, fill::none);
  stats_perm.fill(datum::nan);

  for (int b = 0; b < B; ++b) {
    uvec idx = permute_within_blocks(blocks, n);
    mat Yb = Y.rows(idx);                   // permute rows across all columns
    mat Beta_b = XtX_inv * (F.t() * Yb);    // q x p
    stats_perm.row(b) = (cF * Beta_b);      // 1 x p
  }

  // p-values per column (two-sided, +1 correction)
  vec pval(p, fill::value(datum::nan));
  for (uword j = 0; j < p; ++j) {
    vec col = stats_perm.col(j);
    uvec ok = find_finite(col);
    if (ok.n_elem == 0 || !std::isfinite(stat_obs[j])) continue;
    vec sp = col.elem(ok);
    uword ge = accu(abs(sp) >= std::abs(stat_obs[j]));
    pval[j] = (1.0 + (double)ge) / (1.0 + (double)sp.n_elem);
  }

  return List::create(
    _["stat_obs"]   = stat_obs,    // length p
    _["stats_perm"] = stats_perm,  // B x p
    _["pval"]       = pval
  );
}

// ------------- Freedman–Lane: Xr is n x qx, Yr is n x p (already residualized in R) -------------
// [[Rcpp::export]]
List perm_FL_contrast_mat(const arma::mat& Xr,                 // n x qx
                          const arma::mat& Yr,                 // n x p
                          const arma::vec& contrastX,          // length qx (aligned to cols(Xr))
                          const arma::uvec& blocks,            // length n (0-based) or length 0
                          const int B) {
  const uword n = Xr.n_rows, qx = Xr.n_cols, p = Yr.n_cols;
  if (Yr.n_rows != n) stop("Yr and Xr must have same nrow.");

  mat XtX = Xr.t() * Xr;            // qx x qx
  mat XtX_inv = inv_sympd(XtX);
  // Or: mat XtX_inv = pinv(XtX);

  mat XtY = Xr.t() * Yr;            // qx x p
  mat Beta = XtX_inv * XtY;         // qx x p

  rowvec cX = contrastX.t();        // 1 x qx
  rowvec s_obs_row = cX * Beta;     // 1 x p
  vec stat_obs = s_obs_row.t();     // p

  mat stats_perm(B, p, fill::none);
  stats_perm.fill(datum::nan);

  for (int b = 0; b < B; ++b) {
    uvec idx = permute_within_blocks(blocks, n);
    mat Yb = Yr.rows(idx);
    mat Beta_b = XtX_inv * (Xr.t() * Yb);
    stats_perm.row(b) = (cX * Beta_b);
  }

  vec pval(p, fill::value(datum::nan));
  for (uword j = 0; j < p; ++j) {
    vec col = stats_perm.col(j);
    uvec ok = find_finite(col);
    if (ok.n_elem == 0 || !std::isfinite(stat_obs[j])) continue;
    vec sp = col.elem(ok);
    uword ge = accu(abs(sp) >= std::abs(stat_obs[j]));
    pval[j] = (1.0 + (double)ge) / (1.0 + (double)sp.n_elem);
  }

  return List::create(
    _["stat_obs"]   = stat_obs,
    _["stats_perm"] = stats_perm,
    _["pval"]       = pval
  );
}
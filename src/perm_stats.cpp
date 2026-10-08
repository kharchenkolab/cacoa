// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Per-permutation statistics of the pairwise-effects engine (Track D): the pivotal F of a single contrast
// under a relabeling p of the samples, against a fixed Gower matrix G.
//
//   F_p = [ a_p' G a_p / c'(X'X)^- c ] / [ (tr G - sum(H_p o G)) / (n - q) ],   a_p = a[p],  H_p = H[p, p]
//
// where a = X (X'X)^- c are the regression weights of the contrast and H the hat matrix of X (both row-permuted
// together with the design). The R reference implementation is permutedStats() in R/pairwise_inference.R.

static inline arma::uvec perm_index(const arma::imat& P, arma::uword b) {
  arma::uvec p(P.n_rows);
  for (arma::uword i = 0; i < P.n_rows; ++i) p(i) = static_cast<arma::uword>(P(i, b) - 1);   // 1-based -> 0-based
  return p;
}

// [[Rcpp::export]]
arma::vec permuted_contrast_F(const arma::mat& G, const arma::vec& a, const arma::mat& H, double cXc, int q, const arma::imat& P) {
  const arma::uword n = G.n_rows, B = P.n_cols;
  if (a.n_elem != n || H.n_rows != n || H.n_cols != n || P.n_rows != n) Rcpp::stop("dimension mismatch in permuted_contrast_F");
  const double trG = arma::trace(G);
  const double df_res = static_cast<double>(n) - q;
  arma::vec out(B);
  for (arma::uword b = 0; b < B; ++b) {
    arma::uvec p = perm_index(P, b);
    arma::vec ap = a(p);
    arma::mat Hp = H.submat(p, p);
    double ss = arma::dot(ap, G * ap) / cXc;
    double rss = trG - arma::accu(Hp % G);
    out(b) = ss / (rss / df_res);
  }
  return out;
}

// Freedman-Lane variant: the Gower matrix is permuted (G* = K1 + K2[, p] + K3[p, ] + K4[p, p]) while the
// design stays fixed.
// [[Rcpp::export]]
arma::vec permuted_contrast_F_fl(const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& K4,
                                 const arma::vec& a, const arma::mat& H, double cXc, int q, const arma::imat& P) {
  const arma::uword n = K1.n_rows, B = P.n_cols;
  if (a.n_elem != n || H.n_rows != n || P.n_rows != n) Rcpp::stop("dimension mismatch in permuted_contrast_F_fl");
  const double df_res = static_cast<double>(n) - q;
  arma::uvec all = arma::regspace<arma::uvec>(0, n - 1);
  arma::vec out(B);
  for (arma::uword b = 0; b < B; ++b) {
    arma::uvec p = perm_index(P, b);
    arma::mat Gs = K1 + K2.cols(p) + K3.rows(p) + K4.submat(p, p);
    double ss = arma::dot(a, Gs * a) / cXc;
    double rss = arma::trace(Gs) - arma::accu(H % Gs);
    out(b) = ss / (rss / df_res);
  }
  return out;
}

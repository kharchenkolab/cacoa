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

// ---------------------------------------------------------------------------------------------------------
// Full contrast statistics under a relabeling (step 2 of the engine convergence): F, shift, var, total.
// R reference: permutedStats() in R/pairwise_inference.R. Pieces (q x n matrix A = (X'X)^- X', hat matrix H,
// regression weights a = A'c, Z the dispersion design, znum / zden the dispersion-design rows at the contrast
// endpoints) are precomputed once in R; the relabeling permutes the design (block scheme) or the Gower matrix
// (Freedman-Lane).
// ---------------------------------------------------------------------------------------------------------

// statistics for one configuration: Gs (Gower matrix, possibly permuted), Xp / Ap / Zp / Hp / ap (design pieces,
// possibly row-permuted), dG = diag(Gs), trG = trace(Gs)
static inline arma::rowvec contrast_stats_one(const arma::mat& Gs, const arma::mat& Xp, const arma::mat& Ap, const arma::mat& Zp,
                                              const arma::mat& Hp, const arma::vec& ap, const arma::vec& cvec, double cXc, int q,
                                              const arma::vec& znum, const arma::vec& zden, bool bias_correct, bool need_var) {
  const arma::uword n = Gs.n_rows;
  const double df_res = static_cast<double>(n) - q;
  arma::rowvec out(4); out.fill(arma::datum::nan);
  const double ss = arma::dot(ap, Gs * ap) / cXc;
  const double rss = arma::trace(Gs) - arma::accu(Hp % Gs);
  out(0) = ss / (rss / df_res);
  if (!need_var) return out;
  arma::mat AG = Ap * Gs;                                   // q x n
  arma::mat M = AG * Ap.t();                                // q x q (bias-uncorrected effect matrix)
  arma::vec hG = arma::sum(Xp % AG.t(), 1);                 // diag(H_p G)
  arma::vec hGh = arma::sum((Xp * M) % Xp, 1);              // diag(H_p G H_p)
  arma::vec r = Gs.diag() - 2.0 * hG + hGh;                 // diag(R G R)
  arma::mat Rp = -Hp; Rp.diag() += 1.0;                     // residual projection under the relabeling
  arma::mat W = (Rp % Rp) * Zp;                             // E[r] = (R o R) Z gamma
  arma::vec gamma = arma::pinv(W.t() * W) * (W.t() * r);
  if (!gamma.is_finite()) return out;
  arma::vec s = Zp * gamma;
  if (bias_correct) M -= Ap * (Ap.each_row() % s.t()).t(); // M - A_p diag(s) A_p'
  const double shift = arma::as_scalar(cvec.t() * M * cvec);
  const double var = 2.0 * (arma::dot(znum, gamma) - arma::dot(zden, gamma));
  out(1) = shift; out(2) = var; out(3) = shift + var / 2.0;
  return out;
}

// [[Rcpp::export]]
arma::mat permuted_contrast_stats(const arma::mat& G, const arma::mat& X, const arma::mat& Z, const arma::mat& A, const arma::mat& H,
                                  const arma::vec& a, const arma::vec& cvec, double cXc, int q, const arma::vec& znum, const arma::vec& zden,
                                  const arma::imat& P, bool bias_correct = true, bool need_var = true) {
  const arma::uword n = G.n_rows, B = P.n_cols;
  if (X.n_rows != n || Z.n_rows != n || A.n_cols != n || H.n_rows != n || a.n_elem != n || P.n_rows != n) Rcpp::stop("dimension mismatch in permuted_contrast_stats");
  arma::mat out(B, 4);
  for (arma::uword b = 0; b < B; ++b) {
    arma::uvec p = perm_index(P, b);
    out.row(b) = contrast_stats_one(G, X.rows(p), A.cols(p), Z.rows(p), H.submat(p, p), a(p), cvec, cXc, q, znum, zden, bias_correct, need_var);
  }
  return out;
}

// [[Rcpp::export]]
arma::mat permuted_contrast_stats_fl(const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& K4,
                                     const arma::mat& X, const arma::mat& Z, const arma::mat& A, const arma::mat& H, const arma::vec& a,
                                     const arma::vec& cvec, double cXc, int q, const arma::vec& znum, const arma::vec& zden,
                                     const arma::imat& P, bool bias_correct = true, bool need_var = true) {
  const arma::uword n = K1.n_rows, B = P.n_cols;
  if (X.n_rows != n || Z.n_rows != n || A.n_cols != n || H.n_rows != n || a.n_elem != n || P.n_rows != n) Rcpp::stop("dimension mismatch in permuted_contrast_stats_fl");
  arma::mat out(B, 4);
  for (arma::uword b = 0; b < B; ++b) {
    arma::uvec p = perm_index(P, b);
    arma::mat Gs = K1 + K2.cols(p) + K3.rows(p) + K4.submat(p, p);
    out.row(b) = contrast_stats_one(Gs, X, A, Z, H, a, cvec, cXc, q, znum, zden, bias_correct, need_var);
  }
  return out;
}

// ---------------------------------------------------------------------------------------------------------
// Whole-factor (term) statistics under a relabeling (step 3): location F of the columns of X_full not in
// X_reduced, and the dispersion F of the leverage-corrected residual dispersions regressed on Z_full vs Z_reduced.
// R references: termTestGower() and dispersionTermTest() in R/pairwise_effects.R, looped in
// termPermutationStats() and screenOneMatrix().
// ---------------------------------------------------------------------------------------------------------

// rank-based least squares residual sum of squares of v on Z (Z may be rank deficient)
static inline double rss_ls(const arma::mat& Z, const arma::vec& v) {
  arma::vec beta = arma::pinv(Z) * v;
  arma::vec r = v - Z * beta;
  return arma::dot(r, r);
}

// statistics for one configuration: Gs (possibly permuted), Hf / Hr hat matrices (possibly permuted), Zf / Zr
// dispersion designs (possibly row-permuted); df / nu from the ranks (invariant under relabeling)
static inline arma::rowvec term_stats_one(const arma::mat& Gs, const arma::mat& Hf, const arma::mat& Hr, const arma::mat& Zf, const arma::mat& Zr,
                                          double df, double nu, int qZf, int qZr, bool need_disp) {
  arma::rowvec out(2); out.fill(arma::datum::nan);
  const double trG = arma::trace(Gs);
  const double ss = arma::accu((Hf - Hr) % Gs);
  const double rss = trG - arma::accu(Hf % Gs);
  out(0) = (ss / df) / (rss / nu);
  if (!need_disp) return out;
  arma::mat HG = Hf * Gs;
  arma::vec r = Gs.diag() - 2.0 * arma::sum(Hf % Gs, 1) + arma::sum(HG % Hf, 1);    // diag(R G R)
  arma::vec h = Hf.diag();
  arma::uvec ok = arma::find(h < 0.99);
  const double m = static_cast<double>(ok.n_elem);
  const double df_d = qZf - qZr, nu_d = m - qZf;
  if (df_d <= 0 || nu_d <= 0) return out;
  arma::vec v = arma::sqrt(arma::clamp(r(ok), 0.0, arma::datum::inf) / (1.0 - h(ok)));
  const double rss_f = rss_ls(Zf.rows(ok), v), rss_r = rss_ls(Zr.rows(ok), v);
  out(1) = ((rss_r - rss_f) / df_d) / (rss_f / nu_d);
  return out;
}

// [[Rcpp::export]]
arma::mat permuted_term_stats(const arma::mat& G, const arma::mat& Hf, const arma::mat& Hr, const arma::mat& Zf, const arma::mat& Zr,
                              double df, double nu, int qZf, int qZr, const arma::imat& P, bool need_disp = true) {
  const arma::uword n = G.n_rows, B = P.n_cols;
  if (Hf.n_rows != n || Hr.n_rows != n || Zf.n_rows != n || Zr.n_rows != n || P.n_rows != n) Rcpp::stop("dimension mismatch in permuted_term_stats");
  arma::mat out(B, 2);
  for (arma::uword b = 0; b < B; ++b) {
    arma::uvec p = perm_index(P, b);
    out.row(b) = term_stats_one(G, Hf.submat(p, p), Hr.submat(p, p), Zf.rows(p), Zr.rows(p), df, nu, qZf, qZr, need_disp);
  }
  return out;
}

// [[Rcpp::export]]
arma::mat permuted_term_stats_fl(const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& K4,
                                 const arma::mat& Hf, const arma::mat& Hr, const arma::mat& Zf, const arma::mat& Zr,
                                 double df, double nu, int qZf, int qZr, const arma::imat& P, bool need_disp = true) {
  const arma::uword n = K1.n_rows, B = P.n_cols;
  if (Hf.n_rows != n || Hr.n_rows != n || Zf.n_rows != n || Zr.n_rows != n || P.n_rows != n) Rcpp::stop("dimension mismatch in permuted_term_stats_fl");
  arma::mat out(B, 2);
  for (arma::uword b = 0; b < B; ++b) {
    arma::uvec p = perm_index(P, b);
    arma::mat Gs = K1 + K2.cols(p) + K3.rows(p) + K4.submat(p, p);
    out.row(b) = term_stats_one(Gs, Hf, Hr, Zf, Zr, df, nu, qZf, qZr, need_disp);
  }
  return out;
}

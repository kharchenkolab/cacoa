// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <mutex>
#include "parallel.h"
#include "lm_common.h"

// Batched cluster-free expression-shift kernel (step 4 of the engine convergence).
//
// For every cell (column of Y, the pair-layout distance matrix produced by estimateExpressionShiftsPairsLM), this
// does in C++ what clusterFreeExpressionShifts() used to do in R per cell: select the samples present, check
// estimability of the contrast on the reduced design, Gower-centre the distances, induce the shared permutations
// onto the present samples within the plan's cells, compute the observed shift F and shift estimate, the permuted
// F under the block (relabel the design) or Freedman-Lane (re-index the residual kernel) scheme, the per-cell
// permutation p-value and z-score, and the per-permutation extremes of the z-scores over cells for the
// max-statistic adjustment. R reference: the former oneCell() loop, kept in tests/testthat/helper-reference.R.

namespace {

struct Hat { arma::mat XtXi, A, H; int rank; };
inline Hat hat_of(const arma::mat& X) {
  Hat h; h.XtXi = arma::pinv(X.t() * X); h.A = h.XtXi * X.t(); h.H = X * h.A; h.rank = (int)arma::rank(X); return h;
}
inline bool estimable(const arma::mat& X, const arma::vec& c, double tol = 1e-8) {
  if (c.n_elem == 0 || arma::all(arma::abs(c) < tol)) return false;
  arma::mat XtX = X.t() * X; arma::vec Pc = arma::pinv(XtX) * (XtX * c);
  return arma::norm(Pc - c) <= tol * std::max(1.0, arma::norm(c));
}
inline arma::mat gower(const arma::mat& D2) {
  const arma::uword n = D2.n_rows; arma::vec rm = arma::mean(D2, 1); double gm = arma::mean(rm);
  arma::mat G = -0.5 * (D2 - arma::repmat(rm, 1, n) - arma::repmat(rm.t(), n, 1) + gm);
  return 0.5 * (G + G.t());
}
inline double contrast_F(const arma::mat& Gs, const arma::mat& Hp, const arma::vec& ap, double cXc, int q) {
  const double n = (double)Gs.n_rows;
  return (arma::dot(ap, Gs * ap) / cXc) / ((arma::trace(Gs) - arma::accu(Hp % Gs)) / (n - q));
}
// shift estimate with the intercept-only dispersion model (bias-corrected): c' (M - A diag(s) A') c
inline double shift_estimate(const arma::mat& G, const arma::mat& X, const Hat& h, const arma::vec& c, bool bias) {
  arma::mat AG = h.A * G; arma::mat M = AG * h.A.t();
  if (bias) {
    arma::vec hG = arma::sum(X % AG.t(), 1), hGh = arma::sum((X * M) % X, 1);
    arma::vec r = G.diag() - 2.0 * hG + hGh;
    arma::mat R = -h.H; R.diag() += 1.0;
    arma::vec w = arma::sum(R % R, 1);                       // (R o R) 1
    double gamma = arma::dot(w, r) / arma::dot(w, w);        // least squares on the intercept
    M -= h.A * (h.A.each_row() % arma::rowvec(X.n_rows, arma::fill::value(gamma))).t();
  }
  return arma::as_scalar(c.t() * M * c);
}

} // namespace

// [[Rcpp::export]]
Rcpp::List cluster_free_shift_batch(const arma::mat& Y, const arma::imat& pairs, int n_samples, const arma::mat& X, const arma::vec& cvec,
                                    const arma::ivec& level_code, int min_samp_per_level, const arma::ivec& stratum, const arma::uvec& inset,
                                    const arma::imat& P, bool freedman_lane, bool bias_correct, int n_cores) {
  const arma::uword m_cells = Y.n_cols, B = P.n_cols, n = (arma::uword)n_samples;
  if ((arma::uword)X.n_rows != n || P.n_rows != n || level_code.n_elem != n || stratum.n_elem != n || inset.n_elem != n) Rcpp::stop("dimension mismatch in cluster_free_shift_batch");
  arma::vec stat(m_cells), pval(m_cells), shift(m_cells), z(m_cells); stat.fill(arma::datum::nan); pval.fill(arma::datum::nan); shift.fill(arma::datum::nan); z.fill(arma::datum::nan);
  arma::ivec nsamp(m_cells, arma::fill::zeros);
  arma::vec mx(B); mx.fill(-arma::datum::inf); arma::vec mn(B); mn.fill(arma::datum::inf);
  std::mutex ext_lock;
  const bool use_levels = arma::any(level_code >= 0);
  const arma::uword n_pairs = pairs.n_rows;

  cacoa::parallelFor(0, (int)m_cells, [&](int k) {
    // 1. present samples: every pair of the sample observed
    arma::uvec bad(n, arma::fill::zeros);
    for (arma::uword r = 0; r < n_pairs; ++r) if (!std::isfinite(Y(r, k))) { bad[pairs(r, 0)] = 1; bad[pairs(r, 1)] = 1; }
    arma::uvec s = arma::find(bad == 0); const arma::uword m = s.n_elem;
    if (m < 3) return;
    if (use_levels) {
      int n_den = 0, n_num = 0;
      for (arma::uword j = 0; j < m; ++j) { if (level_code[s[j]] == 1) ++n_den; else if (level_code[s[j]] == 2) ++n_num; }
      if (n_den < min_samp_per_level || n_num < min_samp_per_level) return;
    }
    // 2. reduced design and contrast
    arma::mat Xs = X.rows(s);
    arma::uvec keep = arma::find(arma::sum(arma::abs(Xs), 0).t() > 1e-12);
    for (arma::uword j = 0; j < X.n_cols; ++j) if (!arma::any(keep == j) && std::abs(cvec[j]) > 1e-12) return;
    Xs = Xs.cols(keep); arma::vec c = cvec(keep);
    if (Xs.n_cols >= m || !estimable(Xs, c)) return;
    // 3. Gower matrix of the present samples ("cor" / "cosine" / "js" distances enter as they are)
    arma::mat D(m, m, arma::fill::zeros); arma::ivec pos(n); pos.fill(-1); for (arma::uword j = 0; j < m; ++j) pos[s[j]] = (int)j;
    for (arma::uword r = 0; r < n_pairs; ++r) { int i1 = pos[pairs(r, 0)], i2 = pos[pairs(r, 1)]; if (i1 >= 0 && i2 >= 0) { D(i1, i2) = Y(r, k); D(i2, i1) = Y(r, k); } }
    arma::mat G = gower(D);
    Hat h = hat_of(Xs);
    arma::vec a = h.A.t() * c; const double cXc = arma::as_scalar(c.t() * h.XtXi * c);
    const double Fo = contrast_F(G, h.H, a, cXc, h.rank);
    if (!std::isfinite(Fo)) return;
    // 4. induced permutations within the plan's cells restricted to the present samples
    std::vector<arma::uvec> blocks;
    {
      std::map<int, std::vector<arma::uword>> by_stratum;
      for (arma::uword j = 0; j < m; ++j) if (inset[s[j]]) by_stratum[stratum[s[j]]].push_back(j);
      for (auto& kv : by_stratum) if (kv.second.size() >= 2) blocks.push_back(arma::uvec(kv.second));
    }
    arma::vec Fp(B);
    if (freedman_lane) {
      arma::mat Q, Rq; arma::qr(Q, Rq, c); arma::mat Xr = Xs * Q.cols(1, Q.n_cols - 1);
      Hat hr = hat_of(Xr); arma::mat Rr = -hr.H; Rr.diag() += 1.0;
      arma::mat K1 = hr.H * G * hr.H, K2 = hr.H * G * Rr, K3 = Rr * G * hr.H, K4 = Rr * G * Rr;
      for (arma::uword b = 0; b < B; ++b) {
        arma::uvec q = induced_perm(arma::conv_to<arma::uvec>::from(P.col(b) - 1), s, blocks);
        arma::mat Gs = K1 + K2.cols(q) + K3.rows(q) + K4.submat(q, q);
        Fp[b] = contrast_F(Gs, h.H, a, cXc, h.rank);
      }
    } else {
      for (arma::uword b = 0; b < B; ++b) {
        arma::uvec q = induced_perm(arma::conv_to<arma::uvec>::from(P.col(b) - 1), s, blocks);
        Fp[b] = contrast_F(G, h.H.submat(q, q), a(q), cXc, h.rank);
      }
    }
    // 5. per-cell results and running extremes of the z-scores
    stat[k] = Fo; nsamp[k] = (int)m; shift[k] = shift_estimate(G, Xs, h, c, bias_correct);
    arma::uword ge = 0; for (arma::uword b = 0; b < B; ++b) if (Fp[b] >= Fo - 1e-12) ++ge;
    pval[k] = (ge + 1.0) / (B + 1.0);
    const double mu = arma::mean(Fp), sdv = arma::stddev(Fp);
    if (std::isfinite(sdv) && sdv > 0) {
      z[k] = (Fo - mu) / sdv;
      arma::vec zp = (Fp - mu) / sdv;
      std::lock_guard<std::mutex> lk(ext_lock);
      mx = arma::max(mx, zp); mn = arma::min(mn, zp);
    }
  }, n_cores, false);

  return Rcpp::List::create(Rcpp::_["stat"] = stat, Rcpp::_["p"] = pval, Rcpp::_["shift"] = shift, Rcpp::_["n"] = nsamp, Rcpp::_["z"] = z,
                            Rcpp::_["max"] = mx, Rcpp::_["min"] = mn);
}

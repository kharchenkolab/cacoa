// Weighted Gower-model statistics shared by the permutation kernels (perm_stats.cpp) and the cluster-free kernel.
#ifndef CACOA_GOWER_STATS_H
#define CACOA_GOWER_STATS_H
#include <RcppArmadillo.h>
namespace cacoa_gower {
struct WHat { arma::mat XtWXi, A, H; arma::vec h; int rank; };
inline WHat what_of(const arma::mat& X, const arma::vec& w) {
  WHat o; arma::mat Xw = X.each_col() % w;
  o.XtWXi = arma::pinv(X.t() * Xw); o.A = o.XtWXi * Xw.t(); o.H = X * o.A; o.h = arma::clamp(o.H.diag(), -arma::datum::inf, 1.0); o.rank = (int)arma::rank(X);
  return o;
}
inline double median_of(arma::vec v) { return v.n_elem ? arma::median(v) : arma::datum::nan; }
// robust = 1 huber (iterated by the caller), 2 winsor (one step)
inline arma::vec robust_weights(const arma::vec& score, const arma::vec& base, int robust, double k) {
  arma::uvec ref = arma::find(base >= 0.5 && score == score);        // finite scores of full-weight samples
  if (ref.n_elem < 3) return base;
  arma::vec rs = score(ref); const double med = median_of(rs); const double s = 1.4826 * median_of(arma::abs(rs - med));
  if (!std::isfinite(s) || s <= 1e-12) return base;
  arma::vec w = base;
  for (arma::uword i = 0; i < score.n_elem; ++i) {
    if (!(score[i] == score[i]) || score[i] <= med + k * s) continue;
    w[i] *= (robust == 1) ? (k * s) / (score[i] - med) : std::pow((med + k * s) / score[i], 2.0);
  }
  return w;
}
inline arma::vec disp_gamma_w(const arma::mat& R, const arma::mat& Z, const arma::vec& r, const arma::vec& w) {
  arma::mat Wm = (R % R) * Z; arma::vec sw = arma::sqrt(w);
  arma::mat Ws = Wm.each_col() % sw;
  return arma::pinv(Ws.t() * Ws) * (Ws.t() * (r % sw));
}
struct WFitState { WHat hi; arma::mat R; arma::vec r, score; };
inline WFitState wfit(const arma::mat& G, const arma::mat& X, const arma::vec& w, double h_cut) {
  WFitState st; st.hi = what_of(X, w); st.R = -st.hi.H; st.R.diag() += 1.0;
  st.r = arma::sum((st.R * G) % st.R, 1);
  st.score.set_size(X.n_rows); st.score.fill(arma::datum::nan);
  for (arma::uword i = 0; i < X.n_rows; ++i) if (st.hi.h[i] < h_cut) st.score[i] = std::sqrt(std::max(st.r[i], 0.0) / (1.0 - st.hi.h[i]));
  return st;
}
inline arma::vec iterate_weights(const arma::mat& G, const arma::mat& X, const arma::vec& base, int robust, double k, int maxit, double h_cut, WFitState& st) {
  arma::vec w = base; st = wfit(G, X, w, h_cut);
  if (robust == 0) return w;
  const int its = (robust == 1) ? maxit : 1;
  for (int it = 0; it < its; ++it) {
    arma::vec wn = robust_weights(st.score, base, robust, k);
    if (arma::abs(wn - w).max() < 1e-8) break;
    w = wn; st = wfit(G, X, w, h_cut);
  }
  return w;
}
inline arma::rowvec contrast_stats_w(const arma::mat& G, const arma::mat& X, const arma::mat& Z, const arma::vec& base, const arma::vec& c,
                                     const arma::vec& znum, const arma::vec& zden, bool bias, bool need_var, int robust, double k, int maxit) {
  arma::rowvec out(4); out.fill(arma::datum::nan);
  WFitState st; arma::vec w = iterate_weights(G, X, base, robust, k, maxit, 1.0 - 1e-8, st);
  arma::vec a = st.hi.A.t() * c; const double cXc = arma::as_scalar(c.t() * st.hi.XtWXi * c);
  const double ss = arma::dot(a, G * a) / cXc, rss = arma::dot(w, st.r), df = arma::accu(w) - st.hi.rank;
  out(0) = ss / (rss / df);
  if (!need_var) return out;
  arma::mat M = st.hi.A * G * st.hi.A.t();
  arma::vec gamma = disp_gamma_w(st.R, Z, st.r, w);
  if (!gamma.is_finite()) return out;
  arma::vec s = Z * gamma;
  if (bias) M -= st.hi.A * (st.hi.A.each_row() % s.t()).t();
  const double shift = arma::as_scalar(c.t() * M * c), var = 2.0 * (arma::dot(znum, gamma) - arma::dot(zden, gamma));
  out(1) = shift; out(2) = var; out(3) = shift + var / 2.0;
  return out;
}
inline double rss_wls(const arma::mat& Z, const arma::vec& v, const arma::vec& w) {
  arma::vec sw = arma::sqrt(w); arma::mat Zs = Z.each_col() % sw; arma::vec b = arma::pinv(Zs) * (v % sw);
  arma::vec res = v - Z * b; return arma::dot(w, res % res);
}
inline arma::rowvec term_stats_w(const arma::mat& G, const arma::mat& Xf, const arma::mat& Xr, const arma::mat& Zf, const arma::mat& Zr,
                                 const arma::vec& base, bool need_disp, int robust, double k, int maxit) {
  arma::rowvec out(2); out.fill(arma::datum::nan);
  WFitState st; arma::vec w = iterate_weights(G, Xf, base, robust, k, maxit, 0.99, st);
  WHat hr = what_of(Xr, w);
  arma::mat Hd = st.hi.H - hr.H;
  const double ss = arma::accu((Hd.each_col() % w) % G), rss = arma::dot(w, st.r);
  const double df = st.hi.rank - hr.rank, nu = arma::accu(w) - st.hi.rank;
  if (df > 0 && nu > 0) out(0) = (ss / df) / (rss / nu);
  if (!need_disp) return out;
  arma::uvec ok = arma::find(st.hi.h < 0.99);
  arma::vec v = st.score(ok), ww = w(ok);
  const double dfd = arma::rank(Zf) - arma::rank(Zr), nud = arma::accu(ww) - arma::rank(Zf);
  if (dfd <= 0 || nud <= 0) return out;
  const double rf = rss_wls(Zf.rows(ok), v, ww), rr = rss_wls(Zr.rows(ok), v, ww);
  out(1) = ((rr - rf) / dfd) / (rf / nud);
  return out;
}
} // namespace cacoa_gower
#endif

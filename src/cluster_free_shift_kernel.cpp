// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <mutex>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <memory>
#include "parallel.h"
#include "lm_common.h"
#include "cf_common.h"
#include "gower_stats.h"
#include "cf_accum.h"

// Batched cluster-free expression-shift kernel (step 4 of the engine convergence).
//
// For every cell (column of Y, the pair-layout distance matrix produced by estimateExpressionShiftsPairsLM), this
// does in C++ what clusterFreeExpressionShifts() used to do in R per cell: select the samples present, check
// estimability of the contrast on the reduced design, Gower-centre the distances, apply the shared permutations
// (block scheme: the global relabeling restricted to the present samples; Freedman-Lane: induced onto the present
// samples within the plan's cells), compute the observed shift F and shift estimate, the permuted
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


struct CellResult { double F = arma::datum::nan, p = arma::datum::nan, shift = arma::datum::nan, z = arma::datum::nan; int n = 0, n_perm = 0; bool ok = false; arma::vec zp, Fp; };

// Complete-sample case (every sample present, block scheme, no robust weighting): the hat matrix is the same for
// every such cell and the permuted F is a quadratic form in the cell's Gower matrix, F_b = (vec(G)' vec(a_b a_b') /
// cXc) / ((tr G - vec(G)' vec(H_b)) / (n - q)) with a_b = a[q_b] and H_b = H[q_b, q_b]: the two n^2 x B matrices are
// built once and the permutation loop of a cell becomes two matrix-vector products.
struct CompleteCase { bool ok = false; arma::uword n_cols = 0; arma::mat Aouter, Hperm; double cXc = 0; int rank = 0; };
// Two-group design (intercept + one 0/1 dummy, the contrast on the dummy), block scheme, no robust weighting: with the
// Gower-centred G (zero row sums) the permuted F of a cell with m present samples and n_A / n_B relabeled samples is
// F_b = S_b r_b (m - 2) / (tr G - S_b r_b), S_b = 1_B' G 1_B (the sum of G over the B block), r_b = 1/n_A + 1/n_B:
// no hat matrix per permutation, and all permutations at once as G L with L the m x B matrix of B-membership.
struct TwoGroup { bool ok = false; arma::vec dummy; int min_per_level = 2; arma::uword n_cols = 0; };
struct KernelInputs {
  const arma::imat& pairs; arma::uword n; const arma::mat& X; const arma::vec& cvec; const arma::ivec& level_code; int min_samp_per_level;
  const arma::ivec& stratum; const arma::uvec& inset; const arma::imat& P; bool freedman_lane; bool bias_correct; bool use_levels;
  int robust; double robust_k; const CompleteCase* cc = nullptr; const TwoGroup* tg = nullptr;
};
// Both codings of a two-group design are recognized: intercept + 0/1 dummy with the contrast (0, +-1), and two
// complementary 0/1 indicators (cell means) with the contrast +-(-1, 1). Either way the group of sample i is dummy[i].
static TwoGroup two_group(const KernelInputs& in) {
  TwoGroup tg; if (in.robust > 0 || in.freedman_lane || in.X.n_cols != 2 || in.cvec.n_elem != 2) return tg;
  const arma::vec c0 = in.X.col(0), c1 = in.X.col(1);
  auto binary = [](const arma::vec& v) { for (arma::uword i = 0; i < v.n_elem; ++i) if (std::abs(v[i]) > 1e-12 && std::abs(v[i] - 1.0) > 1e-12) return false; return true; };
  if (!binary(c1)) return tg;
  const bool intercept = arma::all(arma::abs(c0 - 1.0) <= 1e-12) && std::abs(in.cvec[0]) <= 1e-12 && std::abs(std::abs(in.cvec[1]) - 1.0) <= 1e-12;
  const bool one_hot = binary(c0) && arma::all(arma::abs(c0 + c1 - 1.0) <= 1e-12) && std::abs(in.cvec[0] + in.cvec[1]) <= 1e-12 && std::abs(std::abs(in.cvec[1]) - 1.0) <= 1e-12;
  if (!intercept && !one_hot) return tg;
  tg.ok = true; tg.dummy = c1; tg.min_per_level = in.use_levels ? in.min_samp_per_level : 1; tg.n_cols = 2; return tg;
}
static CompleteCase complete_case(const KernelInputs& in) {
  CompleteCase cc; const arma::uword n = in.n, B = in.P.n_cols;
  if (in.robust > 0 || in.freedman_lane || (double)n * n * B * 16.0 > 256.0 * 1048576.0) return cc;
  arma::uvec keep = arma::find(arma::sum(arma::abs(in.X), 0).t() > 1e-12);
  for (arma::uword j = 0; j < in.X.n_cols; ++j) if (!arma::any(keep == j) && std::abs(in.cvec[j]) > 1e-12) return cc;
  arma::mat Xs = in.X.cols(keep); arma::vec c = in.cvec(keep);
  if (Xs.n_cols >= n || !estimable(Xs, c)) return cc;
  Hat h = hat_of(Xs); arma::vec a = h.A.t() * c; cc.cXc = arma::as_scalar(c.t() * h.XtXi * c); cc.rank = h.rank; cc.n_cols = Xs.n_cols;
  std::vector<arma::uvec> blocks; { std::map<int, std::vector<arma::uword>> by; for (arma::uword j = 0; j < n; ++j) if (in.inset[j]) by[in.stratum[j]].push_back(j); for (auto& kv : by) if (kv.second.size() >= 2) blocks.push_back(arma::uvec(kv.second)); }
  const arma::uvec units = arma::regspace<arma::uvec>(0, n - 1);
  cc.Aouter.set_size(n * n, B); cc.Hperm.set_size(n * n, B);
  for (arma::uword b = 0; b < B; ++b) {
    arma::uvec q = induced_perm(arma::conv_to<arma::uvec>::from(in.P.col(b) - 1), units, blocks);
    arma::vec aq = a(q); cc.Aouter.col(b) = arma::vectorise(aq * aq.t()); cc.Hperm.col(b) = arma::vectorise(h.H.submat(q, q));
  }
  cc.ok = true; return cc;
}

// one cell: y holds the pair distances (NaN where a sample has too few cells)
static CellResult test_cell(const arma::vec& y, const KernelInputs& in) {
  CellResult out; const arma::uword n = in.n, n_pairs = in.pairs.n_rows, B = in.P.n_cols;
  // samples with too few cells have every pair distance missing; drop them (and only them) by removing, one at a time,
  // the sample with the most missing pairs until no missing pair is left among the remaining samples
  arma::uvec bad(n, arma::fill::zeros);
  for (;;) {
    arma::uvec miss(n, arma::fill::zeros);
    for (arma::uword r = 0; r < n_pairs; ++r) {
      int a = in.pairs(r, 0), b = in.pairs(r, 1);
      if (bad[a] || bad[b] || std::isfinite(y[r])) continue;
      ++miss[a]; ++miss[b];
    }
    if (miss.max() == 0) break;
    bad[miss.index_max()] = 1;
  }
  arma::uvec s = arma::find(bad == 0); const arma::uword m = s.n_elem;
  if (m < 3) return out;
  if (in.use_levels) {
    int n_den = 0, n_num = 0;
    for (arma::uword j = 0; j < m; ++j) { if (in.level_code[s[j]] == 1) ++n_den; else if (in.level_code[s[j]] == 2) ++n_num; }
    if (n_den < in.min_samp_per_level || n_num < in.min_samp_per_level) return out;
  }
  arma::mat Xs = in.X.rows(s);
  arma::uvec keep = arma::find(arma::sum(arma::abs(Xs), 0).t() > 1e-12);
  for (arma::uword j = 0; j < in.X.n_cols; ++j) if (!arma::any(keep == j) && std::abs(in.cvec[j]) > 1e-12) return out;
  Xs = Xs.cols(keep); arma::vec c = in.cvec(keep);
  if (Xs.n_cols >= m || !estimable(Xs, c)) return out;
  arma::mat D(m, m, arma::fill::zeros); arma::ivec pos(n); pos.fill(-1); for (arma::uword j = 0; j < m; ++j) pos[s[j]] = (int)j;
  for (arma::uword r = 0; r < n_pairs; ++r) { int i1 = pos[in.pairs(r, 0)], i2 = pos[in.pairs(r, 1)]; if (i1 >= 0 && i2 >= 0) { D(i1, i2) = y[r]; D(i2, i1) = y[r]; } }
  arma::mat G = gower(D);
  Hat h = hat_of(Xs);
  arma::vec a = h.A.t() * c; const double cXc = arma::as_scalar(c.t() * h.XtXi * c);
  double Fo = contrast_F(G, h.H, a, cXc, h.rank);
  arma::vec ones(m, arma::fill::ones), w_obs = ones; const arma::mat Z1(m, 1, arma::fill::ones); const arma::vec z1(1, arma::fill::ones);
  double shift_w = arma::datum::nan;
  if (in.robust > 0) {   // robust path: weights re-estimated per relabeling, FL parts from the observed robust weights
    arma::rowvec o = cacoa_gower::contrast_stats_w(G, Xs, Z1, ones, c, z1, z1, in.bias_correct, true, in.robust, in.robust_k, 5);
    cacoa_gower::WFitState st; w_obs = cacoa_gower::iterate_weights(G, Xs, ones, in.robust, in.robust_k, 5, 1.0 - 1e-8, st);
    Fo = o(0); shift_w = o(1);
  }
  if (!std::isfinite(Fo)) return out;
  std::vector<arma::uvec> blocks;
  {
    std::map<int, std::vector<arma::uword>> by_stratum;
    for (arma::uword j = 0; j < m; ++j) if (in.inset[s[j]]) by_stratum[in.stratum[s[j]]].push_back(j);
    for (auto& kv : by_stratum) if (kv.second.size() >= 2) blocks.push_back(arma::uvec(kv.second));
  }
  // Label permutations (the block scheme) are applied as global relabelings: present sample j takes the design row of
  // the sample P(s[j], b) maps it to, which may be absent from this neighbourhood. Inducing the permutation within the
  // present samples instead would keep every cell's level counts fixed and make the relabelings incoherent across
  // cells (each cell's permuted labeling is then conditional on its own sample set): the per-cell p-value stays exact
  // but any statistic pooled over cells (group means, the max statistic) is compared with a null that is too narrow.
  // Under a relabeling that leaves fewer than min_samp_per_level present samples at a compared level, or makes the
  // contrast inestimable, the cell would not have been tested; such permutations are skipped (Fp NaN) and the
  // p-value, z-score and pooled accumulators use the admissible ones. Residual permutations (Freedman-Lane) act on
  // the present samples' residuals and keep the induced permutation.
  const bool complete = (m == n);
  auto permuted_design = [&](arma::uword b, arma::mat& Xp, arma::vec& cp) -> bool {
    arma::uvec ps(m);
    for (arma::uword j = 0; j < m; ++j) ps[j] = (arma::uword)(in.P(s[j], b) - 1);
    if (in.use_levels) {
      int n_den = 0, n_num = 0;
      for (arma::uword j = 0; j < m; ++j) { if (in.level_code[ps[j]] == 1) ++n_den; else if (in.level_code[ps[j]] == 2) ++n_num; }
      if (n_den < in.min_samp_per_level || n_num < in.min_samp_per_level) return false;
    }
    arma::mat Xf = in.X.rows(ps);
    arma::uvec kp = arma::find(arma::sum(arma::abs(Xf), 0).t() > 1e-12);
    for (arma::uword j = 0; j < in.X.n_cols; ++j) if (!arma::any(kp == j) && std::abs(in.cvec[j]) > 1e-12) return false;
    Xp = Xf.cols(kp); cp = in.cvec(kp);
    return Xp.n_cols < m && estimable(Xp, cp);
  };
  arma::vec Fp(B); Fp.fill(arma::datum::nan);
  if (in.robust > 0) {
    if (in.freedman_lane) {
      arma::mat Q, Rq; arma::qr(Q, Rq, c); arma::mat Xr = Xs * Q.cols(1, Q.n_cols - 1);
      cacoa_gower::WHat hr = cacoa_gower::what_of(Xr, w_obs); arma::mat Rr = -hr.H; Rr.diag() += 1.0;
      arma::mat K1 = hr.H * G * hr.H.t(), K2 = hr.H * G * Rr.t(), K3 = Rr * G * hr.H.t(), K4 = Rr * G * Rr.t();
      for (arma::uword b = 0; b < B; ++b) {
        arma::uvec q = induced_perm(arma::conv_to<arma::uvec>::from(in.P.col(b) - 1), s, blocks);
        arma::mat Gs = K1 + K2.cols(q) + K3.rows(q) + K4.submat(q, q);
        Fp[b] = cacoa_gower::contrast_stats_w(Gs, Xs, Z1, ones, c, z1, z1, in.bias_correct, false, in.robust, in.robust_k, 5)(0);
      }
    } else if (complete) {
      for (arma::uword b = 0; b < B; ++b) {
        arma::uvec q = induced_perm(arma::conv_to<arma::uvec>::from(in.P.col(b) - 1), s, blocks);
        Fp[b] = cacoa_gower::contrast_stats_w(G, Xs.rows(q), Z1, ones, c, z1, z1, in.bias_correct, false, in.robust, in.robust_k, 5)(0);
      }
    } else {
      arma::mat Xp; arma::vec cp;
      for (arma::uword b = 0; b < B; ++b) {
        if (!permuted_design(b, Xp, cp)) continue;
        Fp[b] = cacoa_gower::contrast_stats_w(G, Xp, Z1, ones, cp, z1, z1, in.bias_correct, false, in.robust, in.robust_k, 5)(0);
      }
    }
  } else if (in.freedman_lane) {
    arma::mat Q, Rq; arma::qr(Q, Rq, c); arma::mat Xr = Xs * Q.cols(1, Q.n_cols - 1);
    Hat hr = hat_of(Xr); arma::mat Rr = -hr.H; Rr.diag() += 1.0;
    arma::mat K1 = hr.H * G * hr.H, K2 = hr.H * G * Rr, K3 = Rr * G * hr.H, K4 = Rr * G * Rr;
    for (arma::uword b = 0; b < B; ++b) {
      arma::uvec q = induced_perm(arma::conv_to<arma::uvec>::from(in.P.col(b) - 1), s, blocks);
      arma::mat Gs = K1 + K2.cols(q) + K3.rows(q) + K4.submat(q, q);
      Fp[b] = contrast_F(Gs, h.H, a, cXc, h.rank);
    }
  } else if (in.tg && in.tg->ok && Xs.n_cols == 2) {   // two-group design: all permutations as one matrix product, no hat matrix per relabeling
    arma::mat L(m, B);
    for (arma::uword b = 0; b < B; ++b) for (arma::uword j = 0; j < m; ++j) L(j, b) = in.tg->dummy[(arma::uword)(in.P(s[j], b) - 1)];
    const arma::rowvec nB = arma::sum(L, 0), nA = (double)m - nB;
    const arma::rowvec S = arma::sum((G * L) % L, 0), r = 1.0 / nA + 1.0 / nB;
    const double trG = arma::trace(G);
    for (arma::uword b = 0; b < B; ++b) {
      if (nA[b] < in.tg->min_per_level || nB[b] < in.tg->min_per_level) continue;
      const double den = (trG - S[b] * r[b]) / ((double)m - 2.0);
      Fp[b] = (S[b] * r[b]) / den;
    }
  } else if (complete && in.cc && in.cc->ok && in.cc->n_cols == Xs.n_cols) {   // every sample present: vectorized over permutations
    const arma::vec g = arma::vectorise(G);
    const arma::vec num = (in.cc->Aouter.t() * g) / in.cc->cXc;
    const arma::vec den = (arma::trace(G) - in.cc->Hperm.t() * g) / ((double)m - in.cc->rank);
    Fp = num / den;
  } else if (complete) {   // every sample present: the global relabeling is a row permutation of the observed design
    for (arma::uword b = 0; b < B; ++b) {
      arma::uvec q = induced_perm(arma::conv_to<arma::uvec>::from(in.P.col(b) - 1), s, blocks);
      Fp[b] = contrast_F(G, h.H.submat(q, q), a(q), cXc, h.rank);
    }
  } else {
    arma::mat Xp; arma::vec cp;
    for (arma::uword b = 0; b < B; ++b) {
      if (!permuted_design(b, Xp, cp)) continue;
      Hat hp = hat_of(Xp); arma::vec ap = hp.A.t() * cp; const double cXcp = arma::as_scalar(cp.t() * hp.XtXi * cp);
      Fp[b] = contrast_F(G, hp.H, ap, cXcp, hp.rank);
    }
  }
  const arma::uvec adm = arma::find_finite(Fp); const arma::uword B_adm = adm.n_elem;
  if (B_adm == 0) return out;
  out.ok = true; out.F = Fo; out.n = (int)m; out.n_perm = (int)B_adm; out.shift = (in.robust > 0) ? shift_w : shift_estimate(G, Xs, h, c, in.bias_correct);
  arma::uword ge = 0; for (arma::uword i = 0; i < B_adm; ++i) if (Fp[adm[i]] >= Fo - 1e-12) ++ge;
  out.p = (ge + 1.0) / (B_adm + 1.0);
  const arma::vec Fa = Fp(adm); const double mu = arma::mean(Fa), sdv = arma::stddev(Fa);
  if (std::isfinite(sdv) && sdv > 0) { out.z = (Fo - mu) / sdv; out.zp = (Fp - mu) / sdv; }   // NaN under the skipped permutations
  out.Fp = std::move(Fp);
  return out;
}

// The shared accumulator (cf_accum.h) takes the observed and permuted F of every cell; the kernel keeps the per-cell
// shift estimate and sample count next to it.
struct ShiftAccum {
  cacoa_cf::Accum acc; arma::vec shift; arma::ivec nsamp;
  ShiftAccum(arma::uword m, arma::uword B, const arma::vec& z_grid, const arma::ivec& group_in, int n_groups, bool keep_perm, bool average, const std::vector<std::vector<int>>* nn, bool standardize)
    : acc(m, B, z_grid, false, 1, group_in, n_groups, keep_perm, average, nn, standardize), shift(m), nsamp(m, arma::fill::zeros) { shift.fill(arma::datum::nan); }
  void add(arma::uword k, const CellResult& r) { if (!r.ok) return; shift[k] = r.shift; nsamp[k] = r.n; acc.add(k, r.F, r.Fp); }
  Rcpp::List result(int n_cores) {
    acc.finalize(n_cores); Rcpp::List out = acc.result(); out["shift"] = shift; out["n"] = nsamp; return out;
  }
};

} // namespace

// Batched kernel on a precomputed pair-layout distance matrix Y (pairs x cells); used as the reference for the
// streaming kernel and by tests.
// [[Rcpp::export]]
Rcpp::List cluster_free_shift_batch(const arma::mat& Y, const arma::imat& pairs, int n_samples, const arma::mat& X, const arma::vec& cvec,
                                    const arma::ivec& level_code, int min_samp_per_level, const arma::ivec& stratum, const arma::uvec& inset,
                                    const arma::imat& P, bool freedman_lane, bool bias_correct, int n_cores, int robust = 0, double robust_k = 1.345,
                                    Rcpp::NumericVector z_grid = Rcpp::NumericVector::create(), Rcpp::IntegerVector group = Rcpp::IntegerVector::create(), bool keep_perm = false,
                                    bool average = false, Rcpp::Nullable<Rcpp::List> avg_nn = R_NilValue, bool vectorize = true, bool standardize = true) {
  const arma::uword m_cells = Y.n_cols, n = (arma::uword)n_samples;
  if ((arma::uword)X.n_rows != n || P.n_rows != n || level_code.n_elem != n || stratum.n_elem != n || inset.n_elem != n) Rcpp::stop("dimension mismatch in cluster_free_shift_batch");
  KernelInputs in{pairs, n, X, cvec, level_code, min_samp_per_level, stratum, inset, P, freedman_lane, bias_correct, arma::any(level_code >= 0), robust, robust_k};
  CompleteCase cc; if (vectorize) cc = complete_case(in); in.cc = &cc; TwoGroup tg; if (vectorize) tg = two_group(in); in.tg = &tg;
  arma::ivec grp(group.begin(), group.size()); const int n_groups = grp.n_elem ? (int)grp.max() + 1 : 0;
  std::vector<std::vector<int>> nn_avg; if (average) { if (avg_nn.isNull()) Rcpp::stop("average needs avg_nn"); nn_avg = cacoa_cf::neighbour_lists(Rcpp::List(avg_nn), false, (int)m_cells, true); }
  ShiftAccum acc(m_cells, P.n_cols, arma::vec(z_grid.begin(), z_grid.size()), grp, n_groups, keep_perm, average, average ? &nn_avg : nullptr, standardize);
  cacoa::parallelFor(0, (int)m_cells, [&](int k) { acc.add(k, test_cell(Y.col(k), in)); }, n_cores, false);
  return acc.result(n_cores);
}

// Normalized neighbourhood profiles (profile = "normalized" in clusterFreeExpressionShifts()).
// Per unit and sample a: p_ga = Y_ga / N_a (summed counts over the summed full library sizes of the sample's cells in the
// unit), x_ga = log10(scale p_ga + 1). Noise model M2 (label-free; pbar_g pooled over the unit's present samples):
//   v_ga = Var_Poisson(log10(scale Y / N_a + 1)), Y ~ Pois(N_a pbar_g)   (exact, from a lookup table; delta method above it)
//        + t'(pbar_g)^2 phi(mu_g) pbar_g^2 / n_eff,a                        (within-unit heterogeneity; n_eff = N^2 / sum m_i^2)
// with mu = log10(pbar) and t'(p) = scale / ((scale p + 1) ln 10). Gene weights w_g = 1 / (vbar_g + tau2(mu_g)) (tau2: the
// between-sample floor; vbar_g the mean over present samples), sample noise c_a = sum_g w_g v_ga. Distances on
// z_ga = sqrt(w_g) (x_ga - xbar_g): "cor" 1 - <z_a, z_b> / (|z_a| |z_b|), with the known-noise correction |z_a|^2 - c_a
// (floored at disatt_floor |z_a|^2); "euclidean" |z_a - z_b|^2, with the correction - c_a - c_b. Genes absent from the
// unit are dropped when a noise model is used.
namespace cacoa_prof {
struct ProfOpts {
  double scale = 1e3;
  int dist = 1;                // 1 cor (gene-centred), 2 euclidean
  bool noise = false;          // M2 noise model (gene weights 1 / (vbar + tau2))
  bool tau = false;            // tau2(mu) floor in the gene weights
  bool correct = false;        // known-noise correction
  double disatt_floor = 0.05;
  double unit_detection = 0;   // > 0: per unit, only the genes detected in at least this fraction of the unit's cells (label-free)
  std::vector<double> mu_grid, phi, tau2;   // curves on mu = log10(pbar)
  std::vector<double> lam_grid, s_grid;     // Poisson lookup axes: log10(lambda), log10(scale / N)
  arma::mat V1;                             // Var(log10(s Y + 1)), Y ~ Pois(lambda); rows lambda, cols s
};
inline double interp1(const std::vector<double>& x, const std::vector<double>& y, double v) {
  if (x.empty()) return 0.0;
  if (v <= x.front()) return y.front();
  if (v >= x.back()) return y.back();
  size_t j = std::upper_bound(x.begin(), x.end(), v) - x.begin(); double t = (v - x[j - 1]) / (x[j] - x[j - 1]);
  return y[j - 1] + t * (y[j] - y[j - 1]);
}
inline double tprime(double s, double p) { return s / ((s * p + 1.0) * std::log(10.0)); }
// Poisson variance of the transform: bilinear in (log10 lambda, log10 s) inside the table, delta method above it,
// proportional to lambda below it
inline double poisson_var(const ProfOpts& o, double lambda, double N, double p) {
  if (lambda <= 0) return 0.0;
  const double ll = std::log10(lambda), ls = std::log10(o.scale / N);
  if (o.lam_grid.empty() || ll > o.lam_grid.back()) { const double tp = tprime(o.scale, p); return tp * tp * p / N; }
  const auto& X = o.lam_grid; const auto& Y = o.s_grid;
  double lx = std::max(X.front(), ll), ly = std::min(std::max(Y.front(), ls), Y.back());
  size_t i = std::min<size_t>(std::upper_bound(X.begin(), X.end(), lx) - X.begin(), X.size() - 1); if (i == 0) i = 1;
  size_t j = std::min<size_t>(std::upper_bound(Y.begin(), Y.end(), ly) - Y.begin(), Y.size() - 1); if (j == 0) j = 1;
  const double tx = (lx - X[i - 1]) / (X[i] - X[i - 1]), ty = (ly - Y[j - 1]) / (Y[j] - Y[j - 1]);
  double v = (1 - tx) * (1 - ty) * o.V1(i - 1, j - 1) + tx * (1 - ty) * o.V1(i, j - 1) + (1 - tx) * ty * o.V1(i - 1, j) + tx * ty * o.V1(i, j);
  if (ll < X.front()) v *= lambda / std::pow(10.0, X.front());
  return std::max(v, 0.0);
}
// the same with the per-sample part (the s axis) precomputed once per unit and the lambda axis located by direct indexing
// on the uniform grid (corrected to std::upper_bound's answer): bit-identical to poisson_var, without two log10 and two
// binary searches per gene x sample
struct PoissonSample { double ty = 0; size_t j = 1; };
inline PoissonSample poisson_sample(const ProfOpts& o, double N) {
  PoissonSample r; if (o.lam_grid.empty()) return r; const auto& Y = o.s_grid;
  const double ly = std::min(std::max(Y.front(), std::log10(o.scale / N)), Y.back());
  r.j = std::min<size_t>(std::upper_bound(Y.begin(), Y.end(), ly) - Y.begin(), Y.size() - 1); if (r.j == 0) r.j = 1;
  r.ty = (ly - Y[r.j - 1]) / (Y[r.j] - Y[r.j - 1]); return r;
}
inline size_t upper_index(const std::vector<double>& X, double v, double x0, double idx) {   // = upper_bound(X, v) - X.begin() on a (near-)uniform grid
  const size_t n = X.size(); const double f = (v - x0) * idx; size_t i = f <= 0 ? 0 : (f >= (double)n ? n : (size_t)f + 1);
  while (i < n && X[i] <= v) ++i;
  while (i > 0 && X[i - 1] > v) --i;
  return i;
}
struct PoissonGrid { double x0 = 0, idx = 0, low = 1; };
inline PoissonGrid poisson_grid(const ProfOpts& o) {
  PoissonGrid g; if (o.lam_grid.empty()) return g;
  g.x0 = o.lam_grid.front(); g.idx = o.lam_grid.size() > 1 ? (double)(o.lam_grid.size() - 1) / (o.lam_grid.back() - o.lam_grid.front()) : 0.0;
  g.low = std::pow(10.0, o.lam_grid.front()); return g;
}
inline double poisson_var_fast(const ProfOpts& o, double lambda, double N, double p, const PoissonSample& sm, const PoissonGrid& gr) {
  if (lambda <= 0) return 0.0;
  const double ll = std::log10(lambda);
  if (o.lam_grid.empty() || ll > o.lam_grid.back()) { const double tp = tprime(o.scale, p); return tp * tp * p / N; }
  const auto& X = o.lam_grid; const double lx = std::max(X.front(), ll);
  size_t i = std::min<size_t>(upper_index(X, lx, gr.x0, gr.idx), X.size() - 1); if (i == 0) i = 1;
  const size_t j = sm.j; const double tx = (lx - X[i - 1]) / (X[i] - X[i - 1]), ty = sm.ty;
  double v = (1 - tx) * (1 - ty) * o.V1(i - 1, j - 1) + tx * (1 - ty) * o.V1(i, j - 1) + (1 - tx) * ty * o.V1(i - 1, j) + tx * ty * o.V1(i, j);
  if (ll < X.front()) v *= lambda / gr.low;
  return std::max(v, 0.0);
}
struct Unit {
  arma::uvec genes, present;   // rows of cm used; samples with >= min cells
  arma::mat x, v;              // genes x n_samples (columns of absent samples unused)
  arma::vec w, mu, c, t2;      // gene weights, log10 pbar, sample noise c_a (n_samples), tau2(mu_g) (zeros without the floor)
  std::vector<unsigned> nobs;
};
inline Unit unit_profile(const Eigen::SparseMatrix<double>& cm, const std::vector<int>& sample_of_cell, const std::vector<int>& ids,
                         const std::vector<double>& lib, int n_samples, int min_n_obs, const ProfOpts& o) {
  Unit u; const int G = cm.rows();
  u.nobs.assign(n_samples, 0); std::vector<double> N(n_samples, 0.0), S2(n_samples, 0.0);
  arma::mat Y(G, n_samples, arma::fill::zeros);
  std::vector<double> det; if (o.unit_detection > 0) det.assign(G, 0.0);
  for (int id : ids) {
    const int a = sample_of_cell[id]; u.nobs[a]++; N[a] += lib[id]; S2[a] += lib[id] * lib[id];
    for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(cm, id); it; ++it) { Y(it.row(), a) += it.value(); if (!det.empty() && it.value() > 0) det[it.row()] += 1.0; }
  }
  std::vector<arma::uword> pres; for (int a = 0; a < n_samples; ++a) if ((int)u.nobs[a] >= min_n_obs) pres.push_back(a);
  u.present = arma::uvec(pres);
  double Ntot = 0; for (arma::uword a : pres) Ntot += N[a];
  const double min_det = o.unit_detection * (double)ids.size();
  std::vector<arma::uword> keep; std::vector<double> pbar;
  for (int g = 0; g < G; ++g) {
    if (!det.empty() && det[g] < min_det) continue;   // per-unit gene set: detected in >= unit_detection of the unit's cells
    double Ysum = 0; for (arma::uword a : pres) Ysum += Y(g, a);
    const double pb = Ntot > 0 ? Ysum / Ntot : 0;
    if (o.noise && !(pb > 0)) continue;   // genes absent from the unit carry no information
    keep.push_back(g); pbar.push_back(pb);
  }
  u.genes = arma::uvec(keep); const arma::uword g_n = keep.size();
  u.x.zeros(g_n, n_samples); u.v.zeros(g_n, n_samples); u.w.ones(g_n); u.mu.zeros(g_n); u.c.zeros(n_samples); u.t2.zeros(g_n);
  for (arma::uword k = 0; k < g_n; ++k) for (arma::uword a : pres) u.x(k, a) = std::log10(o.scale * (N[a] > 0 ? Y(keep[k], a) / N[a] : 0) + 1.0);
  if (!o.noise) return u;
  std::vector<double> neff(n_samples, 0.0); for (arma::uword a : pres) neff[a] = (S2[a] > 0) ? N[a] * N[a] / S2[a] : 0;
  std::vector<PoissonSample> ps(n_samples); for (arma::uword a : pres) ps[a] = poisson_sample(o, N[a]);   // per-sample half of the lookup, once per unit
  const PoissonGrid gr = poisson_grid(o);
  for (arma::uword k = 0; k < g_n; ++k) {
    const double p = pbar[k], mu = std::log10(p), tp = tprime(o.scale, p), ph = std::max(interp1(o.mu_grid, o.phi, mu), 0.0); u.mu[k] = mu;
    double vbar = 0;
    for (arma::uword a : pres) { const double v = poisson_var_fast(o, N[a] * p, N[a], p, ps[a], gr) + tp * tp * ph * p * p / std::max(neff[a], 1.0); u.v(k, a) = v; vbar += v; }
    vbar /= std::max<size_t>(1, pres.size());
    const double t2 = o.tau ? std::max(interp1(o.mu_grid, o.tau2, mu), 0.0) : 0.0; u.t2[k] = t2;
    const double den = vbar + t2;
    u.w[k] = (den > 0) ? 1.0 / den : 0.0;
  }
  for (arma::uword a : pres) u.c[a] = arma::dot(u.w, u.v.col(a));
  return u;
}
// pair distances of a unit (NaN for pairs with an absent sample)
inline arma::vec unit_distances(const Unit& u, const arma::imat& pairs, const ProfOpts& o) {
  const arma::uword n_pairs = pairs.n_rows, n = u.x.n_cols; arma::vec y(n_pairs); y.fill(arma::datum::nan);
  if (u.genes.n_elem < 2 || u.present.n_elem < 2) return y;
  arma::mat Z = u.x.cols(u.present);
  if (o.dist == 1) { arma::vec xm = arma::mean(Z, 1); Z.each_col() -= xm; }   // gene centring across the present samples
  Z.each_col() %= arma::sqrt(u.w);
  arma::ivec pos(n); pos.fill(-1); for (arma::uword j = 0; j < u.present.n_elem; ++j) pos[u.present[j]] = (int)j;
  arma::vec nrm2 = arma::sum(Z % Z, 0).t();
  for (arma::uword r = 0; r < n_pairs; ++r) {
    const int a = pairs(r, 0), b = pairs(r, 1), ia = pos[a], ib = pos[b]; if (ia < 0 || ib < 0) continue;
    if (o.dist == 2) {
      double d2 = arma::accu(arma::square(Z.col(ia) - Z.col(ib)));
      if (o.correct) d2 -= u.c[a] + u.c[b];
      y[r] = d2;
    } else {
      double na = nrm2[ia], nb = nrm2[ib];
      if (o.correct) { na = std::max(na - u.c[a], o.disatt_floor * na); nb = std::max(nb - u.c[b], o.disatt_floor * nb); }
      if (!(na > 0) || !(nb > 0)) continue;
      const double rr = std::min(1.0, std::max(-1.0, arma::dot(Z.col(ia), Z.col(ib)) / std::sqrt(na * nb)));
      y[r] = 1.0 - rr;
    }
  }
  return y;
}
inline ProfOpts parse_opts(const Rcpp::List& opt) {
  ProfOpts o;
  o.scale = Rcpp::as<double>(opt["scale"]); o.dist = Rcpp::as<int>(opt["dist"]); o.noise = Rcpp::as<bool>(opt["noise"]);
  o.tau = Rcpp::as<bool>(opt["tau"]); o.correct = Rcpp::as<bool>(opt["correct"]);
  if (opt.containsElementNamed("disatt_floor")) o.disatt_floor = Rcpp::as<double>(opt["disatt_floor"]);
  if (opt.containsElementNamed("unit_detection")) o.unit_detection = Rcpp::as<double>(opt["unit_detection"]);
  if (!(o.scale > 0)) Rcpp::stop("scale must be positive");
  if (!(o.unit_detection >= 0 && o.unit_detection < 1)) Rcpp::stop("unit_detection must be in [0, 1)");
  if (o.dist != 1 && o.dist != 2) Rcpp::stop("normalized profiles take dist cor or euclidean");
  if (o.noise) {
    o.mu_grid = Rcpp::as<std::vector<double>>(opt["mu_grid"]); o.phi = Rcpp::as<std::vector<double>>(opt["phi"]); o.tau2 = Rcpp::as<std::vector<double>>(opt["tau2"]);
    o.lam_grid = Rcpp::as<std::vector<double>>(opt["lam_grid"]); o.s_grid = Rcpp::as<std::vector<double>>(opt["s_grid"]); o.V1 = Rcpp::as<arma::mat>(opt["V1"]);
    if (o.phi.size() != o.mu_grid.size() || o.tau2.size() != o.mu_grid.size()) Rcpp::stop("noise curves must match mu_grid");
    if (o.V1.n_rows != o.lam_grid.size() || o.V1.n_cols != o.s_grid.size()) Rcpp::stop("Poisson table must match its axes");
  } else if (o.tau || o.correct) Rcpp::stop("tau and the noise correction need the noise model");
  return o;
}
} // namespace cacoa_prof

// Per-gene combination statistic of a unit (statistic = "gsd" in clusterFreeExpressionShifts()): replaces the distance and
// pseudo-F with a per-gene test combined over genes. For every relabeling b (0 = observed) and gene g of the unit, with the
// unit's normalized profiles x, their noise v (+ tau2) and the contrast estimate beta_gb = alpha_b' x_g:
//   model variance mv = sum_a alpha_ab^2 (v_ga + tau2_g); moderated (mode 2, default): (d0 mv + cXc_b RSS_gb) / (d0 + df_b), the
//   model variance shrunk towards the residual variance of the relabeled fit (two groups: the pooled within-group variance);
//   z^2 = beta^2 / variance; T_gamma = sum_g |z|^gamma (gamma 1, 2, 4, 8) and max_g z^2; the unit statistic is the max over
//   gamma of log T_gamma standardized over the admissible relabelings (aSPU-type; gamma > 0: that T only). It then goes
//   through the same accumulator as the pseudo-F. Two-group designs use subset tables (test_genes); any other single-df
//   contrast (continuous variables, covariates, block relabeling or Freedman-Lane) the general kernel (test_genes_lm).
//   winsor_k > 0: label-free winsorization of each gene's values before the test (robust option).
namespace cacoa_gsd {
struct GsdOpts {
  int mode = 2;          // 1 model variance, 2 moderated
  double gamma = 0;      // 0 adaptive; 1, 2, 4, 8, or >= 99 (max) for one T only
  double d0 = 10;        // prior degrees of freedom of the model variance in the moderated variance
  double winsor_k = 0;   // > 0: winsorize each gene's values at med +- k robust sd (label-free)
  int force_lm = 0;      // 1: the general kernel also for two-group designs (tests)
};
inline GsdOpts parse_gsd(const Rcpp::List& opt) {
  GsdOpts o;
  if (opt.containsElementNamed("mode")) o.mode = Rcpp::as<int>(opt["mode"]);
  if (opt.containsElementNamed("gamma")) o.gamma = Rcpp::as<double>(opt["gamma"]);
  if (opt.containsElementNamed("d0")) o.d0 = Rcpp::as<double>(opt["d0"]);
  if (opt.containsElementNamed("winsor_k")) o.winsor_k = Rcpp::as<double>(opt["winsor_k"]);
  if (opt.containsElementNamed("force_lm")) o.force_lm = Rcpp::as<int>(opt["force_lm"]);
  if (o.mode != 1 && o.mode != 2) Rcpp::stop("gsd mode must be 1 (model variance) or 2 (moderated)");
  if (!(o.d0 >= 0)) Rcpp::stop("gsd d0 must be non-negative");
  if (!(o.gamma == 0 || o.gamma == 1 || o.gamma == 2 || o.gamma == 4 || o.gamma == 8 || o.gamma >= 99)) Rcpp::stop("gsd gamma must be 0 (adaptive), 1, 2, 4, 8 or >= 99 (max)");
  if (!(o.winsor_k >= 0)) Rcpp::stop("winsor_k must be non-negative");
  return o;
}
// per relabeling b: z^2 = d^2 / v (0 when v <= 0) and its sums for gamma 1, 2, 4, 8 and the max (restrict: no aliasing, so it vectorizes)
__attribute__((target_clones("avx2", "default"), optimize("O3", "no-math-errno")))
static void gstat_sums(const double* __restrict xl, const double* __restrict r, const double* __restrict va, double* __restrict t1, double* __restrict t2,
                       double* __restrict t4, double* __restrict t8, double* __restrict tm, arma::uword nc) {
  for (arma::uword b = 0; b < nc; ++b) {
    const double d = xl[b] * r[b], v = va[b], pos = (double)(v > 0.0);
    const double z2 = pos * (d * d) / (v + (1.0 - pos)), z4 = z2 * z2;
    t1[b] += std::sqrt(z2); t2[b] += z2; t4[b] += z4; t8[b] += z4 * z4; tm[b] = tm[b] > z2 ? tm[b] : z2;
  }
}
// label-free winsorization of the per-gene statistic's input (robust option). Per gene (row of X, the m present samples),
// with med the median over the samples and s_mad = 1.4826 median |x - med|, each value is clipped to med +- k sd_a, sd_a = sqrt(max(V_a, s_mad^2))
// (V_a: the sample's modelled noise + tau2). It uses no labels, so the permutation test stays exact; an outlying sample's pull on a gene is
// bounded, a real group difference (which widens the MAD) is mostly kept.
static void winsorize_rows(arma::mat& X, const arma::mat& V, double k) {
  const arma::uword G = X.n_rows, m = X.n_cols; std::vector<double> r(m);
  auto median = [&](std::vector<double>& v) { const size_t h = v.size() / 2; std::nth_element(v.begin(), v.begin() + h, v.end()); double md = v[h];
    if (v.size() % 2 == 0) md = 0.5 * (md + *std::max_element(v.begin(), v.begin() + h)); return md; };
  for (arma::uword g = 0; g < G; ++g) {
    for (arma::uword j = 0; j < m; ++j) r[j] = X(g, j);
    const double med = median(r); for (arma::uword j = 0; j < m; ++j) r[j] = std::abs(X(g, j) - med);
    const double smad2 = std::pow(1.4826 * median(r), 2.0);
    for (arma::uword j = 0; j < m; ++j) { const double sd = std::sqrt(std::max(V(g, j), smad2)), lo = med - k * sd, hi = med + k * sd;
      X(g, j) = std::min(std::max(X(g, j), lo), hi); }
  }
}
// the per-relabeling sums -> T (single gamma, or adaptive: max over gamma of the standardized log T) -> per-cell p / z
static CellResult gstat_finish(const std::vector<double>& t1, const std::vector<double>& t2, const std::vector<double>& t4, const std::vector<double>& t8,
                               const std::vector<double>& tm, const std::vector<char>& adm, arma::uword nc, arma::uword B, arma::uword G, arma::uword m, bool single, int gsel) {
  CellResult out; const int ng = single ? 1 : 5;
  arma::mat T(ng, nc); T.fill(arma::datum::nan); const double mz2_obs = t2[0] / (double)G;
  for (arma::uword b = 0; b < nc; ++b) { if (!adm[b]) continue; const double tv[5] = {t1[b], t2[b], t4[b], t8[b], tm[b]};
    if (single) T(0, b) = tv[gsel]; else for (int k = 0; k < 5; ++k) T(k, b) = tv[k]; }
  arma::vec F(nc); F.fill(arma::datum::nan);
  if (single) F = T.row(0).t();
  else {   // adaptive: per gamma log T standardized by its mean and sd over the admissible permutations; F = max over gamma (continuous,
           // so the pooled null has no floor; a rank-based min p would put about 5 / (B + 1) of every cell's permutations at the ceiling)
    std::vector<arma::uword> ip; for (arma::uword b = 1; b < nc; ++b) if (adm[b]) ip.push_back(b);
    if (ip.size() < 2) return out;
    arma::vec mx(nc); mx.fill(-arma::datum::inf);
    for (int k = 0; k < 5; ++k) {
      arma::vec lt(nc); for (arma::uword b = 0; b < nc; ++b) lt[b] = adm[b] ? std::log(std::max(T(k, b), 1e-300)) : arma::datum::nan;
      double mu = 0; for (arma::uword b : ip) mu += lt[b]; mu /= (double)ip.size();
      double sd = 0; for (arma::uword b : ip) sd += (lt[b] - mu) * (lt[b] - mu); sd = std::sqrt(sd / ((double)ip.size() - 1.0));
      if (!(sd > 0)) continue;
      for (arma::uword b = 0; b < nc; ++b) if (adm[b]) mx[b] = std::max(mx[b], (lt[b] - mu) / sd);
    }
    for (arma::uword b = 0; b < nc; ++b) if (adm[b] && std::isfinite(mx[b])) F[b] = mx[b];
  }
  const double Fo = F[0]; if (!std::isfinite(Fo)) return out;
  arma::vec Fp = F.subvec(1, B); const arma::uvec ad = arma::find_finite(Fp); const arma::uword B_adm = ad.n_elem;
  if (B_adm == 0) return out;
  out.ok = true; out.F = Fo; out.n = (int)m; out.n_perm = (int)B_adm; out.shift = mz2_obs;
  arma::uword ge = 0; for (arma::uword i = 0; i < B_adm; ++i) if (Fp[ad[i]] >= Fo - 1e-12) ++ge;
  out.p = (ge + 1.0) / (B_adm + 1.0);
  const arma::vec Fa = Fp(ad); const double mu = arma::mean(Fa), sdv = arma::stddev(Fa);
  if (std::isfinite(sdv) && sdv > 0) { out.z = (Fo - mu) / sdv; out.zp = (Fp - mu) / sdv; }
  out.Fp = std::move(Fp);
  return out;
}
// O3 + AVX2 clone (selected at run time, "default" otherwise): vectorize the fused loops
__attribute__((target_clones("avx2", "default"), optimize("O3", "no-math-errno"))) CellResult test_genes(const cacoa_prof::Unit& u, const KernelInputs& in, const GsdOpts& o) {
  CellResult out; const arma::uword B = in.P.n_cols; const arma::uvec& s = u.present; const arma::uword m = s.n_elem, G = u.genes.n_elem;
  if (m < 3 || G < 1) return out;
  if (in.use_levels) {
    int n_den = 0, n_num = 0;
    for (arma::uword j = 0; j < m; ++j) { if (in.level_code[s[j]] == 1) ++n_den; else if (in.level_code[s[j]] == 2) ++n_num; }
    if (n_den < in.min_samp_per_level || n_num < in.min_samp_per_level) return out;
  }
  arma::mat X = u.x.cols(s), V = u.v.cols(s); V.each_col() += u.t2;
  if (o.winsor_k > 0) winsorize_rows(X, V, o.winsor_k);
  { arma::vec xm = arma::mean(X, 1); X.each_col() -= xm; }
  const bool mod = o.mode == 2; arma::vec rsX2; if (mod) rsX2 = arma::sum(arma::square(X), 1);   // total sum of squares of the centred profiles
  const arma::vec rsV = arma::sum(V, 1);
  const arma::uword nc = B + 1; arma::mat L(m, nc);
  for (arma::uword j = 0; j < m; ++j) { L(j, 0) = in.tg->dummy[s[j]]; for (arma::uword b = 0; b < B; ++b) L(j, b + 1) = in.tg->dummy[(arma::uword)(in.P(s[j], b) - 1)]; }
  const arma::rowvec nB = arma::sum(L, 0), nA = (double)m - nB;
  std::vector<char> adm(nc, 0); for (arma::uword b = 0; b < nc; ++b) adm[b] = nA[b] >= in.tg->min_per_level && nB[b] >= in.tg->min_per_level;
  if (!adm[0]) return out;
  const bool single = o.gamma > 0;   // (genes with zero variance contribute z = 0)
  const int gsel = !single ? -1 : (o.gamma >= 99 ? 4 : (o.gamma == 1 ? 0 : o.gamma == 2 ? 1 : o.gamma == 4 ? 2 : o.gamma == 8 ? 3 : -2));
  if (gsel == -2) return out;   // (validated in parse_gsd: no R error from a worker thread)
  // one fused pass per gene over all relabelings (contiguous in b, vectorizable), no genes x relabelings matrices:
  // xl_b = sum_j x_gj L_jb (sum over the B group), vl_b = sum_j s_gj L_jb; with x centred the pooled within-group sum of squares is sum x^2 - xl_b^2 r_b
  // (sum_A x = -xl_b), so the moderated variance needs no further pass
  std::vector<double> r(nc), inb(nc), ina(nc), inb2(nc), ina2(nc);
  for (arma::uword b = 0; b < nc; ++b) { const double nb = adm[b] ? nB[b] : 1.0, na = adm[b] ? nA[b] : 1.0; inb[b] = 1.0 / nb; ina[b] = 1.0 / na; inb2[b] = inb[b] * inb[b]; ina2[b] = ina[b] * ina[b]; r[b] = inb[b] + ina[b]; }
  std::vector<double> t1(nc, 0.0), t2(nc, 0.0), t4(nc, 0.0), t8(nc, 0.0), tm(nc, 0.0);
  const double d0 = o.d0, mm2 = (double)m - 2.0, wm = d0 / (d0 + mm2), ws = 1.0 / (d0 + mm2);   // moderated: (d0 s + (m-2) r s2/(m-2)) / (d0 + m - 2)
  std::vector<double> va(nc);
  // group sums by subset tables: the samples are split into chunks of <= 8; per gene and chunk the sums of x (and of s) over all
  // 2^len subsets are built incrementally, and the sum over the B group of relabeling b is one table lookup per chunk
  // (instead of m multiply-adds per relabeling; the label matrix only enters through the per-relabeling chunk masks)
  const arma::uword CW = 8, nq = (m + CW - 1) / CW;
  std::vector<arma::uword> q0(nq), ql(nq), toff(nq + 1, 0);
  for (arma::uword q = 0; q < nq; ++q) { q0[q] = q * CW; ql[q] = std::min(CW, m - q * CW); toff[q + 1] = toff[q] + ((arma::uword)1 << ql[q]); }
  std::vector<uint16_t> mask(nq * nc);
  for (arma::uword q = 0; q < nq; ++q) for (arma::uword b = 0; b < nc; ++b) { unsigned mk = 0; for (arma::uword t = 0; t < ql[q]; ++t) if (L(q0[q] + t, b) > 0.5) mk |= 1u << t; mask[q * nc + b] = (uint16_t)mk; }
  std::vector<double> tx(toff[nq]), tv(toff[nq]), xl(nc), vl(nc);
  for (arma::uword g = 0; g < G; ++g) {
    for (arma::uword q = 0; q < nq; ++q) {
      double* px = &tx[toff[q]]; double* pv = &tv[toff[q]]; px[0] = 0.0; pv[0] = 0.0; const unsigned n2 = 1u << ql[q];
      for (unsigned mk = 1; mk < n2; ++mk) { const unsigned low = mk & (0u - mk); const unsigned t = (unsigned)__builtin_ctz(mk);
        px[mk] = px[mk ^ low] + X(g, q0[q] + t); pv[mk] = pv[mk ^ low] + V(g, q0[q] + t); }
    }
    { const double* px = &tx[0]; const double* pv = &tv[0]; const uint16_t* mq = &mask[0];
      for (arma::uword b = 0; b < nc; ++b) { xl[b] = px[mq[b]]; vl[b] = pv[mq[b]]; } }
    for (arma::uword q = 1; q < nq; ++q) { const double* px = &tx[toff[q]]; const double* pv = &tv[toff[q]]; const uint16_t* mq = &mask[q * nc];
      for (arma::uword b = 0; b < nc; ++b) { xl[b] += px[mq[b]]; vl[b] += pv[mq[b]]; } }
    // branch-free (vectorizable) finish: variance per relabeling, then z^2 and the five sums; a gene with zero variance gives z = 0
    const double rv = rsV[g];
    for (arma::uword b = 0; b < nc; ++b) va[b] = vl[b] * inb2[b] + (rv - vl[b]) * ina2[b];
    if (mod) { const double rx2 = rsX2[g];
      for (arma::uword b = 0; b < nc; ++b) { double s2 = rx2 - xl[b] * xl[b] * r[b]; s2 = s2 > 0.0 ? s2 : 0.0; va[b] = wm * va[b] + ws * r[b] * s2; } }
    gstat_sums(xl.data(), r.data(), va.data(), t1.data(), t2.data(), t4.data(), t8.data(), tm.data(), nc);
  }
  return gstat_finish(t1, t2, t4, t8, tm, adm, nc, B, G, m, single, gsel);
}
// the per-gene statistic for a general single-df contrast c of any design X (continuous variables, covariates,
// level-vs-level contrasts), under the block relabeling or the Freedman-Lane scheme. Per relabeling b and gene g, with y the
// gene's centred profile over the m present samples (Freedman-Lane: its residual on the reduced design):
//   beta_gb = alpha_b' y;  model variance mv_gb = sum_j alpha_bj^2 (v_gj + tau2_g);
//   RSS_gb = |y|^2 - beta_gb^2 / |alpha_b|^2 - |Q_b' y|^2   (Q_b: basis of the rest of the fitted space, besides 1 and alpha_b);
//   moderated variance (d0 mv + cXc_b RSS) / (d0 + df_b), df_b = m - rank;  z^2 = beta^2 / variance.
// Block scheme: alpha_b = X_b (X_b'X_b)^- c on the relabeled design X_b (rows of the samples the global relabeling maps to),
// cXc_b = c'(X_b'X_b)^- c. Freedman-Lane: x* = H0 x + P_b (I - H0) x, so alpha_b = P_b' alpha, Q_b = P_b' Q0 (Q0: the reduced
// design's basis besides 1). For two groups (intercept + dummy) this is the subset-table kernel's statistic: alpha = 1_B / n_B - 1_A / n_A,
// cXc = 1 / n_A + 1 / n_B. The per-relabeling vectors depend only on the set of present samples; they are cached per set.
struct LmPlan {
  bool ok = false, fl = false; arma::uword m = 0, nc = 0, K = 0;
  std::vector<char> adm;
  std::vector<double> al, al2, q;      // al / al2: [j * nc + b]; q: [(k * m + j) * nc + b]
  std::vector<double> iaa, cxc, df;    // 1 / |alpha_b|^2, cXc_b, residual df (inadmissible: 0, 0, 1)
  arma::mat Q0;                         // Freedman-Lane: orthonormal basis of the reduced design (m x r0, 1 included), for y = (I - Q0 Q0') x
};
struct LmPlanCache { std::mutex mu; std::unordered_map<uint64_t, std::shared_ptr<const LmPlan>> map; };
// orthonormal basis of the columns of A (thin SVD, relative tolerance)
static arma::mat orth_basis(const arma::mat& A, double tol = 1e-9) {
  if (A.n_cols == 0) return arma::mat(A.n_rows, 0);
  arma::mat U, V; arma::vec sv; if (!arma::svd_econ(U, sv, V, A)) return arma::mat(A.n_rows, 0);
  const double s0 = sv.n_elem ? sv[0] : 0.0; arma::uword r = 0; for (arma::uword i = 0; i < sv.n_elem; ++i) if (sv[i] > tol * std::max(1.0, s0)) ++r;
  if (r == 0) return arma::mat(A.n_rows, 0);
  return U.cols(0, r - 1);
}
// the fitted space of Xs besides 1 and (optionally) a: basis of (I - W W') col(Xs), W = orth([1, a])
static arma::mat rest_basis(const arma::mat& Qx, const arma::vec* a) {
  const arma::uword m = Qx.n_rows; arma::mat W(m, a ? 2 : 1); W.col(0).fill(1.0); if (a) W.col(1) = *a;
  arma::mat Wo = orth_basis(W); arma::mat R = Qx - Wo * (Wo.t() * Qx);
  return orth_basis(R, 1e-7);
}
static bool has_ones(const arma::mat& Qx) { const arma::vec one(Qx.n_rows, arma::fill::ones); return arma::norm(one - Qx * (Qx.t() * one)) <= 1e-8 * std::sqrt((double)Qx.n_rows); }
static std::shared_ptr<const LmPlan> build_lm_plan(const arma::uvec& s, const KernelInputs& in) {
  auto pl = std::make_shared<LmPlan>(); const arma::uword m = s.n_elem, B = in.P.n_cols, nc = B + 1; pl->m = m; pl->nc = nc; pl->fl = in.freedman_lane;
  pl->adm.assign(nc, 0); pl->iaa.assign(nc, 0.0); pl->cxc.assign(nc, 0.0); pl->df.assign(nc, 1.0);
  std::vector<arma::vec> A(nc); std::vector<arma::mat> Qr(nc);
  auto design = [&](const arma::uvec& ps, arma::mat& Xp, arma::vec& cp) -> bool {
    if (in.use_levels) { int n_den = 0, n_num = 0; for (arma::uword j = 0; j < m; ++j) { if (in.level_code[ps[j]] == 1) ++n_den; else if (in.level_code[ps[j]] == 2) ++n_num; }
      if (n_den < in.min_samp_per_level || n_num < in.min_samp_per_level) return false; }
    arma::mat Xf = in.X.rows(ps); arma::uvec kp = arma::find(arma::sum(arma::abs(Xf), 0).t() > 1e-12);
    for (arma::uword j = 0; j < in.X.n_cols; ++j) if (!arma::any(kp == j) && std::abs(in.cvec[j]) > 1e-12) return false;
    Xp = Xf.cols(kp); cp = in.cvec(kp); return Xp.n_cols < m && estimable(Xp, cp);
  };
  arma::mat Xs; arma::vec cs; if (!design(s, Xs, cs)) return pl;
  Hat h = hat_of(Xs); const arma::vec a0 = h.A.t() * cs; const double cXc0 = arma::as_scalar(cs.t() * h.XtXi * cs);
  const arma::mat Qx0 = orth_basis(Xs); if (!has_ones(Qx0) || std::abs(arma::sum(a0)) > 1e-8 * arma::norm(a0) * std::sqrt((double)m)) return pl;   // needs 1 in the model and a contrast free of it
  if (in.freedman_lane) {
    arma::mat Qc, Rc; arma::qr(Qc, Rc, cs); const arma::mat Xr = Xs * Qc.cols(1, Qc.n_cols - 1);
    pl->Q0 = orth_basis(Xr); if (!has_ones(pl->Q0)) return pl;
    const arma::mat Qn = rest_basis(pl->Q0, nullptr);   // reduced space besides 1
    std::vector<arma::uvec> blocks; { std::map<int, std::vector<arma::uword>> by; for (arma::uword j = 0; j < m; ++j) if (in.inset[s[j]]) by[in.stratum[s[j]]].push_back(j);
      for (auto& kv : by) if (kv.second.size() >= 2) blocks.push_back(arma::uvec(kv.second)); }
    for (arma::uword b = 0; b < nc; ++b) {
      arma::uvec q = b == 0 ? arma::regspace<arma::uvec>(0, m - 1) : induced_perm(arma::conv_to<arma::uvec>::from(in.P.col(b - 1) - 1), s, blocks);
      arma::vec ab(m); arma::mat Qb(m, Qn.n_cols); for (arma::uword j = 0; j < m; ++j) { ab[q[j]] = a0[j]; Qb.row(q[j]) = Qn.row(j); }
      A[b] = ab; Qr[b] = Qb; pl->adm[b] = 1; pl->iaa[b] = 1.0 / arma::dot(a0, a0); pl->cxc[b] = cXc0; pl->df[b] = (double)m - (double)h.rank;
    }
  } else {
    for (arma::uword b = 0; b < nc; ++b) {
      arma::uvec ps(m); for (arma::uword j = 0; j < m; ++j) ps[j] = b == 0 ? s[j] : (arma::uword)(in.P(s[j], b - 1) - 1);
      arma::mat Xp; arma::vec cp; if (!design(ps, Xp, cp)) continue;
      // full-rank fast path (thin QR): X = QR, alpha = Q R^-T c, cXc = |R^-T c|^2; with 2 columns the fitted space is span{1, alpha}
      arma::mat Qx, Rx; bool fast = arma::qr_econ(Qx, Rx, Xp);
      if (fast) { const arma::vec dg = arma::abs(Rx.diag()); fast = dg.n_elem == Xp.n_cols && dg.min() > 1e-9 * std::max(1.0, dg.max()); }
      if (fast) {
        const arma::vec wv = arma::solve(arma::trimatl(Rx.t()), cp); arma::vec ab = Qx * wv;
        if (!has_ones(Qx) || std::abs(arma::sum(ab)) > 1e-8 * arma::norm(ab) * std::sqrt((double)m)) continue;
        A[b] = ab; Qr[b] = Xp.n_cols == 2 ? arma::mat(m, 0) : rest_basis(Qx, &ab); pl->adm[b] = 1; pl->iaa[b] = 1.0 / arma::dot(ab, ab);
        pl->cxc[b] = arma::dot(wv, wv); pl->df[b] = (double)m - (double)Xp.n_cols; continue;
      }
      Hat hp = hat_of(Xp); arma::vec ab = hp.A.t() * cp; Qx = orth_basis(Xp);
      if (!has_ones(Qx) || std::abs(arma::sum(ab)) > 1e-8 * arma::norm(ab) * std::sqrt((double)m)) continue;
      A[b] = ab; Qr[b] = rest_basis(Qx, &ab); pl->adm[b] = 1; pl->iaa[b] = 1.0 / arma::dot(ab, ab);
      pl->cxc[b] = arma::as_scalar(cp.t() * hp.XtXi * cp); pl->df[b] = (double)m - (double)hp.rank;
    }
  }
  if (!pl->adm[0]) return pl;
  arma::uword K = 0; for (arma::uword b = 0; b < nc; ++b) if (pl->adm[b]) K = std::max<arma::uword>(K, Qr[b].n_cols); pl->K = K;
  pl->al.assign(m * nc, 0.0); pl->al2.assign(m * nc, 0.0); pl->q.assign(K * m * nc, 0.0);
  for (arma::uword b = 0; b < nc; ++b) { if (!pl->adm[b]) continue;
    for (arma::uword j = 0; j < m; ++j) { pl->al[j * nc + b] = A[b][j]; pl->al2[j * nc + b] = A[b][j] * A[b][j]; }
    for (arma::uword k = 0; k < Qr[b].n_cols; ++k) for (arma::uword j = 0; j < m; ++j) pl->q[(k * m + j) * nc + b] = Qr[b](j, k); }
  pl->ok = true; return pl;
}
// per gene: beta, model variance and the K residual projections over all relabelings (contiguous in b, vectorizable)
__attribute__((target_clones("avx2", "default"), optimize("O3", "no-math-errno")))
static void lm_gene_sums(const double* __restrict y, const double* __restrict v, const double* __restrict al, const double* __restrict al2, const double* __restrict q,
                         arma::uword m, arma::uword nc, arma::uword K, double* __restrict beta, double* __restrict mv, double* __restrict ek) {
  for (arma::uword b = 0; b < nc; ++b) { beta[b] = 0.0; mv[b] = 0.0; }
  for (arma::uword j = 0; j < m; ++j) { const double yj = y[j], vj = v[j]; const double* a1 = al + j * nc; const double* a2 = al2 + j * nc;
    for (arma::uword b = 0; b < nc; ++b) { beta[b] += yj * a1[b]; mv[b] += vj * a2[b]; } }
  for (arma::uword k = 0; k < K; ++k) { double* e = ek + k * nc; for (arma::uword b = 0; b < nc; ++b) e[b] = 0.0;
    for (arma::uword j = 0; j < m; ++j) { const double yj = y[j]; const double* qq = q + (k * m + j) * nc; for (arma::uword b = 0; b < nc; ++b) e[b] += yj * qq[b]; } }
}
__attribute__((target_clones("avx2", "default"), optimize("O3", "no-math-errno")))
static void lm_gene_var(const double* __restrict beta, const double* __restrict mv, const double* __restrict ek, const double* __restrict iaa, const double* __restrict cxc,
                        const double* __restrict df, arma::uword nc, arma::uword K, double sy2, bool mod, double d0, double* __restrict va) {
  if (!mod) { for (arma::uword b = 0; b < nc; ++b) va[b] = mv[b]; return; }
  for (arma::uword b = 0; b < nc; ++b) va[b] = sy2 - beta[b] * beta[b] * iaa[b];
  for (arma::uword k = 0; k < K; ++k) { const double* e = ek + k * nc; for (arma::uword b = 0; b < nc; ++b) va[b] -= e[b] * e[b]; }
  for (arma::uword b = 0; b < nc; ++b) { const double rss = va[b] > 0.0 ? va[b] : 0.0; va[b] = (d0 * mv[b] + cxc[b] * rss) / (d0 + df[b]); }
}
static CellResult test_genes_lm(const cacoa_prof::Unit& u, const KernelInputs& in, const GsdOpts& o, LmPlanCache* cache) {
  CellResult out; const arma::uword B = in.P.n_cols; const arma::uvec& s = u.present; const arma::uword m = s.n_elem, G = u.genes.n_elem;
  if (m < 3 || G < 1) return out;
  std::shared_ptr<const LmPlan> pl;
  const bool cacheable = cache && in.n <= 64;
  uint64_t key = 0; if (cacheable) for (arma::uword j = 0; j < m; ++j) key |= (uint64_t)1 << s[j];
  if (cacheable) { std::lock_guard<std::mutex> lk(cache->mu); auto it = cache->map.find(key); if (it != cache->map.end()) pl = it->second; }
  if (!pl) { pl = build_lm_plan(s, in); if (cacheable) { std::lock_guard<std::mutex> lk(cache->mu); if (cache->map.size() < 4096) cache->map.emplace(key, pl); } }
  if (!pl->ok) return out;
  const arma::uword nc = pl->nc, K = pl->K;
  arma::mat X = u.x.cols(s), V = u.v.cols(s); V.each_col() += u.t2;
  if (o.winsor_k > 0) winsorize_rows(X, V, o.winsor_k);
  { arma::vec xm = arma::mean(X, 1); X.each_col() -= xm; }
  if (pl->fl) X -= (X * pl->Q0) * pl->Q0.t();   // residual on the reduced design (1 included)
  const arma::mat Xt = X.t(), Vt = V.t();       // m x G: one contiguous column per gene
  const bool mod = o.mode == 2; const bool single = o.gamma > 0;
  const int gsel = !single ? -1 : (o.gamma >= 99 ? 4 : (o.gamma == 1 ? 0 : o.gamma == 2 ? 1 : o.gamma == 4 ? 2 : o.gamma == 8 ? 3 : -2));
  if (gsel == -2) return out;   // (validated in parse_gsd: no R error from a worker thread)
  std::vector<double> t1(nc, 0.0), t2(nc, 0.0), t4(nc, 0.0), t8(nc, 0.0), tm(nc, 0.0), beta(nc), mv(nc), ek(K * nc), va(nc), one(nc, 1.0);
  for (arma::uword g = 0; g < G; ++g) {
    const double* y = Xt.colptr(g); double sy2 = 0; for (arma::uword j = 0; j < m; ++j) sy2 += y[j] * y[j];
    lm_gene_sums(y, Vt.colptr(g), pl->al.data(), pl->al2.data(), pl->q.data(), m, nc, K, beta.data(), mv.data(), ek.data());
    lm_gene_var(beta.data(), mv.data(), ek.data(), pl->iaa.data(), pl->cxc.data(), pl->df.data(), nc, K, sy2, mod, o.d0, va.data());
    gstat_sums(beta.data(), one.data(), va.data(), t1.data(), t2.data(), t4.data(), t8.data(), tm.data(), nc);
  }
  return gstat_finish(t1, t2, t4, t8, tm, pl->adm, nc, B, G, m, single, gsel);
}

} // namespace cacoa_gsd

// Normalized profiles, noise, weights and pair distances of a list of units (0-based cell ids per unit; for the
// noise-curve fit and tests), with one conversion of the count matrix.
// [[Rcpp::export]]
Rcpp::List cf_unit_profiles(const Eigen::SparseMatrix<double>& cm, Rcpp::IntegerVector sample_per_cell, Rcpp::List units, Rcpp::NumericVector lib,
                            int n_samples, int min_n_obs_per_samp, const arma::imat& pairs, Rcpp::List opts) {
  cacoa_prof::ProfOpts o = cacoa_prof::parse_opts(opts);
  if (cm.cols() != sample_per_cell.size() || cm.cols() != lib.size()) Rcpp::stop("cm, sample_per_cell and lib must cover the same cells");
  std::vector<int> soc(sample_per_cell.size());
  for (int i = 0; i < sample_per_cell.size(); ++i) { const int v = sample_per_cell[i]; if (Rcpp::IntegerVector::is_na(v) || v < 1 || v > n_samples) Rcpp::stop("sample_per_cell must be a 1-based factor"); soc[i] = v - 1; }
  std::vector<double> lb(lib.begin(), lib.end());
  Rcpp::List out(units.size());
  for (int k = 0; k < units.size(); ++k) {
    Rcpp::IntegerVector ids = units[k]; std::vector<int> idv(ids.begin(), ids.end());
    for (int id : idv) if (id < 0 || id >= cm.cols()) Rcpp::stop("unit %d has ids outside the cells", k + 1);
    cacoa_prof::Unit u = cacoa_prof::unit_profile(cm, soc, idv, lb, n_samples, min_n_obs_per_samp, o);
    arma::vec y = cacoa_prof::unit_distances(u, pairs, o);
    out[k] = Rcpp::List::create(Rcpp::Named("genes") = u.genes, Rcpp::Named("present") = u.present, Rcpp::Named("x") = u.x, Rcpp::Named("v") = u.v,
                                Rcpp::Named("w") = u.w, Rcpp::Named("mu") = u.mu, Rcpp::Named("c") = u.c, Rcpp::Named("nobs") = u.nobs, Rcpp::Named("y") = y,
                                Rcpp::Named("tau2") = u.t2);
  }
  return out;
}

// Moment estimates of the within-unit heterogeneity phi for the noise-curve fit (fitShiftNoiseCurves()), one pass over each
// unit's cells: per unit, sample a with >= min_n_obs cells (units with >= 3 such samples) and gene g with Y_ga > 0,
// v = [sum_i y_i^2 - 2 p sum_i y_i m_i + p^2 sum_i m_i^2 - Y] / (p^2 sum_i m_i^2), p = Y_ga / N_a, paired with
// mu = log10(pbar_g) (pbar pooled over the unit's present samples). Same estimator as the former R loop.
// [[Rcpp::export]]
Rcpp::List cf_phi_moments(const Eigen::SparseMatrix<double>& cm, Rcpp::IntegerVector sample_per_cell, Rcpp::List units, Rcpp::NumericVector lib,
                          int n_samples, int min_n_obs_per_samp) {
  if (cm.cols() != sample_per_cell.size() || cm.cols() != lib.size()) Rcpp::stop("cm, sample_per_cell and lib must cover the same cells");
  std::vector<int> soc(sample_per_cell.size());
  for (int i = 0; i < sample_per_cell.size(); ++i) { const int v = sample_per_cell[i]; if (Rcpp::IntegerVector::is_na(v) || v < 1 || v > n_samples) Rcpp::stop("sample_per_cell must be a 1-based factor"); soc[i] = v - 1; }
  const int G = cm.rows(); std::vector<double> out_mu, out_v; int n_used = 0;
  arma::mat Y(G, n_samples), Q(G, n_samples), YM(G, n_samples);
  for (int k = 0; k < units.size(); ++k) {
    Rcpp::IntegerVector ids = units[k];
    std::vector<int> nobs(n_samples, 0); for (int id : ids) { if (id < 0 || id >= cm.cols()) Rcpp::stop("unit %d has ids outside the cells", k + 1); nobs[soc[id]]++; }
    std::vector<int> pres; for (int a = 0; a < n_samples; ++a) if (nobs[a] >= min_n_obs_per_samp) pres.push_back(a);
    if (pres.size() < 3) continue;
    ++n_used; Y.zeros(); Q.zeros(); YM.zeros(); std::vector<double> N(n_samples, 0.0), S2(n_samples, 0.0);
    for (int id : ids) {
      const int a = soc[id]; if (nobs[a] < min_n_obs_per_samp) continue; const double m = lib[id]; N[a] += m; S2[a] += m * m;
      for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(cm, id); it; ++it) { const double y = it.value(); Y(it.row(), a) += y; Q(it.row(), a) += y * y; YM(it.row(), a) += y * m; }
    }
    double Ntot = 0; for (int a : pres) Ntot += N[a];
    for (int g = 0; g < G; ++g) {
      double Ysum = 0; for (int a : pres) Ysum += Y(g, a);
      if (!(Ysum > 0)) continue;
      const double mu = std::log10(Ysum / Ntot);
      for (int a : pres) {
        const double yg = Y(g, a); if (!(yg > 0)) continue;
        const double p = yg / N[a], v = (Q(g, a) - 2 * p * YM(g, a) + p * p * S2[a] - yg) / (p * p * S2[a]);
        if (std::isfinite(v)) { out_mu.push_back(mu); out_v.push_back(v); }
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("mu") = out_mu, Rcpp::Named("v") = out_v, Rcpp::Named("n_units") = n_used);
}

// Streaming kernel: neighbourhood profiles, distances and the test per cell, nothing cells-wide kept in memory.
// nn_ids are 0-based cell indices (graph adjacency) unless nn_one_based; pairs are 0-based sample indices.
// [[Rcpp::export]]
Rcpp::List cluster_free_shift_stream(const Eigen::SparseMatrix<double>& cm, Rcpp::IntegerVector sample_per_cell, Rcpp::List nn_ids, bool nn_one_based,
                                     int min_n_obs_per_samp, std::string dist, bool log_vecs,
                                     const arma::imat& pairs, int n_samples, const arma::mat& X, const arma::vec& cvec,
                                     const arma::ivec& level_code, int min_samp_per_level, const arma::ivec& stratum, const arma::uvec& inset,
                                     const arma::imat& P, bool freedman_lane, bool bias_correct, int n_cores, int robust = 0, double robust_k = 1.345,
                                     Rcpp::NumericVector z_grid = Rcpp::NumericVector::create(), Rcpp::IntegerVector group = Rcpp::IntegerVector::create(), bool keep_perm = false,
                                     bool average = false, Rcpp::Nullable<Rcpp::List> avg_nn = R_NilValue, bool vectorize = true, bool standardize = true,
                                     Rcpp::Nullable<Rcpp::List> prof_opts = R_NilValue, Rcpp::NumericVector lib = Rcpp::NumericVector::create(),
                                     Rcpp::Nullable<Rcpp::List> gsd_opts = R_NilValue, bool avg_count = false) {
  const arma::uword n = (arma::uword)n_samples, m_cells = nn_ids.size(), n_pairs = pairs.n_rows;
  if ((arma::uword)X.n_rows != n || P.n_rows != n || level_code.n_elem != n || stratum.n_elem != n || inset.n_elem != n) Rcpp::stop("dimension mismatch in cluster_free_shift_stream");
  if (cm.cols() != sample_per_cell.size()) Rcpp::stop("cm must have one column per cell");
  if (!(dist == "cor" || dist == "cosine" || dist == "js")) Rcpp::stop("dist must be one of cor, cosine, js");
  const bool normalized = prof_opts.isNotNull(); cacoa_prof::ProfOpts po; if (normalized) po = cacoa_prof::parse_opts(Rcpp::List(prof_opts));
  const bool gsd = gsd_opts.isNotNull(); cacoa_gsd::GsdOpts go; if (gsd) go = cacoa_gsd::parse_gsd(Rcpp::List(gsd_opts));
  if (gsd && (!normalized || !po.noise)) Rcpp::stop("the per-gene statistic needs normalized profiles with the noise model");
  if (gsd && robust != 0) Rcpp::stop("the per-gene statistic has no robust re-weighting (use its winsorization)");
  std::vector<double> lb(lib.begin(), lib.end()); if (normalized && (int)lb.size() != cm.cols()) Rcpp::stop("lib must have one entry per cell");
  std::vector<int> sample_of_cell(cm.cols());
  for (int i = 0; i < cm.cols(); ++i) { int v = sample_per_cell[i]; if (Rcpp::IntegerVector::is_na(v) || v < 1 || v > n_samples) Rcpp::stop("sample_per_cell must be a 1-based factor"); sample_of_cell[i] = v - 1; }
  std::vector<std::vector<int>> nn(m_cells);
  for (arma::uword k = 0; k < m_cells; ++k) {
    Rcpp::IntegerVector ids = nn_ids[k]; nn[k].reserve(ids.size());
    for (int t = 0; t < ids.size(); ++t) { if (Rcpp::IntegerVector::is_na(ids[t])) continue; int v = nn_one_based ? ids[t] - 1 : ids[t];
      if (v < 0 || v >= cm.cols()) Rcpp::stop("nn_ids[[%d]] has an index outside the cells", (int)k + 1); nn[k].push_back(v); }
  }
  KernelInputs in{pairs, n, X, cvec, level_code, min_samp_per_level, stratum, inset, P, freedman_lane, bias_correct, arma::any(level_code >= 0), robust, robust_k};
  CompleteCase cc; if (vectorize) cc = complete_case(in); in.cc = &cc; TwoGroup tg; if (vectorize) tg = two_group(in); in.tg = &tg;
  arma::ivec grp(group.begin(), group.size()); const int n_groups = grp.n_elem ? (int)grp.max() + 1 : 0;
  if (grp.n_elem && grp.n_elem != m_cells) Rcpp::stop("group must have one entry per cell");
  std::vector<std::vector<int>> nn_avg; if (average) { if (avg_nn.isNull()) Rcpp::stop("average needs avg_nn"); nn_avg = cacoa_cf::neighbour_lists(Rcpp::List(avg_nn), false, (int)m_cells, true); }
  ShiftAccum acc(m_cells, P.n_cols, arma::vec(z_grid.begin(), z_grid.size()), grp, n_groups, keep_perm, average, average ? &nn_avg : nullptr, standardize);
  if (avg_count) { if (!average) Rcpp::stop("thinned units need the averaging (avg_count without average)"); acc.acc.by_count = true; acc.acc.cnt.assign(m_cells, 0); }
  const bool gs_tables = gsd && tg.ok && go.force_lm == 0;   // two groups: subset tables; otherwise the general kernel
  cacoa_gsd::LmPlanCache lm_cache;
  cacoa::parallelFor(0, (int)m_cells, [&](int k) {
    const auto& ids = nn[k];
    if ((int)ids.size() < min_n_obs_per_samp) return;
    if (normalized) {
      cacoa_prof::Unit u = cacoa_prof::unit_profile(cm, sample_of_cell, ids, lb, n_samples, min_n_obs_per_samp, po);
      if (gsd) { acc.add(k, gs_tables ? cacoa_gsd::test_genes(u, in, go) : cacoa_gsd::test_genes_lm(u, in, go, &lm_cache)); return; }
      acc.add(k, test_cell(cacoa_prof::unit_distances(u, pairs, po), in));
      return;
    }
    std::vector<unsigned> n_obs = count_values(sample_of_cell, ids, n_samples);
    Eigen::MatrixXd prof = collapseMatrixNorm(cm, sample_of_cell, ids, n_obs);
    if (log_vecs) for (int i = 0; i < prof.size(); ++i) prof(i) = std::log10(1e3 * prof(i) + 1.0);
    arma::vec y(n_pairs); y.fill(arma::datum::nan);
    for (arma::uword r = 0; r < n_pairs; ++r) {
      int a = pairs(r, 0), b = pairs(r, 1);
      if ((int)n_obs[a] < min_n_obs_per_samp || (int)n_obs[b] < min_n_obs_per_samp) continue;
      Eigen::VectorXd v1 = prof.col(a), v2 = prof.col(b);
      y[r] = estimateVectorDistance(v1, v2, dist);
    }
    acc.add(k, test_cell(y, in));
  }, n_cores, false);
  return acc.result(n_cores);
}

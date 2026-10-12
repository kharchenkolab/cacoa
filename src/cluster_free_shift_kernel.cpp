// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <mutex>
#include <map>
#include "parallel.h"
#include "lm_common.h"
#include "cf_common.h"
#include "gower_stats.h"

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

struct KernelInputs {
  const arma::imat& pairs; arma::uword n; const arma::mat& X; const arma::vec& cvec; const arma::ivec& level_code; int min_samp_per_level;
  const arma::ivec& stratum; const arma::uvec& inset; const arma::imat& P; bool freedman_lane; bool bias_correct; bool use_levels;
  int robust; double robust_k;
};

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

// Per-permutation extremes of the z-scores over cells (max-statistic adjustment); when a threshold grid is given,
// per-permutation histograms of the permuted z-scores pooled over cells (hist(b, k) counts the cells whose permuted
// z under permutation b falls in [grid[k-1], grid[k]); k = 0: below grid[0]; k = n_grid: at or above the last
// threshold) for the pooled permutation FDR (permutationFDR() in R); and when cell groups are given, the sum and the
// count of the permuted z-scores per group and permutation (group_sum(b, g), group_count(b, g)) next to the
// observed sum per group, for the exact group-level permutation test (the mean z of a group against its own
// permutation distribution). A cell contributes to a permutation only where it is tested under that relabeling.
struct Accum {
  arma::vec stat, pval, shift, z; arma::ivec nsamp, nperm; arma::vec mx, mn; arma::imat hist; arma::vec grid;
  const arma::ivec& group; arma::mat gsum; arma::imat gcount; arma::vec gobs; arma::ivec gn; std::mutex lock;
  arma::mat perm;   // permutations x cells matrix of the permuted statistics (keep_perm only; NaN where skipped)
  Accum(arma::uword m, arma::uword B, const arma::vec& z_grid, const arma::ivec& group_in, int n_groups, bool keep_perm = false)
    : stat(m), pval(m), shift(m), z(m), nsamp(m, arma::fill::zeros), nperm(m, arma::fill::zeros), mx(B), mn(B), grid(z_grid), group(group_in) {
    stat.fill(arma::datum::nan); pval.fill(arma::datum::nan); shift.fill(arma::datum::nan); z.fill(arma::datum::nan); mx.fill(-arma::datum::inf); mn.fill(arma::datum::inf);
    if (grid.n_elem) hist.zeros(B, grid.n_elem + 1);
    if (n_groups > 0) { gsum.zeros(B, n_groups); gcount.zeros(B, n_groups); gobs.zeros(n_groups); gn.zeros(n_groups); }
    if (keep_perm) { perm.set_size(B, m); perm.fill(arma::datum::nan); } }
  void add(arma::uword k, const CellResult& r) {
    if (!r.ok) return;
    stat[k] = r.F; pval[k] = r.p; shift[k] = r.shift; nsamp[k] = r.n; nperm[k] = r.n_perm; z[k] = r.z;
    if (perm.n_elem && r.Fp.n_elem == perm.n_rows) perm.col(k) = r.Fp;   // column k is written by one worker only
    if (r.zp.n_elem) {
      std::vector<arma::uword> bins;
      if (grid.n_elem) { bins.resize(r.zp.n_elem); for (arma::uword b = 0; b < r.zp.n_elem; ++b) bins[b] = std::isfinite(r.zp[b]) ? (arma::uword)(std::upper_bound(grid.begin(), grid.end(), r.zp[b]) - grid.begin()) : 0; }
      const int g = (group.n_elem > k) ? group[k] : -1;
      const bool in_group = g >= 0 && (arma::uword)g < gsum.n_cols && std::isfinite(r.z);
      std::lock_guard<std::mutex> lk(lock);
      for (arma::uword b = 0; b < r.zp.n_elem; ++b) {
        if (!std::isfinite(r.zp[b])) continue;   // the cell is not tested under this relabeling
        if (r.zp[b] > mx[b]) mx[b] = r.zp[b];
        if (r.zp[b] < mn[b]) mn[b] = r.zp[b];
        if (grid.n_elem) ++hist(b, bins[b]);
        if (in_group) { gsum(b, g) += r.zp[b]; ++gcount(b, g); }
      }
      if (in_group) { gobs[g] += r.z; ++gn[g]; }
    }
  }
  Rcpp::List result() const {
    return Rcpp::List::create(Rcpp::_["stat"] = stat, Rcpp::_["p"] = pval, Rcpp::_["shift"] = shift, Rcpp::_["n"] = nsamp, Rcpp::_["n_perm"] = nperm, Rcpp::_["z"] = z,
                              Rcpp::_["max"] = mx, Rcpp::_["min"] = mn, Rcpp::_["hist"] = hist, Rcpp::_["grid"] = grid,
                              Rcpp::_["group_sum"] = gsum, Rcpp::_["group_count"] = gcount, Rcpp::_["group_obs"] = gobs, Rcpp::_["group_n"] = gn, Rcpp::_["perm"] = perm);
  }
};

} // namespace

// Batched kernel on a precomputed pair-layout distance matrix Y (pairs x cells); used as the reference for the
// streaming kernel and by tests.
// [[Rcpp::export]]
Rcpp::List cluster_free_shift_batch(const arma::mat& Y, const arma::imat& pairs, int n_samples, const arma::mat& X, const arma::vec& cvec,
                                    const arma::ivec& level_code, int min_samp_per_level, const arma::ivec& stratum, const arma::uvec& inset,
                                    const arma::imat& P, bool freedman_lane, bool bias_correct, int n_cores, int robust = 0, double robust_k = 1.345,
                                    Rcpp::NumericVector z_grid = Rcpp::NumericVector::create(), Rcpp::IntegerVector group = Rcpp::IntegerVector::create(), bool keep_perm = false) {
  const arma::uword m_cells = Y.n_cols, n = (arma::uword)n_samples;
  if ((arma::uword)X.n_rows != n || P.n_rows != n || level_code.n_elem != n || stratum.n_elem != n || inset.n_elem != n) Rcpp::stop("dimension mismatch in cluster_free_shift_batch");
  KernelInputs in{pairs, n, X, cvec, level_code, min_samp_per_level, stratum, inset, P, freedman_lane, bias_correct, arma::any(level_code >= 0), robust, robust_k};
  arma::ivec grp(group.begin(), group.size()); const int n_groups = grp.n_elem ? (int)grp.max() + 1 : 0;
  Accum acc(m_cells, P.n_cols, arma::vec(z_grid.begin(), z_grid.size()), grp, n_groups, keep_perm);
  cacoa::parallelFor(0, (int)m_cells, [&](int k) { acc.add(k, test_cell(Y.col(k), in)); }, n_cores, false);
  return acc.result();
}

// Streaming kernel: neighbourhood profiles, distances and the test per cell, nothing cells-wide kept in memory.
// nn_ids are 0-based cell indices (graph adjacency) unless nn_one_based; pairs are 0-based sample indices.
// [[Rcpp::export]]
Rcpp::List cluster_free_shift_stream(const Eigen::SparseMatrix<double>& cm, Rcpp::IntegerVector sample_per_cell, Rcpp::List nn_ids, bool nn_one_based,
                                     int min_n_obs_per_samp, std::string dist, bool log_vecs,
                                     const arma::imat& pairs, int n_samples, const arma::mat& X, const arma::vec& cvec,
                                     const arma::ivec& level_code, int min_samp_per_level, const arma::ivec& stratum, const arma::uvec& inset,
                                     const arma::imat& P, bool freedman_lane, bool bias_correct, int n_cores, int robust = 0, double robust_k = 1.345,
                                     Rcpp::NumericVector z_grid = Rcpp::NumericVector::create(), Rcpp::IntegerVector group = Rcpp::IntegerVector::create(), bool keep_perm = false) {
  const arma::uword n = (arma::uword)n_samples, m_cells = nn_ids.size(), n_pairs = pairs.n_rows;
  if ((arma::uword)X.n_rows != n || P.n_rows != n || level_code.n_elem != n || stratum.n_elem != n || inset.n_elem != n) Rcpp::stop("dimension mismatch in cluster_free_shift_stream");
  if (cm.cols() != sample_per_cell.size()) Rcpp::stop("cm must have one column per cell");
  if (!(dist == "cor" || dist == "cosine" || dist == "js")) Rcpp::stop("dist must be one of cor, cosine, js");
  std::vector<int> sample_of_cell(cm.cols());
  for (int i = 0; i < cm.cols(); ++i) { int v = sample_per_cell[i]; if (Rcpp::IntegerVector::is_na(v) || v < 1 || v > n_samples) Rcpp::stop("sample_per_cell must be a 1-based factor"); sample_of_cell[i] = v - 1; }
  std::vector<std::vector<int>> nn(m_cells);
  for (arma::uword k = 0; k < m_cells; ++k) {
    Rcpp::IntegerVector ids = nn_ids[k]; nn[k].reserve(ids.size());
    for (int t = 0; t < ids.size(); ++t) { if (Rcpp::IntegerVector::is_na(ids[t])) continue; int v = nn_one_based ? ids[t] - 1 : ids[t];
      if (v < 0 || v >= cm.cols()) Rcpp::stop("nn_ids[[%d]] has an index outside the cells", (int)k + 1); nn[k].push_back(v); }
  }
  KernelInputs in{pairs, n, X, cvec, level_code, min_samp_per_level, stratum, inset, P, freedman_lane, bias_correct, arma::any(level_code >= 0), robust, robust_k};
  arma::ivec grp(group.begin(), group.size()); const int n_groups = grp.n_elem ? (int)grp.max() + 1 : 0;
  if (grp.n_elem && grp.n_elem != m_cells) Rcpp::stop("group must have one entry per cell");
  Accum acc(m_cells, P.n_cols, arma::vec(z_grid.begin(), z_grid.size()), grp, n_groups, keep_perm);
  cacoa::parallelFor(0, (int)m_cells, [&](int k) {
    const auto& ids = nn[k];
    if ((int)ids.size() < min_n_obs_per_samp) return;
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
  return acc.result();
}

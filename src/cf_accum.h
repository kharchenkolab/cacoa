// Shared accumulator of the cluster-free per-cell permutation tests (density map and expression-shift map).
//
// Every kernel hands over, per cell, the observed statistic and its vector of permuted statistics (NaN where the
// cell is not tested under a relabeling). The accumulator computes the exact permutation p-value, the
// permutation-standardized statistic z (observed minus permutation mean over permutation sd) and its permuted
// values, and accumulates what the pooled procedures need without keeping the permutations x cells matrix:
//   - per-permutation histograms of the permuted z (or |z|) pooled over cells on a fixed grid (the pooled
//     permutation p-values and FDR, pooledPermutationP() / permutationFDR() in R),
//   - per-permutation extremes of z over cells (max-statistic adjustment),
//   - per-group sums and counts of the permuted z and the observed sum (exact group-level permutation test).
// With `average` the statistics are first averaged over each cell and its graph neighbours (self included in
// `nn_graph`), observed and permuted alike, in a float32 cells x permutations buffer filled by scatter-adds
// (NaN counted as 0, divisor the neighbourhood size), and the per-cell quantities above are computed from the
// averaged values in finalize(): the exact test of the one-hop-larger unit. The permutation matrix is kept only
// on request (keep_perm), of the averaged statistics when averaging.
#ifndef CACOA_CF_ACCUM_H
#define CACOA_CF_ACCUM_H

#include <RcppArmadillo.h>
#include <mutex>
#include <vector>
#include <algorithm>
#include <cmath>
#include "parallel.h"

namespace cacoa_cf {

struct Accum {
  arma::uword m, B;
  arma::vec stat, pval, z; arma::ivec nperm;
  arma::vec mx, mn; arma::imat hist; arma::vec grid; bool use_abs; bool negate; int alt_code; bool standardize;   // histogram of the standardized z (standardize) or of the raw statistic; |.| (use_abs), negated (negate); alt_code: 0 two-sided, 1 greater, 2 less
  arma::ivec group; arma::mat gsum; arma::imat gcount; arma::vec gobs; arma::ivec gn;
  arma::mat perm;                                     // B x m (keep_perm)
  bool average; const std::vector<std::vector<int>>* nn; std::vector<float> avg; std::vector<double> avg_obs; std::vector<char> seen;
  std::mutex lock; std::vector<std::mutex> stripes;

  Accum(arma::uword m_, arma::uword B_, const arma::vec& z_grid, bool use_abs_, int alt_code_, const arma::ivec& group_, int n_groups,
        bool keep_perm, bool average_, const std::vector<std::vector<int>>* nn_graph, bool standardize_ = true)
    : m(m_), B(B_), stat(m_), pval(m_), z(m_), nperm(m_, arma::fill::zeros), mx(B_), mn(B_), grid(z_grid), use_abs(use_abs_), negate(alt_code_ == 2), alt_code(alt_code_), standardize(standardize_),
      group(group_), average(average_), nn(nn_graph), stripes(average_ ? 1024 : 0) {
    stat.fill(arma::datum::nan); pval.fill(arma::datum::nan); z.fill(arma::datum::nan); mx.fill(-arma::datum::inf); mn.fill(arma::datum::inf);
    if (grid.n_elem) hist.zeros(B, grid.n_elem + 1);
    if (n_groups > 0) { gsum.zeros(B, n_groups); gcount.zeros(B, n_groups); gobs.zeros(n_groups); gn.zeros(n_groups); }
    if (keep_perm) { perm.set_size(B, m); perm.fill(arma::datum::nan); }
    if (average) {
      if (!nn || nn->size() != m) Rcpp::stop("averaging needs one neighbour list per cell");
      avg.assign((size_t)m * B, 0.0f); avg_obs.assign(m, 0.0); seen.assign(m, 0);
    }
  }

  // per-cell quantities from an observed statistic and its permuted vector; thread-safe
  void cell(arma::uword k, double obs, const arma::vec& fp) {
    if (!std::isfinite(obs) || fp.n_elem != B) return;
    const arma::uvec adm = arma::find_finite(fp); const arma::uword B_adm = adm.n_elem;
    if (B_adm == 0) return;
    arma::uword ge = 0, le = 0, ge_abs = 0;
    const double tol = (average ? 1e-6 : 1e-12) * std::max(1.0, std::abs(obs));   // ties count; float32 accumulation under averaging
    for (arma::uword i = 0; i < B_adm; ++i) { const double v = fp[adm[i]]; if (v >= obs - tol) ++ge; if (v <= obs + tol) ++le; if (std::abs(v) >= std::abs(obs) - tol) ++ge_abs; }
    const double p = ((alt_code == 0 ? ge_abs : alt_code == 1 ? ge : le) + 1.0) / (B_adm + 1.0);
    const arma::vec fa = fp(adm); const double mu = arma::mean(fa), sdv = arma::stddev(fa);
    double zk = arma::datum::nan; arma::vec zp;
    if (std::isfinite(sdv) && sdv > 0) { zk = (obs - mu) / sdv; zp = (fp - mu) / sdv; }
    stat[k] = obs; pval[k] = p; nperm[k] = (int)B_adm; z[k] = zk;
    if (perm.n_elem) perm.col(k) = fp;                       // column k is written by one worker only
    if (!zp.n_elem) return;
    std::vector<arma::uword> bins;
    if (grid.n_elem) { bins.resize(B); for (arma::uword b = 0; b < B; ++b) { const double s = standardize ? zp[b] : fp[b]; const double v = use_abs ? std::abs(s) : (negate ? -s : s);
      bins[b] = std::isfinite(v) ? (arma::uword)(std::upper_bound(grid.begin(), grid.end(), v) - grid.begin()) : 0; } }
    const int g = (group.n_elem > k) ? group[k] : -1;
    const bool in_group = g >= 0 && (arma::uword)g < gsum.n_cols;
    std::lock_guard<std::mutex> lk(lock);
    for (arma::uword b = 0; b < B; ++b) {
      if (!std::isfinite(zp[b])) continue;
      if (zp[b] > mx[b]) mx[b] = zp[b];
      if (zp[b] < mn[b]) mn[b] = zp[b];
      if (grid.n_elem && std::isfinite(standardize ? zp[b] : fp[b])) ++hist(b, bins[b]);
      if (in_group) { gsum(b, g) += zp[b]; ++gcount(b, g); }
    }
    if (in_group) { gobs[g] += zk; ++gn[g]; }
  }

  // entry point of the kernels: direct per-cell accumulation, or scatter-add into the neighbours' averages
  void add(arma::uword k, double obs, const arma::vec& fp) {
    if (!average) { cell(k, obs, fp); return; }
    if (!std::isfinite(obs) || fp.n_elem != B) return;
    seen[k] = 1;
    for (int j : (*nn)[k]) {
      if (j < 0 || (arma::uword)j >= m) continue;
      std::lock_guard<std::mutex> lk(stripes[(size_t)j % stripes.size()]);
      float* row = &avg[(size_t)j * B];
      for (arma::uword b = 0; b < B; ++b) { const double v = fp[b]; if (std::isfinite(v)) row[b] += (float)v; }
      avg_obs[j] += obs;
    }
  }

  // averaging: divide by the neighbourhood size and run the per-cell accumulation on the averaged statistics
  void finalize(int n_cores) {
    if (!average) return;
    cacoa::parallelFor(0, (int)m, [&](int k) {
      const std::vector<int>& ids = (*nn)[k];
      if (ids.empty()) return;
      bool self_tested = seen[k] != 0;
      if (!self_tested) return;                                      // a cell not tested on its own stays NA (its own vector is absent)
      const double inv = 1.0 / (double)ids.size();
      arma::vec fp(B); const float* row = &avg[(size_t)k * B];
      for (arma::uword b = 0; b < B; ++b) fp[b] = (double)row[b] * inv;
      cell((arma::uword)k, avg_obs[k] * inv, fp);
    }, n_cores, false);
    std::vector<float>().swap(avg);
  }

  Rcpp::List result() const {
    return Rcpp::List::create(Rcpp::_["stat"] = stat, Rcpp::_["p"] = pval, Rcpp::_["z"] = z, Rcpp::_["n_perm"] = nperm,
                              Rcpp::_["max"] = mx, Rcpp::_["min"] = mn, Rcpp::_["hist"] = hist, Rcpp::_["grid"] = grid,
                              Rcpp::_["group_sum"] = gsum, Rcpp::_["group_count"] = gcount, Rcpp::_["group_obs"] = gobs, Rcpp::_["group_n"] = gn,
                              Rcpp::_["perm"] = perm, Rcpp::_["averaged"] = average);
  }
};

// 0-based neighbour lists (self included) from an R list of integer vectors
inline std::vector<std::vector<int>> neighbour_lists(Rcpp::List nn_ids, bool one_based, int n_cells, bool add_self) {
  std::vector<std::vector<int>> nn(nn_ids.size());
  for (int k = 0; k < nn_ids.size(); ++k) {
    Rcpp::IntegerVector ids = nn_ids[k]; nn[k].reserve(ids.size() + 1); bool has_self = false;
    for (int t = 0; t < ids.size(); ++t) { if (Rcpp::IntegerVector::is_na(ids[t])) continue; int v = one_based ? ids[t] - 1 : ids[t];
      if (v < 0 || v >= n_cells) Rcpp::stop("neighbour list %d has an index outside the cells", k + 1); if (v == k) has_self = true; nn[k].push_back(v); }
    if (add_self && !has_self) nn[k].push_back(k);
  }
  return nn;
}

} // namespace cacoa_cf
#endif

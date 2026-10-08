// Shared helpers for the cluster-free code (cluster_free.cpp, expression_shifts.cpp).
// Header-only (inline) so that both translation units can use one definition.
#ifndef CACOA_CF_COMMON_H
#define CACOA_CF_COMMON_H

#include <RcppEigen.h>
#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <numeric>
#include <string>
#include <vector>

namespace cf {
// tolerance used when comparing permutation statistics
const double EPS = std::numeric_limits<double>::epsilon();
}

inline void assert_r(bool cond, const std::string &msg) {
  if (!cond) Rcpp::stop("%s", msg.c_str());
}

inline std::vector<unsigned> count_values(const std::vector<int> &values,
                                          const std::vector<int> &sub_ids,
                                          int n_vals = 0) {
  if (n_vals == 0) {
    for (int i : sub_ids) {
      int v = values.at(i);
      if (v < 0) Rcpp::stop("sample_per_cell must contain only positive factors");
      n_vals = std::max(n_vals, v + 1);
    }
  }
  std::vector<unsigned> counts(n_vals, 0);
  for (int id : sub_ids) counts[values[id]]++;
  return counts;
}

inline double median(std::vector<double> &vec) {
  assert_r(!vec.empty(), "vector for median is empty");
  const auto median_it1 = vec.begin() + vec.size() / 2;
  std::nth_element(vec.begin(), median_it1, vec.end());
  if (vec.size() % 2 != 0) return *median_it1;
  const auto median_it2 = vec.begin() + vec.size() / 2 - 1;
  std::nth_element(vec.begin(), median_it2, vec.end());
  return (*median_it1 + *median_it2) / 2.0;
}

inline double mad(const std::vector<double> &vals, double med) {
  std::vector<double> diffs;
  diffs.reserve(vals.size());
  for (double v : vals) diffs.emplace_back(std::abs(v - med));
  return median(diffs) * 1.4826;
}

inline double var(const std::vector<double> &vals, double mean) {
  double res = 0.0;
  for (double v : vals) res += (v - mean) * (v - mean);
  return res / std::max<size_t>(1, vals.size() - 1);
}

inline Eigen::MatrixXd collapseMatrixNorm(const Eigen::SparseMatrix<double> &mtx,
                                          const std::vector<int> &factor,
                                          const std::vector<int> &nn_ids,
                                          const std::vector<unsigned> &n_obs_per_samp,
                                          int max_factor = 0) {
  assert_r(mtx.cols() == (int)factor.size(),
           "Number of columns in matrix (" + std::to_string(mtx.cols()) +
           ") must match the factor size (" + std::to_string(factor.size()) + ")");
  max_factor = std::max(max_factor + 1, int(n_obs_per_samp.size()));
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(mtx.rows(), max_factor);
  for (int id : nn_ids) {
    int fac = factor[id];
    if (fac >= (int)n_obs_per_samp.size() || fac < 0)
      Rcpp::stop("Wrong factor: %d, id: %d", fac, id);
    for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator gene_it(mtx, id); gene_it; ++gene_it) {
      res(gene_it.row(), fac) += gene_it.value() / double(n_obs_per_samp.at(fac));
    }
  }
  return res;
}

inline std::vector<double> applyMedianFilter(const std::vector<double> &signal,
                                             const std::vector<std::vector<int>> &nn_ids,
                                             const std::vector<size_t> &non_zero_ids) {
  std::vector<double> signal_smoothed(signal.size(), 0.0);
  for (size_t si : non_zero_ids) {
    if (std::isnan(signal.at(si))) { signal_smoothed[si] = NAN; continue; }
    std::vector<double> sig_cur;
    for (int nni : nn_ids.at(si)) {
      double val = signal.at(nni);
      if (!std::isnan(val)) sig_cur.emplace_back(val);
    }
    if (sig_cur.empty()) { signal_smoothed[si] = NAN; continue; }
    signal_smoothed[si] = median(sig_cur);
  }
  return signal_smoothed;
}

inline std::vector<double> applyMedianFilter(const std::vector<double> &signal,
                                             const std::vector<std::vector<int>> &nn_ids) {
  std::vector<size_t> non_zero_ids(signal.size());
  std::iota(non_zero_ids.begin(), non_zero_ids.end(), 0);
  return applyMedianFilter(signal, nn_ids, non_zero_ids);
}

inline std::pair<double, double> range(const std::vector<double> &vec) {
  double min_val = std::numeric_limits<double>::max(), max_val = std::numeric_limits<double>::lowest();
  bool all_nans = true;
  for (double v : vec) {
    if (std::isnan(v)) continue;
    all_nans = false;
    min_val = std::min(min_val, v);
    max_val = std::max(max_val, v);
  }
  if (all_nans) return std::make_pair(NAN, NAN);
  return std::make_pair(min_val, max_val);
}

inline std::pair<double, double> range(const std::vector<double> &vec, double wins) {
  assert_r(!vec.empty(), "vector for range is empty");
  if (wins < (2.0 / vec.size())) return range(vec);
  std::vector<double> vec_filt;
  for (double v : vec) if (!std::isnan(v)) vec_filt.emplace_back(v);
  if (vec_filt.empty()) return std::make_pair(NAN, NAN);
  const auto lq_it = vec_filt.begin() + size_t(std::floor(vec_filt.size() * wins));
  const auto uq_it = vec_filt.begin() + size_t(std::ceil(vec_filt.size() * (1 - wins)));
  std::nth_element(vec_filt.begin(), lq_it, vec_filt.end());
  std::nth_element(vec_filt.begin(), uq_it, vec_filt.end());
  return std::make_pair(*lq_it, *uq_it);
}

// ---- distances between sample vectors ----

// defined (and exported to R) in cluster_free.cpp
double estimateCorrelationDistance(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2, bool centered = true);

inline double cf_average(double val1, double val2) { return (val1 + val2) / 2.0; }

inline double estimateKLDivergence(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2) {
  if (v1.size() != v2.size()) Rcpp::stop("Vectors must have the same length");
  double res = 0.0;
  for (Eigen::Index i = 0; i < v1.size(); ++i) {
    const double d1 = v1[i], d2 = v2[i];
    if (std::isnan(d1) || std::isnan(d2)) return NAN;
    if (d1 > 1e-10 && d2 > 1e-10) res += std::log(d1 / d2) * d1;
  }
  return res;
}

inline double estimateJSDivergence(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2) {
  if (v1.size() != v2.size()) Rcpp::stop("Vectors must have the same length");
  Eigen::VectorXd avg = Eigen::VectorXd::Zero(v1.size());
  std::transform(v1.data(), v1.data() + v1.size(), v2.data(), avg.data(), cf_average);
  const double d1 = estimateKLDivergence(v1, avg);
  const double d2 = estimateKLDivergence(v2, avg);
  return std::sqrt(0.5 * (d1 + d2));
}

inline double estimateVectorDistance(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2, const std::string &dist) {
  if (dist == "cosine") return estimateCorrelationDistance(v1, v2, false);
  if (dist == "js")     return estimateJSDivergence(v1, v2);
  if (dist == "cor")    return estimateCorrelationDistance(v1, v2, true);
  Rcpp::stop("Unknown dist: %s", dist.c_str());
}

struct CFShiftResult {
  std::vector<double> dists;
  std::vector<size_t> s1_ids;
  std::vector<size_t> s2_ids;
};

// sample_per_cell must contain ids from 0..n_samples-1
inline CFShiftResult estimateCellExpressionShift(const Eigen::SparseMatrix<double> &cm,
                                                 const std::vector<int> &sample_per_cell,
                                                 const std::vector<int> &nn_ids,
                                                 size_t min_n_obs_per_samp,
                                                 const std::string &dist = "cosine",
                                                 bool log_vecs = true) {
  if (nn_ids.size() < min_n_obs_per_samp) return CFShiftResult();
  auto n_ids_per_samp = count_values(sample_per_cell, nn_ids);
  auto mat_collapsed  = collapseMatrixNorm(cm, sample_per_cell, nn_ids, n_ids_per_samp);
  if (log_vecs) {
    for (Eigen::Index i = 0; i < mat_collapsed.size(); ++i)
      mat_collapsed(i) = std::log10(1e3 * mat_collapsed(i) + 1.0);
  }
  std::vector<double> dists;
  std::vector<size_t> s1_ids, s2_ids;
  for (size_t s1 = 0; s1 < n_ids_per_samp.size(); ++s1) {
    if (n_ids_per_samp.at(s1) < min_n_obs_per_samp) continue;
    Eigen::VectorXd v1 = mat_collapsed.col((Eigen::Index)s1);
    for (size_t s2 = s1 + 1; s2 < n_ids_per_samp.size(); ++s2) {
      if (n_ids_per_samp.at(s2) < min_n_obs_per_samp) continue;
      Eigen::VectorXd v2 = mat_collapsed.col((Eigen::Index)s2);
      dists.push_back(estimateVectorDistance(v1, v2, dist));
      s1_ids.push_back(s1);
      s2_ids.push_back(s2);
    }
  }
  return CFShiftResult{dists, s1_ids, s2_ids};
}

// ---- permutation-based z-score adjustment ----

inline std::vector<double> adjustZScoresWithPermutations(const std::vector<double> &z_scores, double wins,
                                                         const std::vector<double> &max_vals) {
  std::vector<double> z_adj(z_scores.begin(), z_scores.end());
  auto max_val = range(z_adj, wins).second;
  for (double &z : z_adj) {
    if (std::isnan(z)) continue;
    z = std::min(z, max_val);
    size_t n = (max_vals.end() - std::lower_bound(max_vals.begin(), max_vals.end(), z - cf::EPS)); // #elements >= z
    z = std::max(1.0 - (n + 1.0) / (max_vals.size() + 1.0), 0.5);
  }
  for (double &p : z_adj) p = R::qnorm(p, 0.0, 1.0, 1, 0);
  return z_adj;
}

inline std::vector<double> adjustZScoresWithPermutations(const std::vector<double> &z_scores,
                                                         const std::vector<std::vector<int>> &nn_ids,
                                                         const std::vector<size_t> &non_zero_ids,
                                                         double wins, bool smooth,
                                                         const std::vector<double> &min_vals,
                                                         const std::vector<double> &max_vals, std::mutex &r_mut) {
  std::vector<double> z_adj(z_scores.begin(), z_scores.end());
  if (smooth) z_adj = applyMedianFilter(z_adj, nn_ids, non_zero_ids);
  auto rng = range(z_adj, wins);
  for (double &z : z_adj) {
    if (std::isnan(z)) continue;
    z = std::max(std::min(z, rng.second), rng.first);
    size_t n = (z < 0) ?
      (std::upper_bound(min_vals.begin(), min_vals.end(), z + cf::EPS) - min_vals.begin()) : // #elements <= z
      (max_vals.end() - std::lower_bound(max_vals.begin(), max_vals.end(), z - cf::EPS));    // #elements >= z
    z = std::max(1.0 - (n + 1.0) / (max_vals.size() + 1.0), 0.5);
  }
  {
    std::lock_guard<std::mutex> l(r_mut);
    for (double &p : z_adj) p = R::qnorm(p, 0.0, 1.0, 1, 0);
  }
  for (size_t i = 0; i < z_adj.size(); ++i) z_adj[i] = std::copysign(z_adj[i], z_scores[i]);
  return z_adj;
}

#endif

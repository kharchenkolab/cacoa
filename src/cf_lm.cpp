// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#ifndef EPS
#define EPS 1e-12
#endif

#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include <algorithm>
#include <vector>
#include <limits>
#include <numeric>
#include <cmath>
#include <mutex>
#include <unordered_map>
#include <random>
#include <cstdint>

// ---------- minimal helpers ----------

inline void assert_r(bool cond, const std::string& msg) {
  if (!cond) Rcpp::stop("%s", msg.c_str());
}

static inline uint64_t pack_pair32(uint32_t a, uint32_t b) {
  if (a > b) std::swap(a,b);
  return ( (uint64_t)a << 32 ) | (uint64_t)b;
}

std::vector<unsigned> count_values(const std::vector<int> &values,
                                   const std::vector<int> &sub_ids,
                                   int n_vals=0) {
  if (n_vals == 0) {
    for (int i : sub_ids) {
      int v = values.at(i);
      if (v < 0) Rcpp::stop("sample_per_cell must contain only positive factors");
      n_vals = std::max(n_vals, v + 1);
    }
  }

  std::vector<unsigned> counts(n_vals, 0);
  for (int id : sub_ids) {
    counts[values[id]]++;
  }

  return counts;
}

double median(std::vector<double> &vec) {
  assert_r(!vec.empty(), "vector for median is empty");
  const auto median_it1 = vec.begin() + vec.size() / 2;
  std::nth_element(vec.begin(), median_it1 , vec.end());

  if (vec.size() % 2 != 0)
    return *median_it1;

  const auto median_it2 = vec.begin() + vec.size() / 2 - 1;
  std::nth_element(vec.begin(), median_it2 , vec.end());
  return (*median_it1 + *median_it2) / 2.0;
}

double mad(const std::vector<double> &vals, double med) {
  std::vector<double> diffs;
  diffs.reserve(vals.size());
  for (double v : vals) diffs.emplace_back(std::abs(v - med));
  return median(diffs) * 1.4826;
}

double var(const std::vector<double> &vals, double mean) {
  double res = 0.0;
  for (double v : vals) res += (v - mean) * (v - mean);
  return res / std::max<size_t>(1, vals.size() - 1);
}

Eigen::MatrixXd collapseMatrixNorm(const Eigen::SparseMatrix<double> &mtx,
                                   const std::vector<int> &factor,
                                   const std::vector<int> &nn_ids,
                                   const std::vector<unsigned> &n_obs_per_samp,
                                   int max_factor=0) {
  assert_r(mtx.cols() == (int)factor.size(),
           "Number of columns in matrix must match the factor size");
  max_factor = std::max(max_factor + 1, int(n_obs_per_samp.size()));
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(mtx.rows(), max_factor);

  for (int id : nn_ids) {
    int fac = factor[id];
    if (fac >= (int)n_obs_per_samp.size() || fac < 0)
      Rcpp::stop("Wrong factor: %d, id: %d", fac, id);

    for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator gene_it(mtx, id);
         gene_it; ++gene_it) {
      res(gene_it.row(), fac) += gene_it.value() / double(n_obs_per_samp.at(fac));
    }
  }
  return res;
}

struct CFShiftResult{
  std::vector<double> dists;
  std::vector<size_t> s1_ids;
  std::vector<size_t> s2_ids;
};

// [[Rcpp::export]]
double estimateCorrelationDistance(const Eigen::VectorXd &v1,
                                   const Eigen::VectorXd &v2,
                                   bool centered=true) {
  if (v1.size() != v2.size())
    Rcpp::stop("Vectors must have the same length");

  double m1 = 0.0, m2 = 0.0;
  if (centered) {
    m1 = v1.mean();
    m2 = v2.mean();
  }

  double vp = 0.0, v1s = 0.0, v2s = 0.0;
  for (Eigen::Index i = 0; i < v1.size(); ++i) {
    if (std::isnan(v1[i]) || std::isnan(v2[i])) return NAN;
    const double e1 = (v1[i] - m1), e2 = (v2[i] - m2);
    vp  += e1 * e2;
    v1s += e1 * e1;
    v2s += e2 * e2;
  }
  const double denom = std::max(1e-10, std::sqrt(v1s) * std::sqrt(v2s));
  return 1.0 - vp / denom;
}

inline double average(double val1, double val2) { return (val1 + val2) / 2.0; }

double estimateKLDivergence(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2) {
  double res = 0.0;
  if (v1.size() != v2.size())
    Rcpp::stop("Vectors must have the same length");
  for (Eigen::Index i = 0; i < v1.size(); ++i) {
    const double d1 = v1[i], d2 = v2[i];
    if (std::isnan(d1) || std::isnan(d2)) return NAN;
    if (d1 > 1e-10 && d2 > 1e-10) res += std::log(d1 / d2) * d1;
  }
  return res;
}

double estimateJSDivergence(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2) {
  if (v1.size() != v2.size())
    Rcpp::stop("Vectors must have the same length");

  Eigen::VectorXd avg = Eigen::VectorXd::Zero(v1.size());
  std::transform(v1.data(), v1.data() + v1.size(), v2.data(), avg.data(), average);

  const double d1 = estimateKLDivergence(v1, avg);
  const double d2 = estimateKLDivergence(v2, avg);
  return std::sqrt(0.5 * (d1 + d2));
}

double estimateVectorDistance(const Eigen::VectorXd &v1,
                              const Eigen::VectorXd &v2,
                              const std::string &dist) {
  if (dist == "cosine") return estimateCorrelationDistance(v1, v2, false);
  if (dist == "js")     return estimateJSDivergence(v1, v2);
  if (dist == "cor")    return estimateCorrelationDistance(v1, v2, true);
  Rcpp::stop("Unknown dist: %s", dist.c_str());
}

// sample_per_cell must contain ids from 0..n_samples-1
CFShiftResult estimateCellExpressionShift(const Eigen::SparseMatrix<double> &cm,
                                          const std::vector<int> &sample_per_cell,
                                          const std::vector<int> &nn_ids,
                                          size_t min_n_obs_per_samp,
                                          const std::string &dist="cosine",
                                          bool log_vecs=true) {
  if (nn_ids.size() < min_n_obs_per_samp) return CFShiftResult();

  auto n_ids_per_samp = count_values(sample_per_cell, nn_ids);
  auto mat_collapsed  = collapseMatrixNorm(cm, sample_per_cell, nn_ids, n_ids_per_samp);
  if (log_vecs) {
    for (int i = 0; i < mat_collapsed.size(); ++i)
      mat_collapsed(i) = std::log10(1e3 * mat_collapsed(i) + 1.0);
  }

  std::vector<double> dists;
  std::vector<size_t> s1_ids, s2_ids;
  for (size_t s1 = 0; s1 < n_ids_per_samp.size(); ++s1) {
    if (n_ids_per_samp.at(s1) < min_n_obs_per_samp) continue;
    auto v1 = mat_collapsed.col((Eigen::Index)s1);
    for (size_t s2 = s1 + 1; s2 < n_ids_per_samp.size(); ++s2) {
      if (n_ids_per_samp.at(s2) < min_n_obs_per_samp) continue;
      auto v2 = mat_collapsed.col((Eigen::Index)s2);
      double d = estimateVectorDistance(v1, v2, dist);
      dists.push_back(d);
      s1_ids.push_back(s1);
      s2_ids.push_back(s2);
    }
  }
  return CFShiftResult{dists, s1_ids, s2_ids};
}

// ===================== main builder =====================

// Build Y (pairs × neighborhoods) of sample–sample distances computed
// within each neighborhood:
//  - Collapse cm[genes×cells] to sample-level means per neighborhood.
//  - Optional log10(1e3*x+1) transform.
//  - For each row in pairs_mat, compute distance
//    between the two sample vectors; NA if a sample has < min_n_obs_per_samp cells.

// [[Rcpp::export]]
arma::mat estimateExpressionShiftsPairsLM(
    const Eigen::SparseMatrix<double>& cm,            // genes x cells (sparse)
    Rcpp::IntegerVector sample_per_cell,              // length = n_cells, factor (1-based in R)
    Rcpp::List nn_ids,                                // list of int vectors of cell indices (often 1-based)
    const arma::imat& pairs_mat,                      // (n_pairs x 2) sample indices (often 1-based)
    int min_n_obs_per_samp = 1,
    std::string dist = "cor",                         // "cor", "cosine", or "js"
    bool log_vecs = true)
{
  if (cm.cols() != sample_per_cell.size())
    Rcpp::stop("cm.ncols (%d) must equal length(sample_per_cell) (%d).",
               cm.cols(), sample_per_cell.size());
  if (pairs_mat.n_cols != 2)
    Rcpp::stop("pairs_mat must have exactly 2 columns (sample indices).");
  if (!(dist == "cor" || dist == "cosine" || dist == "js"))
    Rcpp::stop("dist must be one of {'cor','cosine','js'}.");

  const int n_cells = sample_per_cell.size();

  // sample_of_cell: 0-based samples per cell
  std::vector<int> sample_of_cell(n_cells);
  int max_fac = -1;
  for (int i = 0; i < n_cells; ++i) {
    int v = sample_per_cell[i];
    if (Rcpp::IntegerVector::is_na(v))
      Rcpp::stop("sample_per_cell contains NA at position %d", i + 1);
    v -= 1;  // convert 1-based -> 0-based
    if (v < 0) Rcpp::stop("sample_per_cell must be a positive (1-based) factor.");
    sample_of_cell[i] = v;
    if (v > max_fac) max_fac = v;
  }
  const int n_samples = max_fac + 1;

  // nn_ids -> 0-based vectors
  const int n_neigh = nn_ids.size();
  std::vector<std::vector<int>> nn_list;
  nn_list.reserve(n_neigh);
  for (int k = 0; k < n_neigh; ++k) {
    Rcpp::IntegerVector ids_r = nn_ids[k];
    std::vector<int> ids; ids.reserve(ids_r.size());
    int min_id = std::numeric_limits<int>::max();
    int max_id = std::numeric_limits<int>::min();
    for (int t = 0; t < ids_r.size(); ++t) {
      int idx = ids_r[t];
      if (Rcpp::IntegerVector::is_na(idx)) continue;
      min_id = std::min(min_id, idx);
      max_id = std::max(max_id, idx);
      ids.push_back(idx);
    }
    if (!ids.empty() && (min_id >= 1 || max_id == cm.cols())) {
      for (int& v : ids) --v; // 1-based -> 0-based
    }
    for (int v : ids) {
      if (v < 0 || v >= cm.cols())
        Rcpp::stop("nn_ids[[%d]] has out-of-range index.", k + 1);
    }
    nn_list.push_back(std::move(ids));
  }

  // pairs_mat -> 0-based
  arma::imat pairs = pairs_mat;
  if (pairs.n_rows > 0) {
    int pmin = pairs.min();
    int pmax = pairs.max();
    if (pmin >= 1 || pmax >= n_samples) pairs -= 1; // looks 1-based
    if (pairs.min() < 0 || pairs.max() >= n_samples)
      Rcpp::stop("pairs_mat has indices outside [0, %d) after normalization.", n_samples);
  }

  // Output: rows = pairs, cols = neighborhoods
  arma::mat Y(pairs.n_rows, n_neigh);
  Y.fill(std::numeric_limits<double>::quiet_NaN());

  for (int k = 0; k < n_neigh; ++k) {
    const auto& ids = nn_list[k];
    if ((int)ids.size() < min_n_obs_per_samp) continue;

    std::vector<unsigned> n_obs_per_samp = count_values(sample_of_cell, ids, n_samples);
    Eigen::MatrixXd collapsed = collapseMatrixNorm(cm, sample_of_cell, ids, n_obs_per_samp);

    if (log_vecs) {
      for (int i = 0; i < collapsed.size(); ++i)
        collapsed(i) = std::log10(1e3 * collapsed(i) + 1.0);
    }

    for (arma::uword r = 0; r < pairs.n_rows; ++r) {
      int a = pairs(r, 0), b = pairs(r, 1);
      if (a < 0 || b < 0 || a >= n_samples || b >= n_samples) continue;
      if ((int)n_obs_per_samp[a] < min_n_obs_per_samp ||
          (int)n_obs_per_samp[b] < min_n_obs_per_samp) continue;

      const Eigen::VectorXd v1 = collapsed.col((Eigen::Index)a);
      const Eigen::VectorXd v2 = collapsed.col((Eigen::Index)b);
      Y(r, k) = estimateVectorDistance(v1, v2, dist);
    }
  }

  return Y;
}


// ===================== Z-score matrix =====================
// [[Rcpp::export]]
std::vector<double> applyMedianFilter(const std::vector<double> &signal, const std::vector<std::vector<int>> &nn_ids, const std::vector<size_t>& non_zero_ids) {
    std::vector<double> signal_smoothed(signal.size(), 0.0);
    for (size_t si : non_zero_ids) {
        if (std::isnan(signal.at(si))) {
            signal_smoothed[si] = NAN;
            continue;
        }

        std::vector<double> sig_cur;
        for (int nni : nn_ids.at(si)) {
            double val = signal.at(nni);
            if (!std::isnan(val)) {
                sig_cur.emplace_back(val);
            }
        }

        if (sig_cur.empty()) {
            signal_smoothed[si] = NAN;
            continue;
        }

        signal_smoothed[si] = median(sig_cur);
    }
    return signal_smoothed;
}

std::vector<double> applyMedianFilter(const std::vector<double> &signal, const std::vector<std::vector<int>> &nn_ids) {
    std::vector<size_t> non_zero_ids(signal.size());
    std::iota(non_zero_ids.begin(), non_zero_ids.end(), 0);
    return applyMedianFilter(signal, nn_ids, non_zero_ids);
}

std::pair<double, double> range(const std::vector<double> &vec) {
    double min_val = std::numeric_limits<double>::max(), max_val = std::numeric_limits<double>::lowest();
    bool all_nans = true;
    for (double v : vec) {
        if (std::isnan(v))
            continue;

        all_nans = false;
        min_val = std::min(min_val, v);
        max_val = std::max(max_val, v);
    }

    if (all_nans)
        return std::make_pair(NAN, NAN);

    return std::make_pair(min_val, max_val);
}

std::pair<double, double> range(const std::vector<double> &vec, double wins) {
    assert_r(!vec.empty(), "vector for range is empty");

    if (wins < (2.0 / vec.size()))
        return range(vec);

    std::vector<double> vec_filt;
    for (double v : vec) {
        if (!std::isnan(v)) {
            vec_filt.emplace_back(v);
        }
    }

    if (vec_filt.empty())
        return std::make_pair(NAN, NAN);

    const auto lq_it = vec_filt.begin() + size_t(std::floor(vec_filt.size() * wins));
    const auto uq_it = vec_filt.begin() + size_t(std::ceil(vec_filt.size() * (1 - wins)));
    std::nth_element(vec_filt.begin(), lq_it, vec_filt.end());
    std::nth_element(vec_filt.begin(), uq_it, vec_filt.end());
    return std::make_pair(*lq_it, *uq_it);
}

// [[Rcpp::export]]
std::vector<double> applyMedianFilterES(
    const std::vector<double>& x,                 // values to smooth (e.g., res$stat.obs)
    Rcpp::List nn_ids,                            // list of integer vectors (neighbors per node)
    Rcpp::Nullable<Rcpp::IntegerVector> non_zero_ids = R_NilValue,
    bool one_based = false                         // set FALSE if neighbors are already 0-based
) {
    const size_t m = x.size();
    if ((size_t)nn_ids.size() != m)
        Rcpp::stop("nn_ids length (%d) must equal length(x) (%d).", (int)nn_ids.size(), (int)m);

    // Build neighbor list (0-based, in-range, not empty)
    std::vector<std::vector<int>> nn_cpp(m);
    for (size_t i = 0; i < m; ++i) {
        Rcpp::IntegerVector v = nn_ids[i];
        nn_cpp[i].reserve(v.size());
        for (int a : v) {
            if (a == NA_INTEGER) continue;
            int idx = one_based ? (a - 1) : a;
            if (idx >= 0 && idx < (int)m) nn_cpp[i].push_back(idx);
        }
        if (nn_cpp[i].empty()) nn_cpp[i].push_back((int)i); // fallback to self
    }

    // Build non-zero mask indices (0-based)
    std::vector<size_t> nz_cpp;
    if (non_zero_ids.isNotNull()) {
        Rcpp::IntegerVector iv(non_zero_ids.get());
        nz_cpp.reserve(iv.size());
        for (int a : iv) {
            if (a == NA_INTEGER) continue;
            int idx = one_based ? (a - 1) : a;
            if (idx >= 0 && idx < (int)m) nz_cpp.push_back((size_t)idx);
        }
    } else {
        nz_cpp.resize(m);
        std::iota(nz_cpp.begin(), nz_cpp.end(), 0);
    }

    return applyMedianFilter(x, nn_cpp, nz_cpp);
}

static inline void tally_perm(double obs, double perm, int alt, int& ge, int& le, int& ge_abs) {
  if (alt == 0) { if (std::fabs(perm) >= std::fabs(obs)) ge_abs++; }
  else if (alt == 1) { if (perm >= obs) ge++; }
  else { if (perm <= obs) le++; }
}

static inline double z_from_p(double p, int alt, double obs, double med_perm) {
  // Handle boundary / invalid p explicitly
  if (!std::isfinite(p)) return NA_REAL;
  
  // If add-one p-value hits the upper boundary (no evidence in the tested tail),
  // return 0 instead of NA/Inf.
  if (p >= 1.0) return 0.0;
  
  // With add-one, p should never be 0, but guard anyway:
  if (p <= 0.0) {
    // Map to a very large finite z with the correct direction.
    // (Users almost never hit this with add-one.)
    const double p_eff = 1e-16;
    if (alt == 0) {
      double zmag = R::qnorm(1.0 - p_eff/2.0, 0.0, 1.0, 1, 0);
      double sgn  = (obs >= med_perm) ? 1.0 : -1.0;
      return sgn * zmag;
    } else if (alt == 1) {
      return R::qnorm(1.0 - p_eff, 0.0, 1.0, 1, 0);
    } else {
      return -R::qnorm(1.0 - p_eff, 0.0, 1.0, 1, 0);
    }
  }
  
  // Regular case: 0 < p < 1
  if (alt == 0) { // two-sided: sign by location vs permutation median
    double zmag = R::qnorm(1.0 - p/2.0, 0.0, 1.0, 1, 0);
    double sgn  = (obs >= med_perm) ? 1.0 : -1.0;
    return sgn * zmag;
  } else if (alt == 1) { // greater: upper-tail
    return R::qnorm(1.0 - p, 0.0, 1.0, 1, 0);
  } else {               // less: lower-tail
    return -R::qnorm(1.0 - p, 0.0, 1.0, 1, 0);
  }
}


std::vector<double> adjustZScoresWithPermutations(const std::vector<double> &z_scores, double wins, const std::vector<double> &max_vals) {
    std::vector<double> z_adj(z_scores.begin(), z_scores.end());

    auto max_val = range(z_adj, wins).second;
    for (double &z : z_adj) {
        if (std::isnan(z))
            continue;

        z = std::min(z, max_val);
        size_t n = (max_vals.end() - std::lower_bound(max_vals.begin(), max_vals.end(), z - EPS)); // Number of elements that are >=z
        z = std::max(1.0 - (n + 1.0) / (max_vals.size() + 1.0), 0.5);
    }

    for (double &p : z_adj) {
      // p currently holds an adjusted *probability* in [0.5, 1).
      p = R::qnorm(p, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
    }
    return z_adj;
}

std::vector<double> adjustZScoresWithPermutations(const std::vector<double> &z_scores, const std::vector<std::vector<int>> &nn_ids,
                                                  const std::vector<size_t>& non_zero_ids,
                                                  double wins, bool smooth, const std::vector<double> &min_vals,
                                                  const std::vector<double> &max_vals, std::mutex &r_mut) {
    std::vector<double> z_adj(z_scores.begin(), z_scores.end());
    if (smooth) {
        z_adj = applyMedianFilter(z_adj, nn_ids, non_zero_ids);
    }

    auto rng = range(z_adj, wins);
    for (double &z : z_adj) {
        if (std::isnan(z))
            continue;

        z = std::max(std::min(z, rng.second), rng.first);
        size_t n = (z < 0) ?
                   (std::upper_bound(min_vals.begin(), min_vals.end(), z + EPS) - min_vals.begin()) : // Number of elements that are <=z
                   (max_vals.end() - std::lower_bound(max_vals.begin(), max_vals.end(), z - EPS)); // Number of elements that are >=z
        z = std::max(1.0 - (n + 1.0) / (max_vals.size() + 1.0), 0.5);
    }

    {
      std::lock_guard<std::mutex> l(r_mut);
      for (double &p : z_adj) {
          p = R::qnorm(p, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
       }
    }

    for (size_t i = 0; i < z_adj.size(); ++i) {
        z_adj[i] = std::copysign(z_adj[i], z_scores[i]);
    }

    return z_adj;
}

// [[Rcpp::export]]
std::vector<double> adjustedZScoresMaxStat(
    const std::vector<double>& T_obs,                 // length m (tests)
    const Rcpp::NumericMatrix& T_perm_mat,            // B x m (rows=perms, cols=tests)
    int alt,                                          // 0=two-sided, 1=greater, 2=less
    double wins = 0.0,
    bool smooth = false,
    Rcpp::Nullable<Rcpp::List> nn_ids = R_NilValue,                // only used if smooth==true
    Rcpp::Nullable<Rcpp::IntegerVector> non_zero_ids = R_NilValue  // only used if smooth==true
) {
    // ---- dims & checks ----
    const size_t B = static_cast<size_t>(T_perm_mat.nrow());
    const size_t m = static_cast<size_t>(T_perm_mat.ncol());
    if (m == 0 || B == 0) Rcpp::stop("T_perm_mat must be non-empty (rows=perms, cols=tests).");
    if (T_obs.size() != m) Rcpp::stop("Length(T_obs) must equal ncol(T_perm_mat).");
    if (!(alt == 0 || alt == 1 || alt == 2)) Rcpp::stop("alt must be 0, 1, or 2.");

    // ---- optional inputs (sanitize when smooth==true) ----
    std::vector<std::vector<int>> nn_cpp;   // will become length m
    std::vector<size_t> nz_cpp;

    if (smooth) {
        // Build/clean neighbor lists
        nn_cpp.assign(m, std::vector<int>{});  // default: empty per test
        if (nn_ids.isNotNull()) {
            Rcpp::List L(nn_ids.get());
            const int Lsz = L.size();
            for (int i = 0; i < Lsz; ++i) {
                const size_t idx = static_cast<size_t>(i);
                if (idx >= m) break; // ignore extra entries
                Rcpp::IntegerVector v = L[i];
                nn_cpp[idx].reserve(v.size());
                for (int a : v) {
                    if (a == NA_INTEGER) continue;
                    // assume incoming are already 0-based; if yours are 1-based, do: a -= 1;
                    if (a >= 0 && a < static_cast<int>(m)) nn_cpp[idx].push_back(a);
                }
                // fallback: ensure at least self as neighbor
                if (nn_cpp[idx].empty()) nn_cpp[idx].push_back(static_cast<int>(idx));
            }
        }
        // If list shorter than m, fill missing entries with self-neighbors
        for (size_t i = 0; i < m; ++i) {
            if (i >= nn_cpp.size()) nn_cpp.push_back(std::vector<int>{static_cast<int>(i)});
            if (nn_cpp[i].empty())  nn_cpp[i].push_back(static_cast<int>(i));
        }

        // Filter/convert non_zero_ids
        if (non_zero_ids.isNotNull()) {
            Rcpp::IntegerVector iv(non_zero_ids.get());
            nz_cpp.reserve(iv.size());
            for (int x : iv) {
                if (x == NA_INTEGER) continue;
                if (x >= 0 && x < static_cast<int>(m)) nz_cpp.push_back(static_cast<size_t>(x));
            }
        } else {
            // default: all tests are eligible
            nz_cpp.resize(m);
            std::iota(nz_cpp.begin(), nz_cpp.end(), 0);
        }
    }

    // ---- Step 1: observed p & z per test ----
    std::vector<double> z_obs(m), med_perm(m);

    for (size_t j = 0; j < m; ++j) {
        // median() mutates its input -> copy column j
        std::vector<double> col(B);
        for (size_t b = 0; b < B; ++b) col[b] = T_perm_mat(b, j);
        med_perm[j] = median(col);

        int ge = 0, le = 0, ge_abs = 0;
        const double tobs = T_obs[j];
        for (size_t b = 0; b < B; ++b) {
            tally_perm(tobs, T_perm_mat(b, j), alt, ge, le, ge_abs);
        }

        double p_obs;
        if (alt == 1)       p_obs = (ge     + 1.0) / (B + 1.0);
        else if (alt == 2)  p_obs = (le     + 1.0) / (B + 1.0);
        else                p_obs = (ge_abs + 1.0) / (B + 1.0);

        z_obs[j] = z_from_p(p_obs, alt, tobs, med_perm[j]);
    }

    // ---- Step 2: permutation z's for each test ----
    std::vector<std::vector<double>> z_perm(m, std::vector<double>(B, 0.0));

    for (size_t j = 0; j < m; ++j) {
        std::vector<double> sorted_vals(B), sorted_abs_vals(B);
        for (size_t b = 0; b < B; ++b) {
            const double v = T_perm_mat(b, j);
            sorted_vals[b]     = v;
            sorted_abs_vals[b] = std::fabs(v);
        }
        std::sort(sorted_vals.begin(), sorted_vals.end());
        std::sort(sorted_abs_vals.begin(), sorted_abs_vals.end());

        for (size_t b = 0; b < B; ++b) {
            const double x = T_perm_mat(b, j);
            double p_b;
            if (alt == 1) {
                auto it = std::lower_bound(sorted_vals.begin(), sorted_vals.end(), x - EPS);
                size_t n_ge = static_cast<size_t>(sorted_vals.end() - it);
                p_b = (n_ge + 1.0) / (B + 1.0);
            } else if (alt == 2) {
                auto it = std::upper_bound(sorted_vals.begin(), sorted_vals.end(), x + EPS);
                size_t n_le = static_cast<size_t>(it - sorted_vals.begin());
                p_b = (n_le + 1.0) / (B + 1.0);
            } else {
                const double ax = std::fabs(x);
                auto it = std::lower_bound(sorted_abs_vals.begin(), sorted_abs_vals.end(), ax - EPS);
                size_t n_ge_abs = static_cast<size_t>(sorted_abs_vals.end() - it);
                p_b = (n_ge_abs + 1.0) / (B + 1.0);
            }
            z_perm[j][b] = z_from_p(p_b, alt, /*obs=*/x, /*med_perm=*/med_perm[j]);
        }
    }

    // ---- Step 3: per-permutation extremes across tests ----
    std::vector<double> max_vals(B, -std::numeric_limits<double>::infinity());
    std::vector<double> min_vals(B,  std::numeric_limits<double>::infinity());
    for (size_t b = 0; b < B; ++b) {
        double mx = -std::numeric_limits<double>::infinity();
        double mn =  std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < m; ++j) {
            const double z = z_perm[j][b];
            if (z > mx) mx = z;
            if (z < mn) mn = z;
        }
        max_vals[b] = mx;
        min_vals[b] = mn;
    }
    std::sort(max_vals.begin(), max_vals.end());
    std::sort(min_vals.begin(), min_vals.end());

    // ---- Step 4: call your adjustors to get adjusted z's ----
    std::vector<double> z_adj;
    if (alt == 1) {
        z_adj = adjustZScoresWithPermutations(z_obs, wins, max_vals);
    } else {
        std::mutex r_mut;
        z_adj = adjustZScoresWithPermutations(
            z_obs,
            nn_cpp, nz_cpp,
            wins, smooth,
            min_vals, max_vals,
            r_mut
        );
    }
    return z_adj;
}
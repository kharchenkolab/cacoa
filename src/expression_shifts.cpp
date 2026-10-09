// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "cf_common.h"

#include <algorithm>
#include <vector>
#include <limits>
#include <numeric>
#include <cmath>
#include <mutex>
#include <unordered_map>
#include <random>
#include <cstdint>

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
    bool log_vecs = true,
    bool nn_one_based = false,                        // nn_ids are 0-based (graph adjacency) unless TRUE
    bool pairs_one_based = true)                      // pairs_mat as produced by combn() in R
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
    for (int t = 0; t < ids_r.size(); ++t) {
      int idx = ids_r[t];
      if (Rcpp::IntegerVector::is_na(idx)) continue;
      ids.push_back(nn_one_based ? idx - 1 : idx);
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
    if (pairs_one_based) pairs -= 1;
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

// [[Rcpp::export]]
std::vector<double> adjustedZScoresMaxStat(
    const std::vector<double>& z_obs,                 // observed test stat (or z), length m
    int alt,                                          // 0=two-sided, 1=greater, 2=less
    const Rcpp::NumericVector& max_vals_in,           // per-permutation max stats, length B; should be on the same scale as z_obs
    const Rcpp::NumericVector& min_vals_in,           // per-permutation min stats, length B; should be on the same scale as z_obs
    double wins = 0.0,
    bool smooth = false,
    Rcpp::Nullable<Rcpp::List> nn_ids = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> non_zero_ids = R_NilValue
) {
    const size_t m = z_obs.size();
    const size_t B = max_vals_in.size();

    if (!(alt == 0 || alt == 1 || alt == 2))
        Rcpp::stop("alt must be 0, 1, or 2.");

    if (min_vals_in.size() != B)
        Rcpp::stop("max_vals and min_vals must have equal lengths.");

    // copy extremes
    std::vector<double> max_vals(max_vals_in.begin(), max_vals_in.end());
    std::vector<double> min_vals(min_vals_in.begin(), min_vals_in.end());

    std::sort(max_vals.begin(), max_vals.end());
    std::sort(min_vals.begin(), min_vals.end());

    // ---- optional neighbor metadata (unchanged) ----
    std::vector<std::vector<int>> nn_cpp;
    std::vector<size_t> nz_cpp;

    if (smooth) {
        nn_cpp.assign(m, std::vector<int>{});

        if (nn_ids.isNotNull()) {
            Rcpp::List L(nn_ids.get());
            const int Lsz = L.size();
            for (int i = 0; i < Lsz; ++i) {
                size_t idx = static_cast<size_t>(i);
                if (idx >= m) break;
                Rcpp::IntegerVector v = L[i];
                for (int a : v) {
                    if (a >= 0 && a < (int)m)
                        nn_cpp[idx].push_back(a);
                }
                if (nn_cpp[idx].empty())
                    nn_cpp[idx].push_back(static_cast<int>(idx));
            }
        }

        for (size_t i = 0; i < m; ++i)
            if (nn_cpp[i].empty())
                nn_cpp[i].push_back(static_cast<int>(i));

        if (non_zero_ids.isNotNull()) {
            Rcpp::IntegerVector iv(non_zero_ids.get());
            for (int x : iv)
                if (x >= 0 && x < (int)m)
                    nz_cpp.push_back((size_t)x);
        }
        if (nz_cpp.empty()) {
            nz_cpp.resize(m);
            std::iota(nz_cpp.begin(), nz_cpp.end(), 0);
        }
    }

    // ---- adjustment ----
    std::vector<double> z_adj;

    if (alt == 1) {
        // one-sided: use max_vals only
        z_adj = adjustZScoresWithPermutations(z_obs, wins, max_vals);
    } else {
        // two-sided or less: use min_vals + max_vals
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

/* P-value processing */

/* PCA */


#include <RcppEigen.h>
#include <iostream>
#include <vector>
#include <cmath>

#include <progress.hpp>

using namespace Rcpp;
using namespace Eigen;

double incrementMean(double cur_mean, double x, size_t n) {
  return cur_mean + (x - cur_mean) / n;
}

double incrementStdAcc(double cur_std, double x, double cur_mean, double old_mean) {
  return cur_std + (x - old_mean) * (x - cur_mean);
}

//' @param adj_mat adjacency matrix with 1 on the position (r,c) if the cell r is adjacent to the cell c
SparseMatrix<double> clusterFreeZScoreMat(const SparseMatrix<bool>& adj_mat, const SparseMatrix<double>& count_mat,
                                          const std::vector<bool> &is_control, bool verbose=true, double min_z=0.01) {
  if (adj_mat.cols() != adj_mat.rows())
    Rcpp::stop("adj_mat must be squared");

  if (count_mat.rows() != adj_mat.rows())
    Rcpp::stop("adj_mat must have the same number of rows as count_mat");

  if (is_control.size() != adj_mat.rows())
    Rcpp::stop("is_control must size equal to the number of rows in count_mat");

  std::vector<Triplet<double>> z_triplets, lfc_triplets;
  Progress p(count_mat.outerSize(), verbose);
  for (size_t gene_id=0; gene_id < count_mat.outerSize(); ++gene_id) {
    auto cur_col = VectorXd(count_mat.col(gene_id));
    for (SparseMatrix<double, ColMajor>::InnerIterator cnt_it(count_mat, gene_id); cnt_it; ++cnt_it) {
      if (Progress::check_abort())
        Rcpp::stop("Aborted");

      if (cnt_it.value() < 1e-20)
        continue;

      auto dst_cell_id = cnt_it.row();
      double mean_val_control = 0.0, mean_val_case = 0.0;
      double std_acc_control = 0.0, std_acc_case = 0.0;
      size_t n_vals_control = 0, n_vals_case = 0;

      // Incremental mean and std (http://datagenetics.com/blog/november22017/index.html)
      for (SparseMatrix<bool, ColMajor>::InnerIterator adj_it(adj_mat, dst_cell_id); adj_it; ++adj_it) {
        if (!adj_it.value())
          continue;

        double cur_val = cur_col(adj_it.row());
        if (is_control.at(adj_it.row())) {
          n_vals_control++;
          double mean_old = mean_val_control;
          mean_val_control = incrementMean(mean_val_control, cur_val, n_vals_control);
          std_acc_control = incrementStdAcc(std_acc_control, cur_val, mean_val_control, mean_old);
        } else {
          n_vals_case++;
          mean_val_case = incrementMean(mean_val_case, cur_val, n_vals_case);
        }
      }

      if ((n_vals_case > 0) && (n_vals_control > 1)) {
        double std_val = std::sqrt(std_acc_control / (n_vals_control - 1));
        double z_score = (mean_val_case - mean_val_control) / std::max(std_val, 1e-20);
        if (std::abs(z_score) > min_z) {
          z_triplets.emplace_back(dst_cell_id, gene_id, z_score);
        }
      }
      else {
        z_triplets.emplace_back(dst_cell_id, gene_id, NA_REAL);
      }
    }

    p.increment();
  }

  SparseMatrix<double> z_mat(count_mat.rows(), count_mat.cols());
  z_mat.setFromTriplets(z_triplets.begin(), z_triplets.end());
  return z_mat;
}

// [[Rcpp::export]]
SEXP clusterFreeZScoreMat(const SEXP adj_mat, const SEXP count_mat, const std::vector<bool> &is_control, bool verbose=true, double min_z=0.01) {
  auto z_mat_eig = clusterFreeZScoreMat(as<SparseMatrix<bool>>(adj_mat), as<SparseMatrix<double>>(count_mat), is_control, verbose, min_z);

  S4 z_mat(wrap(z_mat_eig));
  z_mat.slot("Dimnames") = S4(count_mat).slot("Dimnames");

  return z_mat;
}

////// Expression shifts

double estimateCosineDistance(const std::vector<double> &v1, const std::vector<double> &v2) {
  if (v1.size() != v2.size())
    Rcpp::stop("Vectors must have the same length");

  double vp = 0, v1s = 0, v2s = 0;
  for (size_t i = 0; i < v1.size(); ++i) {
    if (R_IsNA(v1[i]) || R_IsNA(v2[i]))
      return NA_REAL;

    vp += v1[i] * v2[i];
    v1s += v1[i] * v1[i];
    v2s += v2[i] * v2[i];
  }

  return 1 - vp / std::max(std::sqrt(v1s) * std::sqrt(v2s), 1e-10);
}

NumericMatrix collapseMatrixNorm(const SparseMatrix<double> &mat, const std::vector<int> &factor,
                                 const std::vector<int> &nn_ids, const std::vector<int> &n_obs_per_samp) {
  NumericMatrix res(mat.rows(), n_obs_per_samp.size());
  for (int id : nn_ids) {
      int fac = factor[id] - 1;
      if (fac >= res.cols() || fac < 0)
          Rcpp::stop("Wrong factor: " + std::to_string(fac) + ", id: " + std::to_string(id));

      for (SparseMatrix<double, ColMajor>::InnerIterator gene_it(mat, id); gene_it; ++gene_it) {
          res(gene_it.row(), fac) += gene_it.value() / n_obs_per_samp[fac];
      }
  }

  return res;
}

double estimateClusterFreeExpressionShift(const SparseMatrix<double> &mat, const std::vector<int> &factor, const std::vector<int> &nn_ids,
                                          const std::vector<bool> &is_ref, int n_samples) {
  std::vector<int> n_ids_per_samp(n_samples, 0);
  for (int id : nn_ids) {
    n_ids_per_samp[factor[id] - 1]++;
  }

  const auto mat_collapsed = collapseMatrixNorm(mat, factor, nn_ids, n_ids_per_samp);

  double within_dist = 0, between_dist = 0;
  int n_within = 0, n_between = 0;
  for (int s1 = 0; s1 < n_samples; ++s1) {
    if (n_ids_per_samp[s1] == 0)
      continue;

    auto v1 = as<std::vector<double>>(NumericVector(mat_collapsed(_, s1)));
    for (int s2 = s1 + 1; s2 < n_samples; ++s2) {
      if (n_ids_per_samp[s2] == 0)
        continue;

      auto v2 = as<std::vector<double>>(NumericVector(mat_collapsed(_, s2)));
      double dist = estimateCosineDistance(v1, v2);
      if (is_ref[s1] && is_ref[s2]) {
        within_dist = incrementMean(within_dist, dist, ++n_within);
      } else if (is_ref[s1] != is_ref[s2]) {
        between_dist = incrementMean(between_dist, dist, ++n_between);
      }
    }
  }

  if (n_within < 2 || n_between < 2)
    return NA_REAL;

  return between_dist / within_dist;
}

// [[Rcpp::export]]
NumericVector estimateClusterFreeExpressionShifts(SEXP rmat, IntegerVector rfactor, List nn_ids, const std::vector<bool> &is_ref, bool verbose=true) {
  auto mat = as<SparseMatrix<double>>(rmat);
  auto factor = as<std::vector<int>>(rfactor);
  std::vector<double> res_scores;
  Progress p(LENGTH(nn_ids), verbose);
  int n_samples = LENGTH(rfactor.attr("levels"));

  for (auto &ids : nn_ids) {
    if (Progress::check_abort())
      Rcpp::stop("Aborted");

    double score = estimateClusterFreeExpressionShift(mat, factor, as<std::vector<int>>(ids), is_ref, n_samples);
    res_scores.push_back(score);
    p.increment();
  }

  NumericVector res = wrap(res_scores);
  res.attr("names") = nn_ids.names();


  return res;
}

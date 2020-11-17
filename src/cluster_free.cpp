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
                                          const std::vector<bool> is_control, bool verbose=true, double min_z=0.01) {
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
SEXP clusterFreeZScoreMat(const SEXP adj_mat, const SEXP count_mat, const std::vector<bool> is_control, bool verbose=true, double min_z=0.01) {
  auto z_mat_eig = clusterFreeZScoreMat(as<SparseMatrix<bool>>(adj_mat), as<SparseMatrix<double>>(count_mat), is_control, verbose, min_z);

  S4 z_mat(wrap(z_mat_eig));
  z_mat.slot("Dimnames") = S4(count_mat).slot("Dimnames");

  return z_mat;
}

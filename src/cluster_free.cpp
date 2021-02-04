#include <RcppEigen.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <mutex>
#include <functional>
#include <sccore_par.hpp>

#include <progress.hpp>
#include <unistd.h>

using namespace Rcpp;
using namespace Eigen;

double incrementMean(double cur_mean, double x, size_t n) {
  return cur_mean + (x - cur_mean) / n;
}

double incrementStdAcc(double cur_std, double x, double cur_mean, double old_mean) {
  return cur_std + (x - old_mean) * (x - cur_mean);
}

double median(std::vector<double> &vec) {
  assert(!vec.empty());
  const auto median_it1 = vec.begin() + vec.size() / 2;
  std::nth_element(vec.begin(), median_it1 , vec.end());

  if (vec.size() % 2 != 0)
    return *median_it1;

  const auto median_it2 = vec.begin() + vec.size() / 2 - 1;
  std::nth_element(vec.begin(), median_it2 , vec.end());
  return (*median_it1 + *median_it2) / 2;
}

double estimateCellZScore(const SparseMatrix<bool> &adj_mat, const std::vector<bool> &is_control,
                          const VectorXd &cur_col, size_t dst_cell_id, bool normalize_both=false) {
    double mean_val_control = 0.0, mean_val_case = 0.0, mean_val_all = 0.0, std_acc = 0.0;
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

            if (!normalize_both) {
              std_acc = incrementStdAcc(std_acc, cur_val, mean_val_control, mean_old);
            }
        } else {
            n_vals_case++;
            mean_val_case = incrementMean(mean_val_case, cur_val, n_vals_case);
        }

        if (normalize_both) {
          double mean_old = mean_val_all;
          mean_val_all = incrementMean(mean_val_all, cur_val, n_vals_control + n_vals_case);
          std_acc = incrementStdAcc(std_acc, cur_val, mean_val_all, mean_old);
        }
    }

    if ((n_vals_case == 0) || (n_vals_control <= 1))
        return NAN;

    size_t n_vals = normalize_both ? (n_vals_control + n_vals_case) : n_vals_control;
    double std_val = sqrt(std_acc / (n_vals - 1));
    return (mean_val_case - mean_val_control) / std::max(std_val, 1e-20);
}

//' @param adj_mat adjacency matrix with 1 on the position (r,c) if the cell r is adjacent to the cell c
SparseMatrix<double> clusterFreeZScoreMat(const SparseMatrix<bool>& adj_mat, const SparseMatrix<double>& count_mat,
                                          const std::vector<bool> &is_control, const double min_z=0.01,
                                          bool normalize_both=false, bool verbose=true, int n_cores=1) {
  if (adj_mat.cols() != adj_mat.rows())
    stop("adj_mat must be squared");

  if (count_mat.rows() != adj_mat.rows())
    stop("adj_mat must have the same number of rows as count_mat");

  if (is_control.size() != adj_mat.rows())
    stop("is_control must size equal to the number of rows in count_mat");

  std::vector<Triplet<double>> z_triplets;

  std::mutex mut;
  auto task = [&count_mat, &z_triplets, &adj_mat, &is_control, &min_z, &normalize_both, &mut](int gene_id) {
    auto cur_col = VectorXd(count_mat.col(gene_id));
    for (SparseMatrix<double, ColMajor>::InnerIterator cnt_it(count_mat, gene_id); cnt_it; ++cnt_it) {
      if (cnt_it.value() < 1e-20)
        continue;

      auto dst_cell_id = cnt_it.row();
      auto z_cur = estimateCellZScore(adj_mat, is_control, cur_col, dst_cell_id, normalize_both);
      if (std::isnan(z_cur) || std::abs(z_cur) >= min_z) {
        std::lock_guard<std::mutex> l(mut);
        z_triplets.emplace_back(dst_cell_id, gene_id, z_cur);
      }
    }
  };

  sccore::runTaskParallelFor(0, count_mat.cols(), task, n_cores, verbose);
  SparseMatrix<double> z_mat(count_mat.rows(), count_mat.cols());
  z_mat.setFromTriplets(z_triplets.begin(), z_triplets.end());
  return z_mat;
}

// [[Rcpp::export]]
SEXP clusterFreeZScoreMat(const SEXP adj_mat, const SEXP count_mat, const std::vector<bool> &is_control, double min_z=0.01,
                          bool normalize_both=false, bool verbose=true, int n_cores=1) {
  auto z_mat_eig = clusterFreeZScoreMat(as<SparseMatrix<bool>>(adj_mat), as<SparseMatrix<double>>(count_mat),
          is_control, min_z, normalize_both, verbose, n_cores);

  S4 z_mat(wrap(z_mat_eig));
  z_mat.slot("Dimnames") = S4(count_mat).slot("Dimnames");

  return z_mat;
}

////// Expression shifts

double estimateCorrelationDistance(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2, bool centered) {
  if (v1.size() != v2.size())
    stop("Vectors must have the same length");

  double m1 = 0, m2 = 0;
  if (centered) {
    m1 = v1.mean();
    m2 = v2.mean();
  }

  double vp = 0, v1s = 0, v2s = 0;
  for (size_t i = 0; i < v1.size(); ++i) {
    if (std::isnan(v1[i]) || std::isnan(v2[i]))
      return NAN;

    double e1 = (v1[i] - m1), e2 = (v2[i] - m2);
    vp += e1 * e2;
    v1s += e1 * e1;
    v2s += e2 * e2;
  }

  return 1 - vp / std::max(std::sqrt(v1s) * std::sqrt(v2s), 1e-10);
}

inline double average(double val1, double val2) {
  return (val1 + val2) / 2;
}

double estimateKLDivergence(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2) {
  double res = 0;
  if (v1.size() != v2.size())
    stop("Vectors must have the same length");

  for (size_t i = 0; i < v1.size(); ++i) {
    double d1 = v1[i], d2 = v2[i];
    if (std::isnan(d1) || std::isnan(d2))
      return NAN;

    if (d1 > 1e-10 && d2 > 1e-10) {
      res += std::log(d1 / d2) * d1;
    }
  }

  return res;
}

double estimateJSDivergence(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2) {
  if (v1.size() != v2.size())
    stop("Vectors must have the same length");

  VectorXd avg = VectorXd::Zero(v1.size());
  std::transform(v1.data(), v1.data() + v1.size(), v2.data(), avg.data(), average);

  double d1 = estimateKLDivergence(v1, avg);
  double d2 = estimateKLDivergence(v2, avg);

  return std::sqrt(0.5 * (d1 + d2));
}

MatrixXd collapseMatrixNorm(const SparseMatrix<double> &mtx, const std::vector<int> &factor,
                            const std::vector<int> &nn_ids, const std::vector<int> &n_obs_per_samp, bool col_norm=true) {
  MatrixXd res = MatrixXd::Zero(mtx.rows(), n_obs_per_samp.size());
  for (int id : nn_ids) {
    int fac = factor[id];
    if (fac >= n_obs_per_samp.size() || fac < 0)
        stop("Wrong factor: " + std::to_string(fac) + ", id: " + std::to_string(id));

    for (SparseMatrix<double, ColMajor>::InnerIterator gene_it(mtx, id); gene_it; ++gene_it) {
        res(gene_it.row(), fac) += gene_it.value() / n_obs_per_samp.at(fac);
    }
  }

  if (col_norm) {
    for (int j = 0; j < res.cols(); j++){
      res.col(j) /= std::max(res.col(j).sum(), 1e-5);
    }
  }

  return res;
}

//' @param sample_per_cell must contains ids from 0 to n_samples-1
//' @param n_samples must be equal to maximum(sample_per_cell) + 1
double estimateCellExpressionShift(const SparseMatrix<double> &cm, const std::vector<int> &sample_per_cell,
                                   const std::vector<int> &nn_ids, const std::vector<bool> &is_ref, const int n_samples,
                                   const int min_n_between, const int min_n_within, const int min_n_obs_per_samp, bool norm_all,
                                   const std::string &dist = "cosine", bool log_vecs=false) {
  std::vector<int> n_ids_per_samp(n_samples, 0);
  for (int id : nn_ids) {
    n_ids_per_samp[sample_per_cell[id]]++;
  }

  auto mat_collapsed = collapseMatrixNorm(cm, sample_per_cell, nn_ids, n_ids_per_samp);
  if (log_vecs) {
    for (int i = 0; i < mat_collapsed.size(); ++i) {
      mat_collapsed(i) = std::log10(1e3 * mat_collapsed(i) + 1);
    }
  }

  std::vector<double> within_dists, between_dists;
  for (int s1 = 0; s1 < n_samples; ++s1) {
    if (n_ids_per_samp.at(s1) < min_n_obs_per_samp)
      continue;

    auto v1 = mat_collapsed.col(s1);
    for (int s2 = s1 + 1; s2 < n_samples; ++s2) {
      if (n_ids_per_samp.at(s2) < min_n_obs_per_samp)
        continue;

      auto v2 = mat_collapsed.col(s2);
      double d;
      if (dist == "cosine") {
        d = estimateCorrelationDistance(v1, v2, false);
      } else if (dist == "js") {
        d = estimateJSDivergence(v1, v2);
      } else if (dist == "cor") {
        d = estimateCorrelationDistance(v1, v2, true);
      } else {
        stop("Unknown dist: ", dist);
      }

      bool is_within = norm_all ? (is_ref.at(s1) == is_ref.at(s2)) : (is_ref.at(s1) && is_ref.at(s2));
      if (is_within) {
        within_dists.push_back(d);
      } else if (is_ref.at(s1) != is_ref.at(s2)) {
        between_dists.push_back(d);
      }
    }
  }

  if ((within_dists.size() < min_n_within) || (between_dists.size() < min_n_between))
    return NAN;

  return median(between_dists) / median(within_dists);
}

// [[Rcpp::export]]
NumericVector estimateClusterFreeExpressionShiftsC(const Eigen::SparseMatrix<double> &cm, IntegerVector sample_per_cell, List nn_ids, const std::vector<bool> &is_ref,
                                                   const int min_n_between=1, const int min_n_within=1, const int min_n_obs_per_samp=1, bool norm_all=false,
                                                   bool verbose=true, int n_cores=1, const std::string &dist="cosine", bool log_vecs=false) {
  const auto samp_per_cell_c = as<std::vector<int>>(IntegerVector(sample_per_cell - 1));
  int n_samples = 0;
  for (int f : samp_per_cell_c) {
    if (f < 0)
      stop("sample_per_cell must contain only positive factors");

    n_samples = std::max(n_samples, f + 1);
  }

  if (n_samples == 0)
      stop("sample_per_cell must be a samp_per_cell_c with non-empty levels");

  std::vector<std::vector<int>> nn_ids_c;
  for (auto &ids : nn_ids) {
    nn_ids_c.emplace_back(as<std::vector<int>>(ids));
  }

  std::vector<double> res_scores(nn_ids_c.size(), 0);

  auto task = [&cm, &samp_per_cell_c, &nn_ids_c, &is_ref, &n_samples, &res_scores, &min_n_between, &min_n_within, &min_n_obs_per_samp, &norm_all,
               dist, log_vecs](int i) {
    res_scores[i] = estimateCellExpressionShift(cm, samp_per_cell_c, nn_ids_c[i], is_ref, n_samples,
                                                min_n_between, min_n_within, min_n_obs_per_samp, norm_all,
                                                dist, log_vecs);
  };
  sccore::runTaskParallelFor(0, nn_ids_c.size(), task, n_cores, verbose);

  NumericVector res = wrap(res_scores);
  res.attr("names") = nn_ids.names();

  return res;
}

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

/// Utils

inline void assert_r(bool condition, const std::string &message) {
    if (!condition) Rcpp::stop(message);
}

std::vector<unsigned> count_values(const std::vector<int> &values, const std::vector<int> &sub_ids, int n_vals=0) {
    if (n_vals == 0) {
        for (int i : sub_ids) {
            int v = values.at(i);
            if (v < 0) stop("sample_per_cell must contain only positive factors");
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
    return (*median_it1 + *median_it2) / 2;
}

double mad(const std::vector<double> &vals, double med) {
    std::vector<double> diffs;
    for (double v : vals) {
        diffs.emplace_back(std::abs(v - med));
    }
    return median(diffs) * 1.4826;
}

double var(const std::vector<double> &vals, double mean) {
    double res = 0;
    for (double v : vals) {
        res += (v - mean) * (v - mean);
    }
    return res / (vals.size() - 1);
}

Eigen::MatrixXd collapseMatrixNorm(const Eigen::SparseMatrix<double> &mtx, const std::vector<int> &factor,
                                   const std::vector<int> &nn_ids, const std::vector<unsigned> &n_obs_per_samp,
                                   int max_factor=0) {
    assert_r(mtx.cols() == factor.size(),
             "Number of columns in matrix (" + std::to_string(mtx.cols()) +
             ") must match the factor size (" + std::to_string(factor.size()) + ")");
    max_factor = std::max(max_factor + 1, int(n_obs_per_samp.size()));
    MatrixXd res = MatrixXd::Zero(mtx.rows(), max_factor);
    for (int id : nn_ids) {
        int fac = factor[id];
        if (fac >= n_obs_per_samp.size() || fac < 0)
            stop("Wrong factor: " + std::to_string(fac) + ", id: " + std::to_string(id));

        for (SparseMatrix<double, ColMajor>::InnerIterator gene_it(mtx, id); gene_it; ++gene_it) {
            res(gene_it.row(), fac) += gene_it.value() / n_obs_per_samp.at(fac);
        }
    }

    return res;
}

Eigen::MatrixXd buildCellXSampleMatrix(const VectorXd &gene_vec, const std::vector<int> &sample_factor,
                                       const std::vector<std::vector<int>> &nn_ids, const std::vector<std::vector<unsigned>> &n_obs_per_samp) {
    assert_r(gene_vec.size() > 0, "gene_vec is empty");
    assert_r(gene_vec.size() == sample_factor.size(),
             "gene_vec size (" + std::to_string(gene_vec.size()) +
             ") must match the sample_factor size (" + std::to_string(sample_factor.size()) + ")");

    assert_r(n_obs_per_samp.size() == sample_factor.size(),
             "n_obs_per_samp size (" + std::to_string(n_obs_per_samp.size()) +
             ") must match the sample_factor size (" + std::to_string(sample_factor.size()) + ")");

    MatrixXd sample_x_cell_cm = MatrixXd::Zero(n_obs_per_samp.at(0).size(), gene_vec.size());
    for (int ci = 0; ci < gene_vec.size(); ++ci) {
        auto n_obs = n_obs_per_samp.at(ci);
        for (int nni : nn_ids.at(ci)) {
            int fac = sample_factor.at(nni);
            sample_x_cell_cm(fac, ci) += gene_vec(nni) / n_obs.at(fac);
        }
    }

    return sample_x_cell_cm;
}

// [[Rcpp::export]]
NumericVector applyMedianFilter(NumericVector signal, std::vector<std::vector<int>> nn_ids) {
  assert_r(signal.size() == nn_ids.size(), "Input vectors must have the same length");

  std::vector<double> signal_c = as<std::vector<double>>(signal);
  std::vector<double> signal_smoothed(signal_c.size());
  for (size_t si = 0; si < signal_c.size(); ++si) {
    if (std::isnan(signal_c.at(si))) {
      signal_smoothed[si] = NAN;
      continue;
    }

    std::vector<double> sig_cur;
    for (int nni : nn_ids.at(si)) {
      double val = signal_c.at(nni);
      if (!std::isnan(val)) {
        sig_cur.emplace_back(val);
      }
    }

    signal_smoothed[si] = median(sig_cur);
  }

  NumericVector signal_smoothed_r = wrap(signal_smoothed);
  signal_smoothed_r.names() = signal.names();
  return signal_smoothed_r;
}

MatrixXd applyMedianFilterMat(const MatrixXd &signal, const std::vector<std::vector<int>> &nn_ids) {
  assert_r(signal.rows() == nn_ids.size(), "Number of rows in signal must match to the length of nn_ids");
  MatrixXd res = MatrixXd::Zero(signal.rows(), signal.cols());
  for (size_t ci = 0; ci < signal.cols(); ++ci) {
    auto sig_col = signal.col(ci);
    auto res_col = res.col(ci);
    for (size_t ri = 0; ri < signal.rows(); ++ri) {
      if (std::isnan(sig_col(ri))) {
        res_col(ri) = NAN;
        continue;
      }

      std::vector<double> sig_cur;
      for (int nni : nn_ids.at(ri)) {
        double val = sig_col(nni);
        if (!std::isnan(val)) {
          sig_cur.emplace_back(val);
        }
      }

      res_col(ri) = median(sig_cur);
    }
  }

  return res;
}

// [[Rcpp::export]]
SEXP applyMedianFilterMat(SEXP signal, std::vector<std::vector<int>> nn_ids) {
  return wrap(applyMedianFilterMat(as<MatrixXd>(signal), nn_ids));
}

// [[Rcpp::export]]
List applyMedianFilterLM(List signal_lst, std::vector<std::vector<int>> nn_ids, int n_cores=1, bool verbose=true) {
    std::mutex mut;
    std::vector<MatrixXd> signal_mats, res_mats(signal_lst.size());
    for (const auto &exp : signal_lst) {
        MatrixXd sm = as<MatrixXd>(exp);
        assert_r(sm.rows() == nn_ids.size(), "Number of rows in signal must match to the length of nn_ids");
        signal_mats.emplace_back(sm);
    }

    auto task = [&signal_mats, &res_mats, &nn_ids, &mut](int mi) {
        MatrixXd res = applyMedianFilterMat(signal_mats[mi], nn_ids);
        {
            std::lock_guard<std::mutex> l(mut);
            res_mats[mi] = res;
        }
    };

    try {
      sccore::runTaskParallelFor(0, signal_mats.size(), task, n_cores, verbose);
    } catch (std::runtime_error &x) {
      Rcpp::stop(x.what());
    }

    List res_lst = wrap(res_mats);
    res_lst.names() = signal_lst.names();
    return res_lst;
}

/// Z-scores

class ZScoreRes {
public:
    std::vector<double> zs;
    std::vector<double> means_ref;
    std::vector<double> means_target;
};

std::pair<std::vector<double>, std::vector<double>>
splitValuesByCondition(const RowVectorXd &vec, const std::vector<bool> &is_ref, const std::vector<unsigned> &n_ids_per_samp, int min_n_obs_per_samp) {
    std::vector<double> ref_vals, target_vals;
    for (int si = 0; si < vec.size(); ++si) {
        if (n_ids_per_samp[si] < min_n_obs_per_samp)
            continue;

        if (is_ref.at(si)) {
            ref_vals.emplace_back(vec[si]);
        } else {
            target_vals.emplace_back(vec[si]);
        }
    }

    return std::make_pair(ref_vals, target_vals);
}

std::tuple<double, double, double> estimateCellGeneZScore(std::vector<double> &ref_vals, std::vector<double> &target_vals,
                                                          int min_n_samp_per_cond, bool robust, bool norm_both) {
    size_t n_targ = target_vals.size(), n_ref = ref_vals.size();
    if (std::min(n_targ, n_ref) < min_n_samp_per_cond)
        return std::make_tuple(NAN, NAN, NAN);

    double m_ref, m_targ, var_ref, var_targ;
    if (robust) {
        m_ref = median(ref_vals);
        m_targ = median(target_vals);
        var_ref = std::pow(mad(ref_vals, m_ref), 2);
        var_targ = std::pow(mad(target_vals, m_targ), 2);
    } else {
        m_ref = std::accumulate(ref_vals.begin(), ref_vals.end(), 0.0) / n_ref;
        m_targ = std::accumulate(target_vals.begin(), target_vals.end(), 0.0) / n_targ;
        var_ref = var(ref_vals, m_ref);
        var_targ = var(target_vals, m_targ);
    }

    double sd_tot = norm_both ?
                    std::sqrt((n_ref * var_ref + n_targ * var_targ) / (n_ref + n_targ)) :
                    std::sqrt(var_ref);
    double z = (sd_tot < 1e-20) ? NAN : ((m_targ - m_ref) / (sd_tot));
    return std::make_tuple(z, m_ref, m_targ);
}

ZScoreRes estimateGeneZScore(const VectorXd &gene_vec, const std::vector<int> &sample_per_cell,
                             const std::vector<std::vector<int>> &nn_ids, const std::vector<std::vector<unsigned>> &n_obs_per_samp,
                             const std::vector<bool> &is_ref, int min_n_samp_per_cond, int min_n_obs_per_samp, bool robust, bool norm_both) {
    min_n_samp_per_cond = std::max(min_n_samp_per_cond, 2);
    auto sample_x_cell_cm = buildCellXSampleMatrix(gene_vec, sample_per_cell, nn_ids, n_obs_per_samp);

    std::vector<double> res_z, res_ms_ref, res_ms_target;
    for (int ci = 0; ci < sample_x_cell_cm.cols(); ++ci) {
        auto c_vec = sample_x_cell_cm.col(ci);
        auto split_res = splitValuesByCondition(c_vec, is_ref, n_obs_per_samp.at(ci), min_n_obs_per_samp);
        auto z_res = estimateCellGeneZScore(split_res.first, split_res.second, min_n_samp_per_cond, robust, norm_both);

        res_z.emplace_back(std::get<0>(z_res));
        res_ms_ref.emplace_back(std::get<1>(z_res));
        res_ms_target.emplace_back(std::get<2>(z_res));
    }

    return ZScoreRes{res_z, res_ms_ref, res_ms_target};
}

// [[Rcpp::export]]
List clusterFreeZScoreMat2(const SEXP count_mat, IntegerVector sample_per_cell, List nn_ids, const std::vector<bool> &is_ref,
                           int min_n_samp_per_cond=2, int min_n_obs_per_samp=1, bool robust=false, bool norm_both=false,
                           double min_z=0.001, bool verbose=true, int n_cores=1, bool return_demo=false) {
    const auto samp_per_cell_c = as<std::vector<int>>(IntegerVector(sample_per_cell - 1));
    if (sample_per_cell.size() == 0 || (*std::max_element(samp_per_cell_c.begin(), samp_per_cell_c.end()) <= 0))
        stop("sample_per_cell must be a factor vector with non-empty levels");

    std::vector<std::vector<int>> nn_ids_c;
    std::vector<std::vector<unsigned>> n_obs_per_samp;
    int n_samples = (*std::max_element(samp_per_cell_c.begin(), samp_per_cell_c.end())) + 1;
    for (auto &ids : nn_ids) {
        auto c_ids = as<std::vector<int>>(ids);
        nn_ids_c.emplace_back(c_ids);
        n_obs_per_samp.emplace_back(count_values(samp_per_cell_c, c_ids, n_samples));
    }

    auto cm = as<SparseMatrix<double>>(count_mat);
    assert_r(cm.rows() == samp_per_cell_c.size(),
             "Number of rows in matrix (" + std::to_string(cm.rows()) +
             ") must match the sample_per_cell size (" + std::to_string(samp_per_cell_c.size()) + ")");

    if (return_demo) {
        auto tr_mat = buildCellXSampleMatrix(cm.col(0), samp_per_cell_c, nn_ids_c, n_obs_per_samp);
        return List::create(wrap(SparseMatrix<double>(tr_mat.sparseView())));
    }

    std::vector<Triplet<double>> z_triplets, m_ref_triplets, m_targ_triplets;
    std::mutex mut;
    auto task = [&cm, &samp_per_cell_c, &nn_ids_c, &n_obs_per_samp, &is_ref, &min_n_samp_per_cond, &min_n_obs_per_samp,
                 &min_z, &robust, &norm_both, &z_triplets, &m_ref_triplets, &m_targ_triplets, &mut](int gi) {
        VectorXd gene_vec = cm.col(gi);
        auto cell_res = estimateGeneZScore(gene_vec, samp_per_cell_c, nn_ids_c, n_obs_per_samp, is_ref, min_n_samp_per_cond, min_n_obs_per_samp, robust, norm_both);
        std::vector<Triplet<double>> z_trip_c, m_ref_trip_c, m_targ_trip_c;
        for (int ci = 0; ci < cell_res.zs.size(); ++ci) {
            if (gene_vec(ci) < 1e-20)  // TODO: try to remove it or move to buildCellXSampleMatrix
                continue;

            double z_cur = cell_res.zs[ci];
            if (std::isnan(z_cur) || std::abs(z_cur) >= min_z) {
                z_trip_c.emplace_back(ci, gi, z_cur);
                m_ref_trip_c.emplace_back(ci, gi, cell_res.means_ref[ci]);
                m_targ_trip_c.emplace_back(ci, gi, cell_res.means_target[ci]);
            }
        }

        {
            std::lock_guard<std::mutex> l(mut);
            z_triplets.insert(z_triplets.end(), std::make_move_iterator(z_trip_c.begin()), std::make_move_iterator(z_trip_c.end()));
            m_ref_triplets.insert(m_ref_triplets.end(), std::make_move_iterator(m_ref_trip_c.begin()), std::make_move_iterator(m_ref_trip_c.end()));
            m_targ_triplets.insert(m_targ_triplets.end(), std::make_move_iterator(m_targ_trip_c.begin()), std::make_move_iterator(m_targ_trip_c.end()));
        }
    };

    try {
        sccore::runTaskParallelFor(0, cm.cols(), task, n_cores, verbose);
    } catch (std::runtime_error &x) {
        Rcpp::stop(x.what());
    }

    SparseMatrix<double> z_mat(cm.rows(), cm.cols()), m_ref_mat(cm.rows(), cm.cols()), m_targ_mat(cm.rows(), cm.cols());
    z_mat.setFromTriplets(z_triplets.begin(), z_triplets.end());
    m_ref_mat.setFromTriplets(m_ref_triplets.begin(), m_ref_triplets.end());
    m_targ_mat.setFromTriplets(m_targ_triplets.begin(), m_targ_triplets.end());

    S4 z_mat_r(wrap(z_mat)), m_ref_mat_r(wrap(m_ref_mat)), m_targ_mat_r(wrap(m_targ_mat));
    m_targ_mat_r.slot("Dimnames") = m_ref_mat_r.slot("Dimnames") = z_mat_r.slot("Dimnames") =S4(count_mat).slot("Dimnames");

    return List::create(_["z"] = z_mat_r, _["reference"] = m_ref_mat_r, _["target"] = m_targ_mat_r);
}

ZScoreRes estimateCellZScore(const SparseMatrix<double> &cm, const std::vector<int> &sample_per_cell,
                             const std::vector<int> &nn_ids, const std::vector<bool> &is_ref,
                             int min_n_samp_per_cond, int min_n_obs_per_samp, bool robust, bool norm_both) {
    min_n_samp_per_cond = std::max(min_n_samp_per_cond, 2);
    auto n_ids_per_samp = count_values(sample_per_cell, nn_ids);
    auto mat_collapsed = collapseMatrixNorm(cm, sample_per_cell, nn_ids, n_ids_per_samp);

    std::vector<double> res_z, res_ms_ref, res_ms_target;
    for (int gi = 0; gi < mat_collapsed.rows(); ++gi) {
        auto g_vec = mat_collapsed.row(gi);
        auto split_res = splitValuesByCondition(g_vec, is_ref, n_ids_per_samp, min_n_obs_per_samp);
        auto z_res = estimateCellGeneZScore(split_res.first, split_res.second, min_n_samp_per_cond, robust, norm_both);

        res_z.emplace_back(std::get<0>(z_res));
        res_ms_ref.emplace_back(std::get<1>(z_res));
        res_ms_target.emplace_back(std::get<2>(z_res));
    }

    return ZScoreRes{res_z, res_ms_ref, res_ms_target};
}

// [[Rcpp::export]]
List clusterFreeZScoreMat(const SEXP count_mat, IntegerVector sample_per_cell, List nn_ids, const std::vector<bool> &is_ref,
                          int min_n_samp_per_cond=2, int min_n_obs_per_samp=1, bool robust=false, bool norm_both=false,
                          double min_z=0.001, bool verbose=true, int n_cores=1) {
    const auto samp_per_cell_c = as<std::vector<int>>(IntegerVector(sample_per_cell - 1));
    if (sample_per_cell.size() == 0 || (*std::max_element(samp_per_cell_c.begin(), samp_per_cell_c.end()) <= 0))
        stop("sample_per_cell must be a factor vector with non-empty levels");

    std::vector<std::vector<int>> nn_ids_c;
    for (auto &ids : nn_ids) {
        nn_ids_c.emplace_back(as<std::vector<int>>(ids));
    }

    auto cm = as<SparseMatrix<double>>(count_mat);
    assert_r(cm.cols() == samp_per_cell_c.size(),
             "Number of columns in matrix (" + std::to_string(cm.cols()) +
             ") must match the sample_per_cell size (" + std::to_string(samp_per_cell_c.size()) + ")");

    std::vector<Triplet<double>> z_triplets, m_ref_triplets, m_targ_triplets;
    std::mutex mut;
    auto task = [&cm, &samp_per_cell_c, &nn_ids_c, &is_ref, &min_n_samp_per_cond, &min_n_obs_per_samp, &min_z, &robust,
                 &norm_both, &z_triplets, &m_ref_triplets, &m_targ_triplets, &mut](int ci) {
        auto cell_res = estimateCellZScore(cm, samp_per_cell_c, nn_ids_c[ci], is_ref, min_n_samp_per_cond, min_n_obs_per_samp, robust, norm_both);
        VectorXd expr = cm.col(ci);
        for (int gi = 0; gi < cell_res.zs.size(); ++gi) {
            if (expr(gi) < 1e-20)
                continue;

            double z_cur = cell_res.zs[gi];
            if (std::isnan(z_cur) || std::abs(z_cur) >= min_z) {
                std::lock_guard<std::mutex> l(mut);
                z_triplets.emplace_back(ci, gi, z_cur);
                m_ref_triplets.emplace_back(ci, gi, cell_res.means_ref[gi]);
                m_targ_triplets.emplace_back(ci, gi, cell_res.means_target[gi]);
            }
        }

    };

    try {
        sccore::runTaskParallelFor(0, nn_ids_c.size(), task, n_cores, verbose);
    } catch (std::runtime_error &x) {
        Rcpp::stop(x.what());
    }


    SparseMatrix<double> z_mat(cm.cols(), cm.rows()), m_ref_mat(cm.cols(), cm.rows()), m_targ_mat(cm.cols(), cm.rows());
    z_mat.setFromTriplets(z_triplets.begin(), z_triplets.end());
    m_ref_mat.setFromTriplets(m_ref_triplets.begin(), m_ref_triplets.end());
    m_targ_mat.setFromTriplets(m_targ_triplets.begin(), m_targ_triplets.end());

    S4 z_mat_r(wrap(z_mat)), m_ref_mat_r(wrap(m_ref_mat)), m_targ_mat_r(wrap(m_targ_mat));
    auto cm_dimnames = as<List>(S4(count_mat).slot("Dimnames"));
    if (cm_dimnames.size() > 1) {
        m_targ_mat_r.slot("Dimnames") = m_ref_mat_r.slot("Dimnames") = z_mat_r.slot("Dimnames") = List::create(cm_dimnames[1], cm_dimnames[0]);
    }

    return List::create(_["z"] = z_mat_r, _["reference"] = m_ref_mat_r, _["target"] = m_targ_mat_r);
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

//' @param sample_per_cell must contains ids from 0 to n_samples-1
//' @param n_samples must be equal to maximum(sample_per_cell) + 1
double estimateCellExpressionShift(const SparseMatrix<double> &cm, const std::vector<int> &sample_per_cell,
                                   const std::vector<int> &nn_ids, const std::vector<bool> &is_ref,
                                   const int min_n_between, const int min_n_within, const int min_n_obs_per_samp, bool norm_all,
                                   const std::string &dist = "cosine", bool log_vecs=false) {
    auto n_ids_per_samp = count_values(sample_per_cell, nn_ids);
    auto mat_collapsed = collapseMatrixNorm(cm, sample_per_cell, nn_ids, n_ids_per_samp);
    if (log_vecs) {
        for (int i = 0; i < mat_collapsed.size(); ++i) {
            mat_collapsed(i) = std::log10(1e3 * mat_collapsed(i) + 1);
        }
    }

    std::vector<double> within_dists, between_dists;
    for (int s1 = 0; s1 < n_ids_per_samp.size(); ++s1) {
        if (n_ids_per_samp.at(s1) < min_n_obs_per_samp)
            continue;

        auto v1 = mat_collapsed.col(s1);
        for (int s2 = s1 + 1; s2 < n_ids_per_samp.size(); ++s2) {
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
    if (sample_per_cell.size() == 0 || (*std::max_element(samp_per_cell_c.begin(), samp_per_cell_c.end()) <= 0))
        stop("sample_per_cell must be a factor vector with non-empty levels");

    std::vector<std::vector<int>> nn_ids_c;
    for (auto &ids : nn_ids) {
        nn_ids_c.emplace_back(as<std::vector<int>>(ids));
    }

    std::vector<double> res_scores(nn_ids_c.size(), 0);

    auto task = [&cm, &samp_per_cell_c, &nn_ids_c, &is_ref, &res_scores, &min_n_between, &min_n_within, &min_n_obs_per_samp, &norm_all,
            dist, log_vecs](int i) {
        res_scores[i] = estimateCellExpressionShift(cm, samp_per_cell_c, nn_ids_c[i], is_ref,
                                                    min_n_between, min_n_within, min_n_obs_per_samp, norm_all,
                                                    dist, log_vecs);
    };
    sccore::runTaskParallelFor(0, nn_ids_c.size(), task, n_cores, verbose);

    NumericVector res = wrap(res_scores);
    res.attr("names") = nn_ids.names();

    return res;
}

/* P-value processing */

// [[Rcpp::export]]
std::vector<std::vector<int>> mapIds(std::vector<std::vector<int>> ids_vec, std::vector<int> id_map) {
    std::map<int, int> id_map_c;
    for (int i = 0; i < id_map.size(); ++i) {
        id_map_c.emplace(id_map[i], i + 1);
    }

    std::vector<std::vector<int>> res_ids;
    for (auto const &ids : ids_vec) {
        std::vector<int> mapped;
        for (int id : ids) {
            auto iter = id_map_c.find(id);
            if (iter != id_map_c.end()) {
                mapped.emplace_back(iter->second);
            }
        }

        res_ids.emplace_back(mapped);
    }

    return res_ids;
}

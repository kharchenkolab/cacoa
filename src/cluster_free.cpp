#include <RcppEigen.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <mutex>
#include <random>
#include <vector>

#include <sccore_par.hpp>

#include <progress.hpp>
#include <unistd.h>

using namespace Rcpp;
using namespace Eigen;

const double EPS = std::numeric_limits<double>::epsilon();

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

MatrixXd buildCellXSampleMatrix(const VectorXd &gene_vec, const std::vector<int> &sample_factor,
                                const std::vector<std::vector<int>> &nn_ids, const std::vector<std::vector<unsigned>> &n_obs_per_samp,
                                const std::vector<size_t>& non_zero_ids) {
    assert_r(gene_vec.size() > 0, "gene_vec is empty");
    assert_r(gene_vec.size() == sample_factor.size(),
             "gene_vec size (" + std::to_string(gene_vec.size()) +
             ") must match the sample_factor size (" + std::to_string(sample_factor.size()) + ")");

    assert_r(n_obs_per_samp.size() == sample_factor.size(),
             "n_obs_per_samp size (" + std::to_string(n_obs_per_samp.size()) +
             ") must match the sample_factor size (" + std::to_string(sample_factor.size()) + ")");

    MatrixXd sample_x_cell_cm = MatrixXd::Zero(n_obs_per_samp.at(0).size(), gene_vec.size());
    for (size_t ci : non_zero_ids) {
        auto n_obs = n_obs_per_samp.at(ci);
        for (int nni : nn_ids.at(ci)) {
            int fac = sample_factor.at(nni);
            sample_x_cell_cm(fac, ci) += gene_vec(nni) / n_obs.at(fac);
        }
    }

    return sample_x_cell_cm;
}

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

std::vector<size_t> findNonZeroInds(const VectorXd &vec) {
    std::vector<size_t> nz_ids;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (vec(i) > EPS) {
            nz_ids.emplace_back(i);
        }
    }

    return nz_ids;
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
    double z = (sd_tot < EPS) ? NAN : ((m_targ - m_ref) / (sd_tot));
    return std::make_tuple(z, m_ref, m_targ);
}

class ClusterFreeDEParams {
public:
    int min_n_samp_per_cond;
    int min_n_obs_per_samp;
    bool robust;
    bool norm_both;
    double min_z;
    bool verbose;
    int n_cores;
    bool adjust_pvalues;
    int n_permutations;
    bool smooth;
    double wins;
};

ZScoreRes estimateGeneZScore(const VectorXd &gene_vec, const MatrixXd &sample_x_cell_cm, const std::vector<std::vector<unsigned>> &n_obs_per_samp,
                             const std::vector<bool> &is_ref, const ClusterFreeDEParams &pars) {
    std::vector<double> res_z, res_ms_ref, res_ms_target;
    for (int ci = 0; ci < sample_x_cell_cm.cols(); ++ci) {
        auto z_res = std::make_tuple(0.0, 0.0, 0.0);
        if (gene_vec(ci) > EPS) {
            auto c_vec = sample_x_cell_cm.col(ci);
            auto split_res = splitValuesByCondition(c_vec, is_ref, n_obs_per_samp.at(ci), pars.min_n_obs_per_samp);
            z_res = estimateCellGeneZScore(split_res.first, split_res.second, pars.min_n_samp_per_cond, pars.robust, pars.norm_both);
        }

        res_z.emplace_back(std::get<0>(z_res));
        res_ms_ref.emplace_back(std::get<1>(z_res));
        res_ms_target.emplace_back(std::get<2>(z_res));
    }

    return ZScoreRes{res_z, res_ms_ref, res_ms_target};
}

std::pair<std::vector<double>, std::vector<double>>
estimateNullZScoreRanges(const VectorXd &gene_vec, const MatrixXd &sample_x_cell_cm, const std::vector<std::vector<int>> &nn_ids, const std::vector<bool> &is_ref,
                         const std::vector<std::vector<unsigned>> &n_obs_per_samp, const std::vector<size_t>& non_zero_ids,
                         const ClusterFreeDEParams &pars, unsigned seed=0) {
    std::vector<bool> is_ref_shuffled(is_ref.begin(), is_ref.end());
    std::mt19937 g(seed);

    std::vector<double> min_vals, max_vals;
    for (int rep = 0; rep < pars.n_permutations; ++rep) {
        std::shuffle(is_ref_shuffled.begin(), is_ref_shuffled.end(), g);
        auto perm_zs = estimateGeneZScore(gene_vec, sample_x_cell_cm, n_obs_per_samp, is_ref_shuffled, pars).zs;
        if (pars.smooth) { // TODO: estimating z-scores per many permutations at once would allow to run median filter per column, saving a lot of time on random access indexing
            perm_zs = applyMedianFilter(perm_zs, nn_ids, non_zero_ids);
        }

        auto z_range = range(perm_zs, pars.wins);
        if (!std::isnan(z_range.first)) {
            min_vals.emplace_back(z_range.first);
            max_vals.emplace_back(z_range.second);
        }
    }

    std::sort(min_vals.begin(), min_vals.end());
    std::sort(max_vals.begin(), max_vals.end());
    return std::make_pair(min_vals, max_vals);
}

std::vector<double> adjustZScoresWithPermutations(const std::vector<double> &z_scores, const std::vector<std::vector<int>> &nn_ids,
                                                  double wins, bool smooth, const std::vector<double> &max_vals) {
    std::vector<double> z_adj(z_scores.begin(), z_scores.end());
    if (smooth) {
        z_adj = applyMedianFilter(z_adj, nn_ids);
    }

    auto max_val = range(z_adj, wins).second;
    for (double &z : z_adj) {
        if (std::isnan(z))
            continue;

        z = std::min(z, max_val);
        size_t n = (max_vals.end() - std::lower_bound(max_vals.begin(), max_vals.end(), z - EPS)); // Number of elements that are >=z
        z = std::max(1.0 - (n + 1.0) / (max_vals.size() + 1.0), 0.5);
    }

    return as<std::vector<double>>(NumericVector(qnorm(NumericVector(wrap(z_adj)))));
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
        z_adj = as<std::vector<double>>(NumericVector(qnorm(NumericVector(wrap(z_adj)))));
    }

    for (size_t i = 0; i < z_adj.size(); ++i) {
        z_adj[i] = std::copysign(z_adj[i], z_scores[i]);
    }

    return z_adj;
}

std::tuple<SparseMatrix<double>, SparseMatrix<double>, SparseMatrix<double>, SparseMatrix<double>>
clusterFreeZScoreMat(const SparseMatrix<double> &cm, const std::vector<int> &sample_per_cell,
                     const std::vector<std::vector<int>> &nn_ids, const std::vector<bool> &is_ref,
                     const std::vector<std::vector<unsigned>> &n_obs_per_samp, const ClusterFreeDEParams &pars) {
    std::vector<Triplet<double>> z_triplets, z_adj_triplets, m_ref_triplets, m_targ_triplets;
    std::mutex save_mut, r_mut;
    auto task = [&cm, &sample_per_cell, &nn_ids, &n_obs_per_samp, &is_ref, &pars,
                 &z_triplets, &z_adj_triplets, &m_ref_triplets, &m_targ_triplets, &save_mut, &r_mut](int gi) {
        VectorXd gene_vec = cm.col(gi);
        auto nz_ids = findNonZeroInds(gene_vec);
        auto sample_x_cell_cm = buildCellXSampleMatrix(gene_vec, sample_per_cell, nn_ids, n_obs_per_samp, nz_ids);
        auto gene_res = estimateGeneZScore(gene_vec, sample_x_cell_cm, n_obs_per_samp, is_ref, pars);

        std::vector<double> z_adj;
        if (pars.adjust_pvalues) {
            auto min_max_vals = estimateNullZScoreRanges(gene_vec, sample_x_cell_cm, nn_ids, is_ref, n_obs_per_samp, nz_ids, pars);
            z_adj = adjustZScoresWithPermutations(gene_res.zs, nn_ids, nz_ids, pars.wins, pars.smooth,
                                                  min_max_vals.first, min_max_vals.second, r_mut);
        }

        std::vector<Triplet<double>> z_trip_c, z_adj_trip_c, m_ref_trip_c, m_targ_trip_c;
        for (int ci = 0; ci < gene_res.zs.size(); ++ci) {
            double z_raw = gene_res.zs[ci];
            if (std::isnan(z_raw) || std::abs(z_raw) >= pars.min_z) {
                z_trip_c.emplace_back(ci, gi, gene_res.zs[ci]);
                m_ref_trip_c.emplace_back(ci, gi, gene_res.means_ref[ci]);
                m_targ_trip_c.emplace_back(ci, gi, gene_res.means_target[ci]);

                if (pars.adjust_pvalues) {
                    z_adj_trip_c.emplace_back(ci, gi, z_adj[ci]);
                }
            }
        }

        {
            std::lock_guard<std::mutex> l(save_mut);
            z_triplets.insert(z_triplets.end(), std::make_move_iterator(z_trip_c.begin()), std::make_move_iterator(z_trip_c.end()));
            z_adj_triplets.insert(z_adj_triplets.end(), std::make_move_iterator(z_adj_trip_c.begin()), std::make_move_iterator(z_adj_trip_c.end()));
            m_ref_triplets.insert(m_ref_triplets.end(), std::make_move_iterator(m_ref_trip_c.begin()), std::make_move_iterator(m_ref_trip_c.end()));
            m_targ_triplets.insert(m_targ_triplets.end(), std::make_move_iterator(m_targ_trip_c.begin()), std::make_move_iterator(m_targ_trip_c.end()));
        }
    };

    try {
        sccore::runTaskParallelFor(0, cm.cols(), task, pars.n_cores, pars.verbose);
    } catch (std::runtime_error &x) {
        Rcpp::stop(x.what());
    }

    SparseMatrix<double> z_mat(cm.rows(), cm.cols()), z_adj_mat(cm.rows(), cm.cols()),
        m_ref_mat(cm.rows(), cm.cols()), m_targ_mat(cm.rows(), cm.cols());
    z_mat.setFromTriplets(z_triplets.begin(), z_triplets.end());
    m_ref_mat.setFromTriplets(m_ref_triplets.begin(), m_ref_triplets.end());
    m_targ_mat.setFromTriplets(m_targ_triplets.begin(), m_targ_triplets.end());
    z_adj_mat.setFromTriplets(z_adj_triplets.begin(), z_adj_triplets.end());

    return std::make_tuple(z_mat, m_ref_mat, m_targ_mat, z_adj_mat);
}

// [[Rcpp::export]]
List clusterFreeZScoreMat(const SEXP count_mat, IntegerVector sample_per_cell, List nn_ids, const std::vector<bool> &is_ref,
                          int min_n_samp_per_cond=2, int min_n_obs_per_samp=1, bool robust=false, bool norm_both=true,
                          double min_z=0.001, bool verbose=true, int n_cores=1, bool adjust_pvalues=false, int n_permutations=500,
                          bool smooth=true, double wins=0.01) {
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

    const ClusterFreeDEParams pars{std::max(min_n_samp_per_cond, 2), min_n_obs_per_samp, robust, norm_both,
                                   min_z, verbose, n_cores, adjust_pvalues, n_permutations, smooth, wins};
    auto mats = clusterFreeZScoreMat(cm, samp_per_cell_c, nn_ids_c, is_ref, n_obs_per_samp, pars);

    S4 z_mat_r(wrap(std::get<0>(mats))), m_ref_mat_r(wrap(std::get<1>(mats))),
        m_targ_mat_r(wrap(std::get<2>(mats))), z_adj_mat_r(wrap(std::get<3>(mats)));
    m_targ_mat_r.slot("Dimnames") = m_ref_mat_r.slot("Dimnames") = z_mat_r.slot("Dimnames") =
        z_adj_mat_r.slot("Dimnames") = S4(count_mat).slot("Dimnames");

    return List::create(_["z"] = z_mat_r, _["reference"] = m_ref_mat_r,
                        _["target"] = m_targ_mat_r, _["z.adj"] = z_adj_mat_r);
}

////// Expression shifts

struct CFShiftResult{
    std::vector<double> dists;
    std::vector<size_t> s1_ids;
    std::vector<size_t> s2_ids;
};

// [[Rcpp::export]]
double estimateCorrelationDistance(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2, bool centered=true) {
    if (v1.size() != v2.size())
        stop("Vectors must have the same length");

    double m1 = 0, m2 = 0;
    if (centered) { // correlation instead of cosine distance
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

double estimateNormalizedExpressionShift(const std::vector<double> &dists, const std::vector<bool> &is_ref,
                                         const std::vector<size_t> &s1_ids, const std::vector<size_t> &s2_ids,
                                         bool norm_all, size_t min_n_within, size_t min_n_between) {
    std::vector<double> between_dists, ref_dists, non_ref_dists;
    for (size_t i = 0; i < dists.size(); ++i) {
        double d = dists.at(i);
        size_t s1 = s1_ids.at(i), s2 = s2_ids.at(i);
        if (is_ref.at(s1) != is_ref.at(s2)) {
            between_dists.push_back(d);
        } else if (is_ref.at(s1)) {
            ref_dists.push_back(d);
        } else if (norm_all) {
            non_ref_dists.push_back(d);
        }
    }

    bool is_nan = (
        (between_dists.size() < min_n_between) ||
        (ref_dists.size() < min_n_within) ||
        (norm_all && (non_ref_dists.size() < min_n_within))
    );

    if (is_nan)
        return NAN;

    double norm_const = norm_all ? ((median(ref_dists) + median(non_ref_dists)) / 2) : median(ref_dists);

    return median(between_dists) - norm_const;
}

double estimateVectorDistance(const VectorXd &v1, const VectorXd &v2, const std::string &dist) {
    if (dist == "cosine")
        return estimateCorrelationDistance(v1, v2, false);
    if (dist == "js")
        return estimateJSDivergence(v1, v2);
    if (dist == "cor")
        return estimateCorrelationDistance(v1, v2, true);

    stop("Unknown dist: ", dist);
}

//' @param sample_per_cell must contains ids from 0 to n_samples-1
//' @param n_samples must be equal to maximum(sample_per_cell) + 1
CFShiftResult estimateCellExpressionShift(const SparseMatrix<double> &cm, const std::vector<int> &sample_per_cell,
                                          const std::vector<int> &nn_ids,
                                          size_t min_n_obs_per_samp, const std::string &dist="cosine", bool log_vecs=true) {
    auto n_ids_per_samp = count_values(sample_per_cell, nn_ids);
    auto mat_collapsed = collapseMatrixNorm(cm, sample_per_cell, nn_ids, n_ids_per_samp);
    if (log_vecs) {
        for (int i = 0; i < mat_collapsed.size(); ++i) {
            mat_collapsed(i) = std::log10(1e3 * mat_collapsed(i) + 1);
        }
    }

    std::vector<double> dists;
    std::vector<size_t> s1_ids, s2_ids;
    for (size_t s1 = 0; s1 < n_ids_per_samp.size(); ++s1) {
        if (n_ids_per_samp.at(s1) < min_n_obs_per_samp)
            continue;

        auto v1 = mat_collapsed.col(s1);
        for (size_t s2 = s1 + 1; s2 < n_ids_per_samp.size(); ++s2) {
            if (n_ids_per_samp.at(s2) < min_n_obs_per_samp)
                continue;

            auto v2 = mat_collapsed.col(s2);
            double d = estimateVectorDistance(v1, v2, dist);

            dists.push_back(d);
            s1_ids.push_back(s1);
            s2_ids.push_back(s2);
        }
    }

    return CFShiftResult{dists, s1_ids, s2_ids};
}

// [[Rcpp::export]]
List estimateClusterFreeExpressionShiftsC(const Eigen::SparseMatrix<double> &cm, IntegerVector sample_per_cell, List nn_ids,
                                          const std::vector<bool> &is_ref, const int min_n_between=1, const int min_n_within=1,
                                          const int min_n_obs_per_samp=1, bool norm_all=true, bool verbose=true, int n_cores=1,
                                          const std::string &dist="cor", bool log_vecs=true, int n_permutations=100,
                                          double wins=0.01) {
    const auto samp_per_cell_c = as<std::vector<int>>(IntegerVector(sample_per_cell - 1));
    if (sample_per_cell.size() == 0 || (*std::max_element(samp_per_cell_c.begin(), samp_per_cell_c.end()) <= 0))
        stop("sample_per_cell must be a factor vector with non-empty levels");

    std::vector<std::vector<int>> nn_ids_c;
    for (auto &ids : nn_ids) {
        nn_ids_c.emplace_back(as<std::vector<int>>(ids));
    }

    std::vector<double> res_scores(nn_ids_c.size(), 0);
    std::vector<CFShiftResult> res_info(nn_ids_c.size());

    // Prepare info for permutations

    auto task = [&cm, &samp_per_cell_c, &nn_ids_c, &is_ref, &res_scores, &res_info, min_n_obs_per_samp, dist, log_vecs,
                 norm_all, min_n_within, min_n_between](int i) {
        auto res = estimateCellExpressionShift(cm, samp_per_cell_c, nn_ids_c[i], min_n_obs_per_samp, dist, log_vecs);
        double d = estimateNormalizedExpressionShift(res.dists, is_ref, res.s1_ids, res.s2_ids,
            norm_all, min_n_within, min_n_between);
        res_scores[i] = std::max(d, 0.0);
        res_info[i] = res;
    };
    sccore::runTaskParallelFor(0, nn_ids_c.size(), task, n_cores, verbose);

    if (n_permutations == 0) {
        NumericVector scores_r = wrap(res_scores);
        scores_r.attr("names") = nn_ids.names();

        return List::create(_["shifts"] = scores_r);
    }

    // Run permutations

    std::vector<double> max_vals(n_permutations, 0);
    std::vector<double> sum_null_dists(nn_ids_c.size(), 0);
    std::vector<size_t> n_null_dists(nn_ids_c.size(), 0);
    std::mutex mut;

    auto task2 = [&nn_ids_c, &res_info, &max_vals, &is_ref, &mut, &sum_null_dists, &n_null_dists, norm_all,
            min_n_within, min_n_between, wins](int ri) {
        std::mt19937 g(ri);
        std::vector<bool> is_ref_shuffled(is_ref.begin(), is_ref.end());
        std::shuffle(is_ref_shuffled.begin(), is_ref_shuffled.end(), g);
        std::vector<double> shuff_scores(nn_ids_c.size(), 0);
        for (size_t ci = 0; ci < nn_ids_c.size(); ++ci) {
            const auto res = res_info.at(ci);

            double d = estimateNormalizedExpressionShift(res.dists, is_ref_shuffled, res.s1_ids, res.s2_ids,
                                                         norm_all, min_n_within, min_n_between);

            if (std::isnan(d))
                continue;

            shuff_scores.at(ci) = d;

            {
                std::lock_guard<std::mutex> l(mut);
                sum_null_dists.at(ci) += d;
                n_null_dists.at(ci)++;
            }
        }

        shuff_scores = applyMedianFilter(shuff_scores, nn_ids_c);
        max_vals.at(ri) = range(shuff_scores, wins).second;
    };
    sccore::runTaskParallelFor(0, n_permutations, task2, n_cores, verbose);

    auto z_scores = adjustZScoresWithPermutations(res_scores, nn_ids_c, wins, true, max_vals);
    for (size_t ci = 0; ci < res_scores.size(); ++ci) {
        res_scores.at(ci) -= sum_null_dists.at(ci) / n_null_dists.at(ci);
    }

    auto scores_smoothed = applyMedianFilter(res_scores, nn_ids_c);
    NumericVector scores_r = wrap(res_scores), z_scores_r = wrap(z_scores), scores_smoothed_r = wrap(scores_smoothed);
    scores_smoothed_r.attr("names") = z_scores_r.attr("names") = scores_r.attr("names") = nn_ids.names();

    return List::create(_["shifts"]=scores_r, _["shifts_smoothed"]=scores_smoothed_r, _["z_scores"]=z_scores_r);
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

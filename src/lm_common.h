#ifndef STATS_COMMONS_H
#define STATS_COMMONS_H

#include <RcppArmadillo.h>
#include <omp.h>
#include <random>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <ctime>

// Disable Armadillo's internal OpenMP to avoid thread oversubscription.
#define ARMA_DONT_USE_OPENMP 

using namespace Rcpp;
using namespace arma;

// --- SHARED STRUCTS ---
/**
 * @brief Efficiently maps sample pairs (nodes) to row indices (edges).
 * Used during Graph/MRQAP randomization to translate a shuffled sample vector
 * into the corresponding permutation of the distance/data vector.
 */
struct PairLookup {
    std::vector<int> lookup; 
    arma::uword n_samples;
    bool active;

    PairLookup() : active(false) {}

    void init(const arma::umat& pairs, arma::uword n) {
        n_samples = n;
        lookup.assign(n * n, -1); 
        for(arma::uword k=0; k < pairs.n_rows; ++k) {
            arma::uword i = pairs(k, 0);
            arma::uword j = pairs(k, 1);
            if (i < n && j < n) {
                lookup[i * n + j] = k;
                lookup[j * n + i] = k;
            }
        }
        active = true;
    }

    void get_row_perm(arma::uvec& row_perm, const arma::umat& pairs, const arma::uvec& sample_perm) {
        arma::uword N = pairs.n_rows;
        row_perm.set_size(N);
        for(arma::uword k=0; k < N; ++k) {
            int new_idx = lookup[sample_perm[pairs(k,0)] * n_samples + sample_perm[pairs(k,1)]];
            row_perm[k] = (new_idx >= 0) ? (arma::uword)new_idx : k; 
        }
    }
};

struct Config {
  std::string robust, na_mode;
  double huber_k, huber_tol, na_weight, pinv_tol, illcond_rcond;
  int huber_maxit, n_randomizations, alt_code; 
  bool center_mean, ret_res, ret_fits, ret_stats;
};

// --- SHARED RANDOMIZATION HELPERS ---

/**
 * @brief Core Randomization Engine.
 * Generates a permutation vector based on the specified constraints and topology.
 *
 * @param rng Thread-local random number generator.
 * @param blocks List of integer vectors defining exchangeable units (Samples or Rows).
 * @param n_units Total number of units to shuffle.
 * @param is_graph_mode If true, shuffles Samples and maps to Rows. If false, shuffles Rows.
 * @param pair_mapper Lookup table for graph mode.
 * @param pairs_mat Topology matrix for graph mode.
 * @return arma::uvec A vector of indices to reorder the data matrix Y.
 */
template <class URNG>
inline arma::uvec generate_permutation(
    URNG& rng,
    const std::vector<arma::uvec>& blocks, 
    arma::uword n_units,                   
    bool is_graph_mode,                    
    PairLookup& pair_mapper,               
    const arma::umat& pairs_mat            
) {
    arma::uvec p(n_units);
    for(arma::uword i=0; i<n_units; ++i) p[i] = i;

    for(const auto& blk : blocks) {
        arma::uvec shuffled_blk = blk; 
        for (arma::uword i = shuffled_blk.n_elem; i > 1; --i) {
             std::uniform_int_distribution<arma::uword> dist(0, i - 1);
             std::swap(shuffled_blk[i-1], shuffled_blk[dist(rng)]);
        }
        p.elem(blk) = p.elem(shuffled_blk);
    }

    if (is_graph_mode) {
        arma::uvec row_perm;
        pair_mapper.get_row_perm(row_perm, pairs_mat, p);
        return row_perm;
    } 
    return p;
}

/**
 * @brief Maps global randomization blocks to specific subset indices.
 * Required when 'na_mode = drop' creates a subset of valid rows, but user constraints
 * are provided in global indices.
 */
static inline std::vector<arma::uvec> subset_blocks(const std::vector<arma::uvec>& global_blocks, 
                                                    const arma::uvec& subset_indices, 
                                                    arma::uword n_global) {
    std::vector<int> glob_to_sub(n_global, -1);
    for(arma::uword k=0; k < subset_indices.n_elem; ++k) glob_to_sub[subset_indices[k]] = k;
    
    std::vector<arma::uvec> sub_blocks;
    for(const auto& blk : global_blocks) {
        std::vector<arma::uword> sb;
        sb.reserve(blk.n_elem);
        for(arma::uword val : blk) {
            if (val < n_global && glob_to_sub[val] != -1) 
                sb.push_back((arma::uword)glob_to_sub[val]);
        }
        if (sb.size() > 0) sub_blocks.push_back(arma::uvec(sb));
    }
    return sub_blocks;
}

static inline std::mt19937_64 make_rng(std::uint64_t base, arma::uword j) {
  return std::mt19937_64(base ^ (0x9e3779b97f4a7c15ULL + j + (j<<6) + (j>>2)));
}

// --- SHARED MATH HELPERS ---

static inline bool inv_xtx_safe(const arma::mat& X, arma::mat& invXtX, arma::mat& Xt, double pinv_tol) {
  Xt = X.t();
  arma::mat XtX = arma::symmatu(Xt * X); 
  if (!XtX.is_finite()) return false;
  if (inv_sympd(invXtX, XtX)) return true;
  
  if (pinv_tol > 0.0) invXtX = arma::pinv(XtX, pinv_tol);
  else invXtX = arma::pinv(XtX);
  return invXtX.is_finite();
}

// Robust scale estimation (Median Absolute Deviation)
static inline double robust_scale_mad(const arma::vec& r) {
  arma::uvec idx = arma::find_finite(r);
  if (idx.n_elem == 0) return 1e-8;
  arma::vec rf = r.elem(idx);
  double med = arma::median(rf);
  double mad = arma::median(arma::abs(rf - med));
  return (mad > 1e-12) ? 1.4826 * mad : 1e-8;
}

static inline bool qr_basis(const arma::mat& Z, arma::mat& Q, double rank_tol = 1e-12) {
  if (Z.n_cols == 0) { Q.set_size(Z.n_rows, 0); return true; }
  arma::mat R;
  bool ok = arma::qr_econ(Q, R, Z);
  if (!ok || !Q.is_finite() || !R.is_finite()) { Q.reset(); return false; }
  arma::uword r = arma::rank(R, rank_tol);
  if (r == 0) { Q.set_size(Z.n_rows, 0); return true; }
  if (r < Q.n_cols) Q = Q.cols(0, r - 1);
  return true;
}

static inline arma::mat project_out_Q(const arma::mat& Q, const arma::mat& B) {
  if (Q.n_cols == 0) return B;
  return B - Q * (Q.t() * B);
}

static inline double z_from_p(double p, int alt, double obs, double med) {
  if (!std::isfinite(p)) return arma::datum::nan;
  if (p >= 1.0) return 0.0;
  if (p <= 0.0) p = 1e-16; 
  if (alt == 0) { // two-sided
    double z = R::qnorm(1.0 - p/2.0, 0.0, 1.0, 1, 0);
    return (obs >= med) ? z : -z;
  } else if (alt == 1) return R::qnorm(1.0 - p, 0.0, 1.0, 1, 0); // greater
  else return -R::qnorm(1.0 - p, 0.0, 1.0, 1, 0); // less
}

// Common hash for pattern grouping
static inline std::uint64_t hash_vec_mask(const arma::vec& v) {
  std::uint64_t h = 1469598103934665603ULL;
  for (const double& val : v) {
    h ^= (std::uint64_t)(arma::is_finite(val) ? 1u : 0u);
    h *= 1099511628211ULL;
  }
  return h;
}

// Helper to parse R arguments for core rows
static inline arma::uvec parse_core_rows(SEXP core_rows, arma::uword n) {
  if (Rf_isNull(core_rows)) return arma::regspace<arma::uvec>(0, n - 1);
  if (Rf_isLogical(core_rows)) {
    Rcpp::LogicalVector L(core_rows);
    std::vector<arma::uword> v; 
    for(int i=0; i<L.size(); ++i) if(L[i]) v.push_back(i);
    return arma::uvec(v);
  }
  Rcpp::IntegerVector I(core_rows);
  std::vector<arma::uword> v; 
  for(int i : I) if(i >= 1 && i <= (int)n) v.push_back(i-1);
  return arma::uvec(v);
}

#endif

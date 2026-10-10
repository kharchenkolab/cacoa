#ifndef STATS_COMMONS_H
#define STATS_COMMONS_H

#include <RcppArmadillo.h>
#include "parallel.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>


using namespace Rcpp;
using namespace arma;

// --- SHARED STRUCTS ---
struct Config {
  std::string robust, na_mode;
  double huber_k, huber_tol, na_weight, pinv_tol, illcond_rcond;
  int huber_maxit, n_randomizations, alt_code; 
  bool center_mean, ret_res, ret_fits, ret_stats;
};

// --- SHARED RANDOMIZATION HELPERS ---

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

/**
 * Induce a global permutation onto a subset of units (R: inducePermutation()).
 * `p` holds 0-based images of a design-convention permutation over n_global rows (permuted design = X[p, ]);
 * `units` are the global ids of the m units present; `blocks` index into `units` and are the (stratum x in-set)
 * cells. Within each block the members (sorted ascending) are reassigned by the ranks of their images under p;
 * units outside every block keep their place. Returns q over the units (design convention: X_units[q, ]).
 * When no unit is missing and the blocks are the cells p was drawn in, q equals p.
 */
static inline arma::uvec induced_perm(const arma::uvec& p, const arma::uvec& units, const std::vector<arma::uvec>& blocks) {
  const arma::uword m = units.n_elem;
  arma::uvec q = arma::regspace<arma::uvec>(0, (m == 0) ? 0 : m - 1);
  if (m == 0) return arma::uvec();
  std::vector<arma::uword> blk, ord;
  for (const auto& b : blocks) {
    if (b.n_elem < 2) continue;
    blk.assign(b.begin(), b.end()); std::sort(blk.begin(), blk.end());
    ord.resize(blk.size()); for (arma::uword i = 0; i < ord.size(); ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(), [&](arma::uword a, arma::uword c) { return p[units[blk[a]]] < p[units[blk[c]]]; });
    // ord[i] = member with rank i; member k takes the member with rank(k) = position of k in ord
    for (arma::uword i = 0; i < ord.size(); ++i) q[blk[ord[i]]] = blk[i];
  }
  return q;
}

/** Inverse of a permutation over m units: y[inv] is the response aligned with the permuted design X[q, ]. */
static inline arma::uvec inverse_perm(const arma::uvec& q) {
  arma::uvec inv(q.n_elem); for (arma::uword i = 0; i < q.n_elem; ++i) inv[q[i]] = i; return inv;
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

// Common hash for pattern grouping
static inline std::uint64_t hash_vec_mask(const arma::vec& v) {
  std::uint64_t h = 1469598103934665603ULL;
  for (const double& val : v) {
    h ^= (std::uint64_t)(std::isfinite(val) ? 1u : 0u);
    h *= 1099511628211ULL;
  }
  return h;
}

#endif

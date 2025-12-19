// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include "lm_common.h"

// --- LOCAL STRUCTS & HELPERS ---

struct DistDesignGroup {
  bool valid = false;
  arma::uvec valid_samples;
  std::vector<bool> is_valid_sample; 
  arma::uvec valid_pairs; 
  arma::vec sqrt_pair_weights;       
  arma::mat X_pair_sub;
  arma::mat B_pair, invXtX_pair, Xt_pair;
  arma::vec alpha; 
  arma::vec w_sample; 
  std::vector<arma::uvec> perm_blocks; 
};

struct Job { arma::uword col_idx; int group_idx; };

static std::vector<arma::uvec> filter_blocks_global(
    const std::vector<arma::uvec>& global_blocks, 
    const std::vector<bool>& is_valid) {
  
  std::vector<arma::uvec> out;
  for(const auto& blk : global_blocks) {
    std::vector<arma::uword> keep;
    keep.reserve(blk.n_elem);
    for(arma::uword idx : blk) {
      if(idx < is_valid.size() && is_valid[idx]) keep.push_back(idx);
    }
    if(keep.size() > 1) out.push_back(arma::uvec(keep));
  }
  return out;
}

static inline std::uint64_t hash_valid_rows(const arma::mat& M) {
  std::uint64_t h = 1469598103934665603ULL;
  for(uword i=0; i<M.n_rows; ++i) {
    bool fin = arma::is_finite(M(i,0));
    h ^= (std::uint64_t)(fin ? 1u : 0u);
    h *= 1099511628211ULL;
  }
  return h;
}

static arma::vec get_projection_vector(const arma::mat& X, const arma::vec& contrast) {
  arma::mat invXtX, Xt;
  if (!inv_xtx_safe(X, invXtX, Xt, 0.0)) return arma::zeros(X.n_rows);
  return X * (invXtX * contrast);
}

// Calculate distances for EXPLICIT pairs
static arma::vec calc_dist_subset(const arma::mat& M, const arma::uvec& perm, 
                                  const arma::uvec& feat_idx, const arma::umat& pairs, int d_code) {
  arma::mat M_sub = M.submat(perm, feat_idx); 
  
  // We assume M is already Gene-Centered from R if d_code=2 (Cosine) is used.
  
  uword n_pairs = pairs.n_rows;
  arma::vec dists(n_pairs);
  arma::vec norms;
  if (d_code == 2) norms = arma::sqrt(arma::sum(arma::square(M_sub), 1));
  
  for (uword k = 0; k < n_pairs; ++k) {
    uword i = pairs(k, 0);
    uword j = pairs(k, 1);
    double d = 0.0;
    
    if (d_code == 0)      d = arma::norm(M_sub.row(i) - M_sub.row(j), 2);
    else if (d_code == 1) d = arma::norm(M_sub.row(i) - M_sub.row(j), 1);
    else if (d_code == 2) {
      double dot = arma::dot(M_sub.row(i), M_sub.row(j));
      double denom = norms(i) * norms(j);
      d = 1.0 - (dot / (denom + 1e-12));
    }
    dists(k) = d;
  }
  return dists;
}

static void rank_transform_matrix(arma::mat& M) {
  for (uword j = 0; j < M.n_cols; ++j) {
    arma::vec v = M.col(j);
    arma::uvec idx = arma::sort_index(v);
    arma::vec ranks(v.n_elem);
    for(uword i=0; i<idx.n_elem; ++i) ranks(idx(i)) = i + 1.0;
    M.col(j) = ranks;
  }
}

static void impute_matrix_means(arma::mat& M, const arma::uvec& valid_samples, const std::vector<uword>& ghost_samples) {
  if (ghost_samples.empty()) return;
  if (valid_samples.n_elem == 0) return; 
  arma::rowvec means = arma::mean(M.rows(valid_samples), 0);
  for(uword ghost_idx : ghost_samples) M.row(ghost_idx) = means;
}

// [[Rcpp::export]]
Rcpp::List fit_with_focusing(
    Rcpp::List M_list,                 
    const arma::mat& sample_X,
    const arma::vec& sample_contrast,
    const arma::mat& pair_X,
    const arma::umat& pairs,           
    const arma::mat& pair_Z,           
    const arma::vec& pair_contrast,
    Rcpp::Nullable<Rcpp::List> sample_perm_groups = R_NilValue, 
    int n_randomizations = 100,
    int n_top_features = 50,
    std::string test_type = "t-test",  
    std::string dist_type = "l2",      
    std::string robust = "none",       
    double huber_k = 1.345,
    int huber_maxit = 8,
    double huber_tol = 1e-6,
    std::string na_mode = "drop",
    double na_weight = 1e-4,
    int n_cores = 1
) {
  Config cfg; 
  cfg.robust = robust; cfg.huber_k = huber_k; cfg.huber_maxit = huber_maxit; cfg.huber_tol = huber_tol;
  cfg.na_mode = na_mode; cfg.na_weight = na_weight; cfg.alt_code = 0; 
  
  int n_mats = M_list.size();
  if (n_mats == 0) return Rcpp::List::create();
  if (pair_X.n_rows != pairs.n_rows) stop("pair_X rows must match pairs rows.");
  
  std::vector<arma::mat> mats(n_mats);
  for(int i=0; i<n_mats; ++i) mats[i] = Rcpp::as<arma::mat>(M_list[i]);
  
  arma::mat M0 = mats[0];
  arma::uword n = M0.n_rows; 
  
  int d_code = (dist_type == "l1") ? 1 : (dist_type == "cosine") ? 2 : 0; 
  
  bool is_drop = (na_mode == "drop");
  
  std::unordered_map<std::uint64_t, std::vector<int>> pattern_map;
  for (int k = 0; k < n_mats; ++k) {
    if (mats[k].n_rows != n) stop("All matrices in M_list must have same number of rows.");
    pattern_map[hash_valid_rows(mats[k])].push_back(k);
  }
  
  std::vector<DistDesignGroup> groups; groups.reserve(pattern_map.size());
  std::vector<Job> jobs; jobs.reserve(n_mats);
  int grp_idx = 0;
  
  for (auto& kv : pattern_map) {
    DistDesignGroup g;
    int rep_idx = kv.second[0];
    const arma::mat& M_rep = mats[rep_idx]; 
    
    std::vector<uword> valid_s;
    std::vector<bool> is_valid(n, false);
    for(uword i=0; i<n; ++i) {
      if(arma::is_finite(M_rep(i,0))) { valid_s.push_back(i); is_valid[i] = true; }
    }
    g.valid_samples = arma::uvec(valid_s);
    g.is_valid_sample = is_valid;
    
    if (g.valid_samples.n_elem < 2) {
      g.valid = false; groups.push_back(g);
      for(int idx : kv.second) jobs.push_back({(uword)idx, grp_idx});
      grp_idx++; continue;
    }
    
    std::vector<uword> valid_p_idx;
    std::vector<double> pair_weights; 
    if (!is_drop) pair_weights.reserve(pairs.n_rows);
    
    for(uword k=0; k < pairs.n_rows; ++k) {
      bool pair_ok = (is_valid[pairs(k, 0)] && is_valid[pairs(k, 1)]);
      if (is_drop) {
        if (pair_ok) valid_p_idx.push_back(k);
      } else {
        pair_weights.push_back(pair_ok ? 1.0 : na_weight);
      }
    }
    
    if (is_drop) {
      g.valid_pairs = arma::uvec(valid_p_idx);
      if (g.valid_pairs.n_elem == 0) {
        g.valid = false; groups.push_back(g);
        for(int idx : kv.second) jobs.push_back({(uword)idx, grp_idx});
        grp_idx++; continue;
      }
      g.X_pair_sub = pair_X.rows(g.valid_pairs);
      arma::mat SX_sub = sample_X.rows(g.valid_samples);
      arma::vec w_sub = get_projection_vector(SX_sub, sample_contrast);
      g.w_sample = arma::vec(n, arma::fill::zeros);
      g.w_sample.elem(g.valid_samples) = w_sub;
    } else {
      g.X_pair_sub = pair_X; 
      g.valid_pairs = arma::regspace<arma::uvec>(0, pairs.n_rows - 1); 
      arma::vec W_p = arma::vec(pair_weights);
      g.sqrt_pair_weights = arma::sqrt(W_p); 
      g.X_pair_sub = pair_X.each_col() % g.sqrt_pair_weights; 
      arma::vec W_s(n); W_s.fill(na_weight); W_s.elem(g.valid_samples).fill(1.0);
      arma::mat SX_w = sample_X.each_col() % arma::sqrt(W_s);
      arma::mat XtX_s = SX_w.t() * SX_w;
      if(XtX_s.n_rows > 0) XtX_s.diag() += 1e-12;
      g.w_sample = (sample_X.each_col() % W_s) * solve(XtX_s, sample_contrast);
    }
    
    if (inv_xtx_safe(g.X_pair_sub, g.invXtX_pair, g.Xt_pair, 0.0)) {
      g.B_pair = g.invXtX_pair * g.Xt_pair;
      g.alpha = g.B_pair.t() * pair_contrast;
      if (pair_Z.n_elem > 0) {
        arma::mat Z_use;
        if (is_drop) Z_use = pair_Z.rows(g.valid_pairs);
        else {
          arma::vec W_p = arma::vec(pair_weights);
          Z_use = pair_Z.each_col() % arma::sqrt(W_p);
        }
        arma::mat Qz; 
        if (qr_basis(Z_use, Qz, 1e-12)) g.alpha = project_out_Q(Qz, g.alpha);
      }
      std::vector<arma::uvec> raw_blocks;
      if (sample_perm_groups.isNotNull()) {
        Rcpp::List pg(sample_perm_groups);
        for (int i=0; i<pg.size(); ++i) {
          IntegerVector bg = pg[i]; arma::uvec ug = as<arma::uvec>(bg);
          if(ug.max() > 0) ug -= 1; raw_blocks.push_back(ug);
        }
      } else raw_blocks.push_back(arma::regspace<arma::uvec>(0, n - 1));
      if (is_drop) g.perm_blocks = filter_blocks_global(raw_blocks, is_valid);
      else g.perm_blocks = raw_blocks;
      g.valid = true;
    }
    groups.push_back(g);
    for(int idx : kv.second) jobs.push_back({(uword)idx, grp_idx});
    grp_idx++;
  }
  
  arma::vec Stat(n_mats); Stat.fill(datum::nan);
  arma::vec Pval(n_mats); Pval.fill(datum::nan);
  arma::mat PermStats(n_randomizations, n_mats); PermStats.fill(datum::nan);
  arma::mat ObsDist(pairs.n_rows, n_mats); ObsDist.fill(datum::nan);
  uword n_rows_obs = ObsDist.n_rows; 
  
#ifdef _OPENMP
  if (n_cores > 1) omp_set_num_threads(n_cores);
#endif
  
#pragma omp parallel for schedule(dynamic)
  for(size_t k=0; k<jobs.size(); ++k) {
    uword mat_idx = jobs[k].col_idx;
    int g_idx = jobs[k].group_idx;
    DistDesignGroup grp = groups[g_idx]; 
    if (!grp.valid) continue;
    
    arma::mat M_local = mats[mat_idx];
    
    if (!is_drop) {
      std::vector<uword> ghosts; std::vector<uword> valids;
      for(uword i=0; i<n; ++i) {
        if(!arma::is_finite(M_local(i,0))) ghosts.push_back(i); else valids.push_back(i);
      }
      impute_matrix_means(M_local, arma::uvec(valids), ghosts);
    } else {
      M_local.elem(find_nonfinite(M_local)).zeros();
    }
    
    if (test_type == "wilcox") rank_transform_matrix(M_local);
    arma::uvec p_id = arma::regspace<arma::uvec>(0, n - 1);
    
    arma::vec scores_obs = M_local.t() * grp.w_sample;
    arma::uvec sorted_obs = arma::sort_index(arma::abs(scores_obs), "descend");
    arma::uvec top_obs = sorted_obs.head(n_top_features);
    
    arma::vec Y_full = calc_dist_subset(M_local, p_id, top_obs, pairs, d_code);
    
    if (is_drop) {
      arma::uvec global_idx = grp.valid_pairs + (mat_idx * n_rows_obs);
      ObsDist.elem(global_idx) = Y_full.elem(grp.valid_pairs); 
    } else {
      // Enforce NaN for missing pairs to trigger correct imputation in R
      for (uword r = 0; r < pairs.n_rows; ++r) {
        uword i = pairs(r, 0);
        uword j = pairs(r, 1);
        if (!grp.is_valid_sample[i] || !grp.is_valid_sample[j]) {
          ObsDist(r, mat_idx) = datum::nan; 
        } else {
          ObsDist(r, mat_idx) = Y_full(r);
        }
      }
    }
    
    arma::vec Y_use;
    if (is_drop) Y_use = Y_full.elem(grp.valid_pairs);
    else Y_use = Y_full % grp.sqrt_pair_weights; // Weighted Y
    
    double s_obs = arma::dot(grp.alpha, Y_use);
    Stat(mat_idx) = s_obs;
    
    std::mt19937_64 rng = make_rng(0xD1B54A32D192ED03ULL + mat_idx, 0);
    arma::vec perm_vec(n_randomizations);
    PairLookup dummy_map; arma::umat dummy_mat;
    
    for (int r = 0; r < n_randomizations; ++r) {
      arma::uvec p = generate_permutation(rng, grp.perm_blocks, n, false, dummy_map, dummy_mat);
      arma::vec w_perm(n); for(uword i=0; i<n; ++i) w_perm[p[i]] = grp.w_sample[i];
      
      arma::vec scores = M_local.t() * w_perm;
      arma::uvec sorted_perm = arma::sort_index(arma::abs(scores), "descend");
      arma::uvec top_idx = sorted_perm.head(n_top_features);
      
      arma::vec Y_perm_full = calc_dist_subset(M_local, p, top_idx, pairs, d_code);
      arma::vec Y_perm;
      if (is_drop) Y_perm = Y_perm_full.elem(grp.valid_pairs);
      else Y_perm = Y_perm_full % grp.sqrt_pair_weights; // Weighted Y
      
      perm_vec[r] = arma::dot(grp.alpha, Y_perm);
    }
    PermStats.col(mat_idx) = perm_vec;
    double pval = (arma::accu(arma::abs(perm_vec) >= std::abs(s_obs)) + 1.0) / (n_randomizations + 1.0);
    Pval(mat_idx) = pval;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("stat") = Stat, Rcpp::Named("p_value") = Pval,
    Rcpp::Named("perm_stats") = PermStats, Rcpp::Named("Y") = ObsDist
  );
}
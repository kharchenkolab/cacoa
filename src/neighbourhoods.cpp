// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include "parallel.h"

// kNN neighbourhood units on the cell graph (R: knnNeighbourhoodMatrix()). For every cell: itself, every graph
// neighbour (never cut to the strongest edges: those are own-arm edges when expression differs between the arms),
// then, if those are fewer than k, two-hop neighbours ranked by the number of shared neighbours (ties by index).
// `adj` is the symmetric cells x cells adjacency (any non-zero = edge; the diagonal is ignored). Returns 0-based
// index vectors, the cell first.
// [[Rcpp::export]]
Rcpp::List knn_neighbourhoods(const Eigen::SparseMatrix<double>& adj, int k, int n_cores = 1) {
  const int n = adj.rows();
  if (adj.cols() != n) Rcpp::stop("adj must be square");
  if (k < 2) Rcpp::stop("k must be at least 2");
  std::vector<std::vector<int>> nb(n);
  for (int j = 0; j < n; ++j) for (Eigen::SparseMatrix<double>::InnerIterator it(adj, j); it; ++it) { if (it.value() != 0.0 && it.row() != j) nb[it.row()].push_back(j); }
  for (int i = 0; i < n; ++i) std::sort(nb[i].begin(), nb[i].end());
  std::vector<std::vector<int>> out(n);
  cacoa::parallelFor(0, n, [&](int i) {
    std::vector<int>& o = out[i]; o.reserve(std::max((size_t)k, nb[i].size() + 1));
    o.push_back(i); o.insert(o.end(), nb[i].begin(), nb[i].end());
    const int need = k - 1 - (int)nb[i].size();
    if (need <= 0) return;
    std::unordered_map<int, int> shared; shared.reserve(nb[i].size() * 16);
    for (int j : nb[i]) for (int l : nb[j]) ++shared[l];
    std::vector<std::pair<int, int>> cand; cand.reserve(shared.size());
    std::vector<char> is1(0);
    for (auto& kv : shared) { const int l = kv.first; if (l == i || std::binary_search(nb[i].begin(), nb[i].end(), l)) continue; cand.push_back({kv.first, kv.second}); }
    if ((int)cand.size() > need) {
      std::partial_sort(cand.begin(), cand.begin() + need, cand.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) { return a.second != b.second ? a.second > b.second : a.first < b.first; });
      cand.resize(need);
    } else std::sort(cand.begin(), cand.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) { return a.second != b.second ? a.second > b.second : a.first < b.first; });
    for (auto& c : cand) o.push_back(c.first);
  }, n_cores, false);
  Rcpp::List res(n);
  for (int i = 0; i < n; ++i) res[i] = Rcpp::IntegerVector(out[i].begin(), out[i].end());
  return res;
}

// Per-sample counts of the cells in each neighbourhood: cells x samples. `sample_of_cell` is a 1-based factor.
// [[Rcpp::export]]
Rcpp::NumericMatrix neighbourhood_counts(Rcpp::List nn_ids, Rcpp::IntegerVector sample_of_cell, int n_samples, bool one_based = false) {
  const int m = nn_ids.size(), n_cells = sample_of_cell.size();
  Rcpp::NumericMatrix cnt(m, n_samples);
  for (int i = 0; i < m; ++i) {
    Rcpp::IntegerVector ids = nn_ids[i];
    for (int t = 0; t < ids.size(); ++t) { int v = one_based ? ids[t] - 1 : ids[t]; if (v < 0 || v >= n_cells) Rcpp::stop("neighbour index outside the cells");
      int s = sample_of_cell[v]; if (Rcpp::IntegerVector::is_na(s) || s < 1 || s > n_samples) continue; cnt(i, s - 1) += 1.0; }
  }
  return cnt;
}

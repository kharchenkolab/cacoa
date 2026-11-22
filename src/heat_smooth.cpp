// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include "parallel.h"

// Chebyshev polynomial of the graph Laplacian applied to a dense signal (R: chebyshevSmooth(), heatSmooth()).
// The recurrence of chebyshevSmooth(): with L' = L / (l.max / 2),
//   T0 = X, T1 = L' X - X, T_k = 2 (L' T_{k-1} - T_{k-1}) - T_{k-2},  r = 0.5 c_0 T0 + sum_{k >= 1} c_k T_k,
// repeated `n_steps` times (the heat kernel exp(-beta L / l.max) split into steps). The Laplacian is symmetric, so row i of
// L X is computed from column i of the CSC matrix: rows are independent and are split over `n_cores` threads; the signal
// is held row-major so that a row of the product is contiguous. Returns the smoothed signal (vertices x signals).
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RMat;

// [[Rcpp::export]]
Eigen::MatrixXd chebyshev_smooth_cpp(const Eigen::MappedSparseMatrix<double>& lap, const Eigen::Map<Eigen::MatrixXd>& X,
                                     double l_max, const Eigen::Map<Eigen::VectorXd>& coeffs, int n_steps, int n_cores) {
  const int n = lap.cols(), p = X.cols(), m = coeffs.size();
  if (lap.rows() != n || X.rows() != n) Rcpp::stop("chebyshev_smooth_cpp: dimension mismatch");
  if (m < 2) Rcpp::stop("chebyshev_smooth_cpp: at least two coefficients needed");
  const double sc = 2.0 / l_max;
  const int n_chunks = std::max(1, std::min(n, std::max(1, n_cores) * 8));
  std::vector<int> bounds(n_chunks + 1);
  for (int c = 0; c <= n_chunks; ++c) bounds[c] = (int)((long long)n * c / n_chunks);
  RMat cur(n, p), old(n, p), nw(n, p), r(n, p);
  RMat Xs = X;   // row-major copy of the signal
  // y = L' x - x for the rows of chunk c
  auto lapRows = [&](const RMat& src, RMat& dst, int c, double mult, const RMat* sub, double subMult) {
    for (int i = bounds[c]; i < bounds[c + 1]; ++i) {
      auto di = dst.row(i); di.setZero();
      for (Eigen::MappedSparseMatrix<double>::InnerIterator it(lap, i); it; ++it) di.noalias() += (it.value() * sc) * src.row(it.index());
      di -= src.row(i);
      di *= mult;
      if (sub) di -= subMult * sub->row(i);
    }
  };
  for (int s = 0; s < n_steps; ++s) {
    // T0 = Xs, T1 = L' Xs - Xs
    cacoa::parallelFor(0, n_chunks, [&](int c) {
      lapRows(Xs, cur, c, 1.0, nullptr, 0.0);
      for (int i = bounds[c]; i < bounds[c + 1]; ++i) { old.row(i) = Xs.row(i); r.row(i) = 0.5 * coeffs[0] * Xs.row(i) + coeffs[1] * cur.row(i); }
    }, n_cores);
    for (int k = 2; k < m; ++k) {
      const double ck = coeffs[k];
      cacoa::parallelFor(0, n_chunks, [&](int c) {
        lapRows(cur, nw, c, 2.0, &old, 1.0);      // T_k = 2 (L' T_{k-1} - T_{k-1}) - T_{k-2}
        for (int i = bounds[c]; i < bounds[c + 1]; ++i) r.row(i) += ck * nw.row(i);
      }, n_cores);
      std::swap(old, cur); std::swap(cur, nw);   // old <- T_{k-1}, cur <- T_k
    }
    Xs.swap(r);
  }
  return Eigen::MatrixXd(Xs);
}

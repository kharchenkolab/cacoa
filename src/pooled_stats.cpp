// Pooled permutation summaries of the cluster-free maps on threads (R: scanCallCounts(), pooledScanStep(), pooledMinPFDR()).
// Same definitions as the R functions; the pooled values are sorted once (chunks sorted on threads, then merged) and every
// count is a binary search in the sorted pool. Non-finite values are passed as -Inf by the R wrappers.
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "parallel.h"

namespace {
// sort v ascending on n_cores threads: chunks sorted in parallel, then merged pairwise (also in parallel)
void parallelSort(std::vector<double>& v, int n_cores) {
  const size_t n = v.size(); const int T = std::max(1, std::min(n_cores, (int)std::max<size_t>(1, n / 100000)));
  if (T <= 1) { std::sort(v.begin(), v.end()); return; }
  std::vector<size_t> b(T + 1); for (int t = 0; t <= T; ++t) b[t] = n * t / T;
  cacoa::parallelFor(0, T, [&](int t) { std::sort(v.begin() + b[t], v.begin() + b[t + 1]); }, n_cores);
  for (size_t width = 1; width < (size_t)T; width *= 2) {
    std::vector<int> starts; for (size_t t = 0; t + width < (size_t)T; t += 2 * width) starts.push_back((int)t);
    cacoa::parallelFor(0, (int)starts.size(), [&](int k) { size_t t = starts[k];
      std::inplace_merge(v.begin() + b[t], v.begin() + b[t + width], v.begin() + b[std::min(t + 2 * width, (size_t)T)]); }, n_cores);
  }
}
inline double tailValue(double x, int tail) { double y = tail == 1 ? std::fabs(x) : (tail == 2 ? -x : x); return std::isfinite(y) ? y : -INFINITY; }
inline size_t countBelow(const std::vector<double>& s, double x) { return std::lower_bound(s.begin(), s.end(), x) - s.begin(); }        // strictly below

// Exact "number of pool values strictly below x" with a bucket index over the finite range of a sorted pool: the search runs
// only within the buckets around x (one bucket of margin on each side, against rounding at the edges), so a query touches a
// few neighbouring elements instead of the whole array. Gives the same counts as std::lower_bound on the whole pool.
struct PoolIndex {   // the pool holds finite values and -Inf only (non-finite values are mapped to -Inf before sorting)
  const std::vector<double>& s; size_t first_fin, last_fin; double lo, hi, width; long nb; std::vector<size_t> start;
  PoolIndex(const std::vector<double>& pool) : s(pool) {
    first_fin = 0; while (first_fin < s.size() && s[first_fin] == -INFINITY) ++first_fin;
    last_fin = s.size();
    nb = (long)std::max<size_t>(1, std::min<size_t>(s.size() / 8, (size_t)1 << 22)); start.assign(nb + 1, last_fin);
    lo = hi = 0; width = 1;
    if (last_fin <= first_fin) return;
    lo = s[first_fin]; hi = s[last_fin - 1]; width = (hi - lo) / nb; if (!(width > 0)) width = 1;
    size_t j = first_fin;
    for (long b = 0; b <= nb; ++b) { const double e = lo + b * width; while (j < last_fin && s[j] < e) ++j; start[b] = j; }
    start[0] = first_fin; start[nb] = last_fin;
  }
  size_t below(double x) const {   // number of pool values strictly below x (x finite or -Inf)
    if (x == -INFINITY) return 0;
    if (last_fin <= first_fin || x < lo) return first_fin;
    if (x > hi) return last_fin;
    long b = (long)std::floor((x - lo) / width); if (b < 0) b = 0; if (b >= nb) b = nb - 1;
    const size_t a = start[std::max<long>(b - 1, 0)], z = (b + 2 <= nb) ? start[b + 2] : last_fin;
    return std::lower_bound(s.begin() + a, s.begin() + z, x) - s.begin();
  }
};
inline size_t countAtMost(const std::vector<double>& s, double x) { return std::upper_bound(s.begin(), s.end(), x) - s.begin(); }      // at or below
}

// scanCallCounts(): S is the (B + 1) x N matrix of tail-transformed (and, if asked, standardized) statistics, row 1 observed,
// non-finite as -Inf. Pool = all (B + 1) N values. For each data set (row): its finite values in decreasing order m_i, the pool
// count strictly below each, cnt_i, and its own count at or above each, R_i; fdr_i = (Nn - cnt_i - R_i) / B / R_i; the number
// called is R at the last position with fdr < alpha (0 if none).
// [[Rcpp::export]]
Rcpp::IntegerVector scan_call_counts_cpp(const Rcpp::NumericMatrix& S, double alpha, int n_cores) {
  const int nr = S.nrow(), N = S.ncol(), B = nr - 1; const size_t Nn = (size_t)nr * N;
  std::vector<double> pool(S.begin(), S.end()); parallelSort(pool, n_cores); const PoolIndex idx(pool);
  Rcpp::IntegerVector out(nr); std::vector<int> res(nr, 0);
  const double* sp = S.begin();
  cacoa::parallelFor(0, nr, [&](int r) {
    std::vector<double> m; m.reserve(N);
    for (int c = 0; c < N; ++c) { double x = sp[r + (size_t)c * nr]; if (std::isfinite(x)) m.push_back(x); }
    if (m.empty()) { res[r] = 0; return; }
    std::sort(m.begin(), m.end(), std::greater<double>());
    const size_t L = m.size(); int called = 0; size_t i = 0;
    while (i < L) {                                   // tie block [i, j): all share R = j (own values at or above)
      size_t j = i + 1; while (j < L && m[j] == m[i]) ++j;
      const double R = (double)j; const double cnt = (double)idx.below(m[i]);
      const double fdr = ((double)Nn - cnt - R) / B / R;
      if (fdr < alpha) called = (int)j;              // the last position with fdr < alpha (in decreasing order of m) wins
      i = j;
    }
    res[r] = called;
  }, n_cores);
  for (int r = 0; r < nr; ++r) out[r] = res[r];
  return out;
}

// pooledScanStep(): pool = the B x N permuted values after the tail transform (tail 0 identity, 1 abs, 2 negate; non-finite
// as -Inf). Per permuted value pp = (#pool >= x + 1) / (Nn + 1); per observed value s, p = (#pool >= s + 1) / (Nn + 1) and
// ge = #pool >= s, which equals the number of pp at or below p (pp and p are the same decreasing function of the value), the
// count pooledMinPFDR(p, pp) takes from a second sort.
// [[Rcpp::export]]
Rcpp::List pooled_scan_step_cpp(const Rcpp::NumericVector& s_obs, const Rcpp::NumericMatrix& sp, int tail, int n_cores) {
  const size_t Nn = (size_t)sp.nrow() * sp.ncol(); const double* ap = sp.begin();
  std::vector<double> pool(Nn);
  const int nchunk = std::max(1, n_cores * 8); const size_t step = (Nn + nchunk - 1) / nchunk;
  cacoa::parallelFor(0, nchunk, [&](int k) { const size_t lo = k * step, hi = std::min(Nn, lo + step); for (size_t i = lo; i < hi; ++i) pool[i] = tailValue(ap[i], tail); }, n_cores);
  parallelSort(pool, n_cores); const PoolIndex idx(pool);
  Rcpp::NumericMatrix pp(sp.nrow(), sp.ncol()); double* ppp = pp.begin();
  cacoa::parallelFor(0, nchunk, [&](int k) { const size_t lo = k * step, hi = std::min(Nn, lo + step);
    for (size_t i = lo; i < hi; ++i) ppp[i] = ((double)(Nn - idx.below(tailValue(ap[i], tail))) + 1.0) / ((double)Nn + 1.0); }, n_cores);
  const R_xlen_t N = s_obs.size(); Rcpp::NumericVector p(N), ge(N);
  for (R_xlen_t i = 0; i < N; ++i) { const double g = (double)(Nn - idx.below(tailValue(s_obs[i], tail))); ge[i] = g; p[i] = (g + 1.0) / ((double)Nn + 1.0); }
  return Rcpp::List::create(Rcpp::Named("p") = p, Rcpp::Named("pp") = pp, Rcpp::Named("ge") = ge);
}

// scanCallCounts() without standardization: the (B + 1) x N matrix is built here from the observed and permuted statistics (tail
// transform, non-finite as -Inf, a cell whose observed value is not finite set to -Inf in every data set), without R copies.
// [[Rcpp::export]]
Rcpp::IntegerVector scan_call_counts_raw_cpp(const Rcpp::NumericVector& s_obs, const Rcpp::NumericMatrix& sp, int tail, double alpha, int n_cores) {
  const int B = sp.nrow(), N = sp.ncol(), nr = B + 1;
  Rcpp::NumericMatrix S(nr, N); double* Sp = S.begin(); const double* ap = sp.begin();
  cacoa::parallelFor(0, N, [&](int c) {
    const double o = tail == 1 ? std::fabs(s_obs[c]) : (tail == 2 ? -s_obs[c] : s_obs[c]); const bool bad = !std::isfinite(o);
    Sp[(size_t)c * nr] = bad ? -INFINITY : o;
    for (int r = 0; r < B; ++r) Sp[(size_t)c * nr + r + 1] = bad ? -INFINITY : tailValue(ap[r + (size_t)c * B], tail);
  }, n_cores);
  return scan_call_counts_cpp(S, alpha, n_cores);
}

// pooledMinPFDR() counts: per observed value, the number of permuted values at or below it (cnt) and of observed values at or
// below it (R); the R wrapper forms the FDR and its monotone version.
// [[Rcpp::export]]
Rcpp::List pooled_min_counts_cpp(const Rcpp::NumericVector& m_obs, const Rcpp::NumericVector& m_perm, int n_cores) {
  std::vector<double> pool(m_perm.begin(), m_perm.end()); parallelSort(pool, n_cores);
  std::vector<double> obs(m_obs.begin(), m_obs.end()); std::sort(obs.begin(), obs.end());
  const R_xlen_t N = m_obs.size(); Rcpp::NumericVector cnt(N), R(N);
  for (R_xlen_t i = 0; i < N; ++i) { cnt[i] = (double)countAtMost(pool, m_obs[i]); R[i] = (double)countAtMost(obs, m_obs[i]); }
  return Rcpp::List::create(Rcpp::Named("cnt") = cnt, Rcpp::Named("R") = R, Rcpp::Named("Nn") = (double)pool.size());
}

// scanCallCounts() with standardization (the shift map): S = rbind(s_obs, sp) (no tail transform), each cell standardized over
// its B + 1 values with R's arithmetic (colMeans sums in long double; sd = sqrt(max(mean(x^2) - mu^2, 0) (B + 1) / B)); a cell
// whose observed value or sd is not finite or whose sd is not positive set to -Inf in every data set; non-finite to -Inf.
// [[Rcpp::export]]
Rcpp::IntegerVector scan_call_counts_std_cpp(const Rcpp::NumericVector& s_obs, const Rcpp::NumericMatrix& sp, double alpha, int n_cores) {
  const int B = sp.nrow(), N = sp.ncol(), nr = B + 1;
  Rcpp::NumericMatrix S(nr, N); double* Sp = S.begin(); const double* ap = sp.begin();
  cacoa::parallelFor(0, N, [&](int c) {
    const double o = s_obs[c]; long double s1 = o, s2 = (long double)(o * o);
    for (int r = 0; r < B; ++r) { const double x = ap[r + (size_t)c * B]; s1 += x; s2 += (long double)(x * x); }
    const double mu = (double)(s1 / nr), m2 = (double)(s2 / nr);
    const double sdv = std::sqrt(std::max(m2 - mu * mu, 0.0) * (double)nr / (double)B);
    const bool bad = !std::isfinite(o) || !std::isfinite(sdv) || !(sdv > 0);
    double* col = Sp + (size_t)c * nr;
    if (bad) { for (int r = 0; r < nr; ++r) col[r] = -INFINITY; return; }
    const double z0 = (o - mu) / sdv; col[0] = std::isfinite(z0) ? z0 : -INFINITY;
    for (int r = 0; r < B; ++r) { const double z = (ap[r + (size_t)c * B] - mu) / sdv; col[r + 1] = std::isfinite(z) ? z : -INFINITY; }
  }, n_cores);
  return scan_call_counts_cpp(S, alpha, n_cores);
}

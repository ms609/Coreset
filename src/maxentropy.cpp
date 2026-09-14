// maxentropy.cpp
//
// Selectors for the maximum-entropy (maxdet) subset problem: given a similarity
// kernel K (n x n, symmetric, positive-semidefinite) and a target size k, find
// the k-subset S maximising log det K_S -- the spanned volume of the selection,
// the Shewry & Wynn (1987) maximum-entropy sampling criterion and the MAP mode
// of a determinantal point process. The PSD repair is orchestrated in R
// (R/maxentropy.R) over the partial eigendecomposition in src/symeigen.cpp;
// the selectors here take the already-repaired K.
//
// Two selectors, mirroring the rest of the package's heuristic/exact split:
//   * MaxEntropyGreedy_cpp -- greedy pivoted Cholesky (add the point of largest
//     residual conditional variance == largest log-det increment). O(n k^2).
//     The exact argmax is NP-hard, so this is the workhorse.
//   * MaxEntropyExact_cpp  -- exact enumeration of all k-subsets, scoring each
//     by a Cholesky log-determinant; for small instances only (the R wrapper
//     gates on choose(n, k)).
// Both are deterministic and return 1-based indices. The O(n^2) helpers the
// wrapper needs (default bandwidth, RBF kernel, distinct-row count) follow.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cstdint>
using namespace Rcpp;

// Greedy maximum-entropy selection by pivoted Cholesky.
//
// At step t the residual conditional variance of each not-yet-selected point i
// is d[i] = K_ii - sum_{s<t} L[i,s]^2; adding the point of largest d[i] is the
// largest possible log-det increment. The first pivot of a unit-diagonal kernel
// is tied across all points, so a deterministic seed (`seed`, 1-based; pass 0
// or less to fall back to the diagonal argmax) breaks it -- the R wrapper
// supplies the most peripheral point (least total similarity). A redundant
// point has d -> 0 and is never preferred, which is the det -> 0 density-blind
// property. Returns the pick order (length min(k, n), 1-based).
//
// [[Rcpp::export]]
IntegerVector MaxEntropyGreedy_cpp(const NumericMatrix& K, int k, int seed) {
  const int n = K.nrow();
  if (k > n) k = n;
  if (k < 1) return IntegerVector(0);

  // Rows live in physical slots: the still-available points occupy slots
  // 0..nAvail-1 and each pivot is swapped into the last available slot as it
  // is chosen, so every sweep below touches only the rows it can still
  // change. orig[slot] is the original (0-based) point index.
  std::vector<double> d(n);
  std::vector<int> orig(n);
  for (int i = 0; i < n; ++i) { d[i] = K(i, i); orig[i] = i; }
  // L holds the first k Cholesky columns, column-major (column t contiguous)
  // and indexed by slot, so the O(nAvail t) dot products against pivot row j
  // at step t run as vectorisable column sweeps: prev[row] += L[row, s] *
  // L[j, s] with s ascending, the summation order of a per-row dot product,
  // so the factor does not depend on this layout or on the slot order.
  std::vector<double> L(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> prev(n);
  IntegerVector perm(k);
  int nAvail = n;

  for (int t = 0; t < k; ++t) {
    int pj;                               // pivot slot
    if (t == 0 && seed >= 1 && seed <= n) {
      pj = seed - 1;                      // seed is 1-based; slots are still in
    } else {                              // original order at t == 0
      // First argmax of d over available points (ties -> lowest original
      // index, as which.max over the original rows).
      pj = -1;
      double best = R_NegInf;
      int bestOrig = n;
      for (int s = 0; s < nAvail; ++s) {
        if (d[s] > best || (d[s] == best && orig[s] < bestOrig)) {
          best = d[s]; pj = s; bestOrig = orig[s];
        }
      }
    }
    const int j = orig[pj];
    perm[t] = j + 1;                      // store 1-based
    const int last = --nAvail;            // the pivot's slot from now on
    if (pj != last) {
      std::swap(d[pj], d[last]);
      std::swap(orig[pj], orig[last]);
      for (int s = 0; s < t; ++s) {
        double* Ls = L.data() + static_cast<size_t>(s) * n;
        std::swap(Ls[pj], Ls[last]);
      }
    }

    const double Ljt = std::sqrt(d[last] > 0.0 ? d[last] : 0.0);
    double* Lt = L.data() + static_cast<size_t>(t) * n;
    Lt[last] = Ljt;
    if (Ljt > 0.0) {
      std::fill(prev.begin(), prev.begin() + nAvail, 0.0);
      for (int s = 0; s < t; ++s) {
        const double* Ls = L.data() + static_cast<size_t>(s) * n;
        const double Ljs = Ls[last];
        for (int row = 0; row < nAvail; ++row) prev[row] += Ls[row] * Ljs;
      }
      const double* Kj = &K(0, j);
      for (int row = 0; row < nAvail; ++row) {
        const double val = (Kj[orig[row]] - prev[row]) / Ljt;
        Lt[row] = val;
        d[row] -= val * val;
        if (d[row] < 0.0) d[row] = 0.0;
      }
    }
  }
  return perm;
}

// Log-determinant of the k x k submatrix K[idx, idx] by Cholesky, or R_NegInf
// if it is not positive-definite (a near-duplicate / collinear subset, which a
// max-log-det search rejects anyway). Left-looking by column: entry (i, j) is
// K(i, j) minus L[i, s] * L[j, s] for s ascending, then divided by L[j, j] --
// the per-entry operation order of a row-by-row factorisation -- but the
// column's entries advance together, so each s-sweep is a contiguous,
// vectorisable pass over rows j..m-1.
static double SubLogDet(const NumericMatrix& K, const std::vector<int>& idx) {
  const int m = static_cast<int>(idx.size());
  std::vector<double> Lc(static_cast<size_t>(m) * m, 0.0);  // column-major
  std::vector<double> col(m);
  double logdet = 0.0;
  for (int j = 0; j < m; ++j) {
    for (int i = j; i < m; ++i) col[i] = K(idx[i], idx[j]);
    for (int s = 0; s < j; ++s) {
      const double* Ls = Lc.data() + static_cast<size_t>(s) * m;
      const double Ljs = Ls[j];
      for (int i = j; i < m; ++i) col[i] -= Ls[i] * Ljs;
    }
    if (col[j] <= 0.0) return R_NegInf;                // not positive-definite
    const double diag = std::sqrt(col[j]);
    double* Lj = Lc.data() + static_cast<size_t>(j) * m;
    Lj[j] = diag;
    logdet += 2.0 * std::log(diag);
    for (int i = j + 1; i < m; ++i) Lj[i] = col[i] / diag;
  }
  return logdet;
}

// Exact maximum-entropy selection by enumerating all k-subsets in lexicographic
// order and keeping the first attaining the maximum log det (strict `>`, so the
// tie-break matches the greedy/utils::combn order). Returns the optimal subset,
// sorted ascending, 1-based. The R wrapper restricts this to small choose(n, k).
//
// [[Rcpp::export]]
IntegerVector MaxEntropyExact_cpp(const NumericMatrix& K, int k) {
  const int n = K.nrow();
  if (k < 1 || k > n) stop("`k` must satisfy 1 <= k <= nrow(K)");

  std::vector<int> idx(k);
  for (int i = 0; i < k; ++i) idx[i] = i;        // 0-based lexicographic start
  std::vector<int> best(idx);
  double bestVal = R_NegInf;

  while (true) {
    const double v = SubLogDet(K, idx);
    if (v > bestVal) { bestVal = v; best = idx; }
    // Advance to the next combination in lexicographic order.
    int i = k - 1;
    while (i >= 0 && idx[i] == n - k + i) --i;
    if (i < 0) break;
    ++idx[i];
    for (int j = i + 1; j < k; ++j) idx[j] = idx[j - 1] + 1;
  }

  IntegerVector out(k);
  for (int i = 0; i < k; ++i) out[i] = best[i] + 1;   // sorted, 1-based
  return out;
}

// Cholesky log-determinant of the kernel block K[idx, idx] (idx 1-based) -- the
// score the exact selector maximises, exposed so the R wrapper reports the same
// quantity it selects by (one honest notion of log det). Returns R_NegInf for a
// non-positive-definite block: a degenerate selection containing duplicate or
// collinear points, whose Gram matrix is singular.
//
// [[Rcpp::export]]
double MaxEntropyLogDet_cpp(const NumericMatrix& K, const IntegerVector& idx) {
  const int n = K.nrow();
  std::vector<int> id(static_cast<size_t>(idx.size()));
  for (R_xlen_t i = 0; i < idx.size(); ++i) {
    const int v = idx[i];
    if (v < 1 || v > n) stop("`idx` entries must be in 1..nrow(K)");
    id[static_cast<size_t>(i)] = v - 1;               // 1-based -> 0-based
  }
  return SubLogDet(K, id);
}

// ----- kernel construction --------------------------------------------------

// Median of the positive entries of the strict upper triangle of a symmetric
// matrix -- the default RBF bandwidth. Equals stats::median(d[d > 0]): each
// off-diagonal value appears twice in the whole matrix, and the median of a
// doubled multiset is the median of the multiset; the mean of the two middle
// values is accumulated in long double as R's mean() does. A selection, not a
// sort. NA when no entry is positive.
//
// [[Rcpp::export]]
double MedianPositiveUpper_cpp(const NumericMatrix& d) {
  const int n = d.nrow();
  if (d.ncol() != n) stop("`d` must be square");
  std::vector<double> v;
  v.reserve(static_cast<size_t>(n) * (n > 0 ? n - 1 : 0) / 2);
  for (int j = 1; j < n; ++j) {
    for (int i = 0; i < j; ++i) {
      const double x = d(i, j);
      if (x > 0.0) v.push_back(x);
    }
  }
  const size_t m = v.size();
  if (m == 0) return NA_REAL;
  const size_t half = (m + 1) / 2 - 1;               // lower middle, 0-based
  std::nth_element(v.begin(), v.begin() + half, v.end());
  const double lo = v[half];
  if (m % 2 == 1) return lo;
  const double hi = *std::min_element(v.begin() + half + 1, v.end());
  const long double sum = static_cast<long double>(lo) + static_cast<long double>(hi);
  return static_cast<double>(sum / 2.0L);
}

// RBF kernel exp(-d^2 / (2 sigma^2)) of a symmetric distance matrix, computed
// on the upper triangle and mirrored. Bit-identical to the R expression.
//
// [[Rcpp::export]]
NumericMatrix RbfKernel_cpp(const NumericMatrix& d, double sigma) {
  const int n = d.nrow();
  if (d.ncol() != n) stop("`d` must be square");
  NumericMatrix k(n, n);
  const double denom = 2.0 * (sigma * sigma);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i <= j; ++i) {
      const double x = d(i, j);
      const double v = std::exp(-(x * x) / denom);
      k(i, j) = v;
      k(j, i) = v;
    }
  }
  return k;
}

// Number of distinct rows of a numeric matrix, as sum(!duplicated(d)) counts
// them: rows are equal when every entry is == (so 0 and -0 agree; entries are
// finite here, validated upstream). Rows are bucketed by a hash of their bit
// patterns (-0 normalised to 0) and compared exactly within a bucket.
//
// [[Rcpp::export]]
int DistinctRows_cpp(const NumericMatrix& d) {
  const int n = d.nrow(), m = d.ncol();
  if (n == 0) return 0;
  std::vector<uint64_t> h(n, 1469598103934665603ULL);          // FNV-1a
  for (int j = 0; j < m; ++j) {
    const double* col = &d(0, j);
    for (int i = 0; i < n; ++i) {
      double x = col[i];
      if (x == 0.0) x = 0.0;                                      // -0 -> +0
      uint64_t bits;
      std::memcpy(&bits, &x, sizeof bits);
      h[i] = (h[i] ^ bits) * 1099511628211ULL;
    }
  }
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) { return h[a] < h[b]; });
  const auto rowsEqual = [&](int a, int b) {
    for (int j = 0; j < m; ++j) if (d(a, j) != d(b, j)) return false;
    return true;
  };
  int distinct = 0;
  std::vector<int> reps;
  for (int a = 0; a < n; ) {
    int b = a;
    while (b < n && h[ord[b]] == h[ord[a]]) ++b;
    reps.clear();
    for (int q = a; q < b; ++q) {
      bool dup = false;
      for (const int r : reps) if (rowsEqual(ord[q], r)) { dup = true; break; }
      if (!dup) { reps.push_back(ord[q]); ++distinct; }
    }
    a = b;
  }
  return distinct;
}

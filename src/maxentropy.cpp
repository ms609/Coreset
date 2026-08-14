// maxentropy.cpp
//
// Selectors for the maximum-entropy (maxdet) subset problem: given a similarity
// kernel K (n x n, symmetric, positive-semidefinite) and a target size k, find
// the k-subset S maximising log det K_S -- the spanned volume of the selection,
// the Shewry & Wynn (1987) maximum-entropy sampling criterion and the MAP mode
// of a determinantal point process. The kernel construction and its PSD repair
// live in R (R/maxentropy.R); these kernels take the already-repaired K.
//
// Two routines, mirroring the rest of the package's heuristic/exact split:
//   * MaxEntropyGreedy_cpp -- greedy pivoted Cholesky (add the point of largest
//     residual conditional variance == largest log-det increment). O(n k^2).
//     The exact argmax is NP-hard, so this is the workhorse.
//   * MaxEntropyExact_cpp  -- exact enumeration of all k-subsets, scoring each
//     by a Cholesky log-determinant; for small instances only (the R wrapper
//     gates on choose(n, k)).
//
// Both are deterministic and return 1-based indices.

#include <Rcpp.h>
#include <vector>
#include <cmath>
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

  std::vector<double> d(n);
  for (int i = 0; i < n; ++i) d[i] = K(i, i);
  // L[i] holds the first k Cholesky entries for row i (only columns < t used).
  std::vector< std::vector<double> > L(n, std::vector<double>(k, 0.0));
  std::vector<bool> avail(n, true);
  IntegerVector perm(k);

  for (int t = 0; t < k; ++t) {
    int j;
    if (t == 0 && seed >= 1 && seed <= n) {
      j = seed - 1;                       // seed is 1-based; ignored if out of range
    } else {
      // First argmax of d over available points (ties -> first, as which.max).
      j = -1;
      double best = R_NegInf;
      for (int i = 0; i < n; ++i) {
        if (avail[i] && d[i] > best) { best = d[i]; j = i; }
      }
    }
    perm[t] = j + 1;                      // store 1-based
    avail[j] = false;

    const double Ljt = std::sqrt(d[j] > 0.0 ? d[j] : 0.0);
    L[j][t] = Ljt;
    if (Ljt > 0.0) {
      for (int row = 0; row < n; ++row) {
        if (!avail[row]) continue;
        double prev = 0.0;
        for (int s = 0; s < t; ++s) prev += L[row][s] * L[j][s];
        const double val = (K(row, j) - prev) / Ljt;
        L[row][t] = val;
        d[row] -= val * val;
        if (d[row] < 0.0) d[row] = 0.0;
      }
    }
  }
  return perm;
}

// Log-determinant of the k x k submatrix K[idx, idx] by Cholesky, or R_NegInf
// if it is not positive-definite (a near-duplicate / collinear subset, which a
// max-log-det search rejects anyway).
static double SubLogDet(const NumericMatrix& K, const std::vector<int>& idx) {
  const int m = static_cast<int>(idx.size());
  std::vector<double> Lc(static_cast<size_t>(m) * m, 0.0);  // lower triangle
  double logdet = 0.0;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = K(idx[i], idx[j]);
      for (int s = 0; s < j; ++s) {
        sum -= Lc[static_cast<size_t>(i) * m + s] * Lc[static_cast<size_t>(j) * m + s];
      }
      if (i == j) {
        if (sum <= 0.0) return R_NegInf;               // not positive-definite
        const double diag = std::sqrt(sum);
        Lc[static_cast<size_t>(i) * m + j] = diag;
        logdet += 2.0 * std::log(diag);
      } else {
        Lc[static_cast<size_t>(i) * m + j] = sum / Lc[static_cast<size_t>(j) * m + j];
      }
    }
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

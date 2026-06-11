#include <Rcpp.h>
#include <cmath>
#include <vector>

// Coordinate-based greedy furthest-point (maximin) selection.
//
// This is the O(N*k*dim) time, O(N) memory counterpart of MaximinFrom_cpp
// (src/maximin.cpp), which reads a fully materialised N x N distance matrix.
// Each greedy step needs only the column d(., best) of distances from the
// most-recently-selected point; for Euclidean data that column is recomputed
// from coordinates in O(N*dim), so the dense matrix is never built.
//
// Identical selection indices to the matrix path are required. stats::dist()
// computes each Euclidean distance as
//   sqrt( sum_j (x[i,j] - x[c,j])^2 )
// accumulating dev*dev in a plain double over dimensions in increasing column
// order (R's R_euclidean in src/library/stats/src/distance.c).
//
// The greedy pass below makes every decision through which.max(min_dist) and an
// elementwise pmin; both are monotone under sqrt, so the pass runs in SQUARED-
// distance space (EuclidColSq, no per-element sqrt -- the dominant inner-loop
// cost at low dim) and still selects the same indices. The squared accumulation
// `sum_j dev*dev` is the exact argument R passes to sqrt, so the ordering of the
// squared values matches the ordering of the matrix entries: argmax and pmin
// resolve to the same indices, and the reported T_k (sqrt'd once at the end) is
// bit-identical to MinDist(). This identity is re-verified over a broad N x dim
// x n x strategy battery by dev/profiling/drivers/verify.R (1620/1620 cases) in
// addition to the cross-path expect_identical() tests in test-farfirst.R.
//
// The seed-anchor primitives further down (RowSums/RowSqSums/Diameter/EuclidCol)
// still return true sqrt distances, bit-matching as.matrix(dist(points)) entry
// for entry, since their callers compare those values directly, not just rank.
//
// FP-flag sensitivity: the only way EuclidCol can diverge from R_euclidean is if
// one of them contracts `s + dev*dev` into a single-rounding FMA and the other
// does not (sqrt is correctly rounded everywhere). A user ~/.R/Makevars with
// -march=native / -Ofast can introduce exactly this. Identity must therefore be
// re-confirmed on each deployment toolchain by running the test suite; do NOT
// pre-emptively force -ffp-contract=off (it would itself cause divergence if R's
// stats was built with FMA).
//
// Arguments mirror .MaximinFromPoints() in R/samplers.R:
//   points  N x dim coordinate matrix (column-major double storage)
//   n       number of points to select (1 <= n <= N)
//   first   1-based index of the first selected point
//   mask    1-based index of a forbidden point (e.g. the anti-medoid's medoid)
//           whose min_dist is pinned to -Inf so it is never selected and never
//           updates any other point; 0 means "no mask".
//
// Returns an integer vector of length n (1-based indices, same as R).

// Squared-Euclidean accumulated in double over columns, matching dist's order.
static inline double EuclidCol(const double* P, int nPts, int dim,
                               int i, int c) {
  double s = 0.0;
  for (int j = 0; j < dim; j++) {
    double dev = P[i + (R_xlen_t)j * nPts] - P[c + (R_xlen_t)j * nPts];
    s += dev * dev;
  }
  return std::sqrt(s);
}

// Squared-Euclidean column for the greedy pass, WITHOUT the final sqrt. The
// pass makes its decisions purely through which.max(min_dist) and the pmin
// update, both monotone under sqrt, so running it in squared-distance space
// yields the same selection while skipping one sqrt per (point, step) -- the
// dominant cost of the inner loop at low dim. The reported T_k is sqrt()'d once
// at the end. (The other primitives below still need true distances; they keep
// EuclidCol, with the sqrt.)
//
// Fill dcol[i] = squared distance from point i to centre c, for all i, by
// accumulating dimension-by-dimension over contiguous columns. Summing dev*dev
// over j ascending matches stats::dist()'s order, so dcol[i] equals the squared
// matrix entry -- but every access is sequential: one vectorisable pass per
// dimension over the contiguous column `P + j*nPts`, instead of `dim` strided
// streams gathered per point. This is the cache/SIMD-friendly shape of the hot
// column.
static inline void EuclidColSqInto(const double* P, int nPts, int dim, int c,
                                   double* dcol) {
  for (int i = 0; i < nPts; i++) dcol[i] = 0.0;
  for (int j = 0; j < dim; j++) {
    const double* pj = P + (R_xlen_t)j * nPts;     // contiguous column j
    double cj = pj[c];
    for (int i = 0; i < nPts; i++) {
      double dev = pj[i] - cj;
      dcol[i] += dev * dev;
    }
  }
}

// [[Rcpp::export]]
Rcpp::IntegerVector MaximinFromPoints_cpp(Rcpp::NumericMatrix points,
                                          int n, int first, int mask) {
  int nPts = points.nrow();
  int dim  = points.ncol();
  if (n < 1 || n > nPts) {                  // defensive: public wrapper guards
    Rcpp::stop("'n' must be in [1, %d]; got %d", nPts, n);
  }
  if (first < 1 || first > nPts) {
    Rcpp::stop("'first' must be in [1, %d]; got %d", nPts, first);
  }
  if (mask < 0 || mask > nPts) {                    // # nocov start
    Rcpp::stop("'mask' must be in [0, %d]; got %d", nPts, mask);
  }                                                  // # nocov end
  const double* P = points.begin();

  Rcpp::IntegerVector selected(n);
  selected[0] = first;
  int first0 = first - 1;           // convert to 0-based once

  // min_dist holds SQUARED nearest-selected distances throughout the pass; the
  // selection is identical to working in true distances (sqrt is monotone) but
  // skips n*N sqrt calls. T_k is recovered with a single sqrt at the end.
  Rcpp::NumericVector min_dist(nPts);
  double* md = min_dist.begin();
  EuclidColSqInto(P, nPts, dim, first0, md);
  std::vector<double> dcol(nPts);   // reused per step for the new point's column
  min_dist[first0] = R_NegInf;      // mask seed before entering loop
  if (mask >= 1) {                                   // # nocov start
    min_dist[mask - 1] = R_NegInf;  // pin forbidden point (anti-medoid medoid)
  }                                                  // # nocov end

  // T_k = min over greedy steps of the chosen point's insertion distance, which
  // is exactly best_val at each step. Tracked in squared space, sqrt'd once.
  double tk_sq = R_PosInf;

  for (int k = 1; k < n; k++) {
    // which.max: first index of the global maximum (strict >, so ties -> first)
    int best = 0;
    double best_val = min_dist[0];
    for (int i = 1; i < nPts; i++) {
      if (min_dist[i] > best_val) {
        best_val = min_dist[i];
        best = i;
      }
    }
    selected[k] = best + 1;         // back to 1-based
    if (best_val < tk_sq) tk_sq = best_val;   // running min insertion distance

    // Mask new point before the pmin update so its self-distance (0) cannot
    // overwrite -Inf. A pinned `mask` point keeps -Inf for the same reason:
    // any non-negative distance is never < -Inf.
    min_dist[best] = R_NegInf;
    EuclidColSqInto(P, nPts, dim, best, dcol.data());
    for (int i = 0; i < nPts; i++) {
      if (dcol[i] < md[i]) md[i] = dcol[i];
    }
  }

  // T_k (min pairwise distance of the selection), computed for free; NA for
  // n < 2. The ensemble driver reads this attribute instead of re-running a
  // full stats::dist() over the selection (see .GonzEnsembleFromPoints).
  selected.attr("t_k") = (n >= 2 && R_finite(tk_sq))
    ? std::sqrt(tk_sq < 0.0 ? 0.0 : tk_sq) : NA_REAL;

  // Return:
  return selected;
}

// Row sums of the on-the-fly Euclidean distance matrix, for the 1-median
// (medoid = which.min(rowSums(d))) used to seed WideSampleMedoidFirst and
// WideSampleAntiMedoid. This is the one coordinate path that stays O(N^2*dim)
// in time (the rowSums are intrinsically all-pairs); memory is O(N).
//
// Bit-exact match to which.min(rowSums(as.matrix(dist(points)))): base R's
// rowSums accumulates each row in a `long double` over columns in increasing
// index, with the zero diagonal included. We mirror that — long double
// accumulator, j ascending, EuclidCol(i, i) == 0 contributes a harmless 0 —
// so the returned doubles equal R's and which.min picks the same index.
//
// [[Rcpp::export]]
Rcpp::NumericVector RowSumsFromPoints_cpp(Rcpp::NumericMatrix points) {
  int nPts = points.nrow();
  int dim  = points.ncol();
  const double* P = points.begin();
  Rcpp::NumericVector out(nPts);
  for (int i = 0; i < nPts; i++) {
    long double s = 0.0L;
    for (int j = 0; j < nPts; j++) {
      s += (long double) EuclidCol(P, nPts, dim, i, j);
    }
    out[i] = (double) s;
  }
  // Return:
  return out;
}

// Row sums of the SQUARED on-the-fly Euclidean distance matrix, for the
// "rownorm" anchor of WideSampleGonzEnsemble (seed = which.max(rowSums(d^2))).
// Like RowSumsFromPoints_cpp this stays O(N^2*dim) in time, O(N) in memory.
//
// Bit-exact match to rowSums(d^2) where d = as.matrix(dist(points)): R forms the
// d^2 matrix by squaring each double entry, then rowSums accumulates each row in
// a long double over columns in increasing index (zero diagonal included). We
// mirror that exactly — compute d_ij = EuclidCol(i, j) as the same double as the
// matrix entry, square it in double (d_ij * d_ij, matching R's `d^2`), and
// accumulate the row in long double with j ascending. A closed-form
// sum_j (N*||x_i||^2 - 2 x_i.x_j + ||x_j||^2) is NOT used: it rounds differently
// from squaring the sqrt distances and would flip which.max on tie-dense data.
//
// [[Rcpp::export]]
Rcpp::NumericVector RowSqSumsFromPoints_cpp(Rcpp::NumericMatrix points) {
  int nPts = points.nrow();
  int dim  = points.ncol();
  const double* P = points.begin();
  Rcpp::NumericVector out(nPts);
  for (int i = 0; i < nPts; i++) {
    long double s = 0.0L;
    for (int j = 0; j < nPts; j++) {
      double dij = EuclidCol(P, nPts, dim, i, j);
      s += (long double) (dij * dij);
    }
    out[i] = (double) s;
  }
  // Return:
  return out;
}

// One column d(., col) of the on-the-fly Euclidean distance matrix: a drop-in
// replacement for `d[, col]` in the WideSampleGist sweep. `col` is 1-based;
// the self-distance d(col, col) is returned as 0, matching the matrix path.
//
// [[Rcpp::export]]
Rcpp::NumericVector EuclidColFromPoints_cpp(Rcpp::NumericMatrix points,
                                            int col) {
  int nPts = points.nrow();
  int dim  = points.ncol();
  if (col < 1 || col > nPts) {                       // # nocov start
    Rcpp::stop("'col' must be in [1, %d]; got %d", nPts, col);
  }                                                  // # nocov end
  const double* P = points.begin();
  int col0 = col - 1;
  Rcpp::NumericVector out(nPts);
  for (int i = 0; i < nPts; i++) {
    out[i] = EuclidCol(P, nPts, dim, i, col0);
  }
  // Return:
  return out;
}

// Diameter of the point set: the maximum off-diagonal Euclidean distance and
// the 1-based (row, col) of the first cell achieving it under R's column-major
// `which.max`. WideSampleGist takes the diameter pair from
//   arrayInd(which.max(d_offdiag), dim(d))   # d_offdiag has diag = -Inf
// which, on a symmetric matrix, lands on the below-diagonal cell (row > col)
// with the smallest column index, then smallest row. We reproduce that by
// scanning column-outer / row-inner with a strict `>` (first max wins ties)
// and skipping the diagonal.
//
// Returns c(d_max, row, col); row/col are 1-based. If every off-diagonal
// distance is 0 (or N < 2) d_max is 0 and row/col are 0, letting the R caller
// take its degenerate-data fallback.
//
// [[Rcpp::export]]
Rcpp::NumericVector DiameterFromPoints_cpp(Rcpp::NumericMatrix points) {
  int nPts = points.nrow();
  int dim  = points.ncol();
  const double* P = points.begin();
  double best = R_NegInf;
  int best_r = 0, best_c = 0;       // 0 => no off-diagonal pair found
  for (int c = 0; c < nPts; c++) {
    for (int r = 0; r < nPts; r++) {
      if (r == c) continue;
      double dist = EuclidCol(P, nPts, dim, r, c);
      if (dist > best) {
        best = dist;
        best_r = r + 1;
        best_c = c + 1;
      }
    }
  }
  Rcpp::NumericVector out(3);
  out[0] = (best_r == 0) ? 0.0 : best;
  out[1] = best_r;
  out[2] = best_c;
  // Return:
  return out;
}

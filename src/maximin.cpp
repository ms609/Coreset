#include <Rcpp.h>

// Greedy furthest-point (maximin) selection in a single C++ pass.
//
// Arguments mirror .MaximinFrom() in R/samplers.R:
//   d      square numeric distance matrix (N x N)
//   n      number of points to select (1 <= n <= N)
//   first  1-based index of the first selected point
//
// Returns an integer vector of length n (1-based indices, same as R).

// [[Rcpp::export]]
Rcpp::IntegerVector MaximinFrom_cpp(Rcpp::NumericMatrix d, int n, int first) {
  int nPts = d.nrow();
  if (first < 1 || first > nPts) {
    Rcpp::stop("'first' must be in [1, %d]; got %d", nPts, first);
  }
  Rcpp::IntegerVector selected(n);
  selected[0] = first;

  Rcpp::NumericVector min_dist(nPts);
  int first0 = first - 1;           // convert to 0-based once
  for (int i = 0; i < nPts; i++) {
    min_dist[i] = d(i, first0);
  }
  min_dist[first0] = R_NegInf;      // mask seed before entering loop

  for (int k = 1; k < n; k++) {
    // which.max: first index of the global maximum (strict >, so ties → first)
    int best = 0;
    double best_val = min_dist[0];
    for (int i = 1; i < nPts; i++) {
      if (min_dist[i] > best_val) {
        best_val = min_dist[i];
        best = i;
      }
    }
    selected[k] = best + 1;         // back to 1-based

    // Mask new point before pmin so d(best, best) = 0 cannot overwrite -Inf.
    min_dist[best] = R_NegInf;

    // pmin update: element-wise minimum with the new column.
    for (int i = 0; i < nPts; i++) {
      double dval = d(i, best);
      if (dval < min_dist[i]) min_dist[i] = dval;
    }
    // min_dist[best] stays -Inf: pmin(-Inf, d(best,best)=0) = -Inf. ✓
  }

  return selected;
}

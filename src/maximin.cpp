#include <Rcpp.h>
#include <vector>

// Greedy furthest-point (maximin) selection in a single C++ pass.
//
// Arguments mirror .MaximinFrom() in R/samplers.R:
//   d      square numeric distance matrix (N x N)
//   n      number of points to select (1 <= n <= N)
//   first  1-based index of the first selected point
//
// Returns an integer vector of length n (1-based indices, same as R).

// Threads engage past this many points (matches maximin_points.cpp): below
// it the per-step barrier costs more than the parallel sweep saves. Results
// are identical at every thread count (see the chunk note in the kernel).
static const int GONZ_PAR_MIN = 32768;

// One fused greedy-step update over [lo, hi): merge the new point's distance
// column into min_dist (pmin) and track the range's post-update first
// maximum (strict >, ascending i — R's which.max rule). Masked entries
// (-Inf) pass through and can never win.
static inline void PminArgmaxChunk(const double* col, double* md,
                                   int lo, int hi,
                                   int* nb_out, double* nbv_out) {
  int nb = lo;
  double nbv = R_NegInf;
  for (int i = lo; i < hi; i++) {
    double v = md[i];
    if (col[i] < v) { v = col[i]; md[i] = v; }
    if (v > nbv) { nbv = v; nb = i; }
  }
  *nb_out = nb;
  *nbv_out = nbv;
}

// [[Rcpp::export]]
Rcpp::IntegerVector MaximinFrom_cpp(Rcpp::NumericMatrix d, int n, int first,
                                    int n_threads = 1) {
  int nPts = d.nrow();
  int nthr = n_threads < 1 ? 1 : n_threads;
  if (n < 1 || n > nPts) {                  // defensive: public wrapper guards
    Rcpp::stop("'n' must be in [1, %d]; got %d", nPts, n);
  }
  if (first < 1 || first > nPts) {
    Rcpp::stop("'first' must be in [1, %d]; got %d", nPts, first);
  }
  Rcpp::IntegerVector selected(n);
  selected[0] = first;

  Rcpp::NumericVector min_dist(nPts);
  double* md = min_dist.begin();
  const double* dp = d.begin();
  int first0 = first - 1;           // convert to 0-based once
  for (int i = 0; i < nPts; i++) {
    md[i] = dp[i + (R_xlen_t)first0 * nPts];
  }
  md[first0] = R_NegInf;            // mask seed before entering loop
#ifdef _OPENMP
  std::vector<int> tnb(nthr);       // per-chunk argmax candidates (see below)
  std::vector<double> tnbv(nthr);
#else
  (void)nthr;
#endif

  // T_k = min over greedy steps of the chosen point's insertion distance
  // (best_val), which equals the selection's minimum pairwise distance. Tracked
  // here so the ensemble driver need not re-score with a d[idx, idx] subset.
  double tk = R_PosInf;

  // which.max, folded: the standalone scan below seeds (best, best_val); every
  // later step's argmax rides the pmin update pass, which already reads each
  // entry's post-update value in ascending order — strict > keeps the same
  // first-maximum tie rule, and masked entries (-Inf) can never win. The final
  // step's update is skipped outright: nothing reads min_dist after the last
  // pick (t_k comes from best_val at pick time).
  int best = 0;
  double best_val = md[0];
  for (int i = 1; i < nPts; i++) {
    if (md[i] > best_val) {
      best_val = md[i];
      best = i;
    }
  }

  for (int k = 1; k < n; k++) {
    selected[k] = best + 1;         // back to 1-based
    if (best_val < tk) tk = best_val;   // running min insertion distance

    // Mask new point before pmin so d(best, best) = 0 cannot overwrite -Inf.
    md[best] = R_NegInf;

    if (k + 1 < n) {
      // pmin update fused with the next step's argmax (see the note above).
      // Past GONZ_PAR_MIN points the sweep splits into one contiguous chunk
      // per thread (disjoint md ranges) and chunk maxima merge in ascending
      // order with strict > — the global first maximum, identical at every
      // thread count.
      const double* col = dp + (R_xlen_t)best * nPts;
      int nb;
      double nbv;
#ifdef _OPENMP
      if (nthr > 1 && nPts >= GONZ_PAR_MIN) {
        const int C = nthr;
#pragma omp parallel for num_threads(nthr) schedule(static, 1)
        for (int t = 0; t < C; t++) {
          int lo = (int)((R_xlen_t)nPts * t / C);
          int hi = (int)((R_xlen_t)nPts * (t + 1) / C);
          PminArgmaxChunk(col, md, lo, hi, &tnb[t], &tnbv[t]);
        }
        nb = tnb[0]; nbv = tnbv[0];
        for (int t = 1; t < C; t++) {
          if (tnbv[t] > nbv) { nbv = tnbv[t]; nb = tnb[t]; }
        }
      } else
#endif
      {
        PminArgmaxChunk(col, md, 0, nPts, &nb, &nbv);
      }
      best = nb;
      best_val = nbv;
    }
    // min_dist[best] stays -Inf: pmin(-Inf, d(best,best)=0) = -Inf. ✓
  }

  // T_k (min pairwise distance of the selection), free; NA for n < 2.
  selected.attr("t_k") = (n >= 2 && R_finite(tk)) ? tk : NA_REAL;

  return selected;
}

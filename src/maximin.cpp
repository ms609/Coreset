#include <Rcpp.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

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
// For THIS kernel the threshold is effectively unreachable — an N x N
// matrix at N = 32768 is 8.6 GB — and lowering it measured slower at every
// RAM-feasible size (round 9: N = 16384, 8 threads 2x slower than serial;
// the per-step barriers dominate a ~15 us sweep), so the matrix pass runs
// serially in practice. The constant is kept aligned with the points
// kernel, whose large-N regime genuinely engages it.
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

// One complete greedy pass from `first0`, writing the n selected 1-based
// indices into `out` and the pass's T_k (NA for n < 2) into `tk_out`.
// All storage is caller-owned: `md` is length-nPts scratch, and `tnb`/`tnbv`
// are the length-nthr per-chunk argmax buffers, read only when the pass
// threads its own sweep (nthr > 1 and nPts >= GONZ_PAR_MIN). Touches no R
// API, so callers may run several passes concurrently on worker threads.
static void MaximinPass(const double* dp, int nPts, int n, int first0,
                        double* md, int* out, double* tk_out,
                        int nthr, int* tnb, double* tnbv) {
  out[0] = first0 + 1;              // back to 1-based
  for (int i = 0; i < nPts; i++) {
    md[i] = dp[i + (R_xlen_t)first0 * nPts];
  }
  md[first0] = R_NegInf;            // mask seed before entering loop
#ifndef _OPENMP
  (void)nthr; (void)tnb; (void)tnbv;
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
    out[k] = best + 1;              // back to 1-based
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
      // nocov: unreachable in practice on this path -- a square matrix at
      // GONZ_PAR_MIN rows already occupies 8.6 GB (see the note above
      // MaximinMultiFrom_cpp), so no test can supply an input that engages it.
      // The coordinate counterpart in maximin_points.cpp is testable and is
      // exercised at 32k+ points.
      if (nthr > 1 && nPts >= GONZ_PAR_MIN) {            // # nocov start
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
      } else                                             // # nocov end
#endif
      {
        PminArgmaxChunk(col, md, 0, nPts, &nb, &nbv);
      }
      best = nb;
      best_val = nbv;
    }
    // min_dist[best] stays -Inf: pmin(-Inf, d(best,best)=0) = -Inf. ✓
  }

  *tk_out = (n >= 2 && R_finite(tk)) ? tk : NA_REAL;
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
  std::vector<double> md(nPts);
  std::vector<int> tnb(nthr);       // per-chunk argmax candidates
  std::vector<double> tnbv(nthr);
  double tk;
  MaximinPass(d.begin(), nPts, n, first - 1, md.data(), selected.begin(), &tk,
              nthr, tnb.data(), tnbv.data());

  // T_k (min pairwise distance of the selection), free; NA for n < 2.
  selected.attr("t_k") = tk;

  return selected;
}

// Several independent greedy passes over one distance matrix — the ensemble
// path's restarts, run one per thread.
//
// Each pass is a self-contained function of its seed over a read-only `d`, so
// restarts are the coarsest parallel axis available here and the only one the
// matrix path can use: its per-step sweep is serial at every RAM-feasible N
// (a matrix at GONZ_PAR_MIN already occupies 8 GB). Threads are capped at the
// seed count, since a pass cannot be split across them.
//
// Returns `idx`, an n x nSeeds matrix whose column s is the selection seeded
// by firsts[s], and `t_k`, that column's minimum pairwise distance.
// [[Rcpp::export]]
Rcpp::List MaximinMultiFrom_cpp(Rcpp::NumericMatrix d, int n,
                                Rcpp::IntegerVector firsts,
                                int n_threads = 1) {
  const int nPts = d.nrow();
  const int nS = firsts.size();
  const int nthr = n_threads < 1 ? 1 : n_threads;
  if (n < 1 || n > nPts) {                  // defensive: public wrapper guards
    Rcpp::stop("'n' must be in [1, %d]; got %d", nPts, n);
  }
  if (nS < 1) {
    Rcpp::stop("'firsts' must name at least one seed");
  }
  for (int s = 0; s < nS; s++) {
    if (firsts[s] < 1 || firsts[s] > nPts) {
      Rcpp::stop("'firsts' must lie in [1, %d]; got %d", nPts, firsts[s]);
    }
  }

  // Every R allocation happens here, before any parallel region: worker
  // threads run MaximinPass, which touches no R API.
  Rcpp::IntegerMatrix selected(n, nS);
  Rcpp::NumericVector tks(nS);
  int* op = selected.begin();
  double* tp = tks.begin();
  const double* dp = d.begin();
  std::vector<int> f0(nS);
  for (int s = 0; s < nS; s++) f0[s] = firsts[s] - 1;

#ifdef _OPENMP
  const int seedThr = nthr < nS ? nthr : nS;
  // Restart-parallel arm. Skipped once a single pass can use every thread
  // itself: above GONZ_PAR_MIN, splitting one pass nthr ways beats running
  // nS passes concurrently (each then serial), because nS is typically the
  // smaller number.
  if (seedThr > 1 && nPts < GONZ_PAR_MIN) {
    std::vector<double> md((R_xlen_t)nPts * seedThr);   // one buffer per thread
#pragma omp parallel for num_threads(seedThr) schedule(static, 1)
    for (int s = 0; s < nS; s++) {
      double* mine = md.data() + (R_xlen_t)nPts * omp_get_thread_num();
      MaximinPass(dp, nPts, n, f0[s], mine, op + (R_xlen_t)s * n, tp + s,
                  1, NULL, NULL);
    }
    return Rcpp::List::create(Rcpp::_["idx"] = selected,
                              Rcpp::_["t_k"] = tks);
  }
#endif
  {
    std::vector<double> md(nPts);
    std::vector<int> tnb(nthr);
    std::vector<double> tnbv(nthr);
    for (int s = 0; s < nS; s++) {
      MaximinPass(dp, nPts, n, f0[s], md.data(), op + (R_xlen_t)s * n, tp + s,
                  nthr, tnb.data(), tnbv.data());
    }
  }
  return Rcpp::List::create(Rcpp::_["idx"] = selected, Rcpp::_["t_k"] = tks);
}

// One contiguous column range [cLo, cHi) of the off-diagonal first-maximum
// scan: column-outer / row-inner with strict > (first max wins ties),
// diagonal skipped — R's column-major which.max order restricted to the
// range.
static inline void OffDiagMaxChunk(const double* dp, int n, int cLo, int cHi,
                                   double* best_out, int* br_out) {
  double best = R_NegInf;
  int br = 0;
  for (int c = cLo; c < cHi; c++) {
    const double* col = dp + (R_xlen_t)c * n;
    for (int r = 0; r < n; r++) {
      if (r == c) continue;
      if (col[r] > best) { best = col[r]; br = r + 1; }
    }
  }
  *best_out = best;
  *br_out = br;
}

// Off-diagonal maximum of a square distance matrix and the 1-based row of
// the first cell attaining it under R's column-major `which.max`. Replaces
// the seeding idiom
//   dOff <- d; diag(dOff) <- -Inf
//   dMax <- max(dOff); arrayInd(which.max(dOff), dim(dOff))[1L, 1L]
// which copies the full N x N matrix and then scans it twice more. Unlike
// DiameterFromPoints_cpp, the FULL matrix is scanned: `.AsDistMatrix`
// silently accepts asymmetric matrices, where an upper-triangle cell can be
// the winner, so the symmetric lower-triangle shortcut would be wrong here.
//
// Returns c(max, row); row is 0 (and max -Inf) when there is no
// off-diagonal cell (N < 2), letting the R caller take its degenerate-data
// fallback exactly as before.
//
// Parallel over contiguous column chunks (every column holds n - 1
// off-diagonal cells, so plain ranges balance); chunk winners merge in
// ascending chunk order with strict >, preserving the global column-major
// first maximum at every thread count.
//
// [[Rcpp::export]]
Rcpp::NumericVector MatrixOffDiagMax_cpp(Rcpp::NumericMatrix d,
                                         int n_threads = 1) {
  int n = d.nrow();
  const double* dp = d.begin();
  double best = R_NegInf;
  int best_r = 0;
#ifdef _OPENMP
  if (n_threads > 1 && n > 64) {
    const int C = n_threads;
    std::vector<double> tb(C, R_NegInf);
    std::vector<int> tr(C, 0);
#pragma omp parallel for num_threads(n_threads) schedule(static, 1)
    for (int t = 0; t < C; t++) {
      int lo = (int)((R_xlen_t)n * t / C);
      int hi = (int)((R_xlen_t)n * (t + 1) / C);
      OffDiagMaxChunk(dp, n, lo, hi, &tb[t], &tr[t]);
    }
    for (int t = 0; t < C; t++) {
      if (tb[t] > best) { best = tb[t]; best_r = tr[t]; }
    }
  } else
#endif
  {
    OffDiagMaxChunk(dp, n, 0, n, &best, &best_r);
  }
  Rcpp::NumericVector out(2);
  out[0] = best;
  out[1] = best_r;
  // Return:
  return out;
}

// rowSums(d^2) for a square matrix without materialising d^2 (a full N x N
// double allocation in R). Bit-exact: R squares each entry as a double and
// rowSums then accumulates each row in a long double over columns in
// ascending index; here acc[i] receives the identical double squares in the
// identical j-ascending order via one column-sequential pass, so the
// returned doubles equal R's exactly (a store/load round-trip of the
// squared double is exact). (A register-banded row layout was measured
// identical to this accumulator-array pass — see log.md round 9 — so the
// simpler orientation is kept.)
//
// Parallel over contiguous row bands: each row's accumulator lives on one
// thread and still receives its terms in j-ascending order (each thread
// walks all columns, reading its band's contiguous segment), so the result
// is identical at every thread count.
//
// [[Rcpp::export]]
Rcpp::NumericVector RowSqSumsFromMatrix_cpp(Rcpp::NumericMatrix d,
                                            int n_threads = 1) {
  int n = d.nrow();
  const double* dp = d.begin();
  Rcpp::NumericVector out(n);
  double* o = out.begin();
#ifdef _OPENMP
  if (n_threads > 1 && n > 64) {
#pragma omp parallel num_threads(n_threads)
    {
      int nt = omp_get_num_threads();
      int t  = omp_get_thread_num();
      int lo = (int)((R_xlen_t)n * t / nt);
      int hi = (int)((R_xlen_t)n * (t + 1) / nt);
      if (hi > lo) {
        std::vector<long double> acc(hi - lo, 0.0L);
        for (int j = 0; j < n; j++) {
          const double* col = dp + (R_xlen_t)j * n;
          for (int i = lo; i < hi; i++) {
            double v = col[i];
            acc[i - lo] += (long double)(v * v);
          }
        }
        for (int i = lo; i < hi; i++) o[i] = (double)acc[i - lo];
      }
    }
    return out;
  }
#else
  (void)n_threads;
#endif
  std::vector<long double> acc(n, 0.0L);
  for (int j = 0; j < n; j++) {
    const double* col = dp + (R_xlen_t)j * n;
    for (int i = 0; i < n; i++) {
      double v = col[i];
      acc[i] += (long double)(v * v);
    }
  }
  for (int i = 0; i < n; i++) o[i] = (double)acc[i];
  // Return:
  return out;
}

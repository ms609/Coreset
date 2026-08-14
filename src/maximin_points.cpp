#include <Rcpp.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// Coordinate-based greedy furthest-point (maximin) selection.
//
// This is the O(N*k*dim) time, O(N) memory counterpart of MaximinFrom_cpp
// (src/maximin.cpp), which reads a fully materialised N x N distance matrix.
// Each greedy step needs only the column d(., best) of distances from the
// most-recently-selected point; for Euclidean data that column is recomputed
// from coordinates in O(N*dim), so the dense matrix is never built.
//
// Selection indices match the matrix path bit-for-bit -- a testing
// convenience, not a correctness requirement (ties could legitimately
// break either way). stats::dist() computes each Euclidean distance as
//   sqrt( sum_j (x[i,j] - x[c,j])^2 )
// accumulating dev*dev in a plain double over dimensions in increasing column
// order (R's R_euclidean in src/library/stats/src/distance.c).
//
// The squared accumulation
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

// Threads engage on the greedy pass only past this many points: below it a
// step's whole update is a few tens of microseconds and the per-step barrier
// costs more than it buys. Results are identical either way (every element is
// computed by exactly one thread in the same order, and the argmax merge
// preserves the first-maximum rule), so this is a pure tuning constant.
static const int GONZ_PAR_MIN = 32768;

// Squared-Euclidean between two points, accumulated in double over columns
// in ascending order — the exact argument stats::dist() passes to sqrt.
static inline double EuclidSqPair(const double* P, int nPts, int dim,
                                  int i, int c) {
  double s = 0.0;
  for (int j = 0; j < dim; j++) {
    double dev = P[i + (R_xlen_t)j * nPts] - P[c + (R_xlen_t)j * nPts];
    s += dev * dev;
  }
  return s;
}

// True Euclidean distance, bit-matching the stats::dist() matrix entry.
static inline double EuclidCol(const double* P, int nPts, int dim,
                               int i, int c) {
  return std::sqrt(EuclidSqPair(P, nPts, dim, i, c));
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

// One sweep over [lo, hi) covering the NB (1..4) dimensions starting at
// column j0, with the running squared sum held in a register. READ_DC loads
// the running sum from dc (false for a leading block, which starts the sum
// at its first dev^2). MERGE marks the final block: instead of storing the
// sum back to dc it merges into md (pmin) and tracks the range's
// post-update first maximum (strict >, ascending i — R's which.max rule;
// masked -Inf entries pass through untouched and can never win).
//
// Exactness: each element's sum is the same left-associated chain
// (((dev_0^2 + dev_1^2) + ...) + dev_{dim-1}^2) as EuclidColSqInto's
// one-dimension-at-a-time accumulation — a store/load round-trip of a
// double is exact, so batching four adds per sweep instead of one changes
// only where the partial sums live (register vs dc), never their values.
template <int NB, bool READ_DC, bool MERGE>
static inline void SweepBlock(const double* P, int nPts, int c, int j0,
                              double* md, double* dc, int lo, int hi,
                              int* nb_out, double* nbv_out) {
  const double* p0 = P + (R_xlen_t)j0 * nPts;
  const double* p1 = NB > 1 ? p0 + nPts : p0;
  const double* p2 = NB > 2 ? p1 + nPts : p1;
  const double* p3 = NB > 3 ? p2 + nPts : p2;
  const double c0 = p0[c], c1 = p1[c], c2 = p2[c], c3 = p3[c];
  int nb = lo;
  double nbv = R_NegInf;
  for (int i = lo; i < hi; i++) {
    double dev = p0[i] - c0;
    double s = READ_DC ? dc[i] + dev * dev : dev * dev;
    if (NB > 1) { dev = p1[i] - c1; s += dev * dev; }
    if (NB > 2) { dev = p2[i] - c2; s += dev * dev; }
    if (NB > 3) { dev = p3[i] - c3; s += dev * dev; }
    if (MERGE) {
      double m = md[i];
      if (s < m) { m = s; md[i] = m; }
      if (m > nbv) { nbv = m; nb = i; }
    } else {
      dc[i] = s;
    }
  }
  if (MERGE) {
    *nb_out = nb;
    *nbv_out = nbv;
  }
}

// One fused greedy-step update over the contiguous range [lo, hi): the new
// centre's squared column, the pmin merge and the next argmax, processed in
// blocks of up to four dimensions per sweep (SweepBlock). dim <= 4 runs as
// a single sweep touching dc not at all; larger dim writes dc once per
// leading block and reads it once per continuation, cutting dc traffic
// from one read-modify-write per middle dimension to one per block of four.
static inline void FusedUpdateChunk(const double* P, int nPts, int dim, int c,
                                    double* md, double* dc, int lo, int hi,
                                    int* nb_out, double* nbv_out) {
  if (dim == 0) {
    // Degenerate zero-column input: the old EuclidColSqInto produced an
    // all-zero column, i.e. pmin(md, 0). Kept defined.
    int nb = lo;
    double nbv = R_NegInf;
    for (int i = lo; i < hi; i++) {
      double m = md[i];
      if (0.0 < m) { m = 0.0; md[i] = m; }
      if (m > nbv) { nbv = m; nb = i; }
    }
    *nb_out = nb;
    *nbv_out = nbv;
    return;
  }
  int j0 = 0;
  bool leading = true;
  while (dim - j0 > 4) {
    if (leading) {
      SweepBlock<4, false, false>(P, nPts, c, j0, md, dc, lo, hi,
                                  nb_out, nbv_out);
    } else {
      SweepBlock<4, true, false>(P, nPts, c, j0, md, dc, lo, hi,
                                 nb_out, nbv_out);
    }
    leading = false;
    j0 += 4;
  }
  switch ((dim - j0 - 1) * 2 + (leading ? 1 : 0)) {
    case 0: SweepBlock<1, true,  true>(P, nPts, c, j0, md, dc, lo, hi,
                                       nb_out, nbv_out); break;
    case 1: SweepBlock<1, false, true>(P, nPts, c, j0, md, dc, lo, hi,
                                       nb_out, nbv_out); break;
    case 2: SweepBlock<2, true,  true>(P, nPts, c, j0, md, dc, lo, hi,
                                       nb_out, nbv_out); break;
    case 3: SweepBlock<2, false, true>(P, nPts, c, j0, md, dc, lo, hi,
                                       nb_out, nbv_out); break;
    case 4: SweepBlock<3, true,  true>(P, nPts, c, j0, md, dc, lo, hi,
                                       nb_out, nbv_out); break;
    case 5: SweepBlock<3, false, true>(P, nPts, c, j0, md, dc, lo, hi,
                                       nb_out, nbv_out); break;
    case 6: SweepBlock<4, true,  true>(P, nPts, c, j0, md, dc, lo, hi,
                                       nb_out, nbv_out); break;
    default: SweepBlock<4, false, true>(P, nPts, c, j0, md, dc, lo, hi,
                                        nb_out, nbv_out); break;
  }
}

// One complete greedy pass from `first0`, writing the n selected 1-based
// indices into `out` and the pass's T_k (NA for n < 2) into `tk_out`.
// Coordinate counterpart of MaximinPass() in maximin.cpp; all storage is
// caller-owned (`md` length nPts; `dc` length nPts when dim > 4, else unused;
// `tnb`/`tnbv` length nthr, read only when the pass threads its own sweep).
// Touches no R API, so callers may run several passes concurrently on worker
// threads.
static void MaximinPointsPass(const double* P, int nPts, int dim, int n,
                              int first0, int mask,
                              double* md, double* dc, int* out, double* tk_out,
                              int nthr, int* tnb, double* tnbv) {
  out[0] = first0 + 1;              // back to 1-based

  // md holds SQUARED nearest-selected distances throughout the pass; the
  // selection is identical to working in true distances (sqrt is monotone) but
  // skips n*N sqrt calls. T_k is recovered with a single sqrt at the end.
  EuclidColSqInto(P, nPts, dim, first0, md);
#ifndef _OPENMP
  (void)nthr; (void)tnb; (void)tnbv;
#endif
  md[first0] = R_NegInf;            // mask seed before entering loop
  if (mask >= 1) {                                   // # nocov start
    md[mask - 1] = R_NegInf;        // pin forbidden point (anti-medoid medoid)
  }                                                  // # nocov end

  // T_k = min over greedy steps of the chosen point's insertion distance, which
  // is exactly best_val at each step. Tracked in squared space, sqrt'd once.
  double tk_sq = R_PosInf;

  // which.max, folded (matches MaximinFrom_cpp): the standalone scan below
  // seeds (best, best_val); each later step's argmax rides the last pass of
  // the fused update below, which reads every entry's post-update value in
  // ascending order — strict > keeps the first-maximum tie rule, and masked
  // entries (-Inf) can never win. The final step's update is skipped:
  // nothing reads min_dist after the last pick.
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
    if (best_val < tk_sq) tk_sq = best_val;   // running min insertion distance

    // Mask new point before the pmin update so its self-distance (0) cannot
    // overwrite -Inf. A pinned `mask` point keeps -Inf for the same reason:
    // any non-negative distance is never < -Inf.
    md[best] = R_NegInf;

    if (k + 1 < n) {
      // The new point's squared column, the pmin merge and the next argmax in
      // one fused sweep (see FusedUpdateChunk). Past GONZ_PAR_MIN points the
      // sweep splits into one contiguous chunk per thread — disjoint ranges
      // of md/dcol, so no thread ever reads another's writes — and the
      // chunk maxima merge in ascending chunk order with strict >, which is
      // exactly the global first-maximum: identical selection at every
      // thread count.
      const int c = best;
      int nb;
      double nbv;
#ifdef _OPENMP
      if (nthr > 1 && nPts >= GONZ_PAR_MIN) {
        const int C = nthr;
#pragma omp parallel for num_threads(nthr) schedule(static, 1)
        for (int t = 0; t < C; t++) {
          int lo = (int)((R_xlen_t)nPts * t / C);
          int hi = (int)((R_xlen_t)nPts * (t + 1) / C);
          FusedUpdateChunk(P, nPts, dim, c, md, dc, lo, hi,
                           &tnb[t], &tnbv[t]);
        }
        nb = tnb[0]; nbv = tnbv[0];
        for (int t = 1; t < C; t++) {
          if (tnbv[t] > nbv) { nbv = tnbv[t]; nb = tnb[t]; }
        }
      } else
#endif
      {
        FusedUpdateChunk(P, nPts, dim, c, md, dc, 0, nPts, &nb, &nbv);
      }
      best = nb;
      best_val = nbv;
    }
  }

  *tk_out = (n >= 2 && R_finite(tk_sq))
    ? std::sqrt(tk_sq < 0.0 ? 0.0 : tk_sq) : NA_REAL;
}

// [[Rcpp::export]]
Rcpp::IntegerVector MaximinFromPoints_cpp(Rcpp::NumericMatrix points,
                                          int n, int first, int mask,
                                          int n_threads = 1) {
  int nPts = points.nrow();
  int dim  = points.ncol();
  int nthr = n_threads < 1 ? 1 : n_threads;
  if (n < 1 || n > nPts) {                  // defensive: public wrapper guards
    Rcpp::stop("'n' must be in [1, %d]; got %d", nPts, n);
  }
  if (first < 1 || first > nPts) {
    Rcpp::stop("'first' must be in [1, %d]; got %d", nPts, first);
  }
  if (mask < 0 || mask > nPts) {                    // # nocov start
    Rcpp::stop("'mask' must be in [0, %d]; got %d", nPts, mask);
  }                                                  // # nocov end

  Rcpp::IntegerVector selected(n);
  std::vector<double> md(nPts);
  // Per-step running sums between dimension blocks; a single-block dim
  // (<= 4) keeps the whole sum in registers and never touches it.
  std::vector<double> dcol(dim > 4 ? nPts : 0);
  std::vector<int> tnb(nthr);       // per-chunk argmax candidates
  std::vector<double> tnbv(nthr);
  double tk;
  MaximinPointsPass(points.begin(), nPts, dim, n, first - 1, mask,
                    md.data(), dcol.data(), selected.begin(), &tk,
                    nthr, tnb.data(), tnbv.data());

  // T_k (min pairwise distance of the selection), computed for free; NA for
  // n < 2. The ensemble driver reads this attribute instead of re-running a
  // full stats::dist() over the selection (see .GonzEnsembleFromPoints).
  selected.attr("t_k") = tk;

  // Return:
  return selected;
}

// Several independent greedy passes over one coordinate matrix — the ensemble
// path's restarts, run one per thread. Coordinate counterpart of
// MaximinMultiFrom_cpp(); see its note on why restarts are the useful axis
// below GONZ_PAR_MIN and the per-step sweep the useful axis above it.
//
// Returns `idx`, an n x nSeeds matrix whose column s is the selection seeded
// by firsts[s], and `t_k`, that column's minimum pairwise distance.
// [[Rcpp::export]]
Rcpp::List MaximinMultiFromPoints_cpp(Rcpp::NumericMatrix points, int n,
                                      Rcpp::IntegerVector firsts,
                                      int n_threads = 1) {
  const int nPts = points.nrow();
  const int dim  = points.ncol();
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
  // threads run MaximinPointsPass, which touches no R API.
  Rcpp::IntegerMatrix selected(n, nS);
  Rcpp::NumericVector tks(nS);
  int* op = selected.begin();
  double* tp = tks.begin();
  const double* P = points.begin();
  const int dcN = dim > 4 ? nPts : 0;
  std::vector<int> f0(nS);
  for (int s = 0; s < nS; s++) f0[s] = firsts[s] - 1;

#ifdef _OPENMP
  const int seedThr = nthr < nS ? nthr : nS;
  if (seedThr > 1 && nPts < GONZ_PAR_MIN) {
    std::vector<double> md((R_xlen_t)nPts * seedThr);   // one buffer per thread
    std::vector<double> dc((R_xlen_t)dcN * seedThr);
#pragma omp parallel for num_threads(seedThr) schedule(static, 1)
    for (int s = 0; s < nS; s++) {
      const int t = omp_get_thread_num();
      MaximinPointsPass(P, nPts, dim, n, f0[s], 0,
                        md.data() + (R_xlen_t)nPts * t,
                        dc.data() + (R_xlen_t)dcN * t,
                        op + (R_xlen_t)s * n, tp + s, 1, NULL, NULL);
    }
    return Rcpp::List::create(Rcpp::_["idx"] = selected,
                              Rcpp::_["t_k"] = tks);
  }
#endif
  {
    std::vector<double> md(nPts);
    std::vector<double> dcol(dcN);
    std::vector<int> tnb(nthr);
    std::vector<double> tnbv(nthr);
    for (int s = 0; s < nS; s++) {
      MaximinPointsPass(P, nPts, dim, n, f0[s], 0, md.data(), dcol.data(),
                        op + (R_xlen_t)s * n, tp + s,
                        nthr, tnb.data(), tnbv.data());
    }
  }
  return Rcpp::List::create(Rcpp::_["idx"] = selected, Rcpp::_["t_k"] = tks);
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
// Each row's long-double accumulation receives its contributions in the same
// j-ascending order on every path, so the returned doubles equal R's exactly:
//
// - Serial: the scan walks the strict lower triangle once (i < j), adding
//   each pair's distance to both endpoint rows. Row r's accumulator then
//   receives d(r, k) for k = 0..r-1 during the earlier outer iterations and
//   d(r, k) for k = r+1..N-1 during its own — ascending k overall, exactly
//   the reference order, with the diagonal's exact +0 contribution dropped
//   (adding a true zero to a long double never rounds). d(i, j) and d(j, i)
//   are the same double bit-for-bit ((-x)*(-x) == x*x), so halving the
//   EuclidCol/sqrt evaluations changes no summand.
// - Parallel: pair-halving would interleave rows' updates across threads in
//   schedule order, so each thread instead computes whole rows (full N-term
//   sweeps, j ascending) — the identical double per row at every thread
//   count. O(N^2 * dim) with a sqrt per pair — compute-bound, the one
//   FarFirst primitive where threads scale near-linearly.
//
// [[Rcpp::export]]
Rcpp::NumericVector RowSumsFromPoints_cpp(Rcpp::NumericMatrix points,
                                          int n_threads = 1) {
  int nPts = points.nrow();
  int dim  = points.ncol();
  const double* P = points.begin();
  Rcpp::NumericVector out(nPts);
  double* o = out.begin();
#ifdef _OPENMP
  if (n_threads > 1 && nPts > 64) {
#pragma omp parallel for num_threads(n_threads) schedule(static)
    for (int i = 0; i < nPts; i++) {
      long double s = 0.0L;
      for (int j = 0; j < nPts; j++) {
        s += (long double) EuclidCol(P, nPts, dim, i, j);
      }
      o[i] = (double) s;
    }
    return out;
  }
#else
  (void)n_threads;
#endif
  std::vector<long double> acc(nPts, 0.0L);
  for (int i = 0; i < nPts - 1; i++) {
    for (int j = i + 1; j < nPts; j++) {
      long double v = (long double) EuclidCol(P, nPts, dim, i, j);
      acc[i] += v;
      acc[j] += v;
    }
  }
  for (int i = 0; i < nPts; i++) o[i] = (double) acc[i];
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
// Serial pair-halved lower-triangle scan and per-row parallel path exactly
// as RowSumsFromPoints_cpp (see its exactness note); the summand here is
// the squared double (dij * dij, matching R's `d^2`), which is the same
// bit-for-bit from either triangle.
//
// [[Rcpp::export]]
Rcpp::NumericVector RowSqSumsFromPoints_cpp(Rcpp::NumericMatrix points,
                                            int n_threads = 1) {
  int nPts = points.nrow();
  int dim  = points.ncol();
  const double* P = points.begin();
  Rcpp::NumericVector out(nPts);
  double* o = out.begin();
#ifdef _OPENMP
  if (n_threads > 1 && nPts > 64) {
#pragma omp parallel for num_threads(n_threads) schedule(static)
    for (int i = 0; i < nPts; i++) {
      long double s = 0.0L;
      for (int j = 0; j < nPts; j++) {
        double dij = EuclidCol(P, nPts, dim, i, j);
        s += (long double) (dij * dij);
      }
      o[i] = (double) s;
    }
    return out;
  }
#else
  (void)n_threads;
#endif
  std::vector<long double> acc(nPts, 0.0L);
  for (int i = 0; i < nPts - 1; i++) {
    for (int j = i + 1; j < nPts; j++) {
      double dij = EuclidCol(P, nPts, dim, i, j);
      long double v = (long double) (dij * dij);
      acc[i] += v;
      acc[j] += v;
    }
  }
  for (int i = 0; i < nPts; i++) o[i] = (double) acc[i];
  // Return:
  return out;
}

// Fused RowSums + RowSqSums in ONE pair sweep, for ensembles that need both
// aggregate families (anti_medoid/medoid/rowsum read the plain sums;
// rownorm reads the squared sums). Each family's accumulators receive the
// identical summands ((long double)d_ij and (long double)(d_ij * d_ij),
// d_ij the same double as the dedicated kernels') in the identical
// j-ascending per-row order — the two families never mix, so both results
// are bit-identical to RowSumsFromPoints_cpp / RowSqSumsFromPoints_cpp;
// the fusion just computes each pair's distance (and its sqrt) once
// instead of twice. Serial pair-halved triangle scan and per-row parallel
// path exactly as those kernels (see RowSumsFromPoints_cpp's note).
//
// Returns list(sums, sqsums).
//
// [[Rcpp::export]]
Rcpp::List RowSumsSqFromPoints_cpp(Rcpp::NumericMatrix points,
                                   int n_threads = 1) {
  int nPts = points.nrow();
  int dim  = points.ncol();
  const double* P = points.begin();
  Rcpp::NumericVector sums(nPts), sqsums(nPts);
  double* oS = sums.begin();
  double* oQ = sqsums.begin();
#ifdef _OPENMP
  if (n_threads > 1 && nPts > 64) {
#pragma omp parallel for num_threads(n_threads) schedule(static)
    for (int i = 0; i < nPts; i++) {
      long double s = 0.0L, q = 0.0L;
      for (int j = 0; j < nPts; j++) {
        double dij = EuclidCol(P, nPts, dim, i, j);
        s += (long double) dij;
        q += (long double) (dij * dij);
      }
      oS[i] = (double) s;
      oQ[i] = (double) q;
    }
    return Rcpp::List::create(sums, sqsums);
  }
#else
  (void)n_threads;
#endif
  std::vector<long double> accS(nPts, 0.0L), accQ(nPts, 0.0L);
  for (int i = 0; i < nPts - 1; i++) {
    for (int j = i + 1; j < nPts; j++) {
      double dij = EuclidCol(P, nPts, dim, i, j);
      long double v = (long double) dij;
      long double vq = (long double) (dij * dij);
      accS[i] += v;
      accS[j] += v;
      accQ[i] += vq;
      accQ[j] += vq;
    }
  }
  for (int i = 0; i < nPts; i++) {
    oS[i] = (double) accS[i];
    oQ[i] = (double) accQ[i];
  }
  // Return:
  return Rcpp::List::create(sums, sqsums);
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
// which, on a symmetric matrix, always lands on a below-diagonal cell: pair
// {a, b} (a < b) appears at linear indices b + a*N (lower) and a + b*N
// (upper), and b + a*N < a + b*N whenever a < b. Scanning only the strict
// lower triangle in the same column-outer / row-inner order therefore visits
// every possible winner in the reference's relative order — half the pairs,
// tie-identical pick.
//
// The scan itself runs in SQUARED space with a guarded sqrt: `best` always
// equals sqrt(bestSq) (see the loop), so a pair with sq <= bestSq satisfies
// sqrt(sq) <= best and can never win the reference's strict-> comparison —
// skipping its sqrt drops no candidate. A pair with sq > bestSq gets the
// true sqrt and the reference's own strict-> test on distances (two distinct
// squared values can round to the same sqrt, so comparing sq alone could
// pick a later cell the reference would have rejected on the tie rule).
// sqrt therefore runs only on running-max candidates, O(log) of them in
// expectation, and the hot loop is sqrt-free.
static inline void DiameterChunk(const double* P, int nPts, int dim,
                                 int cLo, int cHi,
                                 double* best_out, int* br_out, int* bc_out) {
  double best = R_NegInf;     // true distance of the running winner
  double bestSq = R_NegInf;   // largest squared distance seen so far
  int br = 0, bc = 0;
  for (int c = cLo; c < cHi; c++) {
    for (int r = c + 1; r < nPts; r++) {
      double sq = EuclidSqPair(P, nPts, dim, r, c);
      if (sq > bestSq) {
        bestSq = sq;
        double dist = std::sqrt(sq);
        if (dist > best) { best = dist; br = r + 1; bc = c + 1; }
      }
    }
  }
  *best_out = best;
  *br_out = br;
  *bc_out = bc;
}

// Returns c(d_max, row, col); row/col are 1-based. If every off-diagonal
// distance is 0 (or N < 2) d_max is 0 and row/col are 0, letting the R caller
// take its degenerate-data fallback.
//
// Parallel over column chunks balanced by PAIR count (column c holds
// nPts-1-c pairs, so equal column ranges would leave the last thread nearly
// idle); each chunk scans in the serial order and chunk winners merge in
// ascending chunk order with strict > on the true distance, so the winning
// cell is the global column-major first maximum at every thread count —
// exactly the serial scan's pick.
//
// [[Rcpp::export]]
Rcpp::NumericVector DiameterFromPoints_cpp(Rcpp::NumericMatrix points,
                                           int n_threads = 1) {
  int nPts = points.nrow();
  int dim  = points.ncol();
  const double* P = points.begin();
  double best = R_NegInf;
  int best_r = 0, best_c = 0;       // 0 => no off-diagonal pair found
#ifdef _OPENMP
  if (n_threads > 1 && nPts > 64) {
    const int C = n_threads;
    std::vector<double> tb(C, R_NegInf);
    std::vector<int> tr(C, 0), tc(C, 0), lo(C + 1, 0);
    const double T = (double)nPts * (nPts - 1) / 2.0;
    int cb = 0;
    for (int t = 1; t < C; t++) {
      double target = T * t / C;
      while (cb < nPts &&
             (double)cb * (nPts - 1) - (double)cb * (cb - 1) / 2.0 < target) {
        cb++;
      }
      lo[t] = cb;
    }
    lo[C] = nPts;
#pragma omp parallel for num_threads(n_threads) schedule(static, 1)
    for (int t = 0; t < C; t++) {
      DiameterChunk(P, nPts, dim, lo[t], lo[t + 1], &tb[t], &tr[t], &tc[t]);
    }
    for (int t = 0; t < C; t++) {
      if (tb[t] > best) { best = tb[t]; best_r = tr[t]; best_c = tc[t]; }
    }
  } else
#endif
  {
    DiameterChunk(P, nPts, dim, 0, nPts, &best, &best_r, &best_c);
  }
  Rcpp::NumericVector out(3);
  out[0] = (best_r == 0) ? 0.0 : best;
  out[1] = best_r;
  out[2] = best_c;
  // Return:
  return out;
}

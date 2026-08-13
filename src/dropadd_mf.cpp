#include <Rcpp.h>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cmath>

// Matrix-free (coordinate-based) DropAdd Tabu Search for the MaxMin Diversity
// Problem (Porumbel, Hao & Glover 2011) — C++ inner loop.
//
// This is the O(N) memory counterpart of DropAdd_cpp (src/dropadd.cpp), which
// reads a fully materialised n x n distance matrix. The algorithm is identical
// in every respect — streamlined records (min_dist, sum_dist, min_dist_count),
// FIFO circular-buffer drop, lexicographic (min_dist, sum_dist) add, exclusion
// of the just-dropped x# from the add candidates for that iteration (Porumbel
// et al. 2011, p.281), and the MMDPo objective (min_pairwise + 1e-9*sum_pairwise)
// — the ONLY change is the distance source. Instead of indexing a precomputed
// column dmat[, x], each needed distance column d(., x) is recomputed from the
// `points` matrix on the fly in O(n*dim):
//     d(i, x) = sqrt( sum_c (points[i,c] - points[x,c])^2 ).
// The dense matrix is never built, so n = 58k (matrix ~27 GB, and R's
// as.matrix.dist overflows at n = 46340) becomes feasible.
//
// Per iteration this costs two O(n*dim) column passes (the dropped x# and the
// added x_new) plus O(|need_recompute| * m * dim) for the rare recompute branch
// (a point whose nearest selected peer just vanished). Memory is ~5 length-n
// arrays — strictly O(n).
//
// Construction seed (an O(n) design choice; pass DropAdd(seed=) to override).
// The matrix kernel seeds at the max-row-sum point, argmax_x sum_y d(x, y). From
// coordinates that costs O(n^2 * dim) — the very cost this matrix-free path exists
// to avoid — so here the seed is the O(n * dim) centroid-peripheral point,
//     seed = argmax_x || points[x,] - mean_x ||,
// the point farthest from the coordinate centroid. It closely tracks the max-row-sum
// point (a point far from the centroid has a large total distance to the rest; in
// practice the two rules select the same point on typical data). The seed has little
// bearing on the converged objective regardless: the greedy construction is a
// 2-approximation from ANY start and the drop-add tabu search determines the result.
// Ties break to the smallest index, as in the matrix kernel; the construction and
// search are otherwise identical, so when the two seed rules coincide the entire
// trajectory matches the matrix path bit-for-bit.
//
// Bit-faithful distances: EuclidCol below mirrors stats::dist()'s R_euclidean
// accumulation (plain-double sum of dev*dev over columns in increasing index,
// then sqrt), identical to src/maximin_points.cpp::EuclidCol. Because d is
// symmetric and a-b negates exactly (squaring kills the sign), d(i,x) computed
// here equals d(x,i): the streamlined records stay self-consistent on the fly
// exactly as the symmetric matrix keeps them consistent for the matrix kernel.

using namespace Rcpp;

// Squared-Euclidean accumulated in double over columns, matching dist's order;
// see src/maximin_points.cpp for the bit-equivalence argument and the FP-flag
// caveat (a ~/.R/Makevars with -Ofast / FMA can break bit-identity with R's
// stats; tests re-confirm per toolchain).
static inline double EuclidCol(const double* P, int nPts, int dim,
                               int i, int c) {
  double s = 0.0;
  for (int j = 0; j < dim; ++j) {
    double dev = P[i + (R_xlen_t)j * nPts] - P[c + (R_xlen_t)j * nPts];
    s += dev * dev;
  }
  return std::sqrt(s);
}

// One sweep over all points covering the NB (1..4) dimensions starting at
// column j0, accumulating each point's running squared sum in a register.
// FIRST starts the sum at the block's first squared deviation (equal to
// EuclidCol's `0.0 + dev*dev` bit-for-bit); a continuation resumes from the
// stored sum. Each point's chain is the identical left-associated sequence
// EuclidCol produces — a store/load round-trip of a double is exact — so
// the final squared value, and the sqrt applied to it below, match
// EuclidCol's exactly. The payoff is the access shape: contiguous streams
// down up to four coordinate columns at once, instead of EuclidCol's
// per-point gather across all `dim` columns (the same restructure as
// maximin_points.cpp's SweepBlock, round 9).
template <int NB, bool FIRST>
static inline void FillSqBlock(const double* P, int nPts, int c, int j0,
                               double* col) {
  const double* p0 = P + (R_xlen_t)j0 * nPts;
  const double* p1 = NB > 1 ? p0 + nPts : p0;
  const double* p2 = NB > 2 ? p1 + nPts : p1;
  const double* p3 = NB > 3 ? p2 + nPts : p2;
  const double c0 = p0[c], c1 = p1[c], c2 = p2[c], c3 = p3[c];
  for (int i = 0; i < nPts; ++i) {
    double dev = p0[i] - c0;
    double s = FIRST ? dev * dev : col[i] + dev * dev;
    if (NB > 1) { dev = p1[i] - c1; s += dev * dev; }
    if (NB > 2) { dev = p2[i] - c2; s += dev * dev; }
    if (NB > 3) { dev = p3[i] - c3; s += dev * dev; }
    col[i] = s;
  }
}

// Fill col[i] = SQUARED distance from every point to point c: dimension
// blocks of up to four via FillSqBlock. The consumer applies std::sqrt to
// each element as it reads it — the record-update loops are scalar and
// branchy anyway (the errno-guarded sqrt would only de-vectorise a loop
// that cannot vectorise), while these fills stay pure SIMD streams and no
// extra column round-trip is spent on a separate sqrt pass. Every distance
// is sqrt() of the identical squared double EuclidCol feeds it, so each is
// the same double bit-for-bit. Serves the three whole-column fill sites
// (construction seed and add, search drop and add); the recompute branch
// keeps per-pair EuclidCol (it reads scattered rows, not whole columns).
static void FillSqColumn(const double* P, int nPts, int dim, int c,
                         double* col) {
  if (dim == 0) {
    // Degenerate zero-column input: EuclidCol's sum is 0.
    for (int i = 0; i < nPts; ++i) col[i] = 0.0;
    return;
  }
  int j0 = 0;
  bool first = true;
  while (dim - j0 > 4) {
    if (first) {
      FillSqBlock<4, true>(P, nPts, c, j0, col);
    } else {
      FillSqBlock<4, false>(P, nPts, c, j0, col);
    }
    first = false;
    j0 += 4;
  }
  switch ((dim - j0 - 1) * 2 + (first ? 1 : 0)) {
    case 0: FillSqBlock<1, false>(P, nPts, c, j0, col); break;
    case 1: FillSqBlock<1, true>(P, nPts, c, j0, col); break;
    case 2: FillSqBlock<2, false>(P, nPts, c, j0, col); break;
    case 3: FillSqBlock<2, true>(P, nPts, c, j0, col); break;
    case 4: FillSqBlock<3, false>(P, nPts, c, j0, col); break;
    case 5: FillSqBlock<3, true>(P, nPts, c, j0, col); break;
    case 6: FillSqBlock<4, false>(P, nPts, c, j0, col); break;
    default: FillSqBlock<4, true>(P, nPts, c, j0, col); break;
  }
}

// [[Rcpp::export]]
List DropAdd_points_cpp(NumericMatrix points, int m, double time_budget_s,
                          int max_iter, int max_no_improve, bool want_trace,
                          int seed0 = -1) {
  const int n   = points.nrow();
  const int dim = points.ncol();
  if (m < 2 || m > n) stop("m must satisfy 2 <= m <= nrow(points)");
  const double *P = points.begin();    // column-major double storage

  // The public DropAdd() guards NA/NaN via .AsPointsMatrix(); guard here too so
  // a direct `:::DropAdd_points_cpp()` call cannot degrade the argmax (leaving
  // x_new == -1 and writing in_S[-1] out of bounds; MF-01).
  for (R_xlen_t i = 0; i < points.size(); ++i) {
    if (ISNAN(P[i])) stop("`points` must not contain NA/NaN");
  }

  const double eps = 1e-9;

  std::vector<int> S(m);
  // `int` (not unsigned char) for the membership flag: a same-behaviour 0/1
  // marker that sidesteps a GCC-14 -Wstringop-overflow false positive on the
  // byte write `in_S[x_new] = 1`.
  std::vector<int> in_S(n, 0);
  std::vector<double> min_dist(n);
  std::vector<double> sum_dist(n);
  std::vector<int> min_dist_count(n, 0);

  // Scratch column reused across the construction and the search loop, so the
  // two on-the-fly passes per iteration touch one length-n buffer apiece.
  std::vector<double> col(n);

  // -- Construction (Algorithm 1) -----------------------------------------
  // Seed: farthest point from the coordinate anti_centroid (O(n*dim) proxy for the
  // O(n^2*dim) max-row-sum seed). Ties → smallest index. See header DEVIATION.
  // seed0 >= 0 overrides the default warm-start with a caller-supplied 0-based
  // start index (R/dropadd.R::DropAdd() validates the range); seed0 = -1 (the
  // default) computes the centroid-peripheral seed below. See header note.
  int seed = 0;
  if (seed0 >= 0) {
    seed = seed0;
  } else {
    std::vector<double> anti_centroid(dim, 0.0);
    for (int j = 0; j < dim; ++j) {
      long double s = 0.0L;            // long-double accumulator: stable mean
      for (int i = 0; i < n; ++i) s += P[i + (R_xlen_t)j * n];
      anti_centroid[j] = static_cast<double>(s / static_cast<long double>(n));
    }
    double best_d2 = R_NegInf;
    for (int i = 0; i < n; ++i) {
      double s = 0.0;
      for (int j = 0; j < dim; ++j) {
        double dev = P[i + (R_xlen_t)j * n] - anti_centroid[j];
        s += dev * dev;
      }
      if (s > best_d2) { best_d2 = s; seed = i; }
    }
  }

  S[0] = seed;
  in_S[seed] = 1;
  FillSqColumn(P, n, dim, seed, col.data());
  for (int i = 0; i < n; ++i) {
    double dv = std::sqrt(col[i]);
    min_dist[i] = dv;
    sum_dist[i] = dv;
    min_dist_count[i] = 1;
  }
  min_dist[seed] = R_PosInf;          // mask self
  min_dist_count[seed] = 0;

  for (int h = 1; h < m; ++h) {
    // ADD: argmax (min_dist, sum_dist) over !in_S, ties → smallest idx.
    int x_new = -1;
    double best_md = R_NegInf, best_sd = R_NegInf;
    for (int i = 0; i < n; ++i) {
      if (in_S[i]) continue;
      double md = min_dist[i];
      if (md > best_md) {
        best_md = md;
        best_sd = sum_dist[i];
        x_new = i;
      } else if (md == best_md) {
        double sd = sum_dist[i];
        if (sd > best_sd) {
          best_sd = sd;
          x_new = i;
        }
      }
    }
    S[h] = x_new;
    in_S[x_new] = 1;

    // Update records for ADD. Recompute the SQUARED d(., x_new) column once
    // into `col`; each element's true distance is sqrt'd as it is consumed.
    FillSqColumn(P, n, dim, x_new, col.data());
    for (int i = 0; i < n; ++i) {
      double dv = std::sqrt(col[i]);
      sum_dist[i] += dv;
      if (i == x_new) continue;
      double mdi = min_dist[i];
      if (dv < mdi) {
        min_dist[i] = dv;
        min_dist_count[i] = 1;
      } else if (dv == mdi) {
        ++min_dist_count[i];
      }
    }
    // x_new's own min_dist over S[0..h-1] (reuse the cached column).
    double mn = R_PosInf;
    int cnt = 0;
    for (int j = 0; j < h; ++j) {
      double dv = std::sqrt(col[S[j]]);
      if (dv < mn) { mn = dv; cnt = 1; }
      else if (dv == mn) ++cnt;
    }
    min_dist[x_new] = mn;
    min_dist_count[x_new] = cnt;
  }

  // -- Initial objective ---------------------------------------------------
  double cur_maxmin = R_PosInf;
  long double cur_sumpair_ld = 0.0L;       // match R's sum() long-double accumulator
  for (int j = 0; j < m; ++j) {
    int sj = S[j];
    if (min_dist[sj] < cur_maxmin) cur_maxmin = min_dist[sj];
    cur_sumpair_ld += sum_dist[sj];
  }
  double cur_sumpair = static_cast<double>(cur_sumpair_ld * 0.5L);
  double best_maxmin  = cur_maxmin;
  double best_sumpair = cur_sumpair;
  double best_score   = best_maxmin + eps * best_sumpair;
  std::vector<int> best_S = S;

  // -- Drop-Add tabu search (Algorithm 2) ---------------------------------
  int head = 0;                       // 0-based: drop position
  long long iters_done = 0;
  long long no_improve = 0;                 // consecutive non-improving iterations

  std::vector<double> d_xhash(n);     // cached column for x_hash row recompute
  std::vector<int> need_recompute;
  need_recompute.reserve(32);

  std::vector<int> trace_drops, trace_adds;
  if (want_trace) {                                  // # nocov start
    if (max_iter > 0 && max_iter < (1 << 28)) {
      trace_drops.reserve(max_iter);
      trace_adds.reserve(max_iter);
    }
  }                                                  // # nocov end

  // Anytime trace (want_trace only): best T_k at each improvement, for
  // time-to-quality profiling. The construction value is logged at iteration 0.
  std::vector<int> imp_iter; std::vector<double> imp_tk;   // # nocov start
  if (want_trace) { imp_iter.push_back(0); imp_tk.push_back(best_maxmin); }  // # nocov end

  auto t0 = std::chrono::steady_clock::now();
  // chrono::now() is cheap but not free; at large n each iteration is heavy
  // (two O(n*dim) passes) so a tighter window keeps the budget honoured.
  const int check_every = 256;
  int countdown = check_every;

  // m == n leaves Add X(k) = Z - X(k) empty once x_hash is excluded, so no
  // drop-add move exists; the construction already selected every point.
  while (iters_done < max_iter && m < n) {
    // Primary stopping rule: stagnation (deterministic; matrix-free path is
    // RNG-free, index tie-breaks). With time_budget_s = Inf, reproducible.
    if (no_improve >= max_no_improve) break;
    --countdown;
    if (countdown == 0) {
      Rcpp::checkUserInterrupt();   // honour Ctrl-C / setTimeLimit() even here
      auto now = std::chrono::steady_clock::now();
      double el = std::chrono::duration<double>(now - t0).count();
      if (el >= time_budget_s) break;
      countdown = check_every;
    }

    // 1. DROP: x_hash = S[head] (FIFO via circular buffer).
    int x_hash = S[head];
    in_S[x_hash] = 0;

    // Drop pass: recompute the SQUARED d(., x_hash) column on the fly into
    // d_xhash (each element's true distance is sqrt'd as it is consumed —
    // the same double either way). Update sum_dist, decrement
    // min_dist_count for points where d(., x_hash) <= min_dist[.] (those
    // that had x_hash as a nearest peer). The cached squared column also
    // serves the x_hash self-recompute below.
    need_recompute.clear();
    {
      FillSqColumn(P, n, dim, x_hash, d_xhash.data());
      for (int i = 0; i < n; ++i) {
        double dv = std::sqrt(d_xhash[i]);
        sum_dist[i] -= dv;
        if (i == x_hash) continue;
        if (dv <= min_dist[i]) {
          if (--min_dist_count[i] == 0) need_recompute.push_back(i);
        }
      }
    }

    // Recompute min_dist/count for points whose nearest peer just vanished.
    // Iterate selected columns S[j] in the outer loop so each on-the-fly
    // column is computed once and reused across all need_recompute rows.
    if (!need_recompute.empty()) {
      const int K = static_cast<int>(need_recompute.size());
      std::vector<double> mns(K, R_PosInf);
      std::vector<int>    cnts(K, 0);
      for (int j = 0; j < m; ++j) {
        if (j == head) continue;
        const int sj = S[j];
        for (int r = 0; r < K; ++r) {
          const int xx = need_recompute[r];
          if (xx == sj) continue;     // self mask for in-S need_recompute entries
          const double dv = EuclidCol(P, n, dim, xx, sj);
          if (dv < mns[r]) { mns[r] = dv; cnts[r] = 1; }
          else if (dv == mns[r]) ++cnts[r];
        }
      }
      for (int r = 0; r < K; ++r) {
        const int xx = need_recompute[r];
        if (R_finite(mns[r])) {
          min_dist[xx] = mns[r];
          min_dist_count[xx] = cnts[r];
        } else {
          min_dist[xx] = R_PosInf;
          min_dist_count[xx] = 0;
        }
      }
    }

    // x_hash's own record: distance to surviving peers S[j != head], from
    // the cached squared column.
    {
      double mn = R_PosInf;
      int cnt = 0;
      for (int j = 0; j < m; ++j) {
        if (j == head) continue;
        double dv = std::sqrt(d_xhash[S[j]]);
        if (dv < mn) { mn = dv; cnt = 1; }
        else if (dv == mn) ++cnt;
      }
      if (R_finite(mn)) {
        min_dist[x_hash] = mn;
        min_dist_count[x_hash] = cnt;
      } else {                                       // # nocov start
        min_dist[x_hash] = R_PosInf;                // unreachable: finite d, m>=2
        min_dist_count[x_hash] = 0;
      }                                              // # nocov end
    }

    // 2. ADD: argmax (min_dist, sum_dist) over Add X(k) = Z - X(k), ties →
    // smallest idx. x_hash is excluded for this iteration (Porumbel et al.
    // 2011, p.281): the just-dropped point cannot be re-added immediately, the
    // tabu rule that prevents looping. It is eligible again next iteration,
    // once head has advanced.
    int x_new = -1;
    {
      double best_md = R_NegInf, best_sd = R_NegInf;
      for (int i = 0; i < n; ++i) {
        if (in_S[i] || i == x_hash) continue;
        double md = min_dist[i];
        if (md > best_md) {
          best_md = md;
          best_sd = sum_dist[i];
          x_new = i;
        } else if (md == best_md) {
          double sd = sum_dist[i];
          if (sd > best_sd) {
            best_sd = sd;
            x_new = i;
          }
        }
      }
    }
    in_S[x_new] = 1;

    // Add pass: recompute the SQUARED d(., x_new) column on the fly into
    // `col` (sqrt at each consumption, as in the drop pass). Same case
    // logic as ADD in construction.
    {
      FillSqColumn(P, n, dim, x_new, col.data());
      for (int i = 0; i < n; ++i) {
        double dv = std::sqrt(col[i]);
        sum_dist[i] += dv;
        if (i == x_new) continue;
        double mdi = min_dist[i];
        if (dv < mdi) {
          min_dist[i] = dv;
          min_dist_count[i] = 1;
        } else if (dv == mdi) {
          ++min_dist_count[i];
        }
      }
      // x_new's own min_dist over surviving peers (S \ {x_hash}); S[head] still
      // holds x_hash here. Reuse the cached squared column.
      double mn = R_PosInf;
      int cnt = 0;
      for (int j = 0; j < m; ++j) {
        if (j == head) continue;
        double dv = std::sqrt(col[S[j]]);
        if (dv < mn) { mn = dv; cnt = 1; }
        else if (dv == mn) ++cnt;
      }
      min_dist[x_new] = mn;
      min_dist_count[x_new] = cnt;
    }

    // Write x_new into the head slot and advance FIFO.
    S[head] = x_new;
    head = (head + 1 == m) ? 0 : head + 1;

    // 3. Test improvement of best-known MMDPo solution. long-double
    // accumulator matches R's `sum()` so the secondary objective is computed
    // the same way as the matrix kernel.
    double cm = R_PosInf;
    long double cs_ld = 0.0L;
    for (int j = 0; j < m; ++j) {
      int sj = S[j];
      if (min_dist[sj] < cm) cm = min_dist[sj];
      cs_ld += sum_dist[sj];
    }
    double cs = static_cast<double>(cs_ld * 0.5L);
    double cur_score = cm + eps * cs;
    if (cur_score > best_score) {
      best_S       = S;
      best_maxmin  = cm;
      best_sumpair = cs;
      best_score   = cur_score;
      no_improve   = 0;
      if (want_trace) {                              // # nocov start
        imp_iter.push_back(static_cast<int>(iters_done + 1));
        imp_tk.push_back(cm);
      }                                              // # nocov end
    } else {
      ++no_improve;
    }

    ++iters_done;
    if (want_trace) {                                // # nocov start
      trace_drops.push_back(x_hash + 1);   // 1-based for R
      trace_adds.push_back(x_new + 1);
    }                                                // # nocov end
  }

  // Pack output: 1-based indices for R.
  IntegerVector best_S_R(m);
  for (int j = 0; j < m; ++j) best_S_R[j] = best_S[j] + 1;

  List out = List::create(
    _["indices"]   = best_S_R,
    _["objective"] = best_maxmin,
    _["secondary"] = best_sumpair,
    _["iters"]     = iters_done
  );

  if (want_trace) {                                  // # nocov start
    IntegerVector dR(trace_drops.size()), aR(trace_adds.size());
    std::copy(trace_drops.begin(), trace_drops.end(), dR.begin());
    std::copy(trace_adds.begin(),  trace_adds.end(),  aR.begin());
    out["drops"] = dR;
    out["adds"]  = aR;
    out["imp_iter"] = wrap(imp_iter);
    out["imp_tk"]   = wrap(imp_tk);
  }                                                  // # nocov end

  return out;
}

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
// DEVIATION (documented): seed point.
// Porumbel's construction seed is argmax_x sum_y d(x, y) (max row-sum), which is
// O(n^2 * dim) — the very cost this variant avoids. We substitute a cheap
// O(n * dim) proxy: the point farthest from the coordinate centroid,
//     seed = argmax_x || points[x,] - mean_x ||,
// which approximates the peripheral max-row-sum point (a point far from the
// centroid tends to have a large sum of distances to all others). Ties break to
// the smallest index, as in the matrix kernel. The subsequent greedy max-min
// construction and the drop-add search are faithful; only the single seed
// differs, so on instances where the two seed rules coincide the entire
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

// [[Rcpp::export]]
List DropAdd_points_cpp(NumericMatrix points, int m, double time_budget_s,
                          int max_iter, int max_no_improve, bool want_trace) {
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
  // Seed: farthest point from the coordinate centroid (O(n*dim) proxy for the
  // O(n^2*dim) max-row-sum seed). Ties → smallest index. See header DEVIATION.
  int seed = 0;
  {
    std::vector<double> centroid(dim, 0.0);
    for (int j = 0; j < dim; ++j) {
      long double s = 0.0L;            // long-double accumulator: stable mean
      for (int i = 0; i < n; ++i) s += P[i + (R_xlen_t)j * n];
      centroid[j] = static_cast<double>(s / static_cast<long double>(n));
    }
    double best_d2 = R_NegInf;
    for (int i = 0; i < n; ++i) {
      double s = 0.0;
      for (int j = 0; j < dim; ++j) {
        double dev = P[i + (R_xlen_t)j * n] - centroid[j];
        s += dev * dev;
      }
      if (s > best_d2) { best_d2 = s; seed = i; }
    }
  }

  S[0] = seed;
  in_S[seed] = 1;
  for (int i = 0; i < n; ++i) {
    double dv = EuclidCol(P, n, dim, i, seed);
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

    // Update records for ADD. Recompute the d(., x_new) column once into `col`.
    for (int i = 0; i < n; ++i) col[i] = EuclidCol(P, n, dim, i, x_new);
    for (int i = 0; i < n; ++i) {
      double dv = col[i];
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
      double dv = col[S[j]];
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
  if (want_trace) {                                  // LCOV_EXCL_START
    if (max_iter > 0 && max_iter < (1 << 28)) {
      trace_drops.reserve(max_iter);
      trace_adds.reserve(max_iter);
    }
  }                                                  // LCOV_EXCL_STOP

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
      auto now = std::chrono::steady_clock::now();
      double el = std::chrono::duration<double>(now - t0).count();
      if (el >= time_budget_s) break;
      countdown = check_every;
    }

    // 1. DROP: x_hash = S[head] (FIFO via circular buffer).
    int x_hash = S[head];
    in_S[x_hash] = 0;

    // Drop pass: recompute the d(., x_hash) column on the fly into d_xhash.
    // Update sum_dist, decrement min_dist_count for points where
    // d(., x_hash) <= min_dist[.] (those that had x_hash as a nearest peer).
    // The cached column also serves the x_hash self-recompute below.
    need_recompute.clear();
    {
      for (int i = 0; i < n; ++i) d_xhash[i] = EuclidCol(P, n, dim, i, x_hash);
      for (int i = 0; i < n; ++i) {
        double dv = d_xhash[i];
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

    // x_hash's own record: distance to surviving peers S[j != head] (cached).
    {
      double mn = R_PosInf;
      int cnt = 0;
      for (int j = 0; j < m; ++j) {
        if (j == head) continue;
        double dv = d_xhash[S[j]];     // cached, single read
        if (dv < mn) { mn = dv; cnt = 1; }
        else if (dv == mn) ++cnt;
      }
      if (R_finite(mn)) {
        min_dist[x_hash] = mn;
        min_dist_count[x_hash] = cnt;
      } else {                                       // LCOV_EXCL_START
        min_dist[x_hash] = R_PosInf;                // unreachable: finite d, m>=2
        min_dist_count[x_hash] = 0;
      }                                              // LCOV_EXCL_STOP
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

    // Add pass: recompute the d(., x_new) column on the fly into `col`. Same
    // case logic as ADD in construction.
    {
      for (int i = 0; i < n; ++i) col[i] = EuclidCol(P, n, dim, i, x_new);
      for (int i = 0; i < n; ++i) {
        double dv = col[i];
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
      // holds x_hash here. Reuse the cached column.
      double mn = R_PosInf;
      int cnt = 0;
      for (int j = 0; j < m; ++j) {
        if (j == head) continue;
        double dv = col[S[j]];
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
    } else {
      ++no_improve;
    }

    ++iters_done;
    if (want_trace) {                                // LCOV_EXCL_START
      trace_drops.push_back(x_hash + 1);   // 1-based for R
      trace_adds.push_back(x_new + 1);
    }                                                // LCOV_EXCL_STOP
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

  if (want_trace) {                                  // LCOV_EXCL_START
    IntegerVector dR(trace_drops.size()), aR(trace_adds.size());
    std::copy(trace_drops.begin(), trace_drops.end(), dR.begin());
    std::copy(trace_adds.begin(),  trace_adds.end(),  aR.begin());
    out["drops"] = dR;
    out["adds"]  = aR;
  }                                                  // LCOV_EXCL_STOP

  return out;
}

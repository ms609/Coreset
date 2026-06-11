#include <Rcpp.h>
#include <vector>
#include <chrono>
#include <algorithm>

// DropAdd Tabu Search (Porumbel, Hao & Glover 2011) — C++ inner loop.
//
// Mirrors the R reference at R/competitors_dropadd.R::DropAdd(). The R
// wrapper handles:
//   * argument validation and .AsDistMatrix() coercion,
//   * routing .verify=TRUE or .trace != NULL to the R path (test scaffolding),
// so this file implements only the production fast path.
//
// Streamlined records (per Porumbel, Algs 3-4):
//   min_dist[x]       = min over y in current S, y != x, of d(x, y)
//   sum_dist[x]       = sum over y in S of d(x, y)            (d(x,x)=0 is harmless)
//   min_dist_count[x] = | { y in S, y != x : d(x,y) == min_dist[x] } |
//
// FIFO via circular buffer: S is a length-m queue, `head` indexes the oldest
// member. Each iteration drops S[head], writes x_new into that slot, and
// advances head mod m. This replaces Porumbel's iter_stamp + which.min idiom
// with an O(1) lookup.

using namespace Rcpp;

// [[Rcpp::export]]
List DropAdd_cpp(NumericMatrix dmat, int m, double time_budget_s,
                   int max_iter, int max_no_improve, bool want_trace) {
  const int n = dmat.nrow();
  if (n != dmat.ncol()) stop("dmat must be square");
  if (m < 2 || m > n) stop("m must satisfy 2 <= m <= nrow(dmat)");
  const double *dp = REAL(dmat);          // raw pointer; column-major

  // dmat(i, j) == dp[i + j * n]
  auto D = [&](int i, int j) { return dp[i + j * n]; };

  const double eps = 1e-9;

  std::vector<int> S(m);
  // `int` (not unsigned char) for the membership flag: a same-behaviour 0/1
  // marker that sidesteps a GCC-14 -Wstringop-overflow false positive on the
  // byte write `in_S[x_new] = 1`.
  std::vector<int> in_S(n, 0);
  std::vector<double> min_dist(n);
  std::vector<double> sum_dist(n);
  std::vector<int> min_dist_count(n, 0);

  // -- Construction (Algorithm 1) -----------------------------------------
  // Seed: argmax over rows of row-sum; ties broken by smallest index.
  // Accumulate row-sums column-sequentially. dmat is column-major, so the
  // naive i-outer/j-inner loop reads D(i,j)=dp[i+j*n] with stride n (a cache
  // miss per step at large n); sweeping j outer / i inner over the contiguous
  // column `dp + j*n` is sequential. Each rs[i] still accumulates over j in
  // increasing order, so the sums — and the argmax seed — are bit-identical.
  int seed = 0;
  {
    std::vector<double> rs(n, 0.0);
    for (int j = 0; j < n; ++j) {
      const double *col = dp + (std::size_t)j * n;
      for (int i = 0; i < n; ++i) rs[i] += col[i];
    }
    double best_rs = R_NegInf;
    for (int i = 0; i < n; ++i) {
      if (rs[i] > best_rs) { best_rs = rs[i]; seed = i; }
    }
  }
  S[0] = seed;
  in_S[seed] = 1;
  for (int i = 0; i < n; ++i) {
    double dv = D(i, seed);
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

    // Update records for ADD. dmat(:, x_new) is a contiguous column.
    const double *col = dp + x_new * n;
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
    // x_new's own min_dist over S[0..h-1].
    double mn = R_PosInf;
    int cnt = 0;
    for (int j = 0; j < h; ++j) {
      double dv = D(x_new, S[j]);
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
  if (want_trace) {
    if (max_iter > 0 && max_iter < (1 << 28)) {
      trace_drops.reserve(max_iter);
      trace_adds.reserve(max_iter);
    }
  }

  auto t0 = std::chrono::steady_clock::now();
  const int check_every = 1024;       // chrono::now() is cheap but not free
  int countdown = check_every;

  // m == n leaves Add X(k) = Z - X(k) empty once x_hash is excluded, so no
  // drop-add move exists; the construction already selected every point.
  while (iters_done < max_iter && m < n) {
    // Primary stopping rule: stagnation. Deterministic (no RNG, index
    // tie-breaks), so with time_budget_s = Inf the result is reproducible and
    // machine-independent.
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

    // Drop pass: full column scan. Update sum_dist, decrement min_dist_count
    // for points where d(., x_hash) <= min_dist[.] (those that had x_hash as
    // a nearest peer). Cache the column for the x_hash row recompute below.
    need_recompute.clear();
    {
      const double *col = dp + x_hash * n;
      for (int i = 0; i < n; ++i) {
        double dv = col[i];
        d_xhash[i] = dv;
        sum_dist[i] -= dv;
        if (i == x_hash) continue;
        if (dv <= min_dist[i]) {
          if (--min_dist_count[i] == 0) need_recompute.push_back(i);
        }
      }
    }

    // Recompute min_dist/count for points whose nearest peer just vanished.
    // Iterate columns S[j] (contiguous in the inner over need_recompute) so
    // every dmat read is sequential within the cached column.
    if (!need_recompute.empty()) {
      const int K = static_cast<int>(need_recompute.size());
      std::vector<double> mns(K, R_PosInf);
      std::vector<int>    cnts(K, 0);
      for (int j = 0; j < m; ++j) {
        if (j == head) continue;
        const int sj = S[j];
        const double *col = dp + sj * n;
        for (int r = 0; r < K; ++r) {
          const int xx = need_recompute[r];
          if (xx == sj) continue;     // self mask for in-S need_recompute entries
          const double dv = col[xx];
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

    // x_hash's own record: distance to surviving peers S[j != head].
    {
      double mn = R_PosInf;
      int cnt = 0;
      for (int j = 0; j < m; ++j) {
        if (j == head) continue;
        double dv = d_xhash[S[j]];    // cached, single L1 read
        if (dv < mn) { mn = dv; cnt = 1; }
        else if (dv == mn) ++cnt;
      }
      if (R_finite(mn)) {
        min_dist[x_hash] = mn;
        min_dist_count[x_hash] = cnt;
      } else {                                    // # nocov start
        min_dist[x_hash] = R_PosInf;             // unreachable: finite d, m>=2
        min_dist_count[x_hash] = 0;
      }                                           // # nocov end
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

    // Add pass: full column scan. Same case logic as ADD in construction.
    {
      const double *col = dp + x_new * n;
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
      // x_new's own min_dist over surviving peers (S \ {x_hash}).
      double mn = R_PosInf;
      int cnt = 0;
      for (int j = 0; j < m; ++j) {
        if (j == head) continue;      // S[head] still holds x_hash here
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
    // accumulator matches R's `sum()` so the R fallback and C++ port produce
    // bit-identical `secondary` (R has LDOUBLE > double on Windows x86).
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
    if (want_trace) {
      trace_drops.push_back(x_hash + 1);   // 1-based for R
      trace_adds.push_back(x_new + 1);
    }
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

  if (want_trace) {
    IntegerVector dR(trace_drops.size()), aR(trace_adds.size());
    std::copy(trace_drops.begin(), trace_drops.end(), dR.begin());
    std::copy(trace_adds.begin(),  trace_adds.end(),  aR.begin());
    out["drops"] = dR;
    out["adds"]  = aR;
  }

  return out;
}

// grasp.cpp
//
// Compiled kernel for Grasp(): GRASP with path relinking for the Max-Min
// Diversity Problem (Resende, Marti, Gallego & Duarte 2010, static variant).
//
// This mirrors the pure-R reference .Grasp_R() / .Grasp* helpers step for
// step. Crucially, the randomised construction draws its random indices from
// R's own RNG via R_unif_index(), exactly as R's sample.int(k, 1L) does, so
// that from a common set.seed() the kernel and the R reference consume the
// identical random stream and return bit-identical selections. The parity is
// asserted in tests/testthat/test-grasp.R.
//
// Termination is the deterministic stagnation rule: the refinement loop stops
// after max_no_improve consecutive iterations that do not raise the best
// elite objective. Given a seed, iteration count and result are reproducible
// and machine-independent. An optional finite time_budget_s ceiling is
// honoured but (by reintroducing wall-clock dependence) is off by default.

#include <Rcpp.h>
#include <R_ext/Random.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <chrono>
using namespace Rcpp;

static const double NEG_INF = -std::numeric_limits<double>::infinity();
static const double POS_INF =  std::numeric_limits<double>::infinity();

// d is column-major: d(i, j) = dptr[i + (size_t)j * n].
static inline double D(const double* dptr, int n, int i, int j) {
  return dptr[(size_t)i + (size_t)j * (size_t)n];
}

// Minimum pairwise distance over a (sorted) selection; matches .GraspObjective.
static double objective_of(const double* d, int n, const std::vector<int>& sel) {
  int m = (int)sel.size();
  if (m < 2) return NA_REAL;
  double best = POS_INF;
  for (int a = 0; a < m; ++a) {
    for (int b = a + 1; b < m; ++b) {
      double v = D(d, n, sel[a], sel[b]);
      if (v < best) best = v;
    }
  }
  return best;
}

// The two smallest distances from `x` to `set`, with the VERTEX achieving the
// smallest.
//
// The path-relink walk needs `min over set \ {v} of d(x, .)` for one `x` against
// every candidate excluded `v`. All of those minima are drawn from the same
// x-to-set row, so reading the row once for its two smallest entries answers
// them all in O(1) each (see near_excl), in place of a fresh O(|set|) scan per
// (x, v) pair. The incremental piece of a one-element swap: adding `x` to a base
// selection lowers its min pairwise distance to min(base_min, that minimum).
//
// The local search's candidate scan looks like the same shape but is NOT worth
// converting: it excludes only the (usually two) critical vertices, so it pays
// this row scan back barely twice, and the direct scan it would replace is a
// tighter branchless min. Measured a net loss below m ~ 50 (see T-014).
//
// The result is the exact double a direct scan of `set \ {v}` would return — it
// re-selects the very same matrix cell, with no arithmetic. The tie handling
// makes that hold: `<` on the min1 update means the FIRST minimum wins, so
// arg1 is the vertex a forward scan would settle on; a later value tying with
// min1 falls through to the `else if (v < min2)` branch, leaving min2 == min1,
// — a duplicate of the minimum survives the exclusion.
// Selections hold distinct vertices, so excluding by vertex id is well defined.
//
// arg2 names the vertex holding the min2 slot, so the walk can maintain a
// NearTwo across one-element swaps of `set` (T-017): every element other than
// the two stored witnesses sits at >= min2, so removing a non-witness leaves
// (min1, min2) the two smallest, and only a removal of arg1 or arg2 forces a
// rescan; an insertion is the classic O(1) two-smallest update.
struct NearTwo { double min1, min2; int arg1, arg2; };

static inline NearTwo near_two(const double* d, int n, int x,
                               const std::vector<int>& set) {
  NearTwo t; t.min1 = POS_INF; t.min2 = POS_INF; t.arg1 = -1; t.arg2 = -1;
  for (int k : set) {
    double v = D(d, n, x, k);
    if (v < t.min1) {
      t.min2 = t.min1; t.arg2 = t.arg1; t.min1 = v; t.arg1 = k;
    }
    else if (v < t.min2) { t.min2 = v; t.arg2 = k; }
  }
  return t;
}

// min over `set \ {v}` of d(x, .); POS_INF if that set is empty.
static inline double near_excl(const NearTwo& t, int v) {
  return v == t.arg1 ? t.min2 : t.min1;
}

// The extended-improvement tie-break count — #unordered pairs in (rem ∪ {s})
// with distance <= thr (.GraspMinPairCount on cand = rem ++ s) — split into its
// two halves. The rem×rem half is O(m^2) but does not depend on `s`, so the
// caller memoises it across the candidate scan (see the local-search loop);
// the s×rem half is the only O(m) work left per candidate. The split is an
// integer sum of the same comparisons, so it is exactly the old total
// regardless of accumulation order.
static int count_pairs_within(const double* d, int n,
                              const std::vector<int>& rem, double thr) {
  int c = 0, mm = (int)rem.size();
  for (int a = 0; a < mm; ++a)
    for (int b = a + 1; b < mm; ++b)
      if (D(d, n, rem[a], rem[b]) <= thr) ++c;
  return c;
}

static int count_to_set_le(const double* d, int n, const std::vector<int>& rem,
                           int s, double thr) {
  int c = 0;
  for (int v : rem) if (D(d, n, s, v) <= thr) ++c;
  return c;
}

// |intersection| of two ascending-sorted index vectors.
static int intersect_count(const std::vector<int>& a, const std::vector<int>& b) {
  int i = 0, j = 0, c = 0;
  int na = (int)a.size(), nb = (int)b.size();
  while (i < na && j < nb) {
    if (a[i] == b[j]) { ++c; ++i; ++j; }
    else if (a[i] < b[j]) ++i;
    else ++j;
  }
  return c;
}

// One randomised greedy construction; matches .GraspConstruct. Uses R_unif_index
// so the draw stream matches R's sample.int(k, 1L). Returns ascending sel.
static std::vector<int> grasp_construct(const double* d, int n, int m, double alpha) {
  std::vector<int> sel;
  sel.reserve(m);
  int first = (int)R_unif_index((double)n);
  sel.push_back(first);
  std::vector<double> g(n);
  for (int i = 0; i < n; ++i) g[i] = D(d, n, i, first);
  g[first] = NEG_INF;
  for (int h = 1; h < m; ++h) {
    double gmax = NEG_INF, gmin = POS_INF;
    int gmax_idx = -1;                          // first index achieving gmax
    for (int i = 0; i < n; ++i) {
      if (g[i] > NEG_INF) {
        if (g[i] > gmax) { gmax = g[i]; gmax_idx = i; }
        if (g[i] < gmin) gmin = g[i];
      }
    }
    double thresh = gmin + alpha * (gmax - gmin);
    std::vector<int> rcl;
    for (int i = 0; i < n; ++i)
      if (g[i] > NEG_INF && g[i] >= thresh) rcl.push_back(i);
    // FP rounding at the documented alpha = 1 (and deterministically for
    // alpha > 1) can push thresh just past gmax, emptying the RCL — indexing
    // rcl[0] on an empty vector would segfault. Fall back to the unique
    // greedy-best (argmax-g, first index on ties) WITHOUT an R_unif_index draw,
    // matching .GraspConstruct so the two stay bit-identical.
    if (rcl.empty()) rcl.push_back(gmax_idx);
    int pick;
    if ((int)rcl.size() == 1) {
      pick = rcl[0];
    } else {
      pick = rcl[(int)R_unif_index((double)rcl.size())];
    }
    sel.push_back(pick);
    for (int i = 0; i < n; ++i) {
      double dv = D(d, n, i, pick);
      if (dv < g[i]) g[i] = dv;
    }
    g[pick] = NEG_INF;
  }
  std::sort(sel.begin(), sel.end());
  return sel;
}

// Extended-improvement local search; matches .GraspLocalSearch. `sel` ascending.
//
// T-016: the candidate scan is screened through an incrementally-maintained
// nearest-selected summary instead of an O(m) distance scan per candidate.
// For every out-of-selection point `s`, min1[s] is the min distance from `s`
// to the CURRENT selection and arg1[s] the selected vertex realising it.
// The invariant that carries every branch below: arg1[s] always names a
// currently-selected vertex whose distance to `s` EQUALS min1[s] (a witness).
// The strict `<` in every update preserves it under FP ties — a later value
// tying min1 never displaces the stored witness, which therefore survives
// any exclusion of a DIFFERENT vertex. Consequently, for a candidate drop:
//   arg1[s] != drop  =>  min over sel \ {drop} == min1[s]   (exact, O(1));
//   arg1[s] == drop  =>  the witness is excluded — rescan `rem` for that `s`
//                        (exact, O(m); only ~(n-m)/m candidates per drop).
// Every candidate still gets its exact `nd` in the same ascending order, so
// the extended-improvement tie-break fires on exactly the same candidates
// with the same values: selections, iters and pr_calls are bit-identical to
// the direct scan. The values are exact because `min` merely re-selects a
// D() cell — min over fewer elements cannot fall, and the surviving witness
// pins it. (min1 alone is only a LOWER bound on the post-drop min; the
// arg1 test is what makes reading it exact — do not prune on min1 without it.)
static std::vector<int> grasp_local_search(const double* d, int n,
                                         std::vector<int> sel) {
  int m = (int)sel.size();
  if (m < 2) return sel;
  std::vector<char> in_sel(n, 0);
  for (int v : sel) in_sel[v] = 1;
  // Nearest-selected summary for the out-of-selection points. Initialised
  // column-sequentially (cache-friendly: d(s, v) = dptr[s + v*n] is contiguous
  // in s); processing selected vertices in ascending order with strict `<`
  // makes arg1 the first argmin, exactly as a forward row scan would settle.
  // Entries for selected points are placeholders, never read while selected;
  // a vertex leaving the selection gets a fresh scan before first use.
  std::vector<double> min1(n, POS_INF);
  std::vector<int> arg1(n, -1);
  for (int v : sel) {
    const double* col = d + (size_t)v * (size_t)n;
    for (int s = 0; s < n; ++s) {
      if (!in_sel[s] && col[s] < min1[s]) { min1[s] = col[s]; arg1[s] = v; }
    }
  }
  for (;;) {
    // sel is kept ascending. di[k] = nearest-other-selected distance.
    // (wa, wb) is a witness for the global min edge, reused to hoist base_z.
    // dstar IS gmin (the min of di is the min over all pairs), and because no
    // pair lies below the minimum, the extended-improvement pair count
    // (# pairs <= dstar) is the count of pairs AT the running min, kept by a
    // reset-on-new-min counter: a pair equal to the final minimum is scanned
    // while the running min is either above it (reset to 1) or equal to it
    // (increment), never below, so each is counted exactly once — the same
    // comparisons as a second pass over the pairs, summed in another order.
    std::vector<double> di(m, POS_INF);
    int wa = -1, wb = -1;
    double gmin = POS_INF;
    int pair_count = 0;
    for (int a = 0; a < m; ++a)
      for (int b = a + 1; b < m; ++b) {
        double v = D(d, n, sel[a], sel[b]);
        if (v < di[a]) di[a] = v;
        if (v < di[b]) di[b] = v;
        if (v < gmin) { gmin = v; wa = sel[a]; wb = sel[b]; pair_count = 1; }
        else if (v == gmin) ++pair_count;
      }
    double dstar = gmin;
    // Post-drop rescores for the (<= 2) witness drops, computed lazily at
    // most once per pass: min over pairs excluding an endpoint re-selects
    // exactly the D() cell objective_of(sel \ {endpoint}) would, and one pass
    // serves both endpoints where the old code rescored each separately.
    double excl_wa = POS_INF, excl_wb = POS_INF;
    bool have_excl = false;

    double best_dstar = dstar;
    int best_pc = pair_count, best_drop = -1, best_add = -1;
    std::vector<int> rem;
    rem.reserve(m);
    for (int ci = 0; ci < m; ++ci) {           // critical positions, ascending
      if (di[ci] > dstar) continue;
      int drop = sel[ci];
      // remaining = sel without position ci. Its internal min pairwise distance
      // (base_z) equals the global min dstar unless we dropped a witness vertex,
      // in which case the surviving min must be rescored.
      rem.clear();
      for (int t = 0; t < m; ++t) if (t != ci) rem.push_back(sel[t]);
      double base_z;
      if (rem.size() < 2) base_z = POS_INF;
      else if (drop == wa || drop == wb) {
        if (!have_excl) {
          for (int a = 0; a < m; ++a)
            for (int b = a + 1; b < m; ++b) {
              double v = D(d, n, sel[a], sel[b]);
              if (sel[a] != wa && sel[b] != wa && v < excl_wa) excl_wa = v;
              if (sel[a] != wb && sel[b] != wb && v < excl_wb) excl_wb = v;
            }
          have_excl = true;
        }
        base_z = drop == wa ? excl_wa : excl_wb;
      }
      else base_z = dstar;
      // The rem x rem half of the tie-break count depends only on (ci, thr), so
      // hold it across the scan. `thr` is `best_dstar` on every tie -- the
      // overwhelming majority of the calls that get this far -- and changes only
      // on the rare improvement, so this recomputes a handful of times per
      // critical position instead of once per candidate. The `!=` test is an
      // exact double comparison on purpose: `thr` is the same value flowing down
      // both branches, not an approximation of one, so a tolerance would be
      // wrong here.
      double memoThr = 0;
      int memoWithin = -1;
      for (int s = 0; s < n; ++s) {            // out-of-selection, ascending
        if (in_sel[s]) continue;
        // nd = min(base_z, min dist from s to remaining); identical to
        // eval_cand's dstar on (rem ∪ {s}) but O(1) via the nearest-selected
        // summary — rescan only when the stored witness IS the dropped
        // vertex. The pc tie-break is only needed when this candidate can
        // win (nd >= best_dstar).
        double cross;
        if (arg1[s] != drop) {
          cross = min1[s];
        } else {
          cross = POS_INF;
          for (int v : rem) { double vv = D(d, n, s, v); if (vv < cross) cross = vv; }
        }
        double nd = base_z < cross ? base_z : cross;
        if (nd >= best_dstar) {                // win or tie; NaN falls through
          if (memoWithin < 0 || nd != memoThr) {
            memoWithin = count_pairs_within(d, n, rem, nd);
            memoThr = nd;
          }
          int npc = memoWithin + count_to_set_le(d, n, rem, s, nd);
          if (nd > best_dstar) {
            best_dstar = nd; best_pc = npc;
            best_drop = drop; best_add = s;
          } else if (npc < best_pc) {
            best_pc = npc; best_drop = drop; best_add = s;
          }
        }
      }
    }
    if (best_drop < 0) break;
    in_sel[best_drop] = 0;
    in_sel[best_add]  = 1;
    // sel = sort(sel without best_drop, plus best_add)
    std::vector<int> next;
    next.reserve(m);
    for (int v : sel) if (v != best_drop) next.push_back(v);
    next.push_back(best_add);
    std::sort(next.begin(), next.end());
    sel.swap(next);
    // Maintain the nearest-selected summary under the swap. For the bulk of
    // candidates the stored witness survives, so the new min is
    // min(min1, d(s, best_add)) — one contiguous column read. A fresh scan is
    // needed only where no valid witness remains: the dropped vertex itself
    // (its entry is a stale placeholder from its selected tenure) and the
    // candidates whose witness WAS the dropped vertex (~(n-m)/m of them).
    // Out-of-selection entries other than best_drop are valid by induction,
    // so a stale arg1 can never alias best_drop here.
    {
      const double* col_add = d + (size_t)best_add * (size_t)n;
      for (int s = 0; s < n; ++s) {
        if (in_sel[s]) continue;
        if (s == best_drop || arg1[s] == best_drop) {
          double b = POS_INF; int a = -1;
          for (int v : sel) {
            double vv = D(d, n, s, v);
            if (vv < b) { b = vv; a = v; }
          }
          min1[s] = b; arg1[s] = a;
        } else if (col_add[s] < min1[s]) {
          min1[s] = col_add[s]; arg1[s] = best_add;
        }
      }
    }
  }
  return sel;
}

struct PRResult { std::vector<int> best; double objective; };

// Greedy path relinking from x toward y; matches .GraspPathRelink. x, y ascending.
static PRResult grasp_path_relink(const double* d, int n,
                                const std::vector<int>& x,
                                const std::vector<int>& y) {
  if (x == y) return { x, objective_of(d, n, x) };
  // to_drop = x \ y, to_add = y \ x (both ascending).
  std::vector<int> to_drop, to_add;
  {
    int i = 0, j = 0, nx = (int)x.size(), ny = (int)y.size();
    while (i < nx && j < ny) {
      if (x[i] == y[j]) { ++i; ++j; }
      else if (x[i] < y[j]) { to_drop.push_back(x[i]); ++i; }
      else { to_add.push_back(y[j]); ++j; }
    }
    while (i < nx) { to_drop.push_back(x[i]); ++i; }
    while (j < ny) { to_add.push_back(y[j]); ++j; }
  }
  int r = (int)to_drop.size();
  std::vector<int> pk = x;
  std::vector<int> best_sel = x;
  double best_z = objective_of(d, n, x);
  double z_y = objective_of(d, n, y);
  if (z_y > best_z) { best_sel = y; best_z = z_y; }

  // T-017: the walk's per-step state persists across steps instead of being
  // rebuilt from scratch each step.
  //
  // dropList/addList ARE the per-step rebuilds pk ∩ to_drop and to_add \ pk:
  // pk starts at x (⊇ to_drop, disjoint from to_add) and each step removes
  // one dropList member and adds one addList member, so erasing the chosen
  // pair keeps both lists equal to the rebuilds, in the same ascending order
  // — the pair loop enumerates identically.
  std::vector<int> dropList = to_drop;
  std::vector<int> addList  = to_add;

  // Per-member nearest-other-member summary over pk, indexed by vertex id —
  // the same witness invariant as the local search's min1/arg1: earg[v] is a
  // current member whose distance to v EQUALS edi[v], so the global min edge
  // and a witness for it read off the array in O(m), and a swap disturbs only
  // the entries whose witness left. Only members' entries are ever read;
  // a vertex entering pk is freshly scanned before first use, and dropped
  // vertices never return (adds come from to_add, disjoint from to_drop).
  std::vector<double> edi(n, POS_INF);
  std::vector<int> earg(n, -1);
  {
    int mm = (int)pk.size();
    for (int a = 0; a < mm; ++a)
      for (int b = a + 1; b < mm; ++b) {
        double v = D(d, n, pk[a], pk[b]);
        if (v < edi[pk[a]]) { edi[pk[a]] = v; earg[pk[a]] = pk[b]; }
        if (v < edi[pk[b]]) { edi[pk[b]] = v; earg[pk[b]] = pk[a]; }
      }
  }

  // Two-smallest summary per add candidate over pk (T-014's hoist), now
  // maintained across steps: each step removes one member and adds one, an
  // O(1) update unless the removed member is a stored witness (see NearTwo).
  std::vector<NearTwo> nearList(addList.size());
  for (int jj = 0; jj < (int)addList.size(); ++jj)
    nearList[jj] = near_two(d, n, addList[jj], pk);

  for (int k = 0; k < r; ++k) {
    // Global min edge of pk with a witness, off the maintained summary.
    // Dropping any vertex other than a witness endpoint leaves that edge
    // intact, so base_z == gmin; only witness drops need the post-drop
    // rescore, computed lazily at most once per step for both endpoints
    // together (re-selecting exactly the cells objective_of(pk \ {w}) would).
    int wa = -1, wb = -1;
    double gmin = POS_INF;
    for (int v : pk) {
      if (edi[v] < gmin) { gmin = edi[v]; wa = v; wb = earg[v]; }
    }
    double excl_wa = POS_INF, excl_wb = POS_INF;
    bool have_excl = false;
    double best_pair_z = NEG_INF;
    int bi = -1, bj = -1;
    for (int ii = 0; ii < (int)dropList.size(); ++ii) {
      int di_ = dropList[ii];
      double base_z;
      if ((int)pk.size() < 3) base_z = POS_INF;
      else if (di_ == wa || di_ == wb) {
        if (!have_excl) {
          int mm = (int)pk.size();
          for (int a = 0; a < mm; ++a)
            for (int b = a + 1; b < mm; ++b) {
              double v = D(d, n, pk[a], pk[b]);
              if (pk[a] != wa && pk[b] != wa && v < excl_wa) excl_wa = v;
              if (pk[a] != wb && pk[b] != wb && v < excl_wb) excl_wb = v;
            }
          have_excl = true;
        }
        base_z = di_ == wa ? excl_wa : excl_wb;
      }
      else base_z = gmin;
      for (int jj = 0; jj < (int)addList.size(); ++jj) {
        double cross = near_excl(nearList[jj], di_);
        double zc = base_z < cross ? base_z : cross; // == objective_of(pk\{di_} ∪ {dj_})
        if (zc > best_pair_z) { best_pair_z = zc; bi = di_; bj = addList[jj]; }
      }
    }
    // pk = sort(pk without bi, plus bj)
    std::vector<int> next;
    next.reserve((int)pk.size());
    for (int v : pk) if (v != bi) next.push_back(v);
    next.push_back(bj);
    std::sort(next.begin(), next.end());
    pk.swap(next);
    if (best_pair_z > best_z) { best_sel = pk; best_z = best_pair_z; }
    if (k + 1 >= r) break;                     // final step: nothing to maintain
    dropList.erase(std::find(dropList.begin(), dropList.end(), bi));
    {
      int pos = (int)(std::find(addList.begin(), addList.end(), bj)
                      - addList.begin());
      addList.erase(addList.begin() + pos);
      nearList.erase(nearList.begin() + pos);
    }
    // Maintain edi/earg under the swap: the entering vertex and any member
    // whose witness was just dropped rescan; everyone else at most adopts the
    // entering vertex (min over fewer-plus-one, exact by the witness surviving).
    for (int v : pk) {
      if (v == bj || earg[v] == bi) {
        double b2 = POS_INF; int a2 = -1;
        for (int u : pk) {
          if (u == v) continue;
          double vv = D(d, n, u, v);
          if (vv < b2) { b2 = vv; a2 = u; }
        }
        edi[v] = b2; earg[v] = a2;
      } else {
        double vv = D(d, n, v, bj);
        if (vv < edi[v]) { edi[v] = vv; earg[v] = bj; }
      }
    }
    // Maintain each remaining add candidate's two-smallest over pk: rescan
    // only if a stored witness was dropped, else the removed value sat at
    // >= min2 (or was a min-tie whose stored witness survives) and the pair
    // stands; then the classic O(1) insert of the entering vertex.
    for (int jj = 0; jj < (int)addList.size(); ++jj) {
      NearTwo& t = nearList[jj];
      if (t.arg1 == bi || t.arg2 == bi) {
        t = near_two(d, n, addList[jj], pk);
      } else {
        double w = D(d, n, addList[jj], bj);
        if (w < t.min1) {
          t.min2 = t.min1; t.arg2 = t.arg1; t.min1 = w; t.arg1 = bj;
        }
        else if (w < t.min2) { t.min2 = w; t.arg2 = bj; }
      }
    }
  }
  return { best_sel, best_z };
}

// Elite-set insertion; matches .GraspTryInsert. ES descending by ESz.
static void grasp_try_insert(std::vector<std::vector<int>>& ES,
                           std::vector<double>& ESz,
                           const std::vector<int>& sel, double sel_z, int dth) {
  int B = (int)ES.size();
  double z1 = ESz[0];
  double zb = ESz[B - 1];
  std::vector<int> hamm(B);
  int dmin = INT_MAX;
  int m = (int)sel.size();
  for (int b = 0; b < B; ++b) {
    hamm[b] = m - intersect_count(sel, ES[b]);
    if (hamm[b] < dmin) dmin = hamm[b];
  }
  bool accept = false;
  if (sel_z > z1) accept = true;
  else if (sel_z > zb && dmin >= dth) accept = true;
  if (dmin == 0) accept = false;
  if (!accept) return;
  // Eviction is restricted to members WORSE than the candidate -- Resende et al.
  // (2010) §4.1, "we remove the closest solution to x' in ES among those worse
  // than it in value". (Fig. 4 line 8 states the same step as "closest solution
  // to x' in ES with z(x') > z(x^k)", where the primed symbol is the incoming
  // solution, not the elite member.)
  // That restriction makes ESz[0] monotone: the discarded member is always
  // below the incoming one, so the pool maximum cannot fall. The pool is
  // never empty, since acceptance requires sel_z > z1 or sel_z > zb. `dmin`
  // above remains a distance to the WHOLE elite set: it is the diversity test,
  // not the eviction scan, so `dworse` here may exceed it.
  int dworse = INT_MAX;
  for (int b = 0; b < B; ++b) {
    if (ESz[b] < sel_z && hamm[b] < dworse) dworse = hamm[b];
  }
  // closest member among those worse; ties broken by lowest ESz, and among
  // equal ESz the first (lowest index) -- matches which.min().
  int closest = -1;
  double minz = POS_INF;
  for (int b = 0; b < B; ++b) {
    if (ESz[b] < sel_z && hamm[b] == dworse) {
      if (closest < 0 || ESz[b] < minz) { closest = b; minz = ESz[b]; }
    }
  }
  ES.erase(ES.begin() + closest);
  ESz.erase(ESz.begin() + closest);
  // insert keeping descending order: position = #(ESz >= sel_z)
  int pos = 0;
  for (double z : ESz) if (z >= sel_z) ++pos;
  ES.insert(ES.begin() + pos, sel);
  ESz.insert(ESz.begin() + pos, sel_z);
}

// [[Rcpp::export]]
List Grasp_cpp(NumericMatrix dmat, int m, int max_no_improve, int max_iter,
                 int elite_size, double alpha, double time_budget_s,
                 Rcpp::Nullable<Rcpp::Function> progress_cb = R_NilValue) {
  Rcpp::RNGScope scope;                       // GetRNGstate / PutRNGstate
  int n = dmat.nrow();
  const double* d = dmat.begin();
  const int dth = 5;
  const bool gated = R_finite(time_budget_s);
  auto t0 = std::chrono::steady_clock::now();
  auto elapsed = [&]() {
    return std::chrono::duration<double>(
             std::chrono::steady_clock::now() - t0).count();
  };

  // Phase A: build the initial elite set.
  std::vector<std::vector<int>> ES(elite_size);
  std::vector<double> ESz(elite_size);
  for (int b = 0; b < elite_size; ++b) {
    std::vector<int> x  = grasp_construct(d, n, m, alpha);
    std::vector<int> xp = grasp_local_search(d, n, x);
    ESz[b] = objective_of(d, n, xp);
    ES[b]  = xp;
  }
  // Stable descending sort by objective (ties keep construction order), to
  // match R's order(ESz, decreasing = TRUE) under the radix (stable) method.
  {
    std::vector<int> ord(elite_size);
    for (int i = 0; i < elite_size; ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(),
                     [&](int a, int b) { return ESz[a] > ESz[b]; });
    std::vector<std::vector<int>> ES2(elite_size);
    std::vector<double> ESz2(elite_size);
    for (int i = 0; i < elite_size; ++i) { ES2[i] = ES[ord[i]]; ESz2[i] = ESz[ord[i]]; }
    ES.swap(ES2);
    ESz.swap(ESz2);
  }

  // Phase B: refine until max_no_improve consecutive non-improving iterations
  // (deterministic), an optional iteration cap, or an optional time ceiling.
  long long iters = 0;
  long long no_improve = 0;
  double best_z_B = ESz[0];
  const int check_every = 256;        // interrupt-poll cadence (RNG-neutral)
  int countdown = check_every;
  // Optional progress callback: reports the stall counter no_improve (0 ..
  // max_no_improve) so the caller can render a bar that fills as the search
  // stagnates and snaps back to 0 whenever a better elite objective is found.
  // Read-only -- it draws no random numbers, so the RNG stream and result are
  // unaffected whether or not it is supplied.
  const bool report = progress_cb.isNotNull();
  Rcpp::Function cb = report ? Rcpp::Function(progress_cb.get())
                             : Rcpp::Function("identity");
  for (;;) {
    if (no_improve >= max_no_improve) break;
    if (iters >= max_iter) break;
    if (gated && elapsed() >= time_budget_s) break;
    if (--countdown == 0) {
      Rcpp::checkUserInterrupt();      // honour Ctrl-C / setTimeLimit() even here
      countdown = check_every;
    }
    std::vector<int> x  = grasp_construct(d, n, m, alpha);
    std::vector<int> xp = grasp_local_search(d, n, x);
    double zp = objective_of(d, n, xp);
    grasp_try_insert(ES, ESz, xp, zp, dth);
    ++iters;
    if (ESz[0] > best_z_B) { best_z_B = ESz[0]; no_improve = 0; }
    else ++no_improve;
    if (report) cb((int) no_improve);
  }

  // Phase C: path relinking over all elite pairs (deterministic; no RNG).
  std::vector<int> best_sel = ES[0];
  double best_z = ESz[0];
  int pr_calls = 0;
  int K = (int)ES.size();
  if (K >= 2 && !(gated && elapsed() >= time_budget_s)) {
    bool done = false;
    for (int i = 0; i < K - 1 && !done; ++i) {
      for (int j = i + 1; j < K; ++j) {
        PRResult pr1 = grasp_path_relink(d, n, ES[i], ES[j]);
        PRResult pr2 = grasp_path_relink(d, n, ES[j], ES[i]);
        pr_calls += 2;
        std::vector<int>& y_sel = (pr1.objective >= pr2.objective) ? pr1.best
                                                                   : pr2.best;
        std::vector<int> yp = grasp_local_search(d, n, y_sel);
        double zp = objective_of(d, n, yp);
        if (zp > best_z) { best_z = zp; best_sel = yp; }
        if (gated && elapsed() >= time_budget_s) { done = true; break; }
      }
    }
  }

  std::sort(best_sel.begin(), best_sel.end());
  IntegerVector indices(best_sel.size());
  for (int i = 0; i < (int)best_sel.size(); ++i) indices[i] = best_sel[i] + 1;

  return List::create(
    _["indices"]   = indices,
    _["objective"] = best_z,
    _["time_s"]    = elapsed(),
    _["iters"]     = iters,
    _["pr_calls"]  = pr_calls
  );
}

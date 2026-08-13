// grasp.cpp
//
// Compiled kernel for Grasp(): GRASP with path relinking for the Max-Min
// Diversity Problem (Resende, Marti, Gallego & Duarte 2010, static variant).
//
// This mirrors the pure-R reference .Grasp_R() / .Grasp* helpers step for
// step. The randomised constructions consume uniforms PRE-DRAWN from R's own
// RNG on the main thread — one batch of GRASP_BATCH x m draws per phase-B
// batch (matching R's runif() stream draw for draw) — so that from a common
// set.seed() the kernel and the R reference consume the identical random
// stream and return bit-identical selections. The parity is asserted in
// tests/testthat/test-grasp.R.
//
// T-021 determinism contract: the result depends on the seed alone. Worker
// threads never touch R's RNG (or any R API); every construction is a pure
// function of its pre-drawn uniforms; and batches merge into the elite set
// on the main thread in iteration order — so the returned selection is
// bit-identical at EVERY core count, and n_threads buys wall-clock only.
// GRASP_BATCH is part of the algorithm definition, not a tuning knob: R's
// RNG consumption (whole batches of GRASP_BATCH x m draws) depends on it.
//
// Termination is the deterministic stagnation rule: the refinement loop stops
// after max_no_improve consecutive iterations that do not raise the best
// elite objective, evaluated at merge time in iteration order (at most one
// batch of surplus constructions is computed and discarded past the stopping
// point; discards never touch the elite set). Given a seed, iteration count
// and result are reproducible and machine-independent. An optional finite
// time_budget_s ceiling is honoured but (by reintroducing wall-clock
// dependence) is off by default; when it gates, workers skip slots once the
// budget is spent, so truncation lands at slot granularity.

#include <Rcpp.h>
#include <R_ext/Random.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <algorithm>
#include <limits>
#include <chrono>
#include <utility>
using namespace Rcpp;

// Phase-B batch size: fixed, algorithm-defining (see header note). Sized for
// load balance at realistic core counts without inflating the discarded
// tail past the stopping point.
static const int GRASP_BATCH = 32;

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

// One randomised greedy construction; matches .GraspConstruct. Returns
// ascending sel.
//
// T-021: consumes m PRE-DRAWN uniforms from `u` — u[0] picks the first
// vertex, u[h] step h's RCL member — as index = floor(u * size), clamped to
// size-1 against the (never-observed) u*size == size rounding edge. Steps
// whose RCL is a singleton (or empty; see below) leave their uniform unused:
// consumption is FIXED at m draws either way, which is what lets the caller
// pre-draw whole batches on the main thread and run constructions on worker
// threads that must never touch R's RNG. floor(u * size) departs from
// R_unif_index's rejection-sampled unbiased index: the bias is O(size/2^53),
// unobservable at any real problem size, and the pure-R reference mirrors
// the floor exactly (see .GraspConstruct), so parity is preserved.
//
// T-019: the gmax/gmin/argmax scan folds into the PREVIOUS step's g-update
// pass, which already touches every live entry — setting g[pick] = NEG_INF
// before that pass makes the fold read exactly the values the standalone
// scan would, in the same ascending order (same first-index tie semantics).
// One standalone scan seeds the loop. The final step's update is skipped
// outright: g is dead once the last pick is made, so the old last pass
// computed values nobody read. Three O(n) passes per step become two, and
// the RCL buffer is reused across steps instead of reallocated.
static std::vector<int> grasp_construct(const double* d, int n, int m,
                                        double alpha, const double* u) {
  std::vector<int> sel;
  sel.reserve(m);
  int first = (int)(u[0] * n);
  if (first >= n) first = n - 1;
  sel.push_back(first);
  std::vector<double> g(n);
  for (int i = 0; i < n; ++i) g[i] = D(d, n, i, first);
  g[first] = NEG_INF;
  double gmax = NEG_INF, gmin = POS_INF;
  int gmax_idx = -1;                          // first index achieving gmax
  for (int i = 0; i < n; ++i) {
    if (g[i] > NEG_INF) {
      if (g[i] > gmax) { gmax = g[i]; gmax_idx = i; }
      if (g[i] < gmin) gmin = g[i];
    }
  }
  std::vector<int> rcl;
  rcl.reserve(n);
  for (int h = 1; h < m; ++h) {
    double thresh = gmin + alpha * (gmax - gmin);
    rcl.clear();
    for (int i = 0; i < n; ++i)
      if (g[i] > NEG_INF && g[i] >= thresh) rcl.push_back(i);
    // FP rounding at the documented alpha = 1 (and deterministically for
    // alpha > 1) can push thresh just past gmax, emptying the RCL — indexing
    // rcl[0] on an empty vector would segfault. Fall back to the unique
    // greedy-best (argmax-g, first index on ties); the step's pre-drawn
    // uniform simply goes unused, matching .GraspConstruct exactly.
    if (rcl.empty()) rcl.push_back(gmax_idx);
    int pick;
    if ((int)rcl.size() == 1) {
      pick = rcl[0];
    } else {
      int sz = (int)rcl.size();
      int j = (int)(u[h] * sz);
      if (j >= sz) j = sz - 1;
      pick = rcl[j];
    }
    sel.push_back(pick);
    if (h + 1 < m) {
      g[pick] = NEG_INF;
      const double* col = d + (size_t)pick * (size_t)n;
      gmax = NEG_INF; gmin = POS_INF; gmax_idx = -1;
      for (int i = 0; i < n; ++i) {
        double gi = g[i];
        if (gi > NEG_INF) {
          if (col[i] < gi) { gi = col[i]; g[i] = gi; }
          if (gi > gmax) { gmax = gi; gmax_idx = i; }
          if (gi < gmin) gmin = gi;
        }
      }
    }
  }
  std::sort(sel.begin(), sel.end());
  return sel;
}

// Reusable local-search workspace, allocated once per Grasp_cpp call — the
// local search runs hundreds of times per GRASP run, so per-call vectors are
// measurable. Arrays are indexed by vertex id. `in_sel` must be all-zero
// between calls; grasp_local_search restores that on exit.
struct LSScratch {
  std::vector<char> in_sel;
  std::vector<double> min1;              // candidates: nearest-selected
  std::vector<int> arg1;                 //   ... and its witness
  std::vector<double> mmin1, mmin2;      // members: two nearest other members
  std::vector<int> marg1, marg2;         //   ... with their witnesses
  std::vector<int> mcnt;                 // #partners at exactly mmin1
  explicit LSScratch(int n)
    : in_sel(n, 0), min1(n), arg1(n),
      mmin1(n), mmin2(n), marg1(n), marg2(n), mcnt(n) {}
};

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
//
// T-018: the SELECTED side gets the same treatment. Each member carries its
// two nearest other members (near_two value semantics, witnesses for both
// slots) plus mcnt, the count of partners at exactly mmin1. That one summary
// serves everything the former per-pass m²/2 pair sweep provided:
//   dstar      = min over members of mmin1 (witness: argmin member + marg1);
//   pair count = Σ mcnt over members at dstar, halved — every pair at dstar
//                is counted by both endpoints, so the sum is even and exact;
//   post-drop rescore for witness endpoint w
//              = min over members v != w of near_excl(summary[v], w), each an
//                O(1) read: marg1 and marg2 name DIFFERENT partners, so when
//                marg1[v] == w the mmin2 slot's witness survives w's removal.
// Maintenance under a swap (drop D, add A), for each surviving member v:
// a stored witness == D forces a rescan (rare); otherwise D's distance sat
// at >= mmin2, or tied mmin1 while the stored duplicate witness survives —
// values stand, and only a tied partner (d(v,D) == mmin1[v]) adjusts mcnt.
// The decrement precedes A's insertion, so it compares against the pre-swap
// minimum. Header work per pass falls m²/2 → O(m), swap upkeep is O(m) plus
// rare rescans, and one m²/2 sweep remains at entry to seed the summaries.
static std::vector<int> grasp_local_search(const double* d, int n,
                                         std::vector<int> sel, LSScratch& W) {
  int m = (int)sel.size();
  if (m < 2) return sel;
  std::vector<char>& in_sel = W.in_sel;
  std::vector<double>& min1 = W.min1;
  std::vector<int>& arg1 = W.arg1;
  std::vector<double>& mmin1 = W.mmin1;
  std::vector<double>& mmin2 = W.mmin2;
  std::vector<int>& marg1 = W.marg1;
  std::vector<int>& marg2 = W.marg2;
  std::vector<int>& mcnt = W.mcnt;
  for (int v : sel) in_sel[v] = 1;

  // Insert partner u at distance v into member x's summary: near_two's exact
  // value semantics (a min1 tie lands in the min2 slot; strict < elsewhere)
  // plus the count of partners at mmin1.
  auto member_insert = [&](int x, int u, double v) {
    if (v < mmin1[x]) {
      mmin2[x] = mmin1[x]; marg2[x] = marg1[x];
      mmin1[x] = v; marg1[x] = u; mcnt[x] = 1;
    } else {
      if (v == mmin1[x]) ++mcnt[x];
      if (v < mmin2[x]) { mmin2[x] = v; marg2[x] = u; }
    }
  };
  auto member_rescan = [&](int x) {
    mmin1[x] = POS_INF; mmin2[x] = POS_INF;
    marg1[x] = -1; marg2[x] = -1; mcnt[x] = 0;
    for (int u : sel) {
      if (u != x) member_insert(x, u, D(d, n, x, u));
    }
  };

  // Candidate summary. The first selected column initialises every entry
  // unconditionally — no reset pass on the reused scratch — and the rest
  // min-update. Column-sequential (d(s, v) = dptr[s + v*n] is contiguous in
  // s); ascending vertices with strict `<` make arg1 the first argmin,
  // exactly as a forward row scan would settle. Entries for selected points
  // are placeholders, never read while selected; a vertex leaving the
  // selection gets a fresh scan before first use.
  {
    const double* col0 = d + (size_t)sel[0] * (size_t)n;
    for (int s = 0; s < n; ++s) { min1[s] = col0[s]; arg1[s] = sel[0]; }
    for (int k = 1; k < m; ++k) {
      const double* col = d + (size_t)sel[k] * (size_t)n;
      for (int s = 0; s < n; ++s) {
        if (col[s] < min1[s]) { min1[s] = col[s]; arg1[s] = sel[k]; }
      }
    }
  }
  // Member summaries: seeded by one pair sweep, maintained thereafter.
  for (int v : sel) {
    mmin1[v] = POS_INF; mmin2[v] = POS_INF;
    marg1[v] = -1; marg2[v] = -1; mcnt[v] = 0;
  }
  for (int a = 0; a < m; ++a)
    for (int b = a + 1; b < m; ++b) {
      double v = D(d, n, sel[a], sel[b]);
      member_insert(sel[a], sel[b], v);
      member_insert(sel[b], sel[a], v);
    }

  for (;;) {
    // sel is kept ascending. dstar, its witness edge (wa, wb) and the
    // extended-improvement pair count all read off the member summaries in
    // O(m) (see the T-018 header note).
    double dstar = POS_INF;
    int wa = -1, wb = -1;
    for (int v : sel) {
      if (mmin1[v] < dstar) { dstar = mmin1[v]; wa = v; wb = marg1[v]; }
    }
    int pair_count = 0;
    for (int v : sel) if (mmin1[v] == dstar) pair_count += mcnt[v];
    pair_count /= 2;
    // Post-drop rescores for the two witness endpoints, lazily, O(m) each.
    double excl_wa = POS_INF, excl_wb = POS_INF;
    bool have_wa = false, have_wb = false;
    auto excl_of = [&](int w) {
      double b = POS_INF;
      for (int v : sel) {
        if (v == w) continue;
        double e = marg1[v] == w ? mmin2[v] : mmin1[v];
        if (e < b) b = e;
      }
      return b;
    };

    double best_dstar = dstar;
    int best_pc = pair_count, best_drop = -1, best_add = -1;
    std::vector<int> rem;
    rem.reserve(m);
    for (int ci = 0; ci < m; ++ci) {           // critical positions, ascending
      if (mmin1[sel[ci]] > dstar) continue;
      int drop = sel[ci];
      // remaining = sel without position ci. Its internal min pairwise distance
      // (base_z) equals the global min dstar unless we dropped a witness vertex,
      // in which case the surviving min must be rescored.
      rem.clear();
      for (int t = 0; t < m; ++t) if (t != ci) rem.push_back(sel[t]);
      double base_z;
      if (rem.size() < 2) base_z = POS_INF;
      else if (drop == wa) {
        if (!have_wa) { excl_wa = excl_of(wa); have_wa = true; }
        base_z = excl_wa;
      }
      else if (drop == wb) {
        if (!have_wb) { excl_wb = excl_of(wb); have_wb = true; }
        base_z = excl_wb;
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
    // Maintain the member summaries under the swap (see the T-018 header
    // note): the entering vertex and any member whose stored witness was the
    // dropped vertex rescan; everyone else adjusts mcnt if the departed
    // partner tied their minimum — tested BEFORE inserting the entrant, so
    // the comparison is against the pre-swap minimum — then takes the O(1)
    // insert of the entrant.
    {
      const double* col_drop = d + (size_t)best_drop * (size_t)n;
      const double* col_add  = d + (size_t)best_add * (size_t)n;
      for (int v : sel) {
        if (v == best_add || marg1[v] == best_drop || marg2[v] == best_drop) {
          member_rescan(v);
        } else {
          if (col_drop[v] == mmin1[v]) --mcnt[v];
          member_insert(v, best_add, col_add[v]);
        }
      }
    }
    // Maintain the candidate summary under the swap. For the bulk of
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
  // Restore the workspace's all-zero in_sel invariant for the next call.
  for (int v : sel) in_sel[v] = 0;
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
                 int n_threads = 1,
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

  // One local-search workspace per worker thread; workers touch no R API
  // (no RNG, no allocation through R, no interrupts), so construct + local
  // search + objective are safe off the main thread. Without OpenMP the
  // pragmas vanish and everything runs serially through pool[0].
  int nthr = n_threads < 1 ? 1 : n_threads;
#ifndef _OPENMP
  nthr = 1;
#endif
  std::vector<LSScratch> pool;
  pool.reserve(nthr);
  for (int t = 0; t < nthr; ++t) pool.emplace_back(n);
#ifdef _OPENMP
#define GRASP_TID omp_get_thread_num()
#else
#define GRASP_TID 0
#endif

  // Phase A: build the initial elite set — one batch of elite_size slots,
  // uniforms pre-drawn on the main thread (see the header note).
  std::vector<std::vector<int>> ES(elite_size);
  std::vector<double> ESz(elite_size);
  std::vector<double> ubuf((size_t)elite_size * (size_t)m);
  for (size_t t = 0; t < ubuf.size(); ++t) ubuf[t] = unif_rand();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nthr)
#endif
  for (int b = 0; b < elite_size; ++b) {
    LSScratch& W = pool[GRASP_TID];
    std::vector<int> x  = grasp_construct(d, n, m, alpha,
                                          &ubuf[(size_t)b * (size_t)m]);
    std::vector<int> xp = grasp_local_search(d, n, x, W);
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
  // Constructions and local searches for a whole batch run in parallel —
  // they never read the elite set — and the batch merges on the main thread
  // in iteration order, where every stopping rule is evaluated exactly as
  // the serial loop did. Slots past the stopping point are discarded without
  // touching the elite set, so the merged prefix (and hence the result) is
  // independent of both the batch size and the thread count.
  long long iters = 0;
  long long no_improve = 0;
  double best_z_B = ESz[0];
  // Optional progress callback: reports the stall counter no_improve (0 ..
  // max_no_improve) so the caller can render a bar that fills as the search
  // stagnates and snaps back to 0 whenever a better elite objective is found.
  // Read-only -- it draws no random numbers, so the RNG stream and result are
  // unaffected whether or not it is supplied. Called at merge time, on the
  // main thread only, as is Rcpp::checkUserInterrupt (once per batch).
  const bool report = progress_cb.isNotNull();
  Rcpp::Function cb = report ? Rcpp::Function(progress_cb.get())
                             : Rcpp::Function("identity");
  std::vector<std::vector<int>> bsel(GRASP_BATCH);
  std::vector<double> bz(GRASP_BATCH);
  std::vector<char> bdone(GRASP_BATCH);
  ubuf.resize((size_t)GRASP_BATCH * (size_t)m);
  bool stop_B = false;
  while (!stop_B) {
    if (no_improve >= max_no_improve) break;
    if (iters >= max_iter) break;
    if (gated && elapsed() >= time_budget_s) break;
    Rcpp::checkUserInterrupt();      // honour Ctrl-C / setTimeLimit() even here
    for (size_t t = 0; t < ubuf.size(); ++t) ubuf[t] = unif_rand();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nthr)
#endif
    for (int b = 0; b < GRASP_BATCH; ++b) {
      // Under a finite budget a worker may find the clock already spent;
      // skipping the slot bounds the overshoot at slot granularity. The
      // budget path is wall-clock-defined anyway; with the budget off,
      // every slot completes and determinism is unconditional.
      if (gated && elapsed() >= time_budget_s) { bdone[b] = 0; continue; }
      LSScratch& W = pool[GRASP_TID];
      std::vector<int> x  = grasp_construct(d, n, m, alpha,
                                            &ubuf[(size_t)b * (size_t)m]);
      bsel[b] = grasp_local_search(d, n, x, W);
      bz[b] = objective_of(d, n, bsel[b]);
      bdone[b] = 1;
    }
    for (int b = 0; b < GRASP_BATCH; ++b) {
      if (!bdone[b] || no_improve >= max_no_improve || iters >= max_iter ||
          (gated && elapsed() >= time_budget_s)) {
        stop_B = true;
        break;
      }
      grasp_try_insert(ES, ESz, bsel[b], bz[b], dth);
      ++iters;
      if (ESz[0] > best_z_B) { best_z_B = ESz[0]; no_improve = 0; }
      else ++no_improve;
      if (report) cb((int) no_improve);
    }
  }

  // Phase C: path relinking over all elite pairs (deterministic; no RNG).
  // Pairs are independent given the elite set, so they compute in parallel;
  // the reduce runs on the main thread in pair order, preserving the serial
  // first-best-wins tie behaviour. Under a finite budget workers skip pairs
  // once the clock is spent (the serial loop stopped mid-phase the same way).
  std::vector<int> best_sel = ES[0];
  double best_z = ESz[0];
  int pr_calls = 0;
  int K = (int)ES.size();
  if (K >= 2 && !(gated && elapsed() >= time_budget_s)) {
    std::vector<std::pair<int, int>> prs;
    prs.reserve((size_t)K * (K - 1) / 2);
    for (int i = 0; i < K - 1; ++i)
      for (int j = i + 1; j < K; ++j) prs.emplace_back(i, j);
    int np = (int)prs.size();
    std::vector<std::vector<int>> psel(np);
    std::vector<double> pz(np);
    std::vector<char> pdone(np, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nthr)
#endif
    for (int p = 0; p < np; ++p) {
      if (gated && elapsed() >= time_budget_s) continue;
      LSScratch& W = pool[GRASP_TID];
      PRResult pr1 = grasp_path_relink(d, n, ES[prs[p].first], ES[prs[p].second]);
      PRResult pr2 = grasp_path_relink(d, n, ES[prs[p].second], ES[prs[p].first]);
      std::vector<int>& y_sel = (pr1.objective >= pr2.objective) ? pr1.best
                                                                 : pr2.best;
      psel[p] = grasp_local_search(d, n, y_sel, W);
      pz[p] = objective_of(d, n, psel[p]);
      pdone[p] = 1;
    }
    for (int p = 0; p < np; ++p) {
      if (!pdone[p]) continue;
      pr_calls += 2;
      if (pz[p] > best_z) { best_z = pz[p]; best_sel = psel[p]; }
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

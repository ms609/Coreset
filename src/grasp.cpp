// grasp.cpp
//
// Compiled kernel for GraspPR(): GRASP with path relinking for the Max-Min
// Diversity Problem (Resende, Marti, Gallego & Duarte 2010, static variant).
//
// This mirrors the pure-R reference .GraspPR_R() / .Gpr* helpers step for
// step. Crucially, the randomised construction draws its random indices from
// R's own RNG via R_unif_index(), exactly as R's sample.int(k, 1L) does, so
// that from a common set.seed() the kernel and the R reference consume the
// identical random stream and return bit-identical selections. The parity is
// asserted in tests/testthat/test-gpr.R.
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

// Minimum pairwise distance over a (sorted) selection; matches .GprObjective.
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

// min over k in `set` of d(x, k); POS_INF if set empty. The incremental piece
// of a one-element swap: adding `x` to a base selection lowers its min pairwise
// distance to min(base_min, min_to_set(x, base)).
static inline double min_to_set(const double* d, int n, int x,
                                const std::vector<int>& set) {
  double best = POS_INF;
  for (int k : set) { double v = D(d, n, x, k); if (v < best) best = v; }
  return best;
}

// #unordered pairs in (rem ∪ {s}) with distance <= thr — the extended-
// improvement tie-break count (.GprMinPairCount on cand = rem ++ s). Only paid
// for the rare candidate that can actually win a swap (nd >= best_dstar).
// Counts the rem×rem pairs and the s×rem pairs over the same set.
static int count_pairs_le(const double* d, int n, const std::vector<int>& rem,
                          int s, double thr) {
  int c = 0, mm = (int)rem.size();
  for (int a = 0; a < mm; ++a) {
    if (D(d, n, s, rem[a]) <= thr) ++c;
    for (int b = a + 1; b < mm; ++b)
      if (D(d, n, rem[a], rem[b]) <= thr) ++c;
  }
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

// One randomised greedy construction; matches .GprConstruct. Uses R_unif_index
// so the draw stream matches R's sample.int(k, 1L). Returns ascending sel.
static std::vector<int> gpr_construct(const double* d, int n, int m, double alpha) {
  std::vector<int> sel;
  sel.reserve(m);
  int first = (int)R_unif_index((double)n);
  sel.push_back(first);
  std::vector<double> g(n);
  for (int i = 0; i < n; ++i) g[i] = D(d, n, i, first);
  g[first] = NEG_INF;
  for (int h = 1; h < m; ++h) {
    double gmax = NEG_INF, gmin = POS_INF;
    for (int i = 0; i < n; ++i) {
      if (g[i] > NEG_INF) {
        if (g[i] > gmax) gmax = g[i];
        if (g[i] < gmin) gmin = g[i];
      }
    }
    double thresh = gmin + alpha * (gmax - gmin);
    std::vector<int> rcl;
    for (int i = 0; i < n; ++i)
      if (g[i] > NEG_INF && g[i] >= thresh) rcl.push_back(i);
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

// Extended-improvement local search; matches .GprLocalSearch. `sel` ascending.
static std::vector<int> gpr_local_search(const double* d, int n,
                                         std::vector<int> sel) {
  int m = (int)sel.size();
  if (m < 2) return sel;
  std::vector<char> in_sel(n, 0);
  for (int v : sel) in_sel[v] = 1;
  for (;;) {
    // sel is kept ascending. di[k] = nearest-other-selected distance.
    std::vector<double> di(m, POS_INF);
    for (int a = 0; a < m; ++a)
      for (int b = a + 1; b < m; ++b) {
        double v = D(d, n, sel[a], sel[b]);
        if (v < di[a]) di[a] = v;
        if (v < di[b]) di[b] = v;
      }
    double dstar = POS_INF;
    for (int k = 0; k < m; ++k) if (di[k] < dstar) dstar = di[k];
    int pair_count = 0;
    for (int a = 0; a < m; ++a)
      for (int b = a + 1; b < m; ++b)
        if (D(d, n, sel[a], sel[b]) <= dstar) ++pair_count;

    double best_dstar = dstar;
    int best_pc = pair_count, best_drop = -1, best_add = -1;
    std::vector<int> rem;
    rem.reserve(m);
    for (int ci = 0; ci < m; ++ci) {           // critical positions, ascending
      if (di[ci] > dstar) continue;
      int drop = sel[ci];
      // remaining = sel without position ci. Its internal min pairwise distance
      // (base_z) is constant across all candidate adds, so compute it once.
      rem.clear();
      for (int t = 0; t < m; ++t) if (t != ci) rem.push_back(sel[t]);
      double base_z = (rem.size() >= 2) ? objective_of(d, n, rem) : POS_INF;
      for (int s = 0; s < n; ++s) {            // out-of-selection, ascending
        if (in_sel[s]) continue;
        // nd = min(base_z, min dist from s to remaining); identical to
        // eval_cand's dstar on (rem ∪ {s}) but O(m) not O(m^2). The pc tie-break
        // is only needed when this candidate can win (nd >= best_dstar).
        double cross = min_to_set(d, n, s, rem);
        double nd = base_z < cross ? base_z : cross;
        if (nd > best_dstar) {
          best_dstar = nd; best_pc = count_pairs_le(d, n, rem, s, nd);
          best_drop = drop; best_add = s;
        } else if (nd == best_dstar) {
          int npc = count_pairs_le(d, n, rem, s, nd);
          if (npc < best_pc) { best_pc = npc; best_drop = drop; best_add = s; }
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
  }
  return sel;
}

struct PRResult { std::vector<int> best; double objective; };

// Greedy path relinking from x toward y; matches .GprPathRelink. x, y ascending.
static PRResult gpr_path_relink(const double* d, int n,
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
  for (int k = 0; k < r; ++k) {
    // drop_cands = pk ∩ to_drop ; add_cands = to_add \ pk (both ascending).
    std::vector<int> drop_cands, add_cands;
    {
      int i = 0, j = 0;
      while (i < (int)pk.size() && j < (int)to_drop.size()) {
        if (pk[i] == to_drop[j]) { drop_cands.push_back(pk[i]); ++i; ++j; }
        else if (pk[i] < to_drop[j]) ++i;
        else ++j;
      }
    }
    {
      int i = 0, j = 0;
      while (i < (int)to_add.size() && j < (int)pk.size()) {
        if (to_add[i] == pk[j]) { ++i; ++j; }
        else if (to_add[i] < pk[j]) { add_cands.push_back(to_add[i]); ++i; }
        else ++j;
      }
      while (i < (int)to_add.size()) { add_cands.push_back(to_add[i]); ++i; }
    }
    double best_pair_z = NEG_INF;
    int bi = -1, bj = -1;
    std::vector<int> rem;
    rem.reserve((int)pk.size());
    for (int ii = 0; ii < (int)drop_cands.size(); ++ii) {
      int di_ = drop_cands[ii];
      // remaining = pk without di_; its internal min pairwise (base_z) is
      // constant across the add candidates, so hoist it out of the inner loop.
      rem.clear();
      for (int v : pk) if (v != di_) rem.push_back(v);
      double base_z = (rem.size() >= 2) ? objective_of(d, n, rem) : POS_INF;
      for (int jj = 0; jj < (int)add_cands.size(); ++jj) {
        int dj_ = add_cands[jj];
        double cross = min_to_set(d, n, dj_, rem);
        double zc = base_z < cross ? base_z : cross;   // == objective_of(rem ∪ {dj_})
        if (zc > best_pair_z) { best_pair_z = zc; bi = di_; bj = dj_; }
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
  }
  return { best_sel, best_z };
}

// Elite-set insertion; matches .GprTryInsert. ES descending by ESz.
static void gpr_try_insert(std::vector<std::vector<int>>& ES,
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
  // closest member: smallest Hamming distance; ties broken by lowest ESz, and
  // among equal ESz the first (lowest index) -- matches which.min().
  int closest = -1;
  double minz = POS_INF;
  for (int b = 0; b < B; ++b) {
    if (hamm[b] == dmin) {
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
List GraspPR_cpp(NumericMatrix dmat, int m, int max_no_improve, int max_iter,
                 int elite_size, double alpha, double time_budget_s) {
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
    std::vector<int> x  = gpr_construct(d, n, m, alpha);
    std::vector<int> xp = gpr_local_search(d, n, x);
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
  int iters = 0;
  int no_improve = 0;
  double best_z_B = ESz[0];
  for (;;) {
    if (no_improve >= max_no_improve) break;
    if (iters >= max_iter) break;
    if (gated && elapsed() >= time_budget_s) break;
    std::vector<int> x  = gpr_construct(d, n, m, alpha);
    std::vector<int> xp = gpr_local_search(d, n, x);
    double zp = objective_of(d, n, xp);
    gpr_try_insert(ES, ESz, xp, zp, dth);
    ++iters;
    if (ESz[0] > best_z_B) { best_z_B = ESz[0]; no_improve = 0; }
    else ++no_improve;
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
        PRResult pr1 = gpr_path_relink(d, n, ES[i], ES[j]);
        PRResult pr2 = gpr_path_relink(d, n, ES[j], ES[i]);
        pr_calls += 2;
        std::vector<int>& y_sel = (pr1.objective >= pr2.objective) ? pr1.best
                                                                   : pr2.best;
        std::vector<int> yp = gpr_local_search(d, n, y_sel);
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

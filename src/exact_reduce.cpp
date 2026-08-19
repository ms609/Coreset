// exact_reduce.cpp
//
// Decides one ExactMaxMin feasibility probe combinatorially.
//
// A probe asks whether the threshold graph G(lambda) (edges: pairs closer
// than lambda) has an independent set of size >= k -- equivalently, whether
// the complement graph H (pairs >= lambda apart) contains a k-clique.
//
// Three facts bound the search:
//
//   * Every vertex of a k-clique has H-degree >= k - 1, so iteratively
//     deleting vertices of H-degree < k - 1 -- the (k-1)-core peel -- never
//     removes a witness vertex.
//
//   * A k-clique lies wholly within one H-component.
//
//   * A proper colouring bounds the clique number (omega <= chi): a candidate
//     set greedily coloured with fewer colours than the clique still needs
//     cannot complete it.
//
// The colour bound applies at every node of the depth-first search, not only
// at the root: a node's surviving candidates are greedily coloured and
// visited in descending colour order, so the first candidate whose colour
// cannot lift the current clique to size k prunes every candidate before it.
// The search is exhaustive, so finding no clique proves that none exists.
//
// Candidate sets are 64-bit word bitmaps, making the intersection with a
// vertex's neighbourhood -- the operation the search performs at every node --
// a linear word-AND pass.
//
// All tie-breaks are by vertex index, so the output is deterministic. With
// `threads > 1` the branches of the root node are searched concurrently, and
// determinism survives because of what each verdict needs: an infeasibility
// proof exhausts every branch, so its verdict cannot depend on the order the
// branches were visited in; and a witness found by a worker thread is never
// returned -- it only establishes that one exists, and the probe is then
// re-run serially so the witness reported is the one the serial search finds.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <chrono>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Detect a pending user interrupt without longjmp-ing: R_CheckUserInterrupt
// inside R_ToplevelExec turns the jump into a FALSE return. Main R thread
// only.
static void CheckInterruptFn(void*) {
  R_CheckUserInterrupt();
}
static bool PendingInterrupt() {
  return R_ToplevelExec(CheckInterruptFn, NULL) == FALSE;
}

typedef uint64_t BitWord;
static const int kBits = 64;

// Vertices of one component, relabelled 0 .. nv-1, with bitmap adjacency.
// The search visits candidates in ascending local index within a colour
// class, so the caller's ordering of `vars` is the colouring order.
struct CliqueSearch {
  int nv;
  int nw;
  int k;
  std::vector<BitWord> adjStore;               // nv * nw; empty in a worker
  const BitWord* adj;                          // adjStore's data, or shared
  std::vector<std::vector<BitWord> > cand;     // candidate set per depth
  std::vector<std::vector<int> > order, colour;
  std::vector<BitWord> uncoloured, sameColour; // colouring scratch
  std::vector<int> cur, best;
  bool found;
  bool expired;
  long long nodes;
  std::chrono::steady_clock::time_point deadline;
  // Worker mode (stop != NULL): abandon when *stop is raised, and never touch
  // the R API. The thread that owns the R stack additionally polls for a
  // pending interrupt and raises *stop itself, recording it in *interrupted;
  // the throw happens after the threads join.
  std::atomic<bool>* stop;
  std::atomic<bool>* interrupted;
  bool pollInterrupt;
  // When set, order[0]/colour[0] hold a precomputed root colouring and
  // Expand(0) uses it instead of colouring the root itself.
  bool rootReady;

  inline void SetBit(BitWord* s, int v) const {
    s[v >> 6] |= (BitWord(1) << (v & 63));
  }
  inline void ClearBit(BitWord* s, int v) const {
    s[v >> 6] &= ~(BitWord(1) << (v & 63));
  }
  inline int FirstBit(const BitWord* s, int from) const {
    for (int w = from; w < nw; ++w) {
      if (s[w]) {
        // __builtin_ctzll is GCC/clang; MSVC does not build R packages here.
        return (w << 6) + static_cast<int>(__builtin_ctzll(s[w]));
      }
    }
    return -1;
  }

  CliqueSearch(int nv_, int k_, std::chrono::steady_clock::time_point end)
    : nv(nv_), nw((nv_ + kBits - 1) / kBits), k(k_),
      adjStore(static_cast<size_t>(nv_) * nw, 0), adj(adjStore.data()),
      cand(k_ + 1, std::vector<BitWord>(nw, 0)),
      order(k_ + 1, std::vector<int>(nv_, 0)),
      colour(k_ + 1, std::vector<int>(nv_, 0)),
      uncoloured(nw, 0), sameColour(nw, 0),
      found(false), expired(false), nodes(0), deadline(end),
      stop(NULL), interrupted(NULL), pollInterrupt(false), rootReady(false) {
    cur.reserve(k_ + 1);
  }

  // A worker sharing the master's adjacency read-only. All of its own state
  // is allocated here, before the parallel region, so nothing in the search
  // itself can throw across a thread boundary.
  CliqueSearch(const CliqueSearch& master,
               std::atomic<bool>* stop_, std::atomic<bool>* interrupted_)
    : nv(master.nv), nw(master.nw), k(master.k),
      adjStore(), adj(master.adj),
      cand(master.k + 1, std::vector<BitWord>(master.nw, 0)),
      order(master.k + 1, std::vector<int>(master.nv, 0)),
      colour(master.k + 1, std::vector<int>(master.nv, 0)),
      uncoloured(master.nw, 0), sameColour(master.nw, 0),
      found(false), expired(false), nodes(0), deadline(master.deadline),
      stop(stop_), interrupted(interrupted_), pollInterrupt(false),
      rootReady(false) {
    cur.reserve(master.k + 1);
  }

  // Greedy colouring of `set`, writing its vertices to order/colour sorted by
  // colour ascending. Colour c is a set of pairwise non-adjacent vertices, so
  // at most one of them can join any clique.
  int ColourSort(const BitWord* set, std::vector<int>& order,
                 std::vector<int>& colour) {
    std::copy(set, set + nw, uncoloured.begin());
    int idx = 0;
    int c = 0;
    for (;;) {
      const int seed = FirstBit(uncoloured.data(), 0);
      if (seed < 0) {
        break;
      }
      ++c;
      std::copy(uncoloured.begin(), uncoloured.end(), sameColour.begin());
      for (int v = seed; v >= 0; v = FirstBit(sameColour.data(), v >> 6)) {
        ClearBit(sameColour.data(), v);
        ClearBit(uncoloured.data(), v);
        const BitWord* av = &adj[static_cast<size_t>(v) * nw];
        for (int w = 0; w < nw; ++w) {
          sameColour[w] &= ~av[w];
        }
        order[idx] = v;
        colour[idx] = c;
        ++idx;
      }
    }
    return idx;
  }

  // DSATUR colouring of the full root candidate set, written to
  // order[0]/colour[0] sorted by colour ascending. Choosing each vertex by
  // saturation -- how many distinct colours its neighbours already hold --
  // typically closes a colouring in fewer colours than the one greedy pass
  // ColourSort makes, and at the root a colour saved either refutes the probe
  // before any branch opens (a chi-colouring below k bounds omega below k) or
  // shortens the suffix of branches the root loop must visit. The root pays
  // this once per probe; every deeper node keeps the cheap pass -- Round 14
  // measured re-ordering each node and the ordering cost far outran the
  // pruning it bought.
  void DSaturRoot() {
    std::vector<int> vcol(nv, 0);                // 1-based; 0 = uncoloured
    std::vector<int> sat(nv, 0), deg(nv, 0);
    // Distinct neighbour colours per vertex, one bit per colour. At most nv
    // colours exist, so the mask reuses the bitmap geometry.
    std::vector<BitWord> seen(static_cast<size_t>(nv) * nw, 0);
    for (int v = 0; v < nv; ++v) {
      const BitWord* av = &adj[static_cast<size_t>(v) * nw];
      for (int w = 0; w < nw; ++w) {
        deg[v] += static_cast<int>(__builtin_popcountll(av[w]));
      }
    }
    for (int done = 0; done < nv; ++done) {
      int v = -1;
      for (int u = 0; u < nv; ++u) {
        if (vcol[u]) {
          continue;
        }
        if (v < 0 || sat[u] > sat[v] ||
            (sat[u] == sat[v] && deg[u] > deg[v])) {
          v = u;
        }
      }
      const BitWord* sv = &seen[static_cast<size_t>(v) * nw];
      int c = 0;
      while (sv[(c >> 6)] & (BitWord(1) << (c & 63))) {
        ++c;
      }
      vcol[v] = c + 1;
      const BitWord* av = &adj[static_cast<size_t>(v) * nw];
      for (int u = 0; u < nv; ++u) {
        if (vcol[u] || !(av[u >> 6] & (BitWord(1) << (u & 63)))) {
          continue;
        }
        BitWord& word = seen[static_cast<size_t>(u) * nw + (c >> 6)];
        const BitWord bit = BitWord(1) << (c & 63);
        if (!(word & bit)) {
          word |= bit;
          ++sat[u];
        }
      }
    }
    std::vector<int>& ord = order[0];
    std::vector<int>& col = colour[0];
    for (int v = 0; v < nv; ++v) {
      ord[v] = v;
    }
    std::stable_sort(ord.begin(), ord.end(),
                     [&](int a, int b) { return vcol[a] < vcol[b]; });
    for (int i = 0; i < nv; ++i) {
      col[i] = vcol[ord[i]];
    }
    rootReady = true;
  }

  void Expand(int depth) {
    if (((++nodes) & 1023LL) == 0) {
      if (std::chrono::steady_clock::now() > deadline) {
        expired = true;
      }
      if (stop == NULL) {
        checkUserInterrupt();
      } else {
        if (pollInterrupt && PendingInterrupt()) {  // # nocov start
          interrupted->store(true);
          stop->store(true);
        }                                           // # nocov end
        if (expired || stop->load(std::memory_order_relaxed)) {
          expired = true;
          stop->store(true);
        }
      }
    }
    if (expired) {
      return;
    }
    std::vector<BitWord>& set = cand[depth];
    std::vector<int>& ord = order[depth];
    std::vector<int>& col = colour[depth];
    const int m = (depth == 0 && rootReady) ? nv
                                            : ColourSort(set.data(), ord, col);
    for (int i = m - 1; i >= 0; --i) {
      if (depth + col[i] < k) {
        return;                     // colour bound: no k-clique below here
      }
      const int v = ord[i];
      cur.push_back(v);
      if (depth + 1 >= k) {
        best = cur;
        found = true;
        return;
      }
      const BitWord* av = &adj[static_cast<size_t>(v) * nw];
      std::vector<BitWord>& next = cand[depth + 1];
      bool any = false;
      for (int w = 0; w < nw; ++w) {
        next[w] = set[w] & av[w];
        any = any || next[w];
      }
      if (any) {
        Expand(depth + 1);
        if (found || expired) {
          return;
        }
      }
      cur.pop_back();
      ClearBit(set.data(), v);
    }
  }
};

#ifdef _OPENMP
// Search one component's root branches across `threads` OpenMP threads.
//
// Branch i of the root loop is the search below root vertex ord[i], whose
// candidates are {ord[0..i-1]} & N(ord[i]) -- the serial loop reaches that
// set by clearing each visited root from the candidate bitmap, but it is
// computable directly from a prefix mask, so the branches need no sequential
// prefix and are independent. Workers share the master's adjacency read-only
// and own everything else; the only R-API call in the region is the
// interrupt poll on the thread that owns the R stack.
//
// Returns 0 when every branch was exhausted with no clique (the verdict
// "infeasible", identical to serial because exhaustion has no order), 1 when
// some thread found a witness (the caller re-runs the probe serially, so the
// witness reported is the serial one), 2 when the deadline or an interrupt
// cut the search short. `*interrupted` reports a pending user interrupt; the
// caller must throw for it after the join, and treat the search as expired.
static int RootParallel(CliqueSearch& cs, int threads, bool* interrupted) {
  // The caller has already coloured the root (DSaturRoot), so order[0] and
  // colour[0] stand ready and the serial redo will reuse them unchanged.
  std::vector<int>& ord = cs.order[0];
  std::vector<int>& col = cs.colour[0];
  const int m = cs.nv;
  // Colours ascend along `ord`, so the roots the colour bound admits --
  // those with col[i] >= k -- are the suffix from iLo up. The serial loop
  // visits exactly these before its bound breaks.
  int iLo = 0;
  while (iLo < m && col[iLo] < cs.k) {
    ++iLo;
  }
  const int nBranch = m - iLo;
  if (nBranch <= 0) {
    return 0;
  }

  // prefix[i] holds {ord[0..i-1]}: branch i's candidate pool before the
  // neighbourhood intersection.
  const int nw = cs.nw;
  std::vector<BitWord> prefix(static_cast<size_t>(m) * nw, 0);
  for (int i = 1; i < m; ++i) {
    const BitWord* prev = &prefix[static_cast<size_t>(i - 1) * nw];
    BitWord* here = &prefix[static_cast<size_t>(i) * nw];
    std::copy(prev, prev + nw, here);
    cs.SetBit(here, ord[i - 1]);
  }

  const int nT = threads > nBranch ? nBranch : threads;
  std::atomic<bool> stop(false), witness(false), interruptSeen(false);
  // Workers are built before the region so no allocation can throw inside it.
  std::vector<std::unique_ptr<CliqueSearch> > workers;
  workers.reserve(nT);
  for (int t = 0; t < nT; ++t) {
    workers.emplace_back(new CliqueSearch(cs, &stop, &interruptSeen));
  }

#pragma omp parallel num_threads(nT)
  {
    CliqueSearch& w = *workers[omp_get_thread_num()];
    w.pollInterrupt = omp_get_thread_num() == 0;
#pragma omp for schedule(dynamic)
    for (int j = 0; j < nBranch; ++j) {
      if (stop.load(std::memory_order_relaxed)) {
        continue;
      }
      // Visit in the serial loop's order, colour descending, so the branch
      // likeliest to hold a witness under the colouring heuristic goes first.
      const int i = m - 1 - j;
      const int v = ord[i];
      w.cur.clear();
      w.cur.push_back(v);
      if (1 >= w.k) {                          // # nocov start
        // Unreachable: the caller guards k >= 2.
        witness.store(true);
        stop.store(true);
        continue;
      }                                        // # nocov end
      const BitWord* pv = &prefix[static_cast<size_t>(i) * nw];
      const BitWord* av = &w.adj[static_cast<size_t>(v) * nw];
      std::vector<BitWord>& next = w.cand[1];
      bool any = false;
      for (int t = 0; t < nw; ++t) {
        next[t] = pv[t] & av[t];
        any = any || next[t];
      }
      if (any) {
        w.Expand(1);
        if (w.found) {
          w.found = false;
          witness.store(true);
          stop.store(true);
        }
      }
    }
  }

  for (int t = 0; t < nT; ++t) {
    cs.nodes += workers[t]->nodes;
  }
  if (interruptSeen.load()) {          // # nocov start
    *interrupted = true;
    return 2;
  }                                    // # nocov end
  if (witness.load()) {
    return 1;
  }
  if (std::chrono::steady_clock::now() > cs.deadline) {
    return 2;
  }
  return 0;
}
#endif

// The strict upper triangle of `d`, keeping only entries >= `lowest`: the
// candidate thresholds the search can still reach. Column-major, so a column
// is a contiguous read.
// [[Rcpp::export]]
NumericVector TriangleAtLeast_cpp(NumericMatrix d, double lowest) {
  const int n = d.nrow();
  const double* dp = REAL(d);
  std::vector<double> out;
  for (int j = 1; j < n; ++j) {
    const double* col = dp + static_cast<size_t>(j) * n;
    for (int i = 0; i < j; ++i) {
      if (col[i] >= lowest) {
        out.push_back(col[i]);
      }
    }
  }
  return NumericVector(out.begin(), out.end());
}

// Edges of the complement graph H at threshold `lambda`: the pairs at least
// `lambda` apart, as 1-based endpoint vectors with i < j. The pairs closer
// than `lambda` are the threshold graph G(lambda), whose independent sets
// these edges' cliques are.
// [[Rcpp::export]]
List EdgesAtLeast_cpp(NumericMatrix d, double lambda) {
  const int n = d.nrow();
  const double* dp = REAL(d);
  std::vector<int> hi, hj;
  for (int j = 1; j < n; ++j) {
    const double* col = dp + static_cast<size_t>(j) * n;
    for (int i = 0; i < j; ++i) {
      if (col[i] >= lambda) {
        hi.push_back(i + 1);
        hj.push_back(j + 1);
      }
    }
  }
  return List::create(_["hi"] = IntegerVector(hi.begin(), hi.end()),
                      _["hj"] = IntegerVector(hj.begin(), hj.end()));
}

// Decide one probe on the complement graph H (edge list `hi`/`hj`, 1-based,
// each pair once) over `n` vertices against target clique size `k`.
// Returns list(status, witness):
//   "feasible"     -- witness is a k-clique of H (ascending, 1-based),
//   "infeasible"   -- the search was exhaustive and found none,
//   "inconclusive" -- `maxSeconds` elapsed first.
// With `threads > 1` each component's root branches are searched
// concurrently. The verdict and witness are those of the serial search:
// infeasibility is exhaustive, so it cannot depend on visiting order, and a
// threaded witness is only a signal to re-run the probe serially. What can
// shift is where the deadline falls, and that was never deterministic.
// [[Rcpp::export]]
List ThresholdDecide_cpp(IntegerVector hi, IntegerVector hj,
                         int n, int k, double maxSeconds, int threads = 1) {
  const R_xlen_t nE = hi.size();
  const int need = k - 1;
#ifdef _OPENMP
  const int nT = threads < 1 ? 1 : threads;
#else
  (void)threads;                               // # nocov
#endif
  // One deadline for the whole probe: components share the caller's budget.
  const std::chrono::steady_clock::time_point deadline =
    std::chrono::steady_clock::now() +
    std::chrono::duration_cast<std::chrono::steady_clock::duration>(
      std::chrono::duration<double>(maxSeconds));

  // Adjacency (CSR), both directions per edge.
  std::vector<R_xlen_t> off(n + 1, 0);
  for (R_xlen_t e = 0; e < nE; ++e) {
    ++off[hi[e]];
    ++off[hj[e]];
  }
  for (int v = 0; v < n; ++v) {
    off[v + 1] += off[v];
  }
  std::vector<int> adj(2 * nE);
  {
    std::vector<R_xlen_t> pos(off.begin(), off.end() - 1);
    for (R_xlen_t e = 0; e < nE; ++e) {
      const int a = hi[e] - 1;
      const int b = hj[e] - 1;
      adj[pos[a]++] = b;
      adj[pos[b]++] = a;
    }
  }

  // Peel to the (k-1)-core. After the loop, dg[] of a surviving vertex counts
  // its surviving neighbours.
  std::vector<int> dg(n);
  for (int v = 0; v < n; ++v) {
    dg[v] = static_cast<int>(off[v + 1] - off[v]);
  }
  std::vector<char> dead(n, 0);
  std::vector<int> todo;
  for (int v = 0; v < n; ++v) {
    if (dg[v] < need) {
      dead[v] = 1;
      todo.push_back(v);
    }
  }
  while (!todo.empty()) {
    const int v = todo.back();
    todo.pop_back();
    for (R_xlen_t p = off[v]; p < off[v + 1]; ++p) {
      const int u = adj[p];
      if (!dead[u] && --dg[u] < need) {
        dead[u] = 1;
        todo.push_back(u);
      }
    }
  }

  // Components of the surviving subgraph, numbered by smallest member.
  std::vector<int> comp(n, 0);
  std::vector<int> compSize;
  int nComp = 0;
  for (int s = 0; s < n; ++s) {
    if (dead[s] || comp[s]) {
      continue;
    }
    ++nComp;
    comp[s] = nComp;
    int sz = 1;
    todo.clear();
    todo.push_back(s);
    while (!todo.empty()) {
      const int v = todo.back();
      todo.pop_back();
      for (R_xlen_t p = off[v]; p < off[v + 1]; ++p) {
        const int u = adj[p];
        if (!dead[u] && !comp[u]) {
          comp[u] = nComp;
          ++sz;
          todo.push_back(u);
        }
      }
    }
    compSize.push_back(sz);
  }

  std::vector<int> loc(n, -1);
  for (int c = 1; c <= nComp; ++c) {
    // Unreachable: every surviving vertex has >= need = k - 1 alive
    // neighbours, all within its own component, so a surviving component
    // always has >= k members.
    if (compSize[c - 1] < k) {          // # nocov start
      continue;
    }                                    // # nocov end
    // Highest surviving degree first: the colouring the search builds at each
    // node then follows Welsh-Powell order, which needs fewer colours and so
    // prunes harder.
    std::vector<int> vars;
    vars.reserve(compSize[c - 1]);
    for (int v = 0; v < n; ++v) {
      if (!dead[v] && comp[v] == c) {
        vars.push_back(v);
      }
    }
    std::stable_sort(vars.begin(), vars.end(),
                     [&](int a, int b) { return dg[a] > dg[b]; });
    const int nv = static_cast<int>(vars.size());
    for (int t = 0; t < nv; ++t) {
      loc[vars[t]] = t;
    }

    CliqueSearch cs(nv, k, deadline);
    for (int t = 0; t < nv; ++t) {
      const int u = vars[t];
      BitWord* row = &cs.adjStore[static_cast<size_t>(t) * cs.nw];
      for (R_xlen_t p = off[u]; p < off[u + 1]; ++p) {
        const int w = adj[p];
        if (!dead[w] && comp[w] == c) {
          cs.SetBit(row, loc[w]);
        }
      }
    }
    for (int t = 0; t < nv; ++t) {
      cs.SetBit(cs.cand[0].data(), t);
    }
    cs.DSaturRoot();
#ifdef _OPENMP
    if (nT > 1) {
      bool interrupted = false;
      const int outcome = RootParallel(cs, nT, &interrupted);
      if (interrupted) {
        checkUserInterrupt();                  // # nocov
      }
      if (outcome == 1) {
        // A witness exists; find the serial one. cand[0] is untouched by the
        // root driver, so this is the plain search racing what remains of
        // the deadline.
        cs.Expand(0);
      } else if (outcome == 2) {
        cs.expired = true;
      }
    } else {
      cs.Expand(0);
    }
#else
    cs.Expand(0);                              // # nocov
#endif
    for (int t = 0; t < nv; ++t) {
      loc[vars[t]] = -1;
    }

    if (cs.expired) {
      return List::create(_["status"] = "inconclusive",
                          _["witness"] = IntegerVector(0));
    }
    if (cs.found) {
      std::vector<int> w(cs.best.size());
      for (size_t t = 0; t < cs.best.size(); ++t) {
        w[t] = vars[cs.best[t]] + 1;
      }
      std::sort(w.begin(), w.end());
      return List::create(_["status"] = "feasible",
                          _["witness"] = IntegerVector(w.begin(), w.end()));
    }
  }
  return List::create(_["status"] = "infeasible",
                      _["witness"] = IntegerVector(0));
}

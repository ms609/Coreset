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
#include <tuple>
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
  // Branching-set reduction (see Absorb): 0 = colour bound only, 1 =
  // infra-chromatic, 2 = full MaxSAT propagation, 3 = hybrid of the two. Per-depth flags for the
  // candidates it removes from branching, and colour-class scratch.
  int bound;
  std::vector<std::vector<char> > absorbed;
  std::vector<BitWord> clsStore, curStore;
  std::vector<char> clsUsed, clsProp;
  std::vector<int> clsUnit, unitV, unitCls;
  std::vector<char> needU;
  int fullDepth;                               // hybrid: full propagation to here

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
      stop(NULL), interrupted(NULL), pollInterrupt(false), bound(0),
      absorbed(k_ + 1, std::vector<char>(nv_, 0)),
      clsStore(static_cast<size_t>(k_) * nw, 0),
      curStore(static_cast<size_t>(k_) * nw, 0),
      clsUsed(k_, 0), clsProp(k_, 0), clsUnit(k_, -1), unitV(k_, 0),
      unitCls(k_, 0), needU(k_, 0), fullDepth(k_) {
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
      bound(master.bound),
      absorbed(master.k + 1, std::vector<char>(master.nv, 0)),
      clsStore(static_cast<size_t>(master.k) * master.nw, 0),
      curStore(static_cast<size_t>(master.k) * master.nw, 0),
      clsUsed(master.k, 0), clsProp(master.k, 0), clsUnit(master.k, -1),
      unitV(master.k, 0), unitCls(master.k, 0), needU(master.k, 0),
      fullDepth(master.fullDepth) {
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

  // DSATUR colouring of the full root candidate set, kept as a bound only.
  // Choosing each vertex by saturation -- how many distinct colours its
  // neighbours already hold -- typically closes a colouring in fewer colours
  // than the one greedy pass ColourSort makes, and a colouring below k proves
  // no k-clique before any branch opens. The order is deliberately discarded:
  // descending it in place of the greedy order reshapes the search tree
  // unpredictably (Round 17 measured +40% on the suite's costliest
  // refutation), and re-ordering deeper nodes costs more than it prunes
  // (Round 14). The component pays one O(nv^2) pass per probe.
  int DSaturBound() const {
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
    int chi = 0;
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
      if (c + 1 > chi) {
        chi = c + 1;
      }
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
    return chi;
  }

  // Shrink the branching set by MaxSAT reasoning (Li & Quan 2010; Li, Jiang
  // & Manya 2017). The node's colour classes C_1..C_r, r = k - depth - 1,
  // are soft clauses "the clique meets C_c", and a clique meets each at most
  // once, so a clique drawn from them has at most r members -- one short of
  // what the node needs. A branching candidate v (colour > r) can join that
  // bounded pool without lifting the bound past r in either of two ways:
  //
  //   * v has no neighbour in some class, so it joins that class, which
  //     stays an independent set (re-colouring);
  //   * unit propagation from {v} empties a class: each class left with a
  //     single neighbour of everything chosen so far forces that vertex, and
  //     once some class holds no common neighbour, no clique contains v and
  //     meets every class involved. Those classes plus {v} form an
  //     inconsistent set of soft clauses, which together contribute at most
  //     one fewer member than their count.
  //
  // Inconsistent sets are kept disjoint (a class in one is never reused or
  // enlarged), so the absorbed candidates and the r classes still yield at
  // most r clique members, and the absorbed candidates need not be branched
  // on: every k-clique below this node contains a candidate that remains.
  // Does `x` lie outside N(e)? Such an `e` has eliminated `x` from a class.
  inline bool Misses(int e, int x) const {
    return !(adj[static_cast<size_t>(e) * nw + (x >> 6)] &
             (BitWord(1) << (x & 63)));
  }

  // The first eliminator of `x` among v and the first `limit` forced
  // vertices: -1 for v, t for unitV[t], -2 if none (never, for a vertex
  // that a completed propagation left out of its class).
  int FirstEliminator(int v, int x, int limit) const {
    if (Misses(v, x)) {
      return -1;
    }
    for (int t = 0; t < limit; ++t) {
      if (Misses(unitV[t], x)) {
        return t;
      }
    }
    return -2;                                   // # nocov
  }

  // Try to take branching candidate v into the bounded pool, against the
  // classes not yet used. mode 1: re-colour or infra-chromatic only; mode
  // 2: re-colour or full unit propagation, trimming the conflict to the
  // classes that caused it when `trim` is set. Returns true if v is taken.
  bool TryAbsorb(int v, int r, int mode, bool trim) {
    std::fill(clsProp.begin(), clsProp.begin() + r, 0);
    // Each pass intersects the live classes with one forced vertex's
    // neighbourhood and, in the same sweep, counts what survives (capped at
    // two: empty, unit or more), so a pass costs one read of the classes.
    // The first pass intersects with N(v) itself: a class it empties holds
    // no neighbour of v, so v re-colours into it.
    const BitWord* by = &adj[static_cast<size_t>(v) * nw];
    int conflict = -1;
    int nu = 0;
    bool first = true;
    for (;;) {
      int unitClass = -1;
      int unitVertex = -1;
      for (int c = 0; c < r; ++c) {
        if (clsUsed[c] || clsProp[c]) {
          continue;
        }
        BitWord* cu = &curStore[static_cast<size_t>(c) * nw];
        const BitWord* src = first ? &clsStore[static_cast<size_t>(c) * nw]
                                   : cu;
        int pc = 0;
        int at = -1;
        for (int w = 0; w < nw; ++w) {
          const BitWord x = src[w] & by[w];
          cu[w] = x;
          if (x && pc < 2) {
            if (at < 0) {
              at = (w << 6) + static_cast<int>(__builtin_ctzll(x));
            }
            pc += static_cast<int>(__builtin_popcountll(x));
          }
        }
        clsUnit[c] = pc == 1 ? at : -1;
        if (pc == 0) {
          conflict = c;
          if (first) {
            break;
          }
        } else if (pc == 1 && unitClass < 0) {
          unitClass = c;
          unitVertex = at;
        }
      }
      if (first && conflict >= 0) {
        SetBit(&clsStore[static_cast<size_t>(conflict) * nw], v);
        return true;
      }
      if (first && mode == 1) {
        // Infra-chromatic (San Segundo et al. 2015): one forced vertex u
        // at a time, from each class holding a single neighbour of v; a
        // second class with no common neighbour of u and v makes the
        // three soft clauses {v}, u's class and that class inconsistent.
        for (int i = 0; i < r; ++i) {
          if (clsUsed[i] || clsUnit[i] < 0) {
            continue;
          }
          const BitWord* au = &adj[static_cast<size_t>(clsUnit[i]) * nw];
          for (int c = 0; c < r; ++c) {
            if (c == i || clsUsed[c]) {
              continue;
            }
            const BitWord* cu = &curStore[static_cast<size_t>(c) * nw];
            bool empty = true;
            for (int w = 0; w < nw; ++w) {
              if (cu[w] & au[w]) {
                empty = false;
                break;
              }
            }
            if (empty) {
              clsUsed[i] = 1;
              clsUsed[c] = 1;
              return true;
            }
          }
        }
        return false;
      }
      first = false;
      if (conflict >= 0 || unitClass < 0) {
        break;
      }
      clsProp[unitClass] = 1;
      unitV[nu] = unitVertex;
      unitCls[nu] = unitClass;
      ++nu;
      by = &adj[static_cast<size_t>(unitVertex) * nw];
    }
    if (conflict < 0) {
      return false;
    }
    clsUsed[conflict] = 1;
    if (!trim) {
      for (int t = 0; t < nu; ++t) {
        clsUsed[unitCls[t]] = 1;
      }
      return true;
    }
    // Keep only the forced vertices the conflict depends on: each member
    // of the emptied class was eliminated by v or by some forced vertex,
    // and each needed forced vertex's own class was reduced to it by
    // earlier eliminators. Their classes, the emptied one and {v} are
    // inconsistent on their own; the other classes stay free.
    std::fill(needU.begin(), needU.begin() + nu, 0);
    bool ok = true;
    const BitWord* cs = &clsStore[static_cast<size_t>(conflict) * nw];
    for (int x = FirstBit(cs, 0); x >= 0 && ok; ) {
      const int e = FirstEliminator(v, x, nu);
      if (e >= 0) {
        needU[e] = 1;
      } else if (e == -2) {
        ok = false;                              // # nocov
      }
      const int nx = x + 1;
      x = nx < nv ? FirstBitFrom(cs, nx) : -1;
    }
    for (int t = nu - 1; t >= 0 && ok; --t) {
      if (!needU[t]) {
        continue;
      }
      const BitWord* ct = &clsStore[static_cast<size_t>(unitCls[t]) * nw];
      for (int x = FirstBit(ct, 0); x >= 0 && ok; ) {
        if (x != unitV[t]) {
          const int e = FirstEliminator(v, x, t);
          if (e >= 0) {
            needU[e] = 1;
          } else if (e == -2) {
            ok = false;                          // # nocov
          }
        }
        const int nx = x + 1;
        x = nx < nv ? FirstBitFrom(ct, nx) : -1;
      }
    }
    for (int t = 0; t < nu; ++t) {
      if (needU[t] || !ok) {
        clsUsed[unitCls[t]] = 1;
      }
    }
    return true;
  }

  // First set bit at or after bit `from`.
  inline int FirstBitFrom(const BitWord* s, int from) const {
    int w = from >> 6;
    const BitWord head = s[w] & (~BitWord(0) << (from & 63));
    if (head) {
      return (w << 6) + static_cast<int>(__builtin_ctzll(head));
    }
    return FirstBit(s, w + 1);
  }

  void Absorb(int depth, int m) {
    const std::vector<int>& ord = order[depth];
    const std::vector<int>& col = colour[depth];
    std::vector<char>& abs = absorbed[depth];
    std::fill(abs.begin(), abs.begin() + m, 0);
    const int r = k - depth - 1;
    int iB = 0;
    while (iB < m && col[iB] <= r) {
      ++iB;
    }
    if (r < 1 || iB >= m) {
      return;
    }
    std::fill(clsStore.begin(), clsStore.begin() + static_cast<size_t>(r) * nw,
              BitWord(0));
    for (int i = 0; i < iB; ++i) {
      SetBit(&clsStore[static_cast<size_t>(col[i] - 1) * nw], ord[i]);
    }
    std::fill(clsUsed.begin(), clsUsed.begin() + r, 0);
    if (bound < 3) {
      for (int j = iB; j < m; ++j) {
        abs[j] = TryAbsorb(ord[j], r, bound, false);
      }
      return;
    }
    // Hybrid: every candidate first meets the cheap infra-chromatic test,
    // whose conflicts spend exactly two classes; the full propagation then
    // works through what is left with the classes still free, and trims
    // each conflict to its causes so it spends no more than it needs.
    // Shallow nodes only: the subtree an absorption saves shrinks with depth
    // while the propagation's cost does not.
    for (int j = iB; j < m; ++j) {
      abs[j] = TryAbsorb(ord[j], r, 1, false);
    }
    if (depth > fullDepth) {
      return;
    }
    for (int j = iB; j < m; ++j) {
      if (!abs[j]) {
        abs[j] = TryAbsorb(ord[j], r, 2, true);
      }
    }
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
    const int m = ColourSort(set.data(), ord, col);
    if (bound > 0 && m > 0 && depth + col[m - 1] >= k) {
      Absorb(depth, m);
    }
    for (int i = m - 1; i >= 0; --i) {
      if (depth + col[i] < k) {
        return;                     // colour bound: no k-clique below here
      }
      if (bound > 0 && absorbed[depth][i]) {
        continue;                   // bounded with the classes; not a branch
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
  std::vector<int>& ord = cs.order[0];
  std::vector<int>& col = cs.colour[0];
  const int m = cs.ColourSort(cs.cand[0].data(), ord, col);
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

// Is some live vertex v interchangeable for `u` -- non-adjacent to it, with
// every live neighbour of u among its own? Any clique through u then swaps
// u for v at equal size, so u can be discarded without losing a maximum
// clique. (The closed-neighbourhood form over adjacent pairs, natural for
// independent sets, is unsound here: it eats a triangle.) Candidates come
// from one neighbourhood, not all pairs: a dominator is adjacent to every
// live neighbour of u, so in particular to the first one found.
static bool Dominated(const std::vector<BitWord>& adjBits,
                      const std::vector<BitWord>& alive, int nw, int u) {
  const BitWord* au = &adjBits[static_cast<size_t>(u) * nw];
  int w0 = -1;
  for (int w = 0; w < nw; ++w) {
    const BitWord b = au[w] & alive[w];
    if (b) {
      w0 = (w << 6) + static_cast<int>(__builtin_ctzll(b));
      break;
    }
  }
  // Unreachable: the peel enters the sweep with every live degree >= k - 1
  // >= 1, and a dominance removal cannot take u's last live neighbour --
  // the dominator that removed it is adjacent to everything it was adjacent
  // to, u included, and still live when it dominates.
  if (w0 < 0) {                          // # nocov start
    return false;
  }                                      // # nocov end
  const BitWord* a0 = &adjBits[static_cast<size_t>(w0) * nw];
  for (int w = 0; w < nw; ++w) {
    BitWord cb = a0[w] & alive[w] & ~au[w];
    if (w == (u >> 6)) {
      cb &= ~(BitWord(1) << (u & 63));
    }
    while (cb) {
      const int v = (w << 6) + static_cast<int>(__builtin_ctzll(cb));
      cb &= cb - 1;
      const BitWord* av = &adjBits[static_cast<size_t>(v) * nw];
      bool subset = true;
      for (int x = 0; x < nw; ++x) {
        if (au[x] & alive[x] & ~av[x]) {
          subset = false;
          break;
        }
      }
      if (subset) {
        return true;
      }
    }
  }
  return false;
}

// Shrink a component to a fixpoint of two rules that cannot lose a maximum
// clique: the (k-1)-core peel, and vertex dominance. Each rule feeds the
// other -- a discarded vertex lowers its neighbours' degrees, and a peeled
// vertex shrinks the neighbourhoods the subset test compares -- so they
// alternate until a full dominance sweep removes nothing. The sweep visits
// vertices in component order (degree descending) and removals apply
// immediately, so equal-neighbourhood twins lose exactly one member and the
// result is deterministic. `alive` is trimmed in place.
static void PackReduce(const std::vector<BitWord>& adjBits,
                       std::vector<BitWord>& alive, int nv, int nw, int k) {
  const int need = k - 1;
  for (;;) {
    bool peeled = true;
    while (peeled) {
      peeled = false;
      for (int u = 0; u < nv; ++u) {
        if (!(alive[u >> 6] & (BitWord(1) << (u & 63)))) {
          continue;
        }
        const BitWord* au = &adjBits[static_cast<size_t>(u) * nw];
        int deg = 0;
        for (int w = 0; w < nw; ++w) {
          deg += static_cast<int>(__builtin_popcountll(au[w] & alive[w]));
        }
        if (deg < need) {
          alive[u >> 6] &= ~(BitWord(1) << (u & 63));
          peeled = true;
        }
      }
    }
    bool removed = false;
    for (int u = 0; u < nv; ++u) {
      if (!(alive[u >> 6] & (BitWord(1) << (u & 63)))) {
        continue;
      }
      if (Dominated(adjBits, alive, nw, u)) {
        alive[u >> 6] &= ~(BitWord(1) << (u & 63));
        removed = true;
      }
    }
    if (!removed) {
      return;
    }
  }
}

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
                         int n, int k, double maxSeconds, int threads = 1,
                         int bound = 2, int fullDepth = -1) {
  const R_xlen_t nE = hi.size();
  double nodes = 0;                            // search nodes, all components
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
    std::sort(vars.begin(), vars.end(),
              [&](const int& a, const int& b) {
                return std::tie(dg[a], b) > std::tie(dg[b], a);
              });
    const int nv0 = static_cast<int>(vars.size());
    for (int t = 0; t < nv0; ++t) {
      loc[vars[t]] = t;
    }

    // Peel and dominance to a fixpoint before any search structure is
    // built: on the probes that carry a cell's cost, the peel alone bites
    // nowhere and dominance removes up to seven eighths of the edges
    // (Round 17's audit). Vertices it discards keep the maximum clique
    // reachable through an interchangeable survivor, so the verdict is
    // unchanged; a witness is a clique of the reduced graph, which is a
    // clique of H.
    const int nw0 = (nv0 + kBits - 1) / kBits;
    std::vector<BitWord> tadj(static_cast<size_t>(nv0) * nw0, 0);
    for (int t = 0; t < nv0; ++t) {
      const int u = vars[t];
      BitWord* row = &tadj[static_cast<size_t>(t) * nw0];
      for (R_xlen_t p = off[u]; p < off[u + 1]; ++p) {
        const int w = adj[p];
        if (loc[w] >= 0) {
          row[loc[w] >> 6] |= BitWord(1) << (loc[w] & 63);
        }
      }
    }
    std::vector<BitWord> aliveBits(nw0, 0);
    for (int t = 0; t < nv0; ++t) {
      aliveBits[t >> 6] |= BitWord(1) << (t & 63);
    }
    PackReduce(tadj, aliveBits, nv0, nw0, k);
    {
      int kept = 0;
      for (int t = 0; t < nv0; ++t) {
        if (aliveBits[t >> 6] & (BitWord(1) << (t & 63))) {
          vars[kept++] = vars[t];
        } else {
          loc[vars[t]] = -1;
        }
      }
      vars.resize(kept);
    }
    const int nv = static_cast<int>(vars.size());
    // The reduction ends on a peel, so any surviving vertex has >= k - 1
    // live neighbours and any non-empty remnant has >= k members: `vars` is
    // either big enough to search or empty, with no locs left to reset.
    if (nv < k) {
      continue;
    }
    for (int t = 0; t < nv; ++t) {
      loc[vars[t]] = t;
    }

    CliqueSearch cs(nv, k, deadline);
    cs.bound = bound;
    if (fullDepth >= 0) {
      cs.fullDepth = fullDepth;
    }
    for (int t = 0; t < nv; ++t) {
      const int u = vars[t];
      BitWord* row = &cs.adjStore[static_cast<size_t>(t) * cs.nw];
      for (R_xlen_t p = off[u]; p < off[u + 1]; ++p) {
        const int w = adj[p];
        if (loc[w] >= 0) {
          cs.SetBit(row, loc[w]);
        }
      }
    }
    for (int t = 0; t < nv; ++t) {
      cs.SetBit(cs.cand[0].data(), t);
    }
    // A DSATUR colouring below k refutes the component before any branch
    // opens; when it cannot, the search runs with its own greedy order --
    // the DSATUR order is a bound's by-product, not a descent order
    // (Round 17: descending it grew the costliest refutation tree by 40%).
    if (cs.DSaturBound() < k) {
      for (int t = 0; t < nv; ++t) {
        loc[vars[t]] = -1;
      }
      continue;
    }
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
    nodes += static_cast<double>(cs.nodes);

    if (cs.expired) {
      return List::create(_["status"] = "inconclusive",
                          _["witness"] = IntegerVector(0),
                          _["nodes"] = nodes);
    }
    if (cs.found) {
      std::vector<int> w(cs.best.size());
      for (size_t t = 0; t < cs.best.size(); ++t) {
        w[t] = vars[cs.best[t]] + 1;
      }
      std::sort(w.begin(), w.end());
      return List::create(_["status"] = "feasible",
                          _["witness"] = IntegerVector(w.begin(), w.end()),
                          _["nodes"] = nodes);
    }
  }
  return List::create(_["status"] = "infeasible",
                      _["witness"] = IntegerVector(0),
                      _["nodes"] = nodes);
}

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
// All tie-breaks are by vertex index, so the output is deterministic.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <chrono>
using namespace Rcpp;

typedef uint64_t BitWord;
static const int kBits = 64;

// Vertices of one component, relabelled 0 .. nv-1, with bitmap adjacency.
// The search visits candidates in ascending local index within a colour
// class, so the caller's ordering of `vars` is the colouring order.
struct CliqueSearch {
  int nv;
  int nw;
  int k;
  std::vector<BitWord> adj;                    // nv * nw
  std::vector<std::vector<BitWord> > cand;     // candidate set per depth
  std::vector<std::vector<int> > order, colour;
  std::vector<BitWord> uncoloured, sameColour; // colouring scratch
  std::vector<int> cur, best;
  bool found;
  bool expired;
  long long nodes;
  std::chrono::steady_clock::time_point deadline;

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
      adj(static_cast<size_t>(nv_) * nw, 0),
      cand(k_ + 1, std::vector<BitWord>(nw, 0)),
      order(k_ + 1, std::vector<int>(nv_, 0)),
      colour(k_ + 1, std::vector<int>(nv_, 0)),
      uncoloured(nw, 0), sameColour(nw, 0),
      found(false), expired(false), nodes(0), deadline(end) {}

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

  void Expand(int depth) {
    if (((++nodes) & 1023LL) == 0) {
      if (std::chrono::steady_clock::now() > deadline) {
        expired = true;
      }
      checkUserInterrupt();
    }
    if (expired) {
      return;
    }
    std::vector<BitWord>& set = cand[depth];
    std::vector<int>& ord = order[depth];
    std::vector<int>& col = colour[depth];
    const int m = ColourSort(set.data(), ord, col);
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
// [[Rcpp::export]]
List ThresholdDecide_cpp(IntegerVector hi, IntegerVector hj,
                         int n, int k, double maxSeconds) {
  const R_xlen_t nE = hi.size();
  const int need = k - 1;
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
    if (compSize[c - 1] < k) {
      continue;
    }
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
      BitWord* row = &cs.adj[static_cast<size_t>(t) * cs.nw];
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
    cs.Expand(0);
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

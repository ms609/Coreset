// kcentre_cover.cpp
//
// Decides one ExactKCentre feasibility probe combinatorially -- the covering
// dual of exact_reduce.cpp's node-packing decision, and the reason neither
// solver needs a MILP backend.
//
// A probe asks whether `k` centres can cover every point within radius `r`.
// Since `d(i, j) <= r` is symmetric, the coverage incidence IS the closed
// neighbourhood of the threshold graph G(r): centre i covers point j exactly
// when i and j are adjacent (or equal). The probe is therefore "does G(r)
// admit a dominating set of size <= k", and one array of neighbourhood
// bitmaps serves both roles.
//
// Four facts bound the search:
//
//   * A point whose only available centre is c forces c open. The cascade
//     this starts settles most probes outright.
//
//   * If N[a] is a subset of N[b] then covering a covers b, so point b is
//     redundant. Dually, a centre covering a subset of another centre's
//     points is never needed. Both rules run to a fixpoint, and between them
//     they shrink a several-hundred-point probe to a few dozen.
//
//   * Two points sharing no available centre need separate centres, so a
//     greedy set of pairwise-disjoint neighbourhoods lower-bounds the cover.
//
//   * A cover lies within one component of the coverage graph, and the
//     components' minimum covers add, so each is solved against the budget
//     the earlier ones left.
//
// The branch is complete: some centre covering the least-covered uncovered
// point must be open, so branching over that point's covers reaches every
// cover. The search is exhaustive, so finding none proves that none exists.
//
// All tie-breaks are by index, so the verdict and the witness are
// deterministic.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <chrono>
using namespace Rcpp;

typedef uint64_t BitWord;

namespace {

inline void SetBit(BitWord* s, int v) {
  s[v >> 6] |= (BitWord(1) << (v & 63));
}

inline void ClearBit(BitWord* s, int v) {
  s[v >> 6] &= ~(BitWord(1) << (v & 63));
}

inline bool TestBit(const BitWord* s, int v) {
  return ((s[v >> 6] >> (v & 63)) & BitWord(1)) != 0;
}

inline int PopCount(BitWord x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_popcountll(x);
#else
  int c = 0;
  while (x) {
    x &= x - 1;
    ++c;
  }
  return c;
#endif
}

// Index of the lowest set bit; `x` must be non-zero.
inline int Ctz(BitWord x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_ctzll(x);
#else
  int c = 0;
  while (((x >> c) & BitWord(1)) == 0) {
    ++c;
  }
  return c;
#endif
}

// |a & b| over `nw` words.
inline int AndCount(const BitWord* a, const BitWord* b, int nw) {
  int c = 0;
  for (int w = 0; w < nw; ++w) {
    c += PopCount(a[w] & b[w]);
  }
  return c;
}

// Is (a & mask) contained in b?
inline bool SubsetOf(const BitWord* a, const BitWord* b, const BitWord* mask,
                     int nw) {
  for (int w = 0; w < nw; ++w) {
    if ((a[w] & mask[w] & ~b[w]) != 0) {
      return false;
    }
  }
  return true;
}

inline bool AnyBit(const BitWord* a, int nw) {
  for (int w = 0; w < nw; ++w) {
    if (a[w]) {
      return true;
    }
  }
  return false;
}

// Lowest member of (a & mask), or -1 when that set is empty.
inline int FirstOf(const BitWord* a, const BitWord* mask, int nw) {
  for (int w = 0; w < nw; ++w) {
    const BitWord x = a[w] & mask[w];
    if (x) {
      return (w << 6) + Ctz(x);
    }
  }
  return -1;
}

// The members of (a & mask), ascending.
inline void CollectOf(const BitWord* a, const BitWord* mask, int nw,
                      std::vector<int>* out) {
  out->clear();
  for (int w = 0; w < nw; ++w) {
    BitWord x = a[w] & mask[w];
    while (x) {
      out->push_back((w << 6) + Ctz(x));
      x &= x - 1;
    }
  }
}

// Exhaustive minimum-cover search over one component of the coverage graph.
// Bitmaps are indexed by global point / centre id, so no relabelling is
// needed; the component itself arrives as a point mask and a centre mask.
struct CoverSearch {
  int nw;
  const std::vector<BitWord>* nb;          // n * nw closed neighbourhoods
  std::vector<std::vector<BitWord> > un;   // uncovered points, per depth
  std::vector<std::vector<BitWord> > av;   // available centres, per depth
  std::vector<int> cur;
  std::vector<int> best;
  int bestSize;
  int cap;
  long long nodes;
  bool expired;
  std::chrono::steady_clock::time_point deadline;

  const BitWord* Nb(int v) const {
    return &(*nb)[static_cast<size_t>(v) * nw];
  }

  // Greedy pairwise-disjoint-neighbourhood bound: points chosen so that no
  // available centre covers two of them each need a centre of their own.
  // Returns `cap + 1` if some uncovered point has no available centre at all.
  int LowerBound(int depth) {
    std::vector<BitWord> left(un[depth]);
    const BitWord* avd = av[depth].data();
    int cnt = 0;
    while (true) {
      int j = -1;
      for (int w = 0; w < nw; ++w) {
        if (left[w]) {
          j = (w << 6) + Ctz(left[w]);
          break;
        }
      }
      if (j < 0) {
        return cnt;
      }
      ++cnt;
      const BitWord* nj = Nb(j);
      bool coverable = false;
      for (int w = 0; w < nw; ++w) {
        BitWord x = nj[w] & avd[w];
        while (x) {
          const BitWord* nc = Nb((w << 6) + Ctz(x));
          for (int u = 0; u < nw; ++u) {
            left[u] &= ~nc[u];
          }
          coverable = true;
          x &= x - 1;
        }
      }
      if (!coverable) {
        return cap + 1;
      }
      if (cnt > cap) {
        return cnt;
      }
    }
  }

  bool Expired() {
    if (!expired && ((nodes - 1) & 1023LL) == 0) {   // node 1 included
      if (std::chrono::steady_clock::now() > deadline) {
        expired = true;
      } else {
        R_CheckUserInterrupt();
      }
    }
    return expired;
  }

  void Search(int depth) {
    ++nodes;
    if (Expired()) {
      return;
    }
    if (!AnyBit(un[depth].data(), nw)) {
      bestSize = depth;                    // a complete cover, this small
      best.assign(cur.begin(), cur.end());
      return;
    }
    if (depth + 1 >= bestSize) {
      return;                              // no room to improve
    }
    if (depth + LowerBound(depth) >= bestSize) {
      return;
    }
    // Branch on the uncovered point with the fewest available centres: one of
    // them must be open, and the fewest branches settle first.
    int pick = -1;
    int pickDeg = 0;
    for (int w = 0; w < nw; ++w) {
      BitWord x = un[depth][w];
      while (x) {
        const int j = (w << 6) + Ctz(x);
        const int deg = AndCount(Nb(j), av[depth].data(), nw);
        if (pick < 0 || deg < pickDeg) {
          pick = j;
          pickDeg = deg;
        }
        x &= x - 1;
      }
    }
    if (pickDeg == 0) {
      return;                              // uncoverable: dead branch
    }
    std::vector<int> branch;
    CollectOf(Nb(pick), av[depth].data(), nw, &branch);
    for (size_t t = 0; t < branch.size(); ++t) {
      const int c = branch[t];
      const BitWord* nc = Nb(c);
      for (int w = 0; w < nw; ++w) {
        un[depth + 1][w] = un[depth][w] & ~nc[w];
        av[depth + 1][w] = av[depth][w];
      }
      ClearBit(av[depth + 1].data(), c);
      cur.push_back(c);
      Search(depth + 1);
      cur.pop_back();
      if (expired) {
        return;
      }
      // Covers using `c` have now been enumerated exhaustively, so the
      // remaining branches need not consider it again.
      ClearBit(av[depth].data(), c);
    }
  }

  // The smallest cover of this component, or `cap + 1` if it needs more.
  int Solve(const std::vector<BitWord>& points,
            const std::vector<BitWord>& centres, int capIn) {
    cap = capIn;
    bestSize = capIn + 1;
    best.clear();
    cur.clear();
    un.assign(capIn + 2, std::vector<BitWord>(nw, 0));
    av.assign(capIn + 2, std::vector<BitWord>(nw, 0));
    un[0] = points;
    av[0] = centres;
    Search(0);
    return bestSize;
  }
};

}  // namespace

// Decide one covering probe: can `k` centres cover every point of `d` within
// radius `r`?  Returns list(status, witness, nodes):
//   "feasible"     -- witness is a covering set of size <= k (ascending,
//                     1-based),
//   "infeasible"   -- the search was exhaustive and no such set exists,
//   "inconclusive" -- `maxSeconds` elapsed first.
// [[Rcpp::export]]
List CoverDecide_cpp(NumericMatrix d, double r, int k, double maxSeconds) {
  const int n = d.nrow();
  const int nw = (n + 63) / 64;
  const std::chrono::steady_clock::time_point deadline =
    std::chrono::steady_clock::now() +
    std::chrono::duration_cast<std::chrono::steady_clock::duration>(
      std::chrono::duration<double>(maxSeconds));
  const double* dp = REAL(d);

  // Closed neighbourhoods, one contiguous column read each: bit i of nb[j] is
  // set when centre i covers point j. Symmetric, so nb[j] doubles as the set
  // of points that centre j covers.
  std::vector<BitWord> nb(static_cast<size_t>(n) * nw, 0);
  for (int j = 0; j < n; ++j) {
    const double* col = dp + static_cast<size_t>(j) * n;
    BitWord* row = &nb[static_cast<size_t>(j) * nw];
    for (int i = 0; i < n; ++i) {
      if (col[i] <= r) {
        SetBit(row, i);
      }
    }
  }

  std::vector<BitWord> alive(nw, 0), avail(nw, 0);
  for (int v = 0; v < n; ++v) {
    SetBit(alive.data(), v);
    SetBit(avail.data(), v);
  }
  std::vector<int> opened;
  std::vector<int> buf;

  // ---- reduce to a fixpoint --------------------------------------------
  bool infeasible = false;
  bool timeout = false;
  bool changed = true;
  while (changed && !infeasible && AnyBit(alive.data(), nw)) {
    if (std::chrono::steady_clock::now() > deadline) {
      timeout = true;
      break;
    }
    R_CheckUserInterrupt();
    changed = false;

    // A point with a single available centre forces that centre open.
    for (int j = 0; j < n && !infeasible; ++j) {
      if (!TestBit(alive.data(), j)) {
        continue;
      }
      const BitWord* nj = &nb[static_cast<size_t>(j) * nw];
      // A point never loses every available centre here: a centre is dropped
      // only when another available centre covers a superset, which therefore
      // still covers this point. A point with no cover at all -- which the
      // branching can create -- is refuted by the search's own bound.
      const int deg = AndCount(nj, avail.data(), nw);
      if (deg == 1) {
        const int c = FirstOf(nj, avail.data(), nw);
        opened.push_back(c);
        ClearBit(avail.data(), c);
        const BitWord* nc = &nb[static_cast<size_t>(c) * nw];
        for (int w = 0; w < nw; ++w) {
          alive[w] &= ~nc[w];
        }
        changed = true;
        if (static_cast<int>(opened.size()) > k) {
          infeasible = true;
        }
      }
    }
    if (infeasible || !AnyBit(alive.data(), nw)) {
      break;
    }

    // Point dominance. Any b dominated by a lies in N[c] for the first
    // available c in N[a], so the scan is over one neighbourhood rather than
    // all points.
    for (int a = 0; a < n; ++a) {
      if (!TestBit(alive.data(), a)) {
        continue;
      }
      const BitWord* na = &nb[static_cast<size_t>(a) * nw];
      // An alive point always keeps a cover (see the propagation note above),
      // so `c` is a real centre; the test is structural, not a case.
      const int c = FirstOf(na, avail.data(), nw);
      if (c >= 0) {
        CollectOf(&nb[static_cast<size_t>(c) * nw], alive.data(), nw, &buf);
        for (size_t t = 0; t < buf.size(); ++t) {
          const int b = buf[t];
          if (b == a || !TestBit(alive.data(), b)) {
            continue;
          }
          // Equal neighbourhoods need no tie-break here: `a` ascends, so the
          // lowest member of an equal class is reached while the rest are
          // still alive and drops them, and is itself never reached as `b`
          // afterwards. Some point therefore always survives this pass.
          if (SubsetOf(na, &nb[static_cast<size_t>(b) * nw], avail.data(),
                       nw)) {
            ClearBit(alive.data(), b);
            changed = true;
          }
        }
      }
    }

    // Centre dominance. Any x dominating y covers y's first alive point, so
    // again the scan is over one neighbourhood.
    for (int y = 0; y < n; ++y) {
      if (!TestBit(avail.data(), y)) {
        continue;
      }
      const BitWord* ny = &nb[static_cast<size_t>(y) * nw];
      const int q = FirstOf(ny, alive.data(), nw);
      if (q < 0) {
        ClearBit(avail.data(), y);         // covers nothing still needed
        changed = true;
        continue;
      }
      CollectOf(&nb[static_cast<size_t>(q) * nw], avail.data(), nw, &buf);
      for (size_t t = 0; t < buf.size(); ++t) {
        const int x = buf[t];
        if (x == y || !TestBit(avail.data(), x)) {
          continue;
        }
        const BitWord* nx = &nb[static_cast<size_t>(x) * nw];
        if (!SubsetOf(ny, nx, alive.data(), nw)) {
          continue;
        }
        if (y < x && SubsetOf(nx, ny, alive.data(), nw)) {
          continue;
        }
        ClearBit(avail.data(), y);
        changed = true;
        break;
      }
    }
  }

  double nodes = 0;
  if (timeout) {
    return List::create(_["status"] = "inconclusive",
                        _["witness"] = IntegerVector(0),
                        _["nodes"] = nodes);
  }
  if (infeasible) {
    return List::create(_["status"] = "infeasible",
                        _["witness"] = IntegerVector(0),
                        _["nodes"] = nodes);
  }

  // ---- split into components, solve each against the budget left -------
  // A centre covering two points puts both in one component, so components
  // partition the centres too and their minimum covers add.
  std::vector<int> comp(n, -1);
  std::vector<std::vector<int> > comps;
  std::vector<int> stack;
  for (int s = 0; s < n; ++s) {
    if (!TestBit(alive.data(), s) || comp[s] >= 0) {
      continue;
    }
    const int id = static_cast<int>(comps.size());
    comps.push_back(std::vector<int>());
    comp[s] = id;
    stack.push_back(s);
    while (!stack.empty()) {
      const int j = stack.back();
      stack.pop_back();
      comps[id].push_back(j);
      CollectOf(&nb[static_cast<size_t>(j) * nw], avail.data(), nw, &buf);
      const std::vector<int> centres(buf);
      for (size_t t = 0; t < centres.size(); ++t) {
        std::vector<int> reached;
        CollectOf(&nb[static_cast<size_t>(centres[t]) * nw], alive.data(), nw,
                  &reached);
        for (size_t u = 0; u < reached.size(); ++u) {
          if (comp[reached[u]] < 0) {
            comp[reached[u]] = id;
            stack.push_back(reached[u]);
          }
        }
      }
    }
  }

  CoverSearch search;
  search.nw = nw;
  search.nb = &nb;
  search.deadline = deadline;
  search.expired = false;
  search.nodes = 0;

  int budget = k - static_cast<int>(opened.size());
  std::vector<int> witness(opened);
  bool inconclusive = false;
  for (size_t ci = 0; ci < comps.size() && !infeasible && !inconclusive; ++ci) {
    std::vector<BitWord> pts(nw, 0), ctr(nw, 0);
    for (size_t t = 0; t < comps[ci].size(); ++t) {
      SetBit(pts.data(), comps[ci][t]);
      const BitWord* nj = &nb[static_cast<size_t>(comps[ci][t]) * nw];
      for (int w = 0; w < nw; ++w) {
        ctr[w] |= nj[w] & avail[w];
      }
    }
    const int got = search.Solve(pts, ctr, budget);
    nodes += static_cast<double>(search.nodes);
    search.nodes = 0;
    if (search.expired) {
      inconclusive = true;
    } else if (got > budget) {
      infeasible = true;
    } else {
      budget -= got;
      witness.insert(witness.end(), search.best.begin(), search.best.end());
    }
  }

  if (inconclusive) {
    return List::create(_["status"] = "inconclusive",
                        _["witness"] = IntegerVector(0),
                        _["nodes"] = nodes);
  }
  if (infeasible) {
    return List::create(_["status"] = "infeasible",
                        _["witness"] = IntegerVector(0),
                        _["nodes"] = nodes);
  }
  std::sort(witness.begin(), witness.end());
  IntegerVector out(witness.size());
  for (size_t t = 0; t < witness.size(); ++t) {
    out[t] = witness[t] + 1;
  }
  return List::create(_["status"] = "feasible",
                      _["witness"] = out,
                      _["nodes"] = nodes);
}

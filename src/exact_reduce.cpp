// exact_reduce.cpp
//
// Combinatorial reduction for one ExactMaxMin feasibility probe.
//
// A probe asks whether the threshold graph G(lambda) (edges: pairs closer
// than lambda) has an independent set of size >= k -- equivalently, whether
// the complement graph H (pairs >= lambda apart) contains a k-clique. Every
// vertex of a k-clique has H-degree >= k - 1, so iteratively deleting
// vertices of H-degree < k - 1 (the (k-1)-core peel) never removes a witness
// vertex, and a k-clique lies within a single H-component. A proper
// colouring of H bounds its clique number (omega <= chi), so a component
// greedily coloured with fewer than k colours cannot host a witness.
//
// Same-colour vertices are pairwise H-non-adjacent, hence pairwise
// G-adjacent: each colour class is a G-clique whose single packing row
// sum(x) <= 1 carries all of its pairwise constraints, leaving only
// cross-colour G-edges to emit as pairwise rows. The reduced per-component
// model has the same integer feasible set as the all-pairs formulation, with
// an LP bound of ~chi rather than ~n/2.
//
// All tie-breaks are by vertex index, so the output is deterministic.

#include <Rcpp.h>
using namespace Rcpp;

// Reduce the complement graph H (edge list `hi`/`hj`, 1-based) of one
// threshold probe on `n` vertices against target clique size `k`.
// Returns list(clique, comps):
//   clique -- a greedy H-clique of size >= k (sorted, 1-based) when one is
//             found: a feasibility witness needing no IP; else integer(0).
//   comps  -- one entry per surviving component (size >= k and greedy
//             chi >= k), each list(vars, cls, ei, ej): member vertices
//             (ascending, 1-based), their colour class, and the cross-colour
//             G-edges as local indices into `vars`. Empty when `clique` is
//             non-empty or nothing survives (the latter proves
//             infeasibility).
// [[Rcpp::export]]
List ThresholdReduce_cpp(IntegerVector hi, IntegerVector hj,
                         int n, int k) {
  const R_xlen_t nE = hi.size();
  const int need = k - 1;

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

  // Peel to the (k-1)-core. After the loop, dg[] of a surviving vertex
  // counts its surviving neighbours.
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

  // Greedy colouring, highest surviving degree first. Neighbours are always
  // within the vertex's own component, so per-component chi is the largest
  // colour its members receive.
  std::vector<int> alive;
  alive.reserve(n);
  for (int v = 0; v < n; ++v) {
    if (!dead[v]) {
      alive.push_back(v);
    }
  }
  std::sort(alive.begin(), alive.end(), [&](int a, int b) {
    return dg[a] != dg[b] ? dg[a] > dg[b] : a < b;
  });
  std::vector<int> colour(n, 0);
  std::vector<int> compChi(nComp, 0);
  std::vector<int> mark(n + 2, -1);
  int stamp = 0;
  for (const int v : alive) {
    ++stamp;
    for (R_xlen_t p = off[v]; p < off[v + 1]; ++p) {
      const int u = adj[p];
      if (!dead[u] && colour[u]) {
        mark[colour[u]] = stamp;
      }
    }
    int c = 1;
    while (mark[c] == stamp) {
      ++c;
    }
    colour[v] = c;
    if (c > compChi[comp[v] - 1]) {
      compChi[comp[v] - 1] = c;
    }
  }

  // Greedy clique attempt per surviving component: seed at the highest
  // surviving degree, grow by highest degree among common neighbours.
  // Reaching k proves feasibility with no IP.
  std::vector<int> vmark(n, -1);
  int vstamp = 0;
  IntegerVector clique;
  for (int c = 1; c <= nComp && clique.size() == 0; ++c) {
    if (compSize[c - 1] < k || compChi[c - 1] < k) {
      continue;
    }
    int s = -1;
    for (int v = 0; v < n; ++v) {
      if (!dead[v] && comp[v] == c && (s < 0 || dg[v] > dg[s])) {
        s = v;
      }
    }
    std::vector<int> cand;
    for (R_xlen_t p = off[s]; p < off[s + 1]; ++p) {
      if (!dead[adj[p]]) {
        cand.push_back(adj[p]);
      }
    }
    std::vector<int> cl;
    cl.push_back(s);
    while (static_cast<int>(cl.size()) < k && !cand.empty()) {
      int best = cand[0];
      for (const int u : cand) {
        if (dg[u] > dg[best] || (dg[u] == dg[best] && u < best)) {
          best = u;
        }
      }
      cl.push_back(best);
      ++vstamp;
      for (R_xlen_t p = off[best]; p < off[best + 1]; ++p) {
        vmark[adj[p]] = vstamp;
      }
      std::vector<int> next;
      for (const int u : cand) {
        if (u != best && vmark[u] == vstamp) {
          next.push_back(u);
        }
      }
      cand.swap(next);
    }
    if (static_cast<int>(cl.size()) >= k) {
      std::sort(cl.begin(), cl.end());
      IntegerVector out(cl.size());
      for (R_xlen_t t = 0; t < static_cast<R_xlen_t>(cl.size()); ++t) {
        out[t] = cl[t] + 1;
      }
      clique = out;
    }
  }

  // Emit the reduced model for each surviving component.
  List comps;
  if (clique.size() == 0) {
    for (int c = 1; c <= nComp; ++c) {
      if (compSize[c - 1] < k || compChi[c - 1] < k) {
        continue;
      }
      std::vector<int> vars;
      for (int v = 0; v < n; ++v) {
        if (!dead[v] && comp[v] == c) {
          vars.push_back(v);
        }
      }
      const int nv = static_cast<int>(vars.size());
      std::vector<int> loc(n, -1);
      for (int t = 0; t < nv; ++t) {
        loc[vars[t]] = t;
      }
      std::vector<int> ei;
      std::vector<int> ej;
      for (int a = 0; a < nv; ++a) {
        const int u = vars[a];
        ++vstamp;
        for (R_xlen_t p = off[u]; p < off[u + 1]; ++p) {
          vmark[adj[p]] = vstamp;
        }
        for (int b = a + 1; b < nv; ++b) {
          const int v = vars[b];
          if (vmark[v] == vstamp || colour[u] == colour[v]) {
            continue;
          }
          ei.push_back(a + 1);
          ej.push_back(b + 1);
        }
      }
      IntegerVector vv(nv);
      IntegerVector cc(nv);
      for (int t = 0; t < nv; ++t) {
        vv[t] = vars[t] + 1;
        cc[t] = colour[vars[t]];
      }
      comps.push_back(List::create(_["vars"] = vv, _["cls"] = cc,
                                   _["ei"] = IntegerVector(ei.begin(), ei.end()),
                                   _["ej"] = IntegerVector(ej.begin(), ej.end())));
    }
  }
  return List::create(_["clique"] = clique, _["comps"] = comps);
}

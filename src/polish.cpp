#include <Rcpp.h>
#include <algorithm>
#include <vector>

// Local-search polish for max-min diversity (MMDP) selections.
//
// Given a current selection S of size k, search for a 1-swap (x in S, w not in
// S) that increases T_k = min over (a, b) in S^2, a != b, of d(a, b).
//
// Neighbourhood: critical-edge-anchored. We restrict the candidate removals to
// endpoints of pairs achieving the current minimum, and candidate insertions
// to the i = 1..limit nearest neighbours (in the full data, excluding S) of
// each such endpoint. Cf. Resende, Marti, Gallego & Duarte (2010).
//
// Acceptance: Della Croce tie-break. Accept a swap when new_T > T OR
// (new_T == T AND new_n_critical < n_critical). Both T and -n_critical are
// monotone non-decreasing across accepted swaps, so iteration terminates.
//
// Arguments mirror PolishSelection() in R/samplers.R:
//   d           square numeric distance matrix (N x N)
//   S           integer vector of 1-based selected indices (length k)
//   limit       maximum neighbour rank to scan per endpoint (i = 1..limit)
//   max_passes  cap on outer iterations (safety; convergence is theoretical)
//
// Returns IntegerVector of length k (1-based, same as input) with attributes:
//   "passes" : number of outer iterations that produced an accepting swap
//   "swaps"  : same as passes (alias for clarity at the R boundary)

// [[Rcpp::export]]
Rcpp::IntegerVector PolishMaximin_cpp(Rcpp::NumericMatrix d,
                                      Rcpp::IntegerVector S,
                                      int limit,
                                      int max_passes) {
  int N = d.nrow();
  int k = S.size();

  Rcpp::IntegerVector out = Rcpp::clone(S);
  if (k < 2 || N <= k || limit < 1 || max_passes < 1) {
    out.attr("passes") = 0;
    out.attr("swaps")  = 0;
    return out;
  }

  // 0-based working copy of S; in_S[i] tracks membership.
  std::vector<int> Sv(k);
  std::vector<unsigned char> in_S(N, 0);
  for (int s = 0; s < k; ++s) {
    int idx0 = S[s] - 1;
    if (idx0 < 0 || idx0 >= N) {
      Rcpp::stop("S contains an out-of-range index");
    }
    if (in_S[idx0]) {
      Rcpp::stop("S contains a duplicate index");
    }
    Sv[s]       = idx0;
    in_S[idx0]  = 1;
  }

  // Lambda: compute T and n_critical for the current Sv (k slots).
  // Counts unordered critical pairs (a < b) at distance == T.
  auto compute_T = [&](double &T, int &n_critical) {
    T = R_PosInf;
    for (int a = 0; a < k; ++a) {
      for (int b = a + 1; b < k; ++b) {
        double v = d(Sv[a], Sv[b]);
        if (v < T) T = v;
      }
    }
    n_critical = 0;
    for (int a = 0; a < k; ++a) {
      for (int b = a + 1; b < k; ++b) {
        if (d(Sv[a], Sv[b]) == T) ++n_critical;
      }
    }
  };

  // Lambda: collect distinct slot indices appearing in any critical pair.
  auto critical_slots = [&](double T, std::vector<int> &out_slots) {
    out_slots.clear();
    std::vector<unsigned char> seen(k, 0);
    for (int a = 0; a < k; ++a) {
      for (int b = a + 1; b < k; ++b) {
        if (d(Sv[a], Sv[b]) == T) {
          if (!seen[a]) { seen[a] = 1; out_slots.push_back(a); }
          if (!seen[b]) { seen[b] = 1; out_slots.push_back(b); }
        }
      }
    }
  };

  // Lambda: evaluate the swap (Sv[x] := w) without committing. Returns the
  // pair (new_T, new_n_critical). O(k^2).
  auto evaluate_swap = [&](int x, int w, double &new_T, int &new_n_critical) {
    new_T = R_PosInf;
    // Pairs entirely within S \ {Sv[x]}.
    for (int a = 0; a < k; ++a) {
      if (a == x) continue;
      for (int b = a + 1; b < k; ++b) {
        if (b == x) continue;
        double v = d(Sv[a], Sv[b]);
        if (v < new_T) new_T = v;
      }
    }
    // Pairs (w, s) for s in S \ {Sv[x]}.
    for (int a = 0; a < k; ++a) {
      if (a == x) continue;
      double v = d(Sv[a], w);
      if (v < new_T) new_T = v;
    }
    // Count pairs at new_T.
    new_n_critical = 0;
    for (int a = 0; a < k; ++a) {
      if (a == x) continue;
      for (int b = a + 1; b < k; ++b) {
        if (b == x) continue;
        if (d(Sv[a], Sv[b]) == new_T) ++new_n_critical;
      }
    }
    for (int a = 0; a < k; ++a) {
      if (a == x) continue;
      if (d(Sv[a], w) == new_T) ++new_n_critical;
    }
  };

  // Buffer for top-`limit` neighbours of a query slot.
  std::vector<std::pair<double, int> > nbrs;
  nbrs.reserve(N);

  double T;
  int n_critical;
  compute_T(T, n_critical);

  int swaps = 0;
  for (int pass = 0; pass < max_passes; ++pass) {
    std::vector<int> crit;
    critical_slots(T, crit);
    if (crit.empty()) break;  // defensive: k >= 2 guarantees crit non-empty

    bool improved = false;

    for (size_t ci = 0; ci < crit.size() && !improved; ++ci) {
      int x = crit[ci];
      int x_idx = Sv[x];

      // Build neighbour list of x_idx, ascending by distance, excluding S.
      nbrs.clear();
      for (int j = 0; j < N; ++j) {
        if (j == x_idx) continue;   // self has d = 0; would also be in S anyway
        if (in_S[j])     continue;
        nbrs.push_back(std::make_pair(d(x_idx, j), j));
      }
      int take = std::min<int>(limit, static_cast<int>(nbrs.size()));
      if (take == 0) continue;

      // Partial sort: cheapest top-`take` smallest distances.
      std::partial_sort(
        nbrs.begin(), nbrs.begin() + take, nbrs.end(),
        [](const std::pair<double,int> &a, const std::pair<double,int> &b) {
          return a.first < b.first;
        });

      for (int i = 0; i < take; ++i) {
        int w = nbrs[i].second;

        double new_T;
        int    new_n_critical;
        evaluate_swap(x, w, new_T, new_n_critical);

        bool accept = (new_T > T) ||
                      (new_T == T && new_n_critical < n_critical);
        if (accept) {
          in_S[Sv[x]] = 0;
          Sv[x]       = w;
          in_S[w]     = 1;
          T           = new_T;
          n_critical  = new_n_critical;
          ++swaps;
          improved = true;
          break;
        }
      }
    }

    if (!improved) break;
  }

  for (int s = 0; s < k; ++s) out[s] = Sv[s] + 1;
  out.attr("passes") = swaps;
  out.attr("swaps")  = swaps;
  return out;
}

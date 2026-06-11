#include <Rcpp.h>
#include <vector>
#include <algorithm>

// CDSh heuristic for the discrete (vertex) k-centre problem.
//
// k-centre: choose k centres from the N points so the covering radius
//   R(C) = max_p min_{c in C} d(p, c)
// (the largest distance from any point to its nearest centre) is minimised.
// This is the min-MAX covering objective -- distinct from MaxMin/MMDP, which
// maximises the min PAIRWISE distance within the selection (MinDist / T_k).
//
// Gonzalez farthest-first (FarFirst) is the classical 2-approximation here, but
// empirically lands ~50% above optimum. CDSh (Garcia-Diaz, Menchaca-Mendez et
// al. 2017 J.Heuristics; 2019 IEEE Access) reaches ~1-3.5% of optimum at
// O(n^2 log n). We port the authors' own reference implementation
// (github.com/jesgadiaz/k-center-in-C, main_small.c), which differs from the
// paper's "dominating set on the squared threshold graph" framing: the actual
// construction is a fixed-k farthest-point loop in which each centre is the
// MAX-DEGREE NEIGHBOUR (within radius r) of the worst-covered vertex, and a
// trial radius r is judged feasible by the achieved covering radius max_dist
// <= r. The squared-graph 2-hop idea is realised by "stepping back" from the
// farthest vertex to its best-connected neighbour, never by squaring the graph.
//
// The one source of randomness in the reference -- the i = 0 seed vertex -- is
// here a deterministic parameter (`seed`, 1-based), supplied by the R wrapper
// (a peripheral anchor). All other ties resolve to the first (lowest) index.
//
// Memory access: the distance matrix is column-major and SYMMETRIC, so every
// d(a, b) is read through the pointer to column `a` (`P + a*n`), which makes the
// inner loops stride-1. d(a, b) == d(b, a) makes this bit-identical to direct
// d(a, b) indexing while turning the cache-hostile row scans (score init,
// domination decrement) into sequential ones (cf. the T-005 reorder).

// One CDS construction at a fixed radius r from a fixed seed vertex (0-based).
// `P` is the column-major distance matrix, `n` its dimension. Fills `centres`
// (k of them, selection order) and returns the achieved covering radius.
static double CdsConstruct(const double* P, int n, int k,
                           double r, int seed0,
                           std::vector<int>& centres) {
  std::vector<double> dist(n, R_PosInf);   // distance to nearest chosen centre
  std::vector<char>   dominated(n, 0);      // covered within r by some centre?
  std::vector<int>    score(n, 0);          // undominated degree in G_r

  // score[i] = #{ j != i : d(i, j) <= r } -- the threshold-graph degree. Read
  // column i (= d(., i) = d(i, .) by symmetry), so j scans contiguously.
  for (int i = 0; i < n; i++) {
    const double* col_i = P + (R_xlen_t)i * n;
    int s = 0;
    for (int j = 0; j < n; j++) {
      if (j != i && col_i[j] <= r) s++;
    }
    score[i] = s;
  }

  centres.assign(k, 0);
  for (int i = 0; i < k; i++) {
    int farthest;
    if (i == 0) {
      farthest = seed0;                     // deterministic seed (was rand())
    } else {
      // Fold the previous centre into the nearest-centre distances, then take
      // the worst-covered vertex as the next critical vertex (Gonzalez step).
      const double* col_prev = P + (R_xlen_t)centres[i - 1] * n;
      for (int j = 0; j < n; j++) {
        if (col_prev[j] < dist[j]) dist[j] = col_prev[j];
      }
      double md = R_NegInf;
      farthest = 0;
      for (int j = 0; j < n; j++) {
        if (dist[j] > md) { md = dist[j]; farthest = j; }
      }
    }

    // The centre is the highest-degree neighbour of the farthest vertex within
    // radius r (its closed neighbourhood includes itself, d = 0 <= r). First
    // index wins ties.
    const double* col_f = P + (R_xlen_t)farthest * n;
    int maxScore = -1, centre = farthest;
    for (int j = 0; j < n; j++) {
      if (col_f[j] <= r && score[j] > maxScore) {
        maxScore = score[j];
        centre = j;
      }
    }
    centres[i] = centre;

    // Newly dominated vertices (within r of the new centre) stop contributing
    // to anyone's score, so score[j] keeps counting only undominated neighbours.
    const double* col_c = P + (R_xlen_t)centre * n;
    for (int j = 0; j < n; j++) {
      if (!dominated[j] && col_c[j] <= r) {
        dominated[j] = 1;
        const double* col_j = P + (R_xlen_t)j * n;   // d(j, .) = d(., j)
        for (int b = 0; b < n; b++) {
          if (b != j && col_j[b] <= r) score[b]--;
        }
      }
    }
  }

  // Fold in the last centre and read off the achieved covering radius.
  const double* col_last = P + (R_xlen_t)centres[k - 1] * n;
  for (int j = 0; j < n; j++) {
    if (col_last[j] < dist[j]) dist[j] = col_last[j];
  }
  double maxd = 0.0;
  for (int j = 0; j < n; j++) {
    if (dist[j] > maxd) maxd = dist[j];
  }
  return maxd;
}

// Sorted distinct off-diagonal distances (the candidate radii both k-centre
// solvers search). Replaces the R idiom `sort(unique(d[upper.tri(d)]))`, which
// allocates an n x n logical mask and pushes ~n^2/2 values through R's
// `unique`/`order`. Extracts the upper triangle column-sequentially, then
// std::sort + std::unique -- the same sorted distinct set, bit-identical.
// [[Rcpp::export]]
Rcpp::NumericVector KCentreCandidates_cpp(Rcpp::NumericMatrix d) {
  int n = d.nrow();
  const double* P = d.begin();
  std::vector<double> vals;
  vals.reserve((std::size_t)n * (n - 1) / 2);
  for (int j = 1; j < n; j++) {
    const double* col_j = P + (R_xlen_t)j * n;       // column j; rows i < j
    for (int i = 0; i < j; i++) vals.push_back(col_j[i]);
  }
  std::sort(vals.begin(), vals.end());
  vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
  return Rcpp::NumericVector(vals.begin(), vals.end());
}

// CDSh: binary search over the sorted distinct candidate radii `cand`
// (ascending), running one CDS construction per trial. A trial r is feasible
// iff the construction's achieved covering radius is <= r (push the search to
// smaller radii); otherwise push larger. The best (smallest) achieved covering
// radius over all trials, and its centre set, are returned -- the achieved
// radius is a realised point-to-centre distance and may be strictly below the
// trial r, so we always keep the best seen rather than the final bracket.
//
//   d     N x N distance matrix (assumed symmetric, non-negative, finite)
//   k     number of centres (1 <= k < N; the wrapper handles k >= N)
//   seed  1-based seed vertex for the i = 0 critical vertex
//   cand  ascending vector of distinct candidate radii
//
// Returns the k centre indices (1-based, ascending) with a `radius` attribute.
// [[Rcpp::export]]
Rcpp::IntegerVector KCentreCDSh_cpp(Rcpp::NumericMatrix d, int k, int seed,
                                    Rcpp::NumericVector cand) {
  int n = d.nrow();
  int nCand = cand.size();
  if (k < 1 || k >= n) {
    Rcpp::stop("KCentreCDSh_cpp: expect 1 <= k < n; got k = %d, n = %d", k, n);
  }
  if (seed < 1 || seed > n) {
    Rcpp::stop("KCentreCDSh_cpp: 'seed' must be in [1, %d]; got %d", n, seed);
  }
  const double* P = d.begin();
  int seed0 = seed - 1;

  double bestRadius = R_PosInf;
  std::vector<int> bestCentres;
  std::vector<int> work;

  // Evaluate one candidate index, updating the incumbent. Returns the achieved
  // covering radius so the search can compare it against the trial radius.
  auto evalAt = [&](int idx) -> double {
    double r = cand[idx];
    double achieved = CdsConstruct(P, n, k, r, seed0, work);
    if (achieved < bestRadius) {
      bestRadius = achieved;
      bestCentres = work;
    }
    return achieved;
  };

  // Feasible region is a suffix of `cand` (larger radius -> easier to cover);
  // gallop the bracket via binary search. low = below the smallest candidate
  // (treated infeasible), high = largest candidate (k >= 1 covers all within
  // the diameter, so always feasible). Probe inclusively so the smallest
  // feasible candidate is itself evaluated.
  int low = -1, high = nCand - 1;
  evalAt(high);
  while (high - low > 1) {
    int mid = (high + low) / 2;
    double achieved = evalAt(mid);
    if (achieved <= cand[mid]) high = mid; else low = mid;
  }

  std::sort(bestCentres.begin(), bestCentres.end());
  Rcpp::IntegerVector out(bestCentres.size());
  for (size_t i = 0; i < bestCentres.size(); i++) out[i] = bestCentres[i] + 1;
  out.attr("radius") = bestRadius;
  return out;
}

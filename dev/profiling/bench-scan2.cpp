// Micro-benchmark #2 for the MaxMean inner loop (focus area #6, T-012).
//   (A) scan: V0 current vs V2 monotonicity (track max-p add / min-p remove,
//       2 divisions total instead of n; no per-element sub+div).
//   (B) P-update: branched `if (j!=u)` vs branchless (needs diag==0).
//
// Build: g++ -O2 -mfpmath=sse -msse2 -std=c++17 dev/profiling/bench-scan2.cpp -o dev/profiling/bench-scan2

#include <vector>
#include <cstdio>
#include <chrono>
#include <random>
#include <algorithm>

using clk = std::chrono::steady_clock;

static int scan_v0(const std::vector<double>& p, const std::vector<int>& in_S,
                   const std::vector<long long>& tabu_until,
                   double f, double best_f_local, long long iter, int m, int n,
                   double& out_delta) {
  int best_flip = -1; double best_delta = -1e300;
  for (int i = 0; i < n; ++i) {
    double delta;
    if (!in_S[i]) delta = (p[i] - f) / (double)(m + 1);
    else { if (m <= 2) continue; delta = (f - p[i]) / (double)(m - 1); }
    bool is_tabu = tabu_until[i] > iter;
    if (is_tabu && f + delta <= best_f_local) continue;
    if (delta > best_delta) { best_delta = delta; best_flip = i; }
  }
  out_delta = best_delta; return best_flip;
}

// V2: monotonicity. delta_add increases in p, delta_rem decreases in p.
// Track max-p among admissible adds, min-p among admissible removes; compute
// only 2 deltas at the end. Aspiration becomes a precomputed p-threshold.
static int scan_v2(const std::vector<double>& p, const std::vector<int>& in_S,
                   const std::vector<long long>& tabu_until,
                   double f, double best_f_local, long long iter, int m, int n,
                   double& out_delta) {
  const double gap = best_f_local - f;
  const double addThresh = f + gap * (double)(m + 1);   // tabu add adm. iff p > addThresh
  const double remThresh = f - gap * (double)(m - 1);   // tabu rem adm. iff p < remThresh
  double maxAddP = -1e300; int addIdx = -1;
  double minRemP =  1e300; int remIdx = -1;
  for (int i = 0; i < n; ++i) {
    double pi = p[i];
    bool is_tabu = tabu_until[i] > iter;
    if (!in_S[i]) {
      if ((!is_tabu || pi > addThresh) && pi > maxAddP) { maxAddP = pi; addIdx = i; }
    } else {
      if (m <= 2) continue;
      if ((!is_tabu || pi < remThresh) && pi < minRemP) { minRemP = pi; remIdx = i; }
    }
  }
  int best_flip = -1; double best_delta = -1e300;
  if (addIdx >= 0) { double d = (maxAddP - f) / (double)(m + 1); if (d > best_delta) { best_delta = d; best_flip = addIdx; } }
  if (remIdx >= 0) { double d = (f - minRemP) / (double)(m - 1); if (d > best_delta) { best_delta = d; best_flip = remIdx; } }
  out_delta = best_delta; return best_flip;
}

static void pupd_branched(std::vector<double>& p, const double* col, int u, int n) {
  for (int j = 0; j < n; ++j) if (j != u) p[j] += col[j];
}
static void pupd_branchless(std::vector<double>& p, const double* col, int n) {
  for (int j = 0; j < n; ++j) p[j] += col[j];   // requires col[u]==0 (diag)
}

template <class F>
double timed(F&& fn, long long reps) {
  std::vector<double> meds;
  for (int trial = 0; trial < 5; ++trial) {
    auto t0 = clk::now();
    for (long long r = 0; r < reps; ++r) fn();
    auto t1 = clk::now();
    meds.push_back(std::chrono::duration<double, std::nano>(t1 - t0).count() / (double)reps);
  }
  std::sort(meds.begin(), meds.end());
  return meds[meds.size()/2];
}

int main() {
  const int n = 500, m = 136;
  std::mt19937 rng(5813);
  std::uniform_real_distribution<double> ud(-2000.0, 2000.0);
  std::vector<double> p(n), col(n);
  std::vector<int> in_S(n, 0);
  std::vector<long long> tabu_until(n, 0);
  for (int i = 0; i < n; ++i) { p[i] = ud(rng); col[i] = ud(rng); }
  { std::vector<int> idx(n); for (int i=0;i<n;++i) idx[i]=i;
    std::shuffle(idx.begin(), idx.end(), rng);
    for (int j=0;j<m;++j) in_S[idx[j]] = 1; }
  long long iter = 100000;
  { std::uniform_int_distribution<int> pick(0, n-1);
    for (int t=0;t<60;++t) tabu_until[pick(rng)] = iter + (5 + t % 116); }
  double f = 40.0, best_f_local = 54.5;
  int u = 200; col[u] = 0.0;   // diagonal zero for branchless

  double d0, d2;
  int f0 = scan_v0(p, in_S, tabu_until, f, best_f_local, iter, m, n, d0);
  int f2 = scan_v2(p, in_S, tabu_until, f, best_f_local, iter, m, n, d2);

  volatile double sink = 0;
  const long long reps = 2000000;
  double t0 = timed([&]{ double d; sink += scan_v0(p,in_S,tabu_until,f,best_f_local,iter,m,n,d); }, reps);
  double t2 = timed([&]{ double d; sink += scan_v2(p,in_S,tabu_until,f,best_f_local,iter,m,n,d); }, reps);

  // P-update: copy p each rep is too costly; instead alternate +=/-= to keep bounded.
  std::vector<double> pa = p, pb = p;
  double tb  = timed([&]{ pupd_branched(pa, col.data(), u, n);  pupd_branched(pa, col.data(), u, n);  for(auto&x:pa)x*=0.0; }, 200000);
  double tbl = timed([&]{ pupd_branchless(pb, col.data(), n);   pupd_branchless(pb, col.data(), n);   for(auto&x:pb)x*=0.0; }, 200000);

  printf("scan best_flip: v0=%d v2=%d (same=%s)  delta v0=%.10g v2=%.10g\n",
         f0, f2, f0==f2?"yes":"NO", d0, d2);
  printf("scan ns/call:   v0=%.1f  v2=%.1f   speedup=%.3fx\n", t0, t2, t0/t2);
  printf("pupd ns/2call+clr: branched=%.1f branchless=%.1f speedup=%.3fx (incl clr both)\n", tb, tbl, tb/tbl);
  printf("[sink %.3g]\n", (double)sink);
  return 0;
}

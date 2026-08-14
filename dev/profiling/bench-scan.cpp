// Micro-benchmark for the MaxMean tabu best-flip scan (focus area #6, T-012).
// Verifies the strength-reduction optimisation: replace the per-element
// division by the loop-invariant (m +/- 1) with a multiply by a precomputed
// reciprocal (n divisions -> n multiplies + 2 divisions).
//
// Build (match the package flags: -O2 -mfpmath=sse -msse2):
//   g++ -O2 -mfpmath=sse -msse2 -std=c++17 dev/profiling/bench-scan.cpp -o bench-scan
// Run:
//   ./bench-scan

#include <vector>
#include <cstdio>
#include <chrono>
#include <random>
#include <algorithm>

using clk = std::chrono::steady_clock;

// V0: current kernel — per-element division.
static int scan_v0(const std::vector<double>& p, const std::vector<int>& in_S,
                   const std::vector<long long>& tabu_until,
                   double f, double best_f_local, long long iter, int m, int n,
                   double& out_delta) {
  int best_flip = -1;
  double best_delta = -1e300;
  for (int i = 0; i < n; ++i) {
    double delta;
    if (!in_S[i]) {
      delta = (p[i] - f) / (double)(m + 1);
    } else {
      if (m <= 2) continue;
      delta = (f - p[i]) / (double)(m - 1);
    }
    bool is_tabu = tabu_until[i] > iter;
    if (is_tabu && f + delta <= best_f_local) continue;
    if (delta > best_delta) { best_delta = delta; best_flip = i; }
  }
  out_delta = best_delta;
  return best_flip;
}

// V1: strength reduction — precompute reciprocals, multiply per element.
static int scan_v1(const std::vector<double>& p, const std::vector<int>& in_S,
                   const std::vector<long long>& tabu_until,
                   double f, double best_f_local, long long iter, int m, int n,
                   double& out_delta) {
  int best_flip = -1;
  double best_delta = -1e300;
  const double inv_add = 1.0 / (double)(m + 1);
  const double inv_rem = (m > 2) ? 1.0 / (double)(m - 1) : 0.0;
  for (int i = 0; i < n; ++i) {
    double delta;
    if (!in_S[i]) {
      delta = (p[i] - f) * inv_add;
    } else {
      if (m <= 2) continue;
      delta = (f - p[i]) * inv_rem;
    }
    bool is_tabu = tabu_until[i] > iter;
    if (is_tabu && f + delta <= best_f_local) continue;
    if (delta > best_delta) { best_delta = delta; best_flip = i; }
  }
  out_delta = best_delta;
  return best_flip;
}

template <class F>
double time_scan(F&& fn, long long reps, int& sink_flip, double& sink_delta) {
  // 3 medians of `reps` scans each.
  std::vector<double> meds;
  for (int trial = 0; trial < 5; ++trial) {
    auto t0 = clk::now();
    int bf = 0; double bd = 0;
    for (long long r = 0; r < reps; ++r) {
      double d;
      bf ^= fn(d);          // xor to defeat dead-code elimination
      bd += d;
    }
    auto t1 = clk::now();
    sink_flip ^= bf; sink_delta += bd;
    meds.push_back(std::chrono::duration<double, std::nano>(t1 - t0).count() / (double)reps);
  }
  std::sort(meds.begin(), meds.end());
  return meds[meds.size() / 2];   // median ns/scan
}

int main() {
  const int n = 500;
  const int m = 136;             // |S| occupancy from the n=500 driver
  std::mt19937 rng(5813);
  std::uniform_real_distribution<double> ud(-2000.0, 2000.0);

  std::vector<double> p(n);
  std::vector<int> in_S(n, 0);
  std::vector<long long> tabu_until(n, 0);
  for (int i = 0; i < n; ++i) p[i] = ud(rng);
  // mark m members
  { std::vector<int> idx(n); for (int i=0;i<n;++i) idx[i]=i;
    std::shuffle(idx.begin(), idx.end(), rng);
    for (int j=0;j<m;++j) in_S[idx[j]] = 1; }
  // a realistic tabu sprinkle: ~tenure recent flips still tabu
  long long iter = 100000;
  { std::uniform_int_distribution<int> pick(0, n-1);
    for (int t=0;t<60;++t) tabu_until[pick(rng)] = iter + (5 + t % 116); }
  double f = 40.0, best_f_local = 54.5;

  int sink_flip = 0; double sink_delta = 0;
  double d0, d1; int f0, f1;
  f0 = scan_v0(p, in_S, tabu_until, f, best_f_local, iter, m, n, d0);
  f1 = scan_v1(p, in_S, tabu_until, f, best_f_local, iter, m, n, d1);

  const long long reps = 2000000;
  double t_v0 = time_scan([&](double& d){ return scan_v0(p,in_S,tabu_until,f,best_f_local,iter,m,n,d); }, reps, sink_flip, sink_delta);
  double t_v1 = time_scan([&](double& d){ return scan_v1(p,in_S,tabu_until,f,best_f_local,iter,m,n,d); }, reps, sink_flip, sink_delta);

  printf("best_flip  v0=%d  v1=%d  (same=%s)\n", f0, f1, f0==f1 ? "yes":"NO");
  printf("best_delta v0=%.10g  v1=%.10g\n", d0, d1);
  printf("scan ns/call:  v0=%.1f  v1=%.1f   speedup=%.3fx\n", t_v0, t_v1, t_v0/t_v1);
  printf("[sink %d %.3g]\n", sink_flip, sink_delta);
  return 0;
}

// maxmean.cpp
//
// RLTS algorithm for the Max-Mean Dispersion Problem.
// Dieudonne, Zhao, Gu, Esangbedo & Dominique (2020) doi:10.3934/jimo.2020115
//
// Driven by R/maxmean.R::MaxMean(), which handles argument validation and
// .AsDistMatrix() coercion.
//
// Objective: maximise f(S) = sum_{i<j, i,j in S} d_ij / |S|  (Eq. 1).
// |S| is not fixed a priori; any subset of size >= 2 is feasible (Eqs. 2-3).
//
// Key data structures (Section 2.3):
//   p[i]       = sum_{j in S, j != i} d_ij  ("contribution" of element i)
//   sum_pairs  = sum_{i<j, i,j in S} d_ij   = sum_i p[i] / 2
//   f          = sum_pairs / |S|
//
// Flip deltas (Eq. 9):
//   Add    i: delta = (p[i] - f) / (|S| + 1)
//   Remove i: delta = (f - p[i]) / (|S| - 1)
//
// P-array update after flipping u (Eq. 10):
//   For j != u:  p[j] += d_{u,j}  (add)  or  p[j] -= d_{u,j}  (remove)
//   p[u] is unchanged (it already equals sum_{j in S_before} d_{u,j}).
//
// The RL layer (Sections 2.2, 2.4) guides initial-solution construction across
// restarts.  Q[s*n+a] accumulates state-action experience; R[s*n+a] stores
// the reward for selecting element a after previously selecting element s.
// Both matrices are n x n and persist across restarts within one call.

#include <Rcpp.h>
#include <vector>
#include <chrono>
#include <algorithm>
#include <limits>

using namespace Rcpp;

// [[Rcpp::export]]
List MaxMean_cpp(NumericMatrix dmat,
                 double time_budget_s,
                 double iter_budget,   // total tabu-iteration cap (Inf = off)
                 int    alpha_depth,   // no-improve limit per restart (paper: 50000)
                 int    T_min,         // minimum tabu tenure
                 int    T_max,         // maximum tabu tenure (paper: 120)
                 double epsilon,       // RL greedy probability (paper: 0.70)
                 double alpha_rl,      // Q-learning rate      (paper: 0.50)
                 double gamma_rl,      // Q discount factor    (paper: 0.50)
                 bool   use_rl) {      // enable Q-learning initialization

  const int n = dmat.nrow();
  if (n != dmat.ncol()) stop("dmat must be square");
  if (n < 2) stop("dmat must have at least 2 rows");
  const double *dp = REAL(dmat);
  // Column-major access: dmat(i,j) = dp[i + j*n]; column j = dp + j*n.
  // All d_{u,j} reads below use the contiguous column of u (symmetric).

  // Sync C-level RNG with R's seed (Rcpp does not do this automatically): load
  // .Random.seed on construction, write it back on destruction, so set.seed()
  // makes the search reproducible and the stream is consumed exactly once.
  Rcpp::RNGScope rngScope;

  // ---- RL matrices (n x n, row-major) ------------------------------------
  // Allocated only when use_rl is true to avoid O(n^2) overhead otherwise.
  const std::size_t n2 = use_rl ? (std::size_t)n * n : 1;
  std::vector<double> Q(n2, 0.0);   // experience matrix; Q[s*n + a]
  std::vector<double> R(n2, 0.0);   // reward matrix; Q[s*n + a]. Zero per paper
                                    // Section 2.2.1; the reward *magnitude* (1,
                                    // raised on improvement) is applied in the
                                    // Section 2.4 update, not at initialisation.

  // ---- Global best -------------------------------------------------------
  std::vector<int> best_S_global(n, 0);
  double best_f_global  = -1e300;
  long long total_iters = 0;
  int r_val = 1;   // reward magnitude; incremented on global improvement

  // ---- Reusable per-restart buffers -------------------------------------
  std::vector<int>       in_S(n, 0);
  std::vector<double>    p(n, 0.0);   // contribution array
  std::vector<long long> tabu_until(n, 0);
  std::vector<int>       Sa(n);       // RL optional-action list

  // Iteration budget (Eq.-independent stopping criterion, exposed by the wrapper
  // as `maxIter`). Inf disables it; any astronomically large finite value is
  // treated as unlimited too, so the (long long) cast below never overflows.
  const bool iter_unlimited = !R_finite(iter_budget) || iter_budget > 9.0e18;

  auto t0 = std::chrono::steady_clock::now();
  const int check_every = 256;

  auto elapsed_s = [&]() {
    return std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
  };

  // ========================================================================
  //  Outer restart loop
  // ========================================================================
  for (int restart = 0; ; ++restart) {

    Rcpp::checkUserInterrupt();

    // ---- Build initial solution S0 (Section 2.2) -----------------------
    std::fill(in_S.begin(), in_S.end(), 0);
    int m = 0;

    if (restart == 0 || !use_rl) {
      // First restart: random (each element with probability 0.5, Sect. 2.2.1).
      for (int i = 0; i < n; ++i) {
        if (R::runif(0.0, 1.0) < 0.5) { in_S[i] = 1; ++m; }
      }
      // Guarantee at least 2 elements.
      for (int i = 0; m < 2 && i < n; ++i) {
        if (!in_S[i]) { in_S[i] = 1; ++m; }
      }
    } else {
      // Subsequent restarts: RL-guided construction (Algorithm 2).
      for (int i = 0; i < n; ++i) Sa[i] = i;
      int Sa_num = n;

      // Randomly select initial state.
      int init_pos = (int)(R::runif(0.0, 1.0) * Sa_num);
      if (init_pos >= Sa_num) init_pos = Sa_num - 1;
      int state = Sa[init_pos];
      in_S[state] = 1; ++m;
      // Remove state from Sa (compact shift).
      for (int k = init_pos; k < Sa_num - 1; ++k) Sa[k] = Sa[k + 1];
      --Sa_num;

      while (Sa_num > 0) {
        // Action selection: epsilon-greedy on Q[state, *] over Sa.
        int action_pos;
        if (R::runif(0.0, 1.0) < epsilon) {
          // Greedy: argmax Q[state, a] for a in Sa.
          double best_q = -1e300;
          action_pos = 0;
          for (int k = 0; k < Sa_num; ++k) {
            double q = Q[(std::size_t)state * n + Sa[k]];
            if (q > best_q) { best_q = q; action_pos = k; }
          }
        } else {
          // Random.
          action_pos = (int)(R::runif(0.0, 1.0) * Sa_num);
          if (action_pos >= Sa_num) action_pos = Sa_num - 1;
        }
        int action = Sa[action_pos];

        // Reward for (state, action).
        double r = R[(std::size_t)state * n + action];

        // max Q[action, a'] over remaining actions Sa \ {action_pos}
        // (Eq. 4). Use the true max even when all candidate Q are negative; a
        // terminal next state (no remaining actions) contributes 0.
        double max_next_q = R_NegInf;
        for (int k = 0; k < Sa_num; ++k) {
          if (k == action_pos) continue;
          double q = Q[(std::size_t)action * n + Sa[k]];
          if (q > max_next_q) max_next_q = q;
        }
        if (!R_finite(max_next_q)) max_next_q = 0.0;   // terminal: no successor

        // Q-update (Eq. 4).
        double q_old = Q[(std::size_t)state * n + action];
        double q_new = (1.0 - alpha_rl) * q_old +
                       alpha_rl * (r + gamma_rl * max_next_q);

        // Termination check: stop extending the initial solution once the
        // comprehensive reward stops growing (and the state has been visited).
        // The paper's Eq. 5 and Eq. 7 conflict (the (1-alpha) factor); we follow
        // Eq. 5. Gated on m >= 2 so a feasible solution is always produced.
        if (m >= 2 &&
            q_old != 0.0 &&
            q_new < (1.0 - alpha_rl) * q_old) {
          break;
        }

        // Commit action.
        in_S[action] = 1; ++m;
        Q[(std::size_t)state * n + action] = q_new;
        state = action;

        // Remove action from Sa.
        for (int k = action_pos; k < Sa_num - 1; ++k) Sa[k] = Sa[k + 1];
        --Sa_num;
      }
    }

    // Save S0 for reward update later.
    std::vector<int> S0(in_S);

    // ---- Compute P array and initial f ------------------------------------
    std::fill(p.begin(), p.end(), 0.0);
    for (int i = 0; i < n; ++i) {
      if (!in_S[i]) continue;
      const double *col_i = dp + (std::size_t)i * n;
      // Branchless: the R wrapper zeroed the diagonal, so the j == i term adds
      // 0 and needs no guard. (col_i[i] == d(i,i) == 0.)
      for (int j = 0; j < n; ++j) {
        p[j] += col_i[j];
      }
    }
    double sum_pairs = 0.0;
    for (int i = 0; i < n; ++i) {
      if (in_S[i]) sum_pairs += p[i];
    }
    sum_pairs *= 0.5;
    double f = (m >= 2) ? sum_pairs / m : 0.0;

    // ---- Tabu search (Section 2.3) ----------------------------------------
    std::vector<int> best_S_local(in_S);
    double best_f_local = f;
    long long depth = 0;
    long long iter  = 0;
    int countdown   = check_every;

    // Fresh tabu state per restart: the solution is rebuilt from scratch, so
    // tenures earned in a previous restart must not carry over.
    std::fill(tabu_until.begin(), tabu_until.end(), 0);

    // Per-restart slice of the global iteration budget. Checking `iter <
    // iter_cap` each step makes maxIter an *exact* cap (total_iters never
    // exceeds it), unlike the every-256 time check which may overshoot.
    const long long iter_cap = iter_unlimited
      ? std::numeric_limits<long long>::max()
      : std::max((long long)0, (long long)iter_budget - total_iters);

    while (depth < (long long)alpha_depth && iter < iter_cap) {
      --countdown;
      if (countdown == 0) {
        Rcpp::checkUserInterrupt();
        if (elapsed_s() >= time_budget_s) break;
        countdown = check_every;
      }

      // Best-improvement move selection over the admissible set = {non-tabu
      // moves} U {tabu moves meeting the aspiration criterion (f + delta reaches
      // a new local best)}. Non-improving moves are admissible too; accepting
      // the least-bad move is how tabu search escapes local optima (the depth
      // counter restarts after alpha_depth such moves).
      //
      // delta is monotone in p[i] (add: increasing, since /(m+1) > 0; remove:
      // decreasing, since /(m-1) > 0), so the best add is the admissible
      // non-member of largest p[i] and the best remove the admissible member of
      // smallest p[i]. Tracking those two extrema needs no per-element division;
      // the aspiration test becomes a precomputed threshold on p[i]
      //   add:    f + (p-f)/(m+1) > best  <=>  p > f + gap*(m+1)
      //   remove: f + (f-p)/(m-1) > best  <=>  p < f - gap*(m-1),   gap = best - f >= 0
      // so only the two surviving extrema are divided. Ties break to the
      // smallest index (first encountered), matching the prior scan.
      const double gap      = best_f_local - f;           // >= 0 (f <= best_f_local)
      const double addThresh = f + gap * (double)(m + 1);
      const double remThresh = f - gap * (double)(m - 1);
      double maxAddP = R_NegInf; int addIdx = -1;
      double minRemP = R_PosInf; int remIdx = -1;

      for (int i = 0; i < n; ++i) {
        const double pi = p[i];
        const bool is_tabu = tabu_until[i] > iter;
        if (!in_S[i]) {
          if ((!is_tabu || pi > addThresh) && pi > maxAddP) { maxAddP = pi; addIdx = i; }
        } else {
          if (m <= 2) continue;                            // keep |S| >= 2
          if ((!is_tabu || pi < remThresh) && pi < minRemP) { minRemP = pi; remIdx = i; }
        }
      }

      int    best_flip  = -1;
      double best_delta = -1e300;
      if (addIdx >= 0) {
        double d = (maxAddP - f) / (double)(m + 1);
        if (d > best_delta) { best_delta = d; best_flip = addIdx; }
      }
      if (remIdx >= 0) {
        double d = (f - minRemP) / (double)(m - 1);
        if (d > best_delta) { best_delta = d; best_flip = remIdx; }
      }

      if (best_flip < 0) break;

      // Execute flip.
      int    u   = best_flip;
      double p_u = p[u];   // captured before P-array update
      const double *col_u = dp + (std::size_t)u * n;

      if (!in_S[u]) {
        // Add u. Branchless: col_u[u] == 0 (zeroed diagonal), so p[u] += 0.
        in_S[u] = 1; ++m;
        for (int j = 0; j < n; ++j) {
          p[j] += col_u[j];
        }
        sum_pairs += p_u;
        f = sum_pairs / m;
      } else {
        // Remove u. Branchless: col_u[u] == 0, so p[u] -= 0.
        in_S[u] = 0; --m;
        for (int j = 0; j < n; ++j) {
          p[j] -= col_u[j];
        }
        sum_pairs -= p_u;
        f = (m >= 2) ? sum_pairs / m : 0.0;
      }

      // Dynamic tabu tenure: sawtooth T_min..T_max (Galinier / Lai & Hao).
      int tenure = T_min + (int)(iter % (long long)(T_max - T_min + 1));
      tabu_until[u] = iter + (long long)tenure;

      // Track improvement.
      if (f > best_f_local) {
        best_f_local = f;
        best_S_local = in_S;
        depth = 0;
      } else {
        ++depth;
      }
      ++iter;
    }
    total_iters += iter;

    // ---- Update global best -----------------------------------------------
    bool improved_global = (best_f_local > best_f_global);
    if (improved_global) {
      best_f_global = best_f_local;
      best_S_global = best_S_local;
    }

    // ---- RL reward-matrix update (Section 2.4) ----------------------------
    // Elements in S0 that survive tabu get positive reward; dropped elements
    // get negative reward.  Magnitude is the base 1, raised on improvement.
    if (use_rl) {
      double cur_r = improved_global ? (double)r_val : 1.0;
      for (int i = 0; i < n; ++i) {
        if (!S0[i]) continue;
        for (int j = 0; j < n; ++j) {
          if (!S0[j] || j == i) continue;
          R[(std::size_t)i * n + j] = best_S_local[j] ? cur_r : -cur_r;
        }
      }
    }
    // Raise the reward magnitude *after* use, so the first improvement rewards
    // at the base magnitude 1 (paper Section 2.4) rather than 2.
    if (improved_global) ++r_val;

    // Budget check at the END of the cycle (the paper's Algorithm 1 is a
    // do-while), so at least one full restart always completes and the returned
    // selection satisfies |S| >= 2 even for a vanishingly small budget. Either
    // budget stops the search; the iteration check is load-bearing when the
    // time budget is Inf (without it, once iter_cap hits 0 the loop would spin
    // forever building fresh S0s that run zero tabu iterations).
    if (elapsed_s() >= time_budget_s) break;
    if (!iter_unlimited && total_iters >= iter_budget) break;
  } // restart loop

  // ---- Pack output (1-based indices for R) --------------------------------
  int final_m = 0;
  for (int i = 0; i < n; ++i) final_m += best_S_global[i];

  IntegerVector idx(final_m);
  for (int i = 0, pos = 0; i < n; ++i) {
    if (best_S_global[i]) idx[pos++] = i + 1;
  }

  return List::create(
    _["indices"]   = idx,
    _["objective"] = best_f_global,
    _["iters"]     = total_iters
  );
}

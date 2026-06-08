# Tests for DropAddTSPoints (matrix-free DropAdd-TS; Porumbel, Hao & Glover
# 2011). The matrix-free kernel recomputes each distance column from
# coordinates on the fly, so it carries no dense n x n matrix and runs at n far
# beyond the matrix path's ceiling.
#
# The reference is the matrix path DropAddTS(as.matrix(dist(points))). The two
# differ only in the construction seed: the matrix path uses Porumbel's
# max-row-sum seed (O(n^2*dim)), the matrix-free path uses the O(n*dim)
# farthest-from-centroid proxy. Empirically the proxy picks the SAME seed as
# Porumbel on ~97-98% of random Gaussian instances (so the whole trajectory is
# then bit-identical); on the minority where the seeds differ, the final
# objective is within ~3.7% (worst case observed over ~800 random instances at
# a 1-3 s budget). The comparability tests below therefore use a 5% tolerance.

# Brute-force MaxMin objective for a selection, from coordinates.
.maxmin_pts <- function(S, pts) {
  if (length(S) < 2L) return(NA_real_)
  min(stats::dist(pts[S, , drop = FALSE]))
}

# The two seed rules, for the coincidence/divergence assertions. unname() drops
# the dimname which.max inherits from a named rowSums, so the two integer seeds
# compare on value alone.
.rowsum_seed <- function(pts) {
  unname(which.max(rowSums(as.matrix(stats::dist(pts)))))
}
.centroid_seed <- function(pts) {
  unname(which.max(rowSums(sweep(pts, 2L, colMeans(pts))^2)))
}

# ---------------------------------------------------------------------------
# 1. Returned objective equals the true min-pairwise distance over the indices
# ---------------------------------------------------------------------------
test_that("DropAddTSPoints objective equals recomputed min-pairwise distance", {
  set.seed(7)
  pts <- matrix(rnorm(80 * 5), ncol = 5)
  res <- DropAddTSPoints(pts, m = 8L, max_no_improve = 500L)
  expect_length(res$indices, 8L)
  expect_equal(length(unique(res$indices)), 8L)
  expect_true(all(res$indices %in% seq_len(80L)))
  # The reported objective is exactly the MaxMin over the returned selection.
  expect_equal(res$objective, .maxmin_pts(res$indices, pts), tolerance = 1e-12)
  expect_gt(res$objective, 0)
  expect_gte(res$iters, 1L)
})

# ---------------------------------------------------------------------------
# 2a. Comparable quality to the matrix path on a COINCIDENT-seed instance.
#
# On instances where the centroid proxy picks Porumbel's max-row-sum seed
# (the common case), the matrix-free trajectory matches the matrix path
# exactly: the on-the-fly Euclidean distances reproduce stats::dist()'s bits,
# so every argmax resolves to the same index. We assert seed coincidence, then
# the identical final objective.
# ---------------------------------------------------------------------------
test_that("DropAddTSPoints matches matrix path when seeds coincide", {
  set.seed(1)
  pts <- matrix(rnorm(120 * 5), ncol = 5)       # seed 1 -> coincident seed (15)
  expect_equal(.centroid_seed(pts), .rowsum_seed(pts))   # documents the regime

  dmat <- as.matrix(stats::dist(pts))
  # Deterministic convergence (stagnation) on both paths: with a coincident
  # seed the trajectories are bit-identical and stop at the same iteration.
  rm <- DropAddTSPoints(pts, m = 12L, max_no_improve = 1000L)
  rx <- DropAddTS(dmat, m = 12L, max_no_improve = 1000L)
  # Identical objective (bit-identical under a matched-FP toolchain; a hair of
  # tolerance guards against an aggressive-FP build where dist contracts FMAs).
  expect_equal(rm$objective, rx$objective, tolerance = 1e-9)
  expect_equal(rm$objective, .maxmin_pts(rm$indices, pts), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 2b. Comparable quality on a confirmed-DIVERGENT-seed instance.
#
# Seed 75 (n=150, dim=6) is a confirmed case where the proxy seed differs from
# Porumbel's (centroid_seed=121, rowsum_seed=45). This exercises the
# non-identical path the proxy deviation creates. The final objective stays
# within 5% of the matrix path (observed ratio 0.988 at this budget; the worst
# divergent case over ~800 random instances was 3.7%).
# ---------------------------------------------------------------------------
test_that("DropAddTSPoints is within tolerance when seeds diverge", {
  set.seed(75)
  pts <- matrix(rnorm(150 * 6), ncol = 6)
  dmat <- as.matrix(stats::dist(pts))
  expect_true(.centroid_seed(pts) != .rowsum_seed(pts))  # confirms divergence

  tol <- 0.05    # documented comparability tolerance (5%)
  rm <- DropAddTSPoints(pts, m = 15L, max_no_improve = 1000L)
  rx <- DropAddTS(dmat, m = 15L, max_no_improve = 1000L)
  expect_gte(rm$objective, rx$objective * (1 - tol))
  expect_equal(rm$objective, .maxmin_pts(rm$indices, pts), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 3. Time budget honoured.
# ---------------------------------------------------------------------------
test_that("DropAddTSPoints respects time_budget_s within reasonable slack", {
  set.seed(99)
  pts <- matrix(runif(2000 * 8), ncol = 8)
  t0 <- Sys.time()
  # Disable stagnation so the wall-clock ceiling is the binding criterion.
  res <- DropAddTSPoints(pts, m = 20L, time_budget_s = 1, max_no_improve = 100000000L)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_lte(res$time_s, 1.5)
  expect_lte(elapsed, 2.0)
  expect_length(res$indices, 20L)
  expect_gte(res$iters, 1L)
})

# ---------------------------------------------------------------------------
# 4. max_iter = 0 returns the constructive solution, scoring it correctly.
# ---------------------------------------------------------------------------
test_that("DropAddTSPoints construction-only result matches its MaxMin score", {
  set.seed(11)
  pts <- matrix(rnorm(60 * 4), ncol = 4)
  res <- DropAddTSPoints(pts, m = 6L, max_iter = 0L)
  expect_length(res$indices, 6L)
  expect_equal(res$iters, 0L)
  expect_equal(res$objective, .maxmin_pts(res$indices, pts), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 5. Scales past the matrix path's memory ceiling.
#
# n = 5000 already exceeds a comfortable in-test dense matrix; the matrix-free
# path stays O(n) memory and returns a scored selection inside the budget. (The
# manuscript-scale targets are n = 20k / 58k; this is a fast CI proxy.)
# ---------------------------------------------------------------------------
test_that("DropAddTSPoints returns a meaningful result at n = 5000", {
  set.seed(5000)
  pts <- matrix(rnorm(5000 * 12), ncol = 12)
  # Scale smoke: wall-clock-bounded run (stagnation disabled) confirming the
  # matrix-free path stays responsive at n = 5000 within a short budget.
  res <- DropAddTSPoints(pts, m = 10L, time_budget_s = 2, max_no_improve = 100000000L)
  expect_length(res$indices, 10L)
  expect_equal(length(unique(res$indices)), 10L)
  expect_gt(res$objective, 0)
  expect_equal(res$objective, .maxmin_pts(res$indices, pts), tolerance = 1e-12)
  expect_gte(res$iters, 100L)   # not iteration-starved
})

# ---------------------------------------------------------------------------
# 6. seed parameter, input validation, and progress output
# ---------------------------------------------------------------------------
test_that("DropAddTSPoints seed parameter and input validation", {
  pts <- matrix(rnorm(30 * 3), ncol = 3)

  # seed: set.seed(1) path is executed
  r1 <- DropAddTSPoints(pts, m = 4L, max_iter = 0L, seed = 1L)
  r2 <- DropAddTSPoints(pts, m = 4L, max_iter = 0L, seed = 1L)
  expect_identical(r1$indices, r2$indices)

  # m validation
  expect_error(DropAddTSPoints(pts, m = 1L), "2 <= m")
  expect_error(DropAddTSPoints(pts, m = 100L), "2 <= m")

  # max_no_improve validation
  expect_error(DropAddTSPoints(pts, m = 4L, max_no_improve = 0L),
               "max_no_improve")

  # max_iter validation
  expect_error(DropAddTSPoints(pts, m = 4L, max_iter = -1L), "max_iter")

  # time_budget_s validation
  expect_error(DropAddTSPoints(pts, m = 4L, time_budget_s = 0),
               "time_budget_s")
  expect_error(DropAddTSPoints(pts, m = 4L, time_budget_s = -1),
               "time_budget_s")
  expect_error(DropAddTSPoints(pts, m = 4L, time_budget_s = NA_real_),
               "time_budget_s")
})

test_that("DropAddTSPoints progress = TRUE fires the cli hooks", {
  pts <- matrix(rnorm(20 * 2), ncol = 2)
  expect_no_error(
    DropAddTSPoints(pts, m = 3L, max_iter = 2L, progress = TRUE)
  )
})

# ---------------------------------------------------------------------------
# 7. C++ path coverage (dropadd_mf.cpp)
# ---------------------------------------------------------------------------

test_that("DropAddTSPoints C++ construction covers sum_dist tie-break (lines 131-134)", {
  # 5-point geometry: A=(0,0), B=(3,0), C=(0,2), P=(1,1), Q=(2,1).
  # Centroid=(1.2,0.8); B is farthest -> seed=B.
  # Construction: B->C->A, then P and Q are tied on min_dist=sqrt(2)
  # but sum_dist[Q]=2*sqrt(5)+sqrt(2) > sum_dist[P]=2*sqrt(2)+sqrt(5) -> Q wins.
  # When A is added: d(P,A)=sqrt(2)==min_dist[P] -> line 152 (equality) fires.
  pts5 <- rbind(c(0, 0), c(3, 0), c(0, 2), c(1, 1), c(2, 1))
  res <- DropAddTSPoints(pts5, m = 4L, max_iter = 0L)
  expect_length(res$indices, 4L)
  expect_equal(res$iters, 0L)
})

test_that("DropAddTSPoints C++ m=2 covers DROP else-branch (lines 263-264)", {
  # Rhombus: (0,0),(1,1),(2,0),(1,-1). All from-centroid distances equal ->
  # seed = index 1 (ties -> smallest). Construction: {(0,0),(2,0)}.
  # Drop (0,0): (2,0) loses its sole selected peer; need_recompute fires;
  # self-mask inside recompute leaves mns=Inf -> else branch (lines 263-264).
  pts_rh <- rbind(c(0, 0), c(1, 1), c(2, 0), c(1, -1))
  res <- DropAddTSPoints(pts_rh, m = 2L, max_iter = 4L)
  expect_length(res$indices, 2L)
})

test_that("DropAddTSPoints C++ main-loop ADD covers tie-break (304-307) and equality (327)", {
  # 7-point: P=(0,0),Q=(4,0),R=(0,4),ANC=(10,10),X=(1,0),Y=(3,0),W=(3.5,0).
  # ANC is farthest from centroid -> seed. Construction: ANC->P->Q->R.
  # Drop ANC: X and Y tied on min_dist=1; sum_dist[Y]=9>sum_dist[X]~8.12 ->
  # lines 304-307 fire. When Y added: d(Y,W)=0.5==min_dist[W]=0.5 -> line 327.
  pts7 <- rbind(c(0, 0), c(4, 0), c(0, 4), c(10, 10), c(1, 0), c(3, 0), c(3.5, 0))
  res <- DropAddTSPoints(pts7, m = 4L, max_iter = 1L)
  expect_length(res$indices, 4L)
  expect_equal(res$iters, 1L)
})

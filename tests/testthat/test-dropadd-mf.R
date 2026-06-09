# Tests for DropAdd points path (matrix-free DropAdd-TS; Porumbel, Hao &
# Glover 2011). The matrix-free kernel recomputes each distance column from
# coordinates on the fly, so it carries no dense n x n matrix and runs at n far
# beyond the matrix path's ceiling.
#
# The reference is the matrix path DropAdd(as.matrix(dist(points))). The two
# differ only in the construction seed: the matrix path uses Porumbel's
# max-row-sum seed (O(n^2*dim)), the matrix-free path uses the O(n*dim)
# farthest-from-centroid proxy. Empirically the proxy picks the SAME seed as
# Porumbel on ~97-98% of random Gaussian instances (so the whole trajectory is
# then bit-identical); on the minority where the seeds differ, the final
# objective is within ~3.7% (worst case observed over ~800 random instances at
# a 1-3 s budget). The comparability tests below therefore use a 5% tolerance.

# Brute-force MaxMin objective for a selection, from coordinates.
.MaxminPts <- function(S, pts) {
  if (length(S) < 2L) return(NA_real_)
  min(stats::dist(pts[S, , drop = FALSE]))
}

# The two seed rules, for the coincidence/divergence assertions. unname() drops
# the dimname which.max inherits from a named rowSums, so the two integer seeds
# compare on value alone.
.RowsumSeed <- function(pts) {
  unname(which.max(rowSums(as.matrix(stats::dist(pts)))))
}
.CentroidSeed <- function(pts) {
  unname(which.max(rowSums(sweep(pts, 2L, colMeans(pts))^2)))
}

# ---------------------------------------------------------------------------
# 1. Returned objective equals the true min-pairwise distance over the indices
# ---------------------------------------------------------------------------
test_that("DropAdd points path: objective equals recomputed min-pairwise distance", {
  set.seed(7)
  pts <- matrix(rnorm(80 * 5), ncol = 5)
  res <- DropAdd(points = pts, m = 8L, plateau = 500L)
  expect_length(res, 8L)
  expect_equal(length(unique(res)), 8L)
  expect_true(all(res %in% seq_len(80L)))
  # The reported objective is exactly the MaxMin over the returned selection.
  expect_equal(attr(res, "score"), .MaxminPts(res, pts), tolerance = 1e-12)
  expect_gt(attr(res, "score"), 0)
  expect_gte(attr(res, "iters"), 1L)
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
test_that("DropAdd points path matches matrix path when seeds coincide", {
  set.seed(1)
  pts <- matrix(rnorm(120 * 5), ncol = 5)       # seed 1 -> coincident seed (15)
  expect_equal(.CentroidSeed(pts), .RowsumSeed(pts))   # documents the regime

  dmat <- as.matrix(stats::dist(pts))
  # Deterministic convergence (stagnation) on both paths: with a coincident
  # seed the trajectories are bit-identical and stop at the same iteration.
  rm <- DropAdd(points = pts, m = 12L, plateau = 1000L)
  rx <- DropAdd(dmat, m = 12L, plateau = 1000L)
  # Identical objective (bit-identical under a matched-FP toolchain; a hair of
  # tolerance guards against an aggressive-FP build where dist contracts FMAs).
  expect_equal(attr(rm, "score"), attr(rx, "score"), tolerance = 1e-9)
  expect_equal(attr(rm, "score"), .MaxminPts(rm, pts), tolerance = 1e-12)
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
test_that("DropAdd points path is within tolerance when seeds diverge", {
  set.seed(75)
  pts <- matrix(rnorm(150 * 6), ncol = 6)
  dmat <- as.matrix(stats::dist(pts))
  expect_true(.CentroidSeed(pts) != .RowsumSeed(pts))  # confirms divergence

  tol <- 0.05    # documented comparability tolerance (5%)
  rm <- DropAdd(points = pts, m = 15L, plateau = 1000L)
  rx <- DropAdd(dmat, m = 15L, plateau = 1000L)
  expect_gte(attr(rm, "score"), attr(rx, "score") * (1 - tol))
  expect_equal(attr(rm, "score"), .MaxminPts(rm, pts), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 3. Time budget honoured.
# ---------------------------------------------------------------------------
test_that("DropAdd points path respects timeBudgetS within reasonable slack", {
  set.seed(99)
  pts <- matrix(runif(2000 * 8), ncol = 8)
  t0 <- Sys.time()
  # Disable stagnation so the wall-clock ceiling is the binding criterion.
  res <- DropAdd(points = pts, m = 20L, timeBudgetS = 1, plateau = 100000000L)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_lte(attr(res, "time_s"), 1.5)
  expect_lte(elapsed, 2.0)
  expect_length(res, 20L)
  expect_gte(attr(res, "iters"), 1L)
})

# ---------------------------------------------------------------------------
# 4. maxIter = 0 returns the constructive solution, scoring it correctly.
# ---------------------------------------------------------------------------
test_that("DropAdd points path construction-only result matches its MaxMin score", {
  set.seed(11)
  pts <- matrix(rnorm(60 * 4), ncol = 4)
  res <- DropAdd(points = pts, m = 6L, maxIter = 0L)
  expect_length(res, 6L)
  expect_equal(attr(res, "iters"), 0L)
  expect_equal(attr(res, "score"), .MaxminPts(res, pts), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 5. Scales past the matrix path's memory ceiling.
#
# n = 5000 already exceeds a comfortable in-test dense matrix; the matrix-free
# path stays O(n) memory and returns a scored selection inside the budget. (The
# manuscript-scale targets are n = 20k / 58k; this is a fast CI proxy.)
# ---------------------------------------------------------------------------
test_that("DropAdd points path returns a meaningful result at n = 5000", {
  set.seed(5000)
  pts <- matrix(rnorm(5000 * 12), ncol = 12)
  # Scale smoke: wall-clock-bounded run (stagnation disabled) confirming the
  # matrix-free path stays responsive at n = 5000 within a short budget.
  res <- DropAdd(points = pts, m = 10L, timeBudgetS = 2, plateau = 100000000L)
  expect_length(res, 10L)
  expect_equal(length(unique(res)), 10L)
  expect_gt(attr(res, "score"), 0)
  expect_equal(attr(res, "score"), .MaxminPts(res, pts), tolerance = 1e-12)
  expect_gte(attr(res, "iters"), 100L)   # not iteration-starved
})

# ---------------------------------------------------------------------------
# 6. seed parameter, input validation, and progress output
# ---------------------------------------------------------------------------
test_that("DropAdd points path: seed parameter and input validation", {
  pts <- matrix(rnorm(30 * 3), ncol = 3)

  # seed: set.seed(1) path is executed
  r1 <- DropAdd(points = pts, m = 4L, maxIter = 0L, seed = 1L)
  r2 <- DropAdd(points = pts, m = 4L, maxIter = 0L, seed = 1L)
  expect_identical(r1, r2)

  # m validation
  expect_error(DropAdd(points = pts, m = 1L), "2 <= m")
  expect_error(DropAdd(points = pts, m = 100L), "2 <= m")

  # plateau validation
  expect_error(DropAdd(points = pts, m = 4L, plateau = 0L),
               "plateau")

  # maxIter validation
  expect_error(DropAdd(points = pts, m = 4L, maxIter = -1L), "maxIter")

  # timeBudgetS validation
  expect_error(DropAdd(points = pts, m = 4L, timeBudgetS = 0),
               "timeBudgetS")
  expect_error(DropAdd(points = pts, m = 4L, timeBudgetS = -1),
               "timeBudgetS")
  expect_error(DropAdd(points = pts, m = 4L, timeBudgetS = NA_real_),
               "timeBudgetS")

  # both d and points supplied
  expect_error(DropAdd(d = as.matrix(dist(pts)), points = pts, m = 4L),
               "supply")
})

test_that("DropAdd points path: progress = TRUE fires the cli hooks", {
  pts <- matrix(rnorm(20 * 2), ncol = 2)
  expect_no_error(
    DropAdd(points = pts, m = 3L, maxIter = 2L, progress = TRUE)
  )
})

# ---------------------------------------------------------------------------
# 7. C++ path coverage (dropadd_mf.cpp)
# ---------------------------------------------------------------------------

test_that("DropAdd points path C++: construction covers sum_dist tie-break (lines 131-134)", {
  # 5-point geometry: A=(0,0), B=(3,0), C=(0,2), P=(1,1), Q=(2,1).
  # Centroid=(1.2,0.8); B is farthest -> seed=B.
  # Construction: B->C->A, then P and Q are tied on min_dist=sqrt(2)
  # but sum_dist[Q]=2*sqrt(5)+sqrt(2) > sum_dist[P]=2*sqrt(2)+sqrt(5) -> Q wins.
  # When A is added: d(P,A)=sqrt(2)==min_dist[P] -> line 152 (equality) fires.
  pts5 <- rbind(c(0, 0), c(3, 0), c(0, 2), c(1, 1), c(2, 1))
  res <- DropAdd(points = pts5, m = 4L, maxIter = 0L)
  expect_length(res, 4L)
  expect_equal(attr(res, "iters"), 0L)
})

test_that("DropAdd points path C++: m=2 covers DROP else-branch (lines 263-264)", {
  # Rhombus: (0,0),(1,1),(2,0),(1,-1). All from-centroid distances equal ->
  # seed = index 1 (ties -> smallest). Construction: {(0,0),(2,0)}.
  # Drop (0,0): (2,0) loses its sole selected peer; need_recompute fires;
  # self-mask inside recompute leaves mns=Inf -> else branch (lines 263-264).
  pts_rh <- rbind(c(0, 0), c(1, 1), c(2, 0), c(1, -1))
  res <- DropAdd(points = pts_rh, m = 2L, maxIter = 4L)
  expect_length(res, 2L)
})

test_that("DropAdd points path C++: main-loop ADD covers tie-break (304-307) and equality (327)", {
  # 7-point: P=(0,0),Q=(4,0),R=(0,4),ANC=(10,10),X=(1,0),Y=(3,0),W=(3.5,0).
  # ANC is farthest from centroid -> seed. Construction: ANC->P->Q->R.
  # Drop ANC: X and Y tied on min_dist=1; sum_dist[Y]=9>sum_dist[X]~8.12 ->
  # lines 304-307 fire. When Y added: d(Y,W)=0.5==min_dist[W]=0.5 -> line 327.
  pts7 <- rbind(c(0, 0), c(4, 0), c(0, 4), c(10, 10), c(1, 0), c(3, 0), c(3.5, 0))
  res <- DropAdd(points = pts7, m = 4L, maxIter = 1L)
  expect_length(res, 4L)
  expect_equal(attr(res, "iters"), 1L)
})

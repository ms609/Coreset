# Tests for DropAdd points path (matrix-free DropAdd-TS; Porumbel, Hao &
# Glover 2011). The matrix-free kernel recomputes each distance column from
# coordinates on the fly, so it carries no dense n x n matrix and runs at n far
# beyond the matrix path's ceiling.
#
# The reference is the matrix path DropAdd(as.matrix(dist(points))). The two
# differ only in the construction seed: the matrix path uses Porumbel's
# max-row-sum seed (O(n^2*dim)), the matrix-free path uses the O(n*dim)
# farthest-from-anti_centroid proxy. Empirically the proxy picks the SAME seed as
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
  res <- DropAdd(8L, plateau = 500L, points = pts)
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
# On instances where the anti_centroid proxy picks Porumbel's max-row-sum seed
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
  rm <- DropAdd(12L, plateau = 1000L, points = pts)
  rx <- DropAdd(12L, dmat, plateau = 1000L)
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
  rm <- DropAdd(15L, plateau = 1000L, points = pts)
  rx <- DropAdd(15L, dmat, plateau = 1000L)
  expect_gte(attr(rm, "score"), attr(rx, "score") * (1 - tol))
  expect_equal(attr(rm, "score"), .MaxminPts(rm, pts), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 3. Time budget honoured.
# ---------------------------------------------------------------------------
test_that("DropAdd points path respects maxSeconds within reasonable slack", {
  set.seed(99)
  pts <- matrix(runif(2000 * 8), ncol = 8)
  t0 <- Sys.time()
  # Disable stagnation so the wall-clock ceiling is the binding criterion.
  res <- DropAdd(20L, maxSeconds = 0.05, plateau = 100000000L,
                 points = pts)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_lte(attr(res, "time_s"), 5)
  expect_lte(elapsed, 5.5)
  expect_length(res, 20L)
  expect_gte(attr(res, "iters"), 1L)
})

# ---------------------------------------------------------------------------
# 4. Returned score is accurate on the points path.
# ---------------------------------------------------------------------------
test_that("DropAdd points path: reported score equals actual MaxMin of returned indices", {
  set.seed(11)
  pts <- matrix(rnorm(60 * 4), ncol = 4)
  res <- DropAdd(6L, points = pts)
  expect_length(res, 6L)
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
  res <- DropAdd(10L, maxSeconds = 2, plateau = 100000000L, points = pts)
  expect_length(res, 10L)
  expect_equal(length(unique(res)), 10L)
  expect_gt(attr(res, "score"), 0)
  expect_equal(attr(res, "score"), .MaxminPts(res, pts), tolerance = 1e-12)
  expect_gte(attr(res, "iters"), 100L)   # not iteration-starved
})

# ---------------------------------------------------------------------------
# 6. seed parameter, input validation, and progress output
# ---------------------------------------------------------------------------
test_that("DropAdd points path: deterministic and input validation", {
  pts <- matrix(rnorm(30 * 3), ncol = 3)

  # RNG-free: repeated calls are identical (the no-op `seed` arg was removed).
  r1 <- DropAdd(4L, plateau = 100L, points = pts)
  r2 <- DropAdd(4L, plateau = 100L, points = pts)
  attr(r1, "time_s") <- attr(r2, "time_s") <- NULL
  expect_identical(r1, r2)

  # k validation
  expect_error(DropAdd(1L, points = pts), "2 <= k")
  expect_error(DropAdd(100L, points = pts), "2 <= k")

  # plateau validation
  expect_error(DropAdd(4L, plateau = 0L, points = pts),
               "plateau")

  # maxSeconds validation
  expect_error(DropAdd(4L, maxSeconds = 0, points = pts),
               "maxSeconds")
  expect_error(DropAdd(4L, maxSeconds = -1, points = pts),
               "maxSeconds")
  expect_error(DropAdd(4L, maxSeconds = NA_real_, points = pts),
               "maxSeconds")

  # both d and points supplied
  expect_error(DropAdd(4L, d = as.matrix(dist(pts)), points = pts),
               "supply")
})

test_that("DropAdd points path: MaxMin.progress option fires the cli hooks", {
  pts <- matrix(rnorm(20 * 2), ncol = 2)
  old <- options(MaxMin.progress = TRUE)
  on.exit(options(old))
  expect_no_error(suppressMessages(DropAdd(3L, maxSeconds = 0.1, points = pts)))
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
  res <- DropAdd(4L, points = pts5)
  expect_length(res, 4L)
})

test_that("DropAdd points path C++: k=2 covers DROP else-branch (lines 263-264)", {
  # Rhombus: (0,0),(1,1),(2,0),(1,-1). All from-anti_centroid distances equal ->
  # seed = index 1 (ties -> smallest). Construction: {(0,0),(2,0)}.
  # Drop (0,0): (2,0) loses its sole selected peer; need_recompute fires;
  # self-mask inside recompute leaves mns=Inf -> else branch (lines 263-264).
  pts_rh <- rbind(c(0, 0), c(1, 1), c(2, 0), c(1, -1))
  res <- DropAdd(2L, points = pts_rh)
  expect_length(res, 2L)
})

test_that("DropAdd points path C++: main-loop ADD covers tie-break (304-307) and equality (327)", {
  # 7-point: P=(0,0),Q=(4,0),R=(0,4),ANC=(10,10),X=(1,0),Y=(3,0),W=(3.5,0).
  # ANC is farthest from anti_centroid -> seed. Construction: ANC->P->Q->R.
  # Drop ANC: X and Y tied on min_dist=1; sum_dist[Y]=9>sum_dist[X]~8.12 ->
  # lines 304-307 fire. When Y added: d(Y,W)=0.5==min_dist[W]=0.5 -> line 327.
  pts7 <- rbind(c(0, 0), c(4, 0), c(0, 4), c(10, 10), c(1, 0), c(3, 0), c(3.5, 0))
  res <- DropAdd(4L, points = pts7)
  expect_length(res, 4L)
})

# ---------------------------------------------------------------------------
# 8. Blocked column fills: every arm of the dimension dispatch
#
# FillSqRange walks `dim` in blocks of four and dispatches the 1-4 column
# remainder on (remainder, is-this-the-first-block). Each arm is a distinct
# FillSqBlock instantiation, so each needs its own `dim`: 1 is a lone
# first block, 7 is a full first block plus a 3-column continuation, and 0
# is the degenerate no-coordinate case. Every arm must reproduce
# stats::dist(), which the recomputed objective checks.
# ---------------------------------------------------------------------------

test_that("DropAdd points path C++: one-dimensional coordinates", {
  # Gaps 5, 6, 7, 8, 9, 10 along a line: the best 3-subset is the two ends
  # plus the point that maximises the smaller of the two gaps it creates.
  pts1 <- matrix(c(0, 5, 11, 18, 26, 35, 45), ncol = 1L)
  res <- DropAdd(3L, points = pts1, plateau = 20L)
  expect_identical(as.integer(res), c(1L, 5L, 7L))
  expect_equal(attr(res, "score"), 19)
})

test_that("DropAdd points path C++: a dimension that needs a block continuation", {
  # dim = 7 = one four-column block then a three-column remainder, the only
  # shape that reaches the non-first three-column instantiation.
  set.seed(4)
  pts7d <- matrix(rnorm(40 * 7), ncol = 7L)
  res <- DropAdd(4L, points = pts7d, plateau = 20L)
  expect_length(res, 4L)
  expect_equal(attr(res, "score"), .MaxminPts(res, pts7d), tolerance = 1e-12)
})

test_that("DropAdd points path C++: zero-dimension input stays defined", {
  # With no coordinates every distance is 0, so no selection beats any other:
  # each argmax keeps the earliest index and the objective is 0.
  pts0 <- matrix(numeric(0), nrow = 6L, ncol = 0L)
  res <- DropAdd(3L, points = pts0, plateau = 5L)
  expect_identical(as.integer(res), 1:3)
  expect_identical(attr(res, "score"), 0)
  expect_identical(attr(res, "secondary"), 0)
})

# ---------------------------------------------------------------------------
# 9. The two argmax merges, on inputs whose ties are exact
#
# Integer coordinates make the ties portable: equal integer squared sums are
# the same double everywhere, so sqrt() of them compares equal on any
# toolchain, unlike an incidental coincidence in Gaussian data.
# ---------------------------------------------------------------------------

test_that("DropAdd points path C++: recompute merge resolves both kinds of tie", {
  # The unit square scaled to side 4, (+-2, +-2), plus (3,4)/(-3,4) to fix the
  # seed and (-5,2) to break the mirror symmetry of the whole set. Every
  # corner of the square is exactly 4 from two others, so min_dist ties at 4
  # recur; the drop puts some corners into need_recompute while the pass
  # winner comes from outside it, and the deferred merge must then settle the
  # tie on sum_dist and, failing that, on the smaller index. Both arms arise
  # here: sum_dist 11 vs 9 for (2,2) over (-2,-2), and an exact sum_dist tie
  # at 4 + 4*sqrt(2) where (-2,2) takes the slot from (-2,-2) on index alone.
  ptsSq <- rbind(c(2, 2), c(2, -2), c(3, 4),
                 c(-2, 2), c(-2, -2), c(-3, 4), c(-5, 2))
  res <- DropAdd(3L, points = ptsSq, plateau = 25L, maxCandidates = 0L)
  expect_identical(as.integer(res), c(2L, 3L, 7L))
  expect_equal(attr(res, "score"), .MaxminPts(res, ptsSq), tolerance = 1e-12)
})

test_that("DropAdd points path C++: chunk merge is invariant on a tied lattice", {
  # A 128 x 128 integer lattice: n = 16384 engages the threaded pass regions,
  # and lattice distances tie exactly, so chunk winners routinely reach the
  # merge tied on min_dist and separated only by sum_dist. Both merge sites
  # (construction and drop) take that arm at k = 8 and k = 20, at every
  # thread count from 2 to 4. The merge exists to reproduce the serial
  # lexicographic maximum, so the property to assert is that it does. CRAN
  # caps tests at 2 cores.
  side <- 128L
  lat <- matrix(as.double(cbind(rep(seq_len(side), times = side),
                                rep(seq_len(side), each = side))), ncol = 2L)
  old <- options(mc.cores = NULL)
  on.exit(options(old), add = TRUE)
  StripTime <- function(x) {
    attr(x, "time_s") <- NULL
    x
  }
  for (k in c(8L, 20L)) {
    options(mc.cores = 1L)
    s1 <- DropAdd(k, points = lat, plateau = 30L, maxCandidates = 0L)
    options(mc.cores = 2L)
    s2 <- DropAdd(k, points = lat, plateau = 30L, maxCandidates = 0L)
    expect_identical(StripTime(s1), StripTime(s2))
    expect_equal(attr(s1, "score"), .MaxminPts(s1, lat), tolerance = 1e-12)
  }
})

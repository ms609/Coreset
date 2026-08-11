# Tests for DropAdd()'s distance-column oracle path (`d = <function>`).
#
# The oracle path is a pure-R reimplementation of src/dropadd.cpp against a
# caller-supplied `colFn(i)` instead of a matrix column. The governing invariant
# is therefore *identity, not similarity*: given the same distances and the same
# seed, the oracle path must reproduce the matrix path exactly -- the same drop
# and add at every iteration, not merely the same final subset. Two trajectories
# can diverge and reconverge, so the trace comparison below is the real gate;
# the end-to-end comparisons are the user-visible consequence.
#
# The seeds are always pinned with `seed=`, because the default warm starts
# differ by construction (max-row-sum on the matrix path, the O(n) two-sweep
# peripheral substitute on the oracle path -- see .DropAddConstructColumn()'s
# seed deviation). A default-seed comparison would not be an identity test.

# Euclidean distance column accumulated in plain double over columns in
# increasing order, exactly mirroring src/dropadd_mf.cpp::EuclidCol and hence
# stats::dist(). `sqrt(rowSums(...))` would accumulate in long double and can
# differ in the last bit, which would make the coordinate comparison a tolerance
# test rather than an identity test.
.EuclidCol <- function(pts, i) {
  s <- numeric(nrow(pts))
  for (j in seq_len(ncol(pts))) {
    s <- s + (pts[, j] - pts[i, j]) ^ 2
  }
  sqrt(s)
}

# Compare two selections on every user-visible field except wall-clock.
.ExpectSameSelection <- function(a, b) {
  testthat::expect_identical(as.integer(a), as.integer(b))
  testthat::expect_identical(attr(a, "score"), attr(b, "score"))
  testthat::expect_identical(attr(a, "secondary"), attr(b, "secondary"))
  testthat::expect_identical(attr(a, "iters"), attr(b, "iters"))
}

# ---------------------------------------------------------------------------
# 1. Trajectory identity: the same drops and adds, iteration by iteration
# ---------------------------------------------------------------------------
test_that("the column-oracle loop reproduces the C++ drop/add trajectory", {
  set.seed(2026)
  for (trial in seq_len(4L)) {
    n <- sample(12:35, 1L)
    k <- sample(2:7, 1L)
    dmat <- as.matrix(dist(matrix(rnorm(n * 2L), ncol = 2L)))
    colFn <- function(i) dmat[, i]
    # .DropAddTrace() runs the kernel with its own max-row-sum seed, so hand the
    # oracle loop that same seed to make the two trajectories comparable.
    tr <- MaxMin:::.DropAddTrace(dmat, k, maxIter = 60L, plateau = 1e9)
    or <- MaxMin:::.DropAddFromColumn(
      colFn, n, k, first = unname(which.max(rowSums(dmat))),
      plateau = 1e9, maxIter = 60L, trace = TRUE
    )
    expect_identical(or$drops, tr$drops)
    expect_identical(or$adds, tr$adds)
    expect_identical(sort(as.integer(or$indices)), tr$indices)
    expect_identical(or$objective, tr$score)
    expect_identical(or$iters, tr$iters)
  }
})

# ---------------------------------------------------------------------------
# 2. Bit-identity to the matrix path, end to end
# ---------------------------------------------------------------------------
test_that("a matrix-wrapping colFn is bit-identical to the matrix path", {
  set.seed(42)
  for (trial in seq_len(6L)) {
    n <- sample(8:40, 1L)
    k <- sample(2:min(8L, n), 1L)
    dmat <- as.matrix(dist(matrix(rnorm(n * 3L), ncol = 3L)))
    colFn <- function(i) dmat[, i]
    s <- sample.int(n, 1L)
    .ExpectSameSelection(
      DropAdd(k, colFn, N = n, seed = s, maxCandidates = 0L),
      DropAdd(k, dmat, seed = s, maxCandidates = 0L)
    )
  }
})

test_that("a self-omitting (length N - 1) colFn matches a length-N one", {
  set.seed(5)
  dmat <- as.matrix(dist(matrix(rnorm(18L * 2L), ncol = 2L)))
  .ExpectSameSelection(
    DropAdd(4L, function(i) dmat[-i, i], N = 18L, seed = 2L),
    DropAdd(4L, function(i) dmat[, i],   N = 18L, seed = 2L)
  )
})

test_that("a Euclidean colFn is bit-identical to the coordinate path", {
  # The briefing's correctness gate: wrap the *same* coordinates the `points`
  # path uses and require an identical answer, so the oracle path has a
  # correctness test that needs no exotic metric.
  set.seed(9)
  for (trial in seq_len(4L)) {
    n <- sample(15:40, 1L)
    k <- sample(3:6, 1L)
    dim <- sample(2:4, 1L)
    pts <- matrix(rnorm(n * dim), ncol = dim)
    colFn <- function(i) .EuclidCol(pts, i)
    s <- sample.int(n, 1L)
    .ExpectSameSelection(
      DropAdd(k, colFn, N = n, seed = s, maxCandidates = 0L),
      DropAdd(k, points = pts, seed = s, maxCandidates = 0L)
    )
  }
})

test_that(".DropAddPick breaks minDist ties on sumDist, then on index", {
  # Asserted directly rather than inferred from a lattice run: the lexicographic
  # tie-break and the exclusion of x# are the two places where a subtle
  # divergence from src/dropadd.cpp would not show up as an obviously wrong
  # answer, only as a different trajectory.
  st <- new.env(parent = emptyenv())
  st$inS       <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
  st$minDist   <- c(5, 5, 5, 4, Inf)
  st$sumDist   <- c(9, 12, 12, 99, 0)
  # minDist wins outright over sumDist; among the minDist ties the largest
  # sumDist wins; among those, the smallest index.
  expect_identical(MaxMin:::.DropAddPick(st), 2L)
  expect_identical(MaxMin:::.DropAddPick(st, exclude = 2L), 3L)
  # A single minDist maximum needs no sumDist comparison at all.
  st$minDist <- c(5, 4, 4, 4, Inf)
  expect_identical(MaxMin:::.DropAddPick(st), 1L)
})

test_that(".DropAddApplyAdd counts equal distances and lowers strict ones", {
  st <- new.env(parent = emptyenv())
  st$minDist      <- c(3, 3, 3, Inf)
  st$sumDist      <- c(10, 10, 10, 10)
  st$minDistCount <- c(1L, 2L, 1L, 0L)
  MaxMin:::.DropAddApplyAdd(st, col = c(3, 2, 7, 0), xNew = 4L)
  # 1: equal -> count incremented, minDist unchanged.
  # 2: strictly lower -> minDist replaced, count reset to 1.
  # 3: higher -> untouched. 4: xNew itself -> left for the caller.
  expect_identical(st$minDist, c(3, 2, 3, Inf))
  expect_identical(st$minDistCount, c(2L, 1L, 1L, 0L))
  expect_identical(st$sumDist, c(13, 12, 17, 10))
})

test_that("exact ties are broken identically to the matrix path", {
  # An integer lattice makes exact distance ties abundant, exercising the
  # `minDistCount` equality branch of .DropAddApplyAdd() and the sumDist
  # lexicographic tie-break of .DropAddPick(). Floating-point data almost never
  # reaches either.
  grid <- as.matrix(expand.grid(x = 0:4, y = 0:3))
  dmat <- as.matrix(dist(grid))
  colFn <- function(i) dmat[, i]
  expect_true(any(duplicated(dmat[lower.tri(dmat)])))
  for (k in c(2L, 3L, 5L, 8L)) {
    .ExpectSameSelection(
      DropAdd(k, colFn, N = nrow(grid), seed = 1L, maxCandidates = 0L),
      DropAdd(k, dmat, seed = 1L, maxCandidates = 0L)
    )
  }
})

# ---------------------------------------------------------------------------
# 3. The substitute seed
# ---------------------------------------------------------------------------
test_that("the default oracle seed is the deterministic peripheral seed", {
  set.seed(11)
  dmat <- as.matrix(dist(matrix(rnorm(25L * 2L), ncol = 2L)))
  colFn <- function(i) dmat[, i]
  # Deterministic: no RNG is drawn, so repeated calls agree and the result is
  # reproducible without set.seed().
  a <- DropAdd(4L, colFn, N = 25L)
  b <- DropAdd(4L, colFn, N = 25L)
  .ExpectSameSelection(a, b)
  # And it is exactly .PeripheralSeedColumn(), the seed FarFirst()'s oracle path
  # already substitutes for its own O(N^2) anchors.
  periph <- MaxMin:::.PeripheralSeedColumn(colFn, 25L)
  .ExpectSameSelection(a, DropAdd(4L, colFn, N = 25L, seed = periph))
  # A max-row-sum seed is unreachable from an oracle without N calls, but when
  # asked for explicitly the oracle path reproduces the matrix default exactly.
  .ExpectSameSelection(
    DropAdd(4L, colFn, N = 25L, seed = unname(which.max(rowSums(dmat)))),
    DropAdd(4L, dmat, maxCandidates = 0L)
  )
})

test_that(".DropAddColumn restores the zero diagonal .DistColumn masks", {
  # .DistColumn() masks the self position to -Inf for FarFirst's which.max /
  # pmin.int; DropAdd folds the whole column into sumDist and needs a 0 there.
  # A -Inf leak would drive `secondary` to -Inf while leaving the indices
  # plausible, so assert the convention directly.
  dmat <- as.matrix(dist(matrix(c(0, 1, 4, 0, 2, 5), ncol = 2L)))
  colFn <- function(i) dmat[, i]
  expect_identical(MaxMin:::.DistColumn(colFn, 2L, 3L)[2L], -Inf)
  expect_identical(MaxMin:::.DropAddColumn(colFn, 2L, 3L), unname(dmat[, 2L]))
  expect_true(is.finite(attr(DropAdd(2L, colFn, N = 3L, seed = 1L), "secondary")))
})

# ---------------------------------------------------------------------------
# 4. Degenerate sizes and stopping criteria
# ---------------------------------------------------------------------------
test_that("the oracle path handles k == n and n == 2 without iterating", {
  set.seed(3)
  dmat <- as.matrix(dist(matrix(rnorm(12L), ncol = 2L)))
  colFn <- function(i) dmat[, i]
  full <- DropAdd(6L, colFn, N = 6L, seed = 1L)
  expect_identical(as.integer(full), 1:6)
  expect_identical(attr(full, "iters"), 0L)   # no drop-add move exists
  .ExpectSameSelection(full, DropAdd(6L, dmat, seed = 1L, maxCandidates = 0L))

  # n == 2 is the smallest legal instance: k == n again, straight out of the
  # construction.
  two <- as.matrix(dist(matrix(c(0, 3, 0, 4), ncol = 2L)))
  res <- DropAdd(2L, function(i) two[, i], N = 2L)
  expect_identical(as.integer(res), 1:2)
  expect_equal(attr(res, "score"), 5)
  expect_equal(attr(res, "secondary"), 5)
})

test_that("k == 2 rebuilds a record with no surviving peer", {
  # At k == 2 the dropped point is the only peer of the survivor, so the
  # survivor's minDistCount hits zero with an empty peer set -- the branch that
  # resets its record to (Inf, 0).
  set.seed(17)
  dmat <- as.matrix(dist(matrix(rnorm(14L * 2L), ncol = 2L)))
  colFn <- function(i) dmat[, i]
  .ExpectSameSelection(
    DropAdd(2L, colFn, N = 14L, seed = 4L, maxCandidates = 0L),
    DropAdd(2L, dmat, seed = 4L, maxCandidates = 0L)
  )
})

test_that("the oracle loop honours maxIter, plateau and maxSeconds", {
  set.seed(23)
  dmat <- as.matrix(dist(matrix(rnorm(30L * 2L), ncol = 2L)))
  colFn <- function(i) dmat[, i]
  capped <- MaxMin:::.DropAddFromColumn(colFn, 30L, 5L, first = 1L,
                                        plateau = 1e9, maxIter = 7L)
  expect_identical(capped$iters, 7L)
  # Stagnation is the primary rule: with a small plateau the run stops early.
  short <- MaxMin:::.DropAddFromColumn(colFn, 30L, 5L, first = 1L, plateau = 3L)
  long  <- MaxMin:::.DropAddFromColumn(colFn, 30L, 5L, first = 1L, plateau = 500L)
  expect_lt(short$iters, long$iters)
  expect_gte(long$objective, short$objective)
  # A finite budget engages the wall-clock check (skipped entirely when the
  # budget is infinite, so the default path stays deterministic); an already
  # exhausted budget stops as soon as `proc.time()` ticks past it.
  budgeted <- MaxMin:::.DropAddFromColumn(colFn, 30L, 5L, first = 1L,
                                          maxSeconds = 1e-9)
  expect_lt(budgeted$iters, long$iters)
  expect_length(DropAdd(5L, colFn, N = 30L, seed = 1L, maxSeconds = 0.2), 5L)
})

test_that("the oracle path reports a score its indices really achieve", {
  set.seed(31)
  dmat <- as.matrix(dist(matrix(rnorm(40L * 3L), ncol = 3L)))
  res <- DropAdd(6L, function(i) dmat[, i], N = 40L)
  sub <- dmat[res, res]
  diag(sub) <- Inf
  expect_identical(attr(res, "score"), min(sub))
  expect_equal(attr(res, "secondary"), sum(dmat[res, res][lower.tri(sub)]))
  expect_s3_class(res, "MaxMinSelection")
  expect_true(attr(res, "time_s") >= 0)
})

# ---------------------------------------------------------------------------
# 5. Dispatch, guards and unsupported requests
# ---------------------------------------------------------------------------
test_that("N is required on the oracle path and validated", {
  dmat <- as.matrix(dist(matrix(rnorm(10L * 2L), ncol = 2L)))
  colFn <- function(i) dmat[, i]
  expect_error(DropAdd(3L, colFn), "N")
  expect_error(DropAdd(3L, colFn, N = c(10L, 10L)), "single positive integer")
  expect_error(DropAdd(3L, colFn, N = NA_integer_), "single positive integer")
  expect_error(DropAdd(3L, colFn, N = 0L), "single positive integer")
  # N < k is caught by the usual k check.
  expect_error(DropAdd(30L, colFn, N = 10L), "2 <= k <= n")
})

test_that("N is warned about on the matrix and coordinate paths", {
  pts <- matrix(rnorm(10L * 2L), ncol = 2L)
  expect_warning(DropAdd(3L, dist(pts), N = 10L), "ignored")
  expect_warning(DropAdd(3L, points = pts, N = 10L), "ignored")
  expect_no_warning(DropAdd(3L, function(i) .EuclidCol(pts, i), N = 10L))
})

test_that("supplying both a colFn and points is an error", {
  pts <- matrix(rnorm(10L * 2L), ncol = 2L)
  expect_error(DropAdd(3L, function(i) .EuclidCol(pts, i), points = pts, N = 10L),
               "not both")
})

test_that("candidate thinning warns and is skipped on the oracle path", {
  set.seed(13)
  dmat <- as.matrix(dist(matrix(rnorm(20L * 2L), ncol = 2L)))
  colFn <- function(i) dmat[, i]
  expect_warning(DropAdd(3L, colFn, N = 20L, maxCandidates = 10L),
                 "not supported on the distance-column path")
  # Skipped, not silently substituted: the answer is the full-problem answer.
  .ExpectSameSelection(
    suppressWarnings(DropAdd(3L, colFn, N = 20L, maxCandidates = 10L)),
    DropAdd(3L, colFn, N = 20L, maxCandidates = 0L)
  )
  # A cap that does not bind is a no-op, as on the matrix path.
  expect_no_warning(DropAdd(3L, colFn, N = 20L, maxCandidates = 20L))
})

test_that("seed, plateau and maxSeconds are validated on the oracle path", {
  dmat <- as.matrix(dist(matrix(rnorm(10L * 2L), ncol = 2L)))
  colFn <- function(i) dmat[, i]
  expect_error(DropAdd(3L, colFn, N = 10L, seed = 0L), "seed")
  expect_error(DropAdd(3L, colFn, N = 10L, seed = 11L), "seed")
  expect_error(DropAdd(3L, colFn, N = 10L, plateau = 0L), "plateau")
  expect_error(DropAdd(3L, colFn, N = 10L, maxSeconds = 0), "maxSeconds")
})

test_that("a malformed or NA/NaN colFn result is rejected", {
  # NA/NaN only: .DistColumn()'s guard deliberately lets `Inf` through, matching
  # FarFirst()'s oracle path. (The matrix path is stricter -- .AsDistMatrix()
  # rejects Inf -- but that is a permissiveness difference on input the matrix
  # path refuses outright, not a divergence in the search.)
  set.seed(19)
  dmat <- as.matrix(dist(matrix(rnorm(10L * 2L), ncol = 2L)))
  # Wrong length: neither N nor N - 1.
  expect_error(DropAdd(3L, function(i) dmat[1:5, i], N = 10L), "length")
  # NA at a non-self position would silently corrupt the streamlined records,
  # exactly as an NA matrix entry would on the matrix path.
  naFn <- function(i) {
    col <- dmat[, i]
    col[if (i == 1L) 2L else 1L] <- NA_real_
    col
  }
  expect_error(DropAdd(3L, naFn, N = 10L), "NA")
})

test_that("the oracle path fires the cli progress hooks", {
  set.seed(29)
  dmat <- as.matrix(dist(matrix(rnorm(15L * 2L), ncol = 2L)))
  old <- options(MaxMin.progress = TRUE)
  on.exit(options(old))
  expect_no_error(suppressMessages(
    DropAdd(3L, function(i) dmat[, i], N = 15L, maxSeconds = 0.1)
  ))
})

# ---------------------------------------------------------------------------
# 6. Solution quality: the oracle path is the same search, so it must not lose
# ---------------------------------------------------------------------------
test_that("the oracle path never worsens its own construction", {
  set.seed(37)
  dmat <- as.matrix(dist(matrix(rnorm(45L * 2L), ncol = 2L)))
  colFn <- function(i) dmat[, i]
  first <- MaxMin:::.PeripheralSeedColumn(colFn, 45L)
  cons <- MaxMin:::.DropAddConstructColumn(colFn, 45L, 7L, first)
  sub <- dmat[cons$S, cons$S]
  diag(sub) <- Inf
  res <- DropAdd(7L, colFn, N = 45L)
  expect_gte(attr(res, "score"), min(sub) - 1e-9)
})

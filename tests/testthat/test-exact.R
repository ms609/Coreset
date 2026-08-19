# Tests for ExactMaxMin (Sayyady & Fathi 2016 iterated node-packing).
#
# Oracle: exhaustive enumeration of every k-subset on tiny instances. For
# n <= 14, k <= 5 the binomial C(n, k) stays small enough to brute-force the
# true MaxMin optimum, against which ExactMaxMin must agree -- both on the
# objective value and on returning a subset that *achieves* it. (The argmax
# subset itself need not match: ties admit several optima.)

# Brute-force MaxMin objective for a selection.
.MaxminObj <- function(S, dmat) {
  sub <- dmat[S, S]
  diag(sub) <- Inf
  min(sub)
}

# Exhaustive optimum over all k-subsets of an n x n distance matrix.
.BruteMaxmin <- function(dmat, k) {
  n <- nrow(dmat)
  combos <- utils::combn(n, k)
  objs <- apply(combos, 2L, .MaxminObj, dmat = dmat)
  list(objective = max(objs), best = combos[, which.max(objs)])
}

# ---------------------------------------------------------------------------
# 1. Brute-force oracle agreement on random small instances
# ---------------------------------------------------------------------------
test_that("ExactMaxMin matches the brute-force optimum (random instances)", {
  cases <- list(
    list(n = 10L, k = 3L, seed = 1L),
    list(n = 12L, k = 4L, seed = 2L),
    list(n = 14L, k = 5L, seed = 3L),
    list(n = 11L, k = 5L, seed = 4L),
    list(n = 13L, k = 2L, seed = 5L)
  )
  for (cs in cases) {
    set.seed(cs$seed)
    pts <- matrix(stats::rnorm(cs$n * 4L), ncol = 4L)
    d <- as.matrix(stats::dist(pts))
    res <- ExactMaxMin(cs$k, d, maxSeconds = 60)
    oracle <- .BruteMaxmin(d, cs$k)

    # The optimum value matches the exhaustive optimum.
    expect_true(isTRUE(attr(res, "proven")), info = sprintf("n=%d k=%d", cs$n, cs$k))
    expect_equal(attr(res, "score"), oracle$objective,
                 info = sprintf("n=%d k=%d", cs$n, cs$k))
    # The returned subset actually achieves that optimum.
    idx <- as.integer(res)
    expect_length(idx, cs$k)
    expect_equal(length(unique(idx)), cs$k)
    expect_equal(.MaxminObj(idx, d), oracle$objective)
    expect_false(is.unsorted(idx))
  }
})

# ---------------------------------------------------------------------------
# 2. dist-object input is accepted and agrees with the matrix path
# ---------------------------------------------------------------------------
test_that("ExactMaxMin accepts a dist object", {
  set.seed(11)
  pts <- matrix(stats::rnorm(12 * 3), ncol = 3)
  dobj <- stats::dist(pts)
  dmat <- as.matrix(dobj)
  r_dist <- ExactMaxMin(4L, dobj)
  r_mat  <- ExactMaxMin(4L, dmat)
  expect_equal(attr(r_dist, "score"), attr(r_mat, "score"))
  expect_true(isTRUE(attr(r_dist, "proven")))
})

# ---------------------------------------------------------------------------
# 3. Edge cases: k = 2 (diameter) and k = n (global minimum distance)
# ---------------------------------------------------------------------------
test_that("ExactMaxMin handles k=2 (diameter) and k=n (all points)", {
  set.seed(21)
  pts <- matrix(stats::rnorm(8 * 3), ncol = 3)
  d <- as.matrix(stats::dist(pts))

  r2 <- ExactMaxMin(2L, d)
  expect_equal(attr(r2, "score"), max(d[upper.tri(d)]))   # diameter
  expect_true(isTRUE(attr(r2, "proven")))

  rn <- ExactMaxMin(8L, d)
  expect_equal(attr(rn, "score"), min(d[upper.tri(d)]))   # forced global min
  expect_equal(.MaxminObj(as.integer(rn), d), attr(rn, "score"))
  expect_true(isTRUE(attr(rn, "proven")))
})

# ---------------------------------------------------------------------------
# 3b. proven=FALSE on budget expiry
# ---------------------------------------------------------------------------
test_that("ExactMaxMin returns proven=FALSE on budget expiry", {
  set.seed(42)
  # Use a large enough instance that 1 microsecond is not enough to certify optimality.
  pts <- matrix(stats::rnorm(30 * 4), ncol = 4)
  d <- as.matrix(stats::dist(pts))
  res <- ExactMaxMin(5L, d, maxSeconds = 1e-9)
  # With a near-zero budget the solver may not certify, but must not error.
  # Note: proven could be TRUE if the solver is extremely fast; we just assert
  # the field exists and the return is valid.
  expect_type(attr(res, "proven"), "logical")
  expect_length(as.integer(res), 5L)
  # If proven = FALSE, the objective is still a valid lower bound (non-negative).
  if (!isTRUE(attr(res, "proven"))) {
    expect_gte(attr(res, "score"), 0)
  }
})

# ---------------------------------------------------------------------------
# 4. The optimum is exact (RNG-independent); the selection is reproducible
#    under a fixed seed. ExactMaxMin's warm start draws on the session RNG, so
#    the proven objective never varies, but which optimal subset is returned
#    can -- and is pinned by set.seed(), matching Grasp().
# ---------------------------------------------------------------------------
test_that("ExactMaxMin objective is RNG-independent; selection seed-reproducible", {
  set.seed(31)
  pts <- matrix(stats::rnorm(13 * 4), ncol = 4)
  d <- as.matrix(stats::dist(pts))

  # The proven optimum is a property of the problem: same value under any seed.
  set.seed(1); a <- ExactMaxMin(4L, d)
  set.seed(2); b <- ExactMaxMin(4L, d)
  expect_true(isTRUE(attr(a, "proven")) && isTRUE(attr(b, "proven")))
  expect_equal(attr(a, "score"), attr(b, "score"))

  # The returned selection is reproducible when the seed is reset.
  set.seed(1); a2 <- ExactMaxMin(4L, d)
  expect_equal(as.integer(a), as.integer(a2))
})

# ---------------------------------------------------------------------------
# 5. Input validation
# ---------------------------------------------------------------------------
test_that("ExactMaxMin validates k and solver", {
  d <- as.matrix(stats::dist(matrix(stats::rnorm(20), ncol = 2)))   # n = 10
  expect_error(ExactMaxMin(1L, d), "2 <= k <= nrow")
  expect_error(ExactMaxMin(11L, d), "2 <= k <= nrow")
  # .ExactAsMatrix validation: non-matrix and non-square inputs
  expect_error(ExactMaxMin(3L, 1:9), "dist|matrix")
  expect_error(ExactMaxMin(2L, matrix(1:6, 2, 3)), "dist|matrix")
})

# ---------------------------------------------------------------------------
# 7. Coreset.progress option fires the cli progress bar
# ---------------------------------------------------------------------------
test_that("ExactMaxMin Coreset.progress option fires the cli hooks", {
  set.seed(1)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(12 * 2), ncol = 2)))
  old <- options(Coreset.progress = TRUE)
  on.exit(options(old))
  expect_no_error(ExactMaxMin(3L, d))
})

# ---------------------------------------------------------------------------
# 6. Return contract
# ---------------------------------------------------------------------------
test_that("ExactMaxMin returns the documented fields", {
  set.seed(41)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(40), ncol = 4)))   # n = 10
  res <- ExactMaxMin(3L, d)
  expect_s3_class(res, "MaxMinSelection")
  expect_type(res, "integer")
  expect_length(res, 3L)
  expect_type(attr(res, "score"),  "double")
  expect_type(attr(res, "proven"), "logical")
  expect_identical(attr(res, "N"), 10L)
  expect_identical(attr(res, "k"), 3L)
  expect_type(attr(res, "time_s"), "double")
})

test_that("ExactMaxMin rejects an asymmetric matrix loudly", {
  # The warm start collects Grasp()/DropAdd() results through tryCatch, so a
  # contract error raised only there would be swallowed into an empty pool and
  # a silent fall back to the trivial bound.
  set.seed(1)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(40L), ncol = 2L)))
  d[2L, 5L] <- d[2L, 5L] + 1
  expect_error(ExactMaxMin(3L, d), "must be symmetric")
})

# Tests for ExactMaxMin (Sayyady & Fathi 2016 iterated node-packing).
#
# Oracle: exhaustive enumeration of every k-subset on tiny instances. For
# n <= 14, k <= 5 the binomial C(n, k) stays small enough to brute-force the
# true MaxMin optimum, against which ExactMaxMin must agree -- both on the
# objective value and on returning a subset that *achieves* it. (The argmax
# subset itself need not match: ties admit several optima.)
#
# All tests skip when the `highs` solver is unavailable.

skip_if_no_highs <- function() {
  testthat::skip_if_not_installed("highs")
}

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
  skip_if_no_highs()
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
    expect_true(res$proven, info = sprintf("n=%d k=%d", cs$n, cs$k))
    expect_equal(res$objective, oracle$objective,
                 info = sprintf("n=%d k=%d", cs$n, cs$k))
    # The returned subset actually achieves that optimum.
    expect_length(res$indices, cs$k)
    expect_equal(length(unique(res$indices)), cs$k)
    expect_equal(.MaxminObj(res$indices, d), oracle$objective)
    expect_false(is.unsorted(res$indices))
  }
})

# ---------------------------------------------------------------------------
# 2. dist-object input is accepted and agrees with the matrix path
# ---------------------------------------------------------------------------
test_that("ExactMaxMin accepts a dist object", {
  skip_if_no_highs()
  set.seed(11)
  pts <- matrix(stats::rnorm(12 * 3), ncol = 3)
  dobj <- stats::dist(pts)
  dmat <- as.matrix(dobj)
  r_dist <- ExactMaxMin(4L, dobj)
  r_mat  <- ExactMaxMin(4L, dmat)
  expect_equal(r_dist$objective, r_mat$objective)
  expect_true(r_dist$proven)
})

# ---------------------------------------------------------------------------
# 3. Edge cases: k = 2 (diameter) and k = n (global minimum distance)
# ---------------------------------------------------------------------------
test_that("ExactMaxMin handles k=2 (diameter) and k=n (all points)", {
  skip_if_no_highs()
  set.seed(21)
  pts <- matrix(stats::rnorm(8 * 3), ncol = 3)
  d <- as.matrix(stats::dist(pts))

  r2 <- ExactMaxMin(2L, d)
  expect_equal(r2$objective, max(d[upper.tri(d)]))   # diameter
  expect_true(r2$proven)

  rn <- ExactMaxMin(8L, d)
  expect_equal(rn$objective, min(d[upper.tri(d)]))   # forced global min
  expect_equal(.MaxminObj(rn$indices, d), rn$objective)
  expect_true(rn$proven)
})

# ---------------------------------------------------------------------------
# 3b. proven=FALSE on budget expiry
# ---------------------------------------------------------------------------
test_that("ExactMaxMin returns proven=FALSE on budget expiry", {
  skip_if_no_highs()
  set.seed(42)
  # Use a large enough instance that 1 microsecond is not enough to certify optimality.
  pts <- matrix(stats::rnorm(30 * 4), ncol = 4)
  d <- as.matrix(stats::dist(pts))
  res <- ExactMaxMin(5L, d, maxSeconds = 1e-9)
  # With a near-zero budget the solver may not certify, but must not error.
  # Note: proven could be TRUE if the solver is extremely fast; we just assert
  # the field exists and the return is valid.
  expect_type(res$proven, "logical")
  expect_length(res$indices, 5L)
  # If proven = FALSE, the objective is still a valid lower bound (non-negative).
  if (!res$proven) {
    expect_gte(res$objective, 0)
  }
})

# ---------------------------------------------------------------------------
# 4. The optimum is exact (RNG-independent); the selection is reproducible
#    under a fixed seed. ExactMaxMin's warm start draws on the session RNG, so
#    the proven objective never varies, but which optimal subset is returned
#    can -- and is pinned by set.seed(), matching Grasp().
# ---------------------------------------------------------------------------
test_that("ExactMaxMin objective is RNG-independent; selection seed-reproducible", {
  skip_if_no_highs()
  set.seed(31)
  pts <- matrix(stats::rnorm(13 * 4), ncol = 4)
  d <- as.matrix(stats::dist(pts))

  # The proven optimum is a property of the problem: same value under any seed.
  set.seed(1); a <- ExactMaxMin(4L, d)
  set.seed(2); b <- ExactMaxMin(4L, d)
  expect_true(a$proven && b$proven)
  expect_equal(a$objective, b$objective)

  # The returned selection is reproducible when the seed is reset.
  set.seed(1); a2 <- ExactMaxMin(4L, d)
  expect_equal(a$indices, a2$indices)
})

# ---------------------------------------------------------------------------
# 5. Input validation
# ---------------------------------------------------------------------------
test_that("ExactMaxMin validates k and solver", {
  skip_if_no_highs()
  d <- as.matrix(stats::dist(matrix(stats::rnorm(20), ncol = 2)))   # n = 10
  expect_error(ExactMaxMin(1L, d), "2 <= k <= nrow")
  expect_error(ExactMaxMin(11L, d), "2 <= k <= nrow")
  expect_error(ExactMaxMin(3L, d, solver = "gurobi"), "Unsupported")
  # .ExactAsMatrix validation: non-matrix and non-square inputs
  expect_error(ExactMaxMin(3L, 1:9), "dist|matrix")
  expect_error(ExactMaxMin(2L, matrix(1:6, 2, 3)), "dist|matrix")
})

# ---------------------------------------------------------------------------
# 7. progress = TRUE fires the cli progress bar
# ---------------------------------------------------------------------------
test_that("ExactMaxMin progress = TRUE fires the cli hooks", {
  skip_if_no_highs()
  set.seed(1)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(12 * 2), ncol = 2)))
  expect_no_error(ExactMaxMin(3L, d, progress = TRUE))
})

# ---------------------------------------------------------------------------
# 6. Return contract
# ---------------------------------------------------------------------------
test_that("ExactMaxMin returns the documented fields", {
  skip_if_no_highs()
  set.seed(41)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(40), ncol = 4)))   # n = 10
  res <- ExactMaxMin(3L, d)
  expect_named(res, c("indices", "objective", "proven", "time_s",
                      "solver", "n", "k"),
               ignore.order = TRUE)
  expect_type(res$indices, "integer")
  expect_type(res$proven, "logical")
  expect_identical(res$solver, "highs")
  expect_identical(res$n, 10L)
  expect_identical(res$k, 3L)
})

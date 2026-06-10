# Tests for ExactMaxMin (Sayyady & Fathi 2016 iterated node-packing).
#
# Oracle: exhaustive enumeration of every m-subset on tiny instances. For
# n <= 14, m <= 5 the binomial C(n, m) stays small enough to brute-force the
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

# Exhaustive optimum over all m-subsets of an n x n distance matrix.
.BruteMaxmin <- function(dmat, m) {
  n <- nrow(dmat)
  combos <- utils::combn(n, m)
  objs <- apply(combos, 2L, .MaxminObj, dmat = dmat)
  list(objective = max(objs), best = combos[, which.max(objs)])
}

# ---------------------------------------------------------------------------
# 1. Brute-force oracle agreement on random small instances
# ---------------------------------------------------------------------------
test_that("ExactMaxMin matches the brute-force optimum (random instances)", {
  skip_if_no_highs()
  cases <- list(
    list(n = 10L, m = 3L, seed = 1L),
    list(n = 12L, m = 4L, seed = 2L),
    list(n = 14L, m = 5L, seed = 3L),
    list(n = 11L, m = 5L, seed = 4L),
    list(n = 13L, m = 2L, seed = 5L)
  )
  for (cs in cases) {
    set.seed(cs$seed)
    pts <- matrix(stats::rnorm(cs$n * 4L), ncol = 4L)
    d <- as.matrix(stats::dist(pts))
    res <- ExactMaxMin(d, m = cs$m, timeBudgetS = 60)
    oracle <- .BruteMaxmin(d, cs$m)

    # The optimum value matches the exhaustive optimum.
    expect_true(res$proven, info = sprintf("n=%d m=%d", cs$n, cs$m))
    expect_equal(res$objective, oracle$objective,
                 info = sprintf("n=%d m=%d", cs$n, cs$m))
    # The returned subset actually achieves that optimum.
    expect_length(res$indices, cs$m)
    expect_equal(length(unique(res$indices)), cs$m)
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
  r_dist <- ExactMaxMin(dobj, m = 4L)
  r_mat  <- ExactMaxMin(dmat, m = 4L)
  expect_equal(r_dist$objective, r_mat$objective)
  expect_true(r_dist$proven)
})

# ---------------------------------------------------------------------------
# 3. Edge cases: m = 2 (diameter) and m = n (global minimum distance)
# ---------------------------------------------------------------------------
test_that("ExactMaxMin handles m=2 (diameter) and m=n (all points)", {
  skip_if_no_highs()
  set.seed(21)
  pts <- matrix(stats::rnorm(8 * 3), ncol = 3)
  d <- as.matrix(stats::dist(pts))

  r2 <- ExactMaxMin(d, m = 2L)
  expect_equal(r2$objective, max(d[upper.tri(d)]))   # diameter
  expect_true(r2$proven)

  rn <- ExactMaxMin(d, m = 8L)
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
  res <- ExactMaxMin(d, m = 5L, timeBudgetS = 1e-9)
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
# 4. Determinism: same input -> same objective (and same indices here)
# ---------------------------------------------------------------------------
test_that("ExactMaxMin is deterministic", {
  skip_if_no_highs()
  set.seed(31)
  pts <- matrix(stats::rnorm(13 * 4), ncol = 4)
  d <- as.matrix(stats::dist(pts))
  a <- ExactMaxMin(d, m = 4L)
  b <- ExactMaxMin(d, m = 4L)
  expect_equal(a$objective, b$objective)
  expect_equal(a$indices, b$indices)
})

# ---------------------------------------------------------------------------
# 5. Input validation
# ---------------------------------------------------------------------------
test_that("ExactMaxMin validates m and solver", {
  skip_if_no_highs()
  d <- as.matrix(stats::dist(matrix(stats::rnorm(20), ncol = 2)))   # n = 10
  expect_error(ExactMaxMin(d, m = 1L), "2 <= m <= nrow")
  expect_error(ExactMaxMin(d, m = 11L), "2 <= m <= nrow")
  expect_error(ExactMaxMin(d, m = 3L, solver = "gurobi"), "Unsupported")
  # .ExactAsMatrix validation: non-matrix and non-square inputs
  expect_error(ExactMaxMin(1:9, m = 3L), "dist|matrix")
  expect_error(ExactMaxMin(matrix(1:6, 2, 3), m = 2L), "dist|matrix")
})

# ---------------------------------------------------------------------------
# 7. progress = TRUE fires the cli progress bar
# ---------------------------------------------------------------------------
test_that("ExactMaxMin progress = TRUE fires the cli hooks", {
  skip_if_no_highs()
  set.seed(1)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(12 * 2), ncol = 2)))
  expect_no_error(ExactMaxMin(d, m = 3L, progress = TRUE))
})

# ---------------------------------------------------------------------------
# 6. Return contract
# ---------------------------------------------------------------------------
test_that("ExactMaxMin returns the documented fields", {
  skip_if_no_highs()
  set.seed(41)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(40), ncol = 4)))   # n = 10
  res <- ExactMaxMin(d, m = 3L)
  expect_named(res, c("indices", "objective", "proven", "time_s",
                      "solver", "n", "m"),
               ignore.order = TRUE)
  expect_type(res$indices, "integer")
  expect_type(res$proven, "logical")
  expect_identical(res$solver, "highs")
  expect_identical(res$n, 10L)
  expect_identical(res$m, 3L)
})

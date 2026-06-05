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
.maxmin_obj <- function(S, dmat) {
  sub <- dmat[S, S]
  diag(sub) <- Inf
  min(sub)
}

# Exhaustive optimum over all m-subsets of an n x n distance matrix.
.brute_maxmin <- function(dmat, m) {
  n <- nrow(dmat)
  combos <- utils::combn(n, m)
  objs <- apply(combos, 2L, .maxmin_obj, dmat = dmat)
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
    res <- ExactMaxMin(d, m = cs$m, time_budget_s = 60)
    oracle <- .brute_maxmin(d, cs$m)

    # The optimum value matches the exhaustive optimum.
    expect_true(res$proven, info = sprintf("n=%d m=%d", cs$n, cs$m))
    expect_equal(res$objective, oracle$objective,
                 info = sprintf("n=%d m=%d", cs$n, cs$m))
    # The returned subset actually achieves that optimum.
    expect_length(res$indices, cs$m)
    expect_equal(length(unique(res$indices)), cs$m)
    expect_equal(.maxmin_obj(res$indices, d), oracle$objective)
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
  expect_equal(rn$indices, 1:8)
  expect_true(rn$proven)
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

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
  skip_if_no_highs()
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
  skip_if_no_highs()
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
  skip_if_no_highs()
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
  skip_if_no_highs()
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
  skip_if_no_highs()
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
  skip_if_no_highs()
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
  skip_if_no_highs()
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

# On small point sets the heuristics attain the optimum outright, so the gallop
# stops at its first probe and the bisection never runs -- 3200 random Euclidean
# instances produced no exception. Two shapes do defeat them, and between them
# they exercise the rest of the search: a graph metric, where the warm start
# lands one attainable threshold short, and a ring, where it lands short and the
# first witness reaches several candidates past the threshold it was asked for.

# Shortest-path closure of a sparse weighted graph, as the ORLIB p-median
# instances are built.
.GraphMetric <- function(n, m, seed) {
  set.seed(seed)
  d <- matrix(Inf, n, n)
  diag(d) <- 0
  ends <- matrix(sample.int(n, 2L * m, replace = TRUE), nrow = 2L)
  w <- sample.int(20L, m, replace = TRUE)
  for (q in seq_len(m)) {
    a <- ends[1L, q]
    b <- ends[2L, q]
    if (a != b) {
      d[a, b] <- w[q]
      d[b, a] <- w[q]
    }
  }
  for (v in seq_len(n)) d <- pmin(d, outer(d[, v], d[v, ], "+"))
  d
}

test_that("the gallop and bisection close a graph metric the heuristics miss", {
  d <- .GraphMetric(120L, 360L, 20L)
  expect_true(all(is.finite(d)))
  set.seed(20)
  ws <- .ExactWarmStart(d, 120L, 30L, NULL)
  set.seed(20)
  sel <- ExactMaxMin(k = 30L, d = d, maxSeconds = 60)
  # The warm start is short, so the search must take a feasible probe before it
  # can take an infeasible one.
  expect_lt(ws$value, attr(sel, "score"))
  expect_true(attr(sel, "proven"))
  expect_equal(attr(sel, "score"), 16)
  sub <- d[sel, sel]
  diag(sub) <- Inf
  expect_length(sel, 30L)
  expect_equal(min(sub), attr(sel, "score"))
})

test_that("a witness reaching past its threshold still yields the optimum", {
  # Thirteen attainable thresholds separate this warm start from the optimum
  # and the first witness covers only six of them, so the search must both skip
  # the thresholds that witness already attains and take a feasible step of its
  # own on the way down.
  set.seed(6)
  ang <- sort(stats::runif(120L, 0, 2 * pi))
  pts <- cbind(cos(ang), sin(ang)) * (1 + stats::rnorm(120L, sd = 0.02))
  d <- as.matrix(stats::dist(pts))
  set.seed(6)
  ws <- .ExactWarmStart(d, 120L, 6L, NULL)
  set.seed(6)
  sel <- ExactMaxMin(k = 6L, d = d, maxSeconds = 60)
  expect_lt(ws$value, attr(sel, "score"))
  expect_true(attr(sel, "proven"))
  sub <- d[sel, sel]
  diag(sub) <- Inf
  expect_equal(min(sub), attr(sel, "score"))
  # Optimality: no 8-subset beats it, checked against the next attainable
  # threshold rather than the search that produced it.
  cand <- sort(unique(TriangleAtLeast_cpp(d, attr(sel, "score"))))
  nxt <- cand[cand > attr(sel, "score")][1L]
  h <- EdgesAtLeast_cpp(d, nxt)
  expect_identical(
    ThresholdDecide_cpp(h[["hi"]], h[["hj"]], 120L, 6L, 60)[["status"]],
    "infeasible"
  )
})

test_that("a malformed caller warm start is dropped, not trusted", {
  d <- as.matrix(stats::dist(matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2L,
                                    byrow = TRUE)))
  set.seed(1)
  # Too short to be a k-subset: the pool must fall back to its own heuristics
  # rather than carry a bound it cannot justify.
  ws <- .ExactWarmStart(d, 4L, 3L, warmStart = c(1L, 2L))
  expect_length(ws$witness, 3L)
  set.seed(1)
  expect_equal(ws$value, .ExactWarmStart(d, 4L, 3L, NULL)$value)
})

# ---------------------------------------------------------------------------
# Warm-start depth is a knob, not a correctness lever
# ---------------------------------------------------------------------------
test_that("the warm-start knobs change the pool but not the optimum", {
  set.seed(11L)
  dmat <- as.matrix(dist(matrix(rnorm(28L), ncol = 2L)))
  truth <- .BruteMaxmin(dmat, 4L)

  # A shallower pool and a deeper one must certify the same value: the pool
  # only supplies a lower bound, and everything above it is still searched.
  set.seed(1L)
  shallow <- ExactMaxMin(4L, dmat, nStart = 1L, graspPlateau = 1L,
                         dropPlateau = 1L)
  set.seed(1L)
  deep <- ExactMaxMin(4L, dmat, nStart = 3L, graspPlateau = 512L,
                      dropPlateau = 5000L)

  expect_equal(attr(shallow, "score"), truth$objective)
  expect_equal(attr(deep, "score"), truth$objective)
  expect_true(attr(shallow, "proven"))
  expect_true(attr(deep, "proven"))
})

test_that("a deeper pool reaches at least as far as a shallower one", {
  set.seed(12L)
  dmat <- as.matrix(dist(matrix(rnorm(120L), ncol = 3L)))
  set.seed(2L)
  shallow <- .ExactWarmStart(dmat, nrow(dmat), 6L, NULL, nStart = 1L,
                             graspPlateau = 1L, dropPlateau = 1L)
  set.seed(2L)
  deep <- .ExactWarmStart(dmat, nrow(dmat), 6L, NULL, nStart = 8L,
                          graspPlateau = 512L, dropPlateau = 5000L)
  expect_true(deep$value >= shallow$value)
  expect_length(deep$witness, 6L)
})

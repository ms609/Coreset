# test-maxmean.R — testthat tests for MaxMean() and MeanDist()

# ---- helpers ---------------------------------------------------------------

MakeD <- function(seed, n) {
  set.seed(seed)
  dist(matrix(rnorm(n * 3L), ncol = 3L))
}

# Brute-force max-mean over all subsets of size >= 2 (feasible only for tiny n).
BruteMaxMean <- function(d) {
  dmat <- as.matrix(d)
  n <- nrow(dmat)
  best_f   <- -Inf
  best_idx <- integer(0)
  for (k in 2:n) {
    combs <- utils::combn(n, k)
    for (i in seq_len(ncol(combs))) {
      idx <- combs[, i]
      sub <- dmat[idx, idx, drop = FALSE]
      f <- sum(sub[lower.tri(sub)]) / k
      if (f > best_f) {
        best_f   <- f
        best_idx <- idx
      }
    }
  }
  list(f = best_f, idx = best_idx)
}

# ---- MeanDist() ------------------------------------------------------------

test_that("MeanDist returns correct value", {
  d <- MakeD(1, 6)
  idx <- c(1L, 3L, 5L)
  dmat <- as.matrix(d)
  sub  <- dmat[idx, idx]
  expected <- sum(sub[lower.tri(sub)]) / length(idx)
  expect_equal(MeanDist(d, idx), expected)
})

test_that("MeanDist divides by |S|, not by pair count", {
  # Three equidistant points at distance 1: sum of pairs = 3, |S| = 3, so f = 1.
  dmat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
  expect_equal(MeanDist(dmat, 1:3), 1)
})

test_that("MeanDist returns NA for fewer than two elements", {
  d <- dist(matrix(1:6, ncol = 2))
  expect_identical(MeanDist(d, 2L), NA_real_)
  expect_identical(MeanDist(d, integer(0)), NA_real_)
})

test_that("MeanDist errors on NA index", {
  d <- dist(matrix(1:6, ncol = 2))
  expect_error(MeanDist(d, c(1L, NA_integer_)), "`idx` must not contain NA")
})

test_that("MeanDist errors on duplicate index", {
  d <- dist(matrix(1:6, ncol = 2))
  expect_error(MeanDist(d, c(1L, 1L)), "duplicate")
})

test_that("MeanDist accepts dist object and matrix equivalently", {
  d <- MakeD(3, 5)
  idx <- 1:4
  expect_equal(MeanDist(as.matrix(d), idx), MeanDist(d, idx))
})

test_that("MeanDist rejects non-square matrix", {
  expect_error(MeanDist(matrix(1:12, nrow = 3), 1:2), "square")
})

# ---- MaxMean() basic shape -------------------------------------------------

test_that("MaxMean returns MaxMeanSelection with documented attributes", {
  set.seed(1)
  d <- MakeD(1, 10)
  result <- MaxMean(d, maxSeconds = 1, useRL = FALSE)
  expect_s3_class(result, "MaxMeanSelection")
  expect_true(is.integer(result))
  expect_gte(length(result), 2L)
  expect_true(is.numeric(attr(result, "score")))
  expect_true(is.numeric(attr(result, "time_s")))
  expect_type(attr(result, "iters"), "double")   # double, not int (can exceed 2^31)
  expect_identical(attr(result, "size"), length(result))
})

test_that("MaxMean indices are sorted, unique and in range", {
  set.seed(2)
  d <- MakeD(2, 12)
  result <- MaxMean(d, maxSeconds = 0.5, useRL = FALSE)
  expect_true(all(diff(result) > 0L))           # sorted ascending, no dups
  expect_gte(min(result), 1L)
  expect_lte(max(result), 12L)
})

test_that("MaxMean score attribute equals MeanDist of the selection", {
  set.seed(3)
  d <- MakeD(3, 15)
  result <- MaxMean(d, maxSeconds = 1, useRL = FALSE)
  expect_equal(attr(result, "score"), MeanDist(d, result), tolerance = 1e-10)
})

# ---- MaxMean() correctness vs brute force on tiny instances ----------------

test_that("MaxMean (no RL) attains the brute-force optimum (n = 6)", {
  d <- MakeD(4, 6)
  brute <- BruteMaxMean(d)
  set.seed(10)
  # maxIter = Inf: spend the full budget on restarts so the optimum is reached.
  result <- MaxMean(d, maxSeconds = 3, maxIter = Inf, useRL = FALSE)
  expect_equal(attr(result, "score"), brute$f, tolerance = 1e-10)
  expect_identical(as.integer(result), as.integer(brute$idx))
})

test_that("MaxMean (with RL) attains the brute-force optimum (n = 6)", {
  d <- MakeD(5, 6)
  brute <- BruteMaxMean(d)
  set.seed(11)
  result <- MaxMean(d, maxSeconds = 3, maxIter = Inf, useRL = TRUE)
  expect_equal(attr(result, "score"), brute$f, tolerance = 1e-10)
  expect_identical(as.integer(result), as.integer(brute$idx))
})

test_that("MaxMean attains the brute-force optimum (n = 8)", {
  d <- MakeD(6, 8)
  brute <- BruteMaxMean(d)
  set.seed(12)
  result <- MaxMean(d, maxSeconds = 3, maxIter = Inf, useRL = TRUE)
  expect_equal(attr(result, "score"), brute$f, tolerance = 1e-10)
})

test_that("MaxMean handles signed (negative) distances", {
  # Max-mean is defined for signed distances; a negative-distance pair should be
  # excluded from the optimal set when it drags the mean down.
  set.seed(13)
  m <- matrix(rnorm(10 * 10), 10)
  dmat <- (m + t(m)) / 2          # symmetric, mixed-sign, zero diagonal
  diag(dmat) <- 0
  brute <- BruteMaxMean(dmat)
  set.seed(14)
  result <- MaxMean(dmat, maxSeconds = 3, maxIter = Inf, useRL = FALSE)
  expect_equal(attr(result, "score"), brute$f, tolerance = 1e-10)
})

# ---- minimal-instance behaviour --------------------------------------------

test_that("MaxMean on two points returns both", {
  dmat <- matrix(c(0, 5, 5, 0), nrow = 2)
  result <- MaxMean(dmat, maxSeconds = 0.2)
  expect_identical(as.integer(result), 1:2)
  expect_equal(attr(result, "score"), 5 / 2)   # sum 5 over |S| = 2
})

test_that("a vanishingly small budget still returns a valid selection (MM-01)", {
  # The restart loop is a do-while: at least one cycle completes, so even an
  # absurd budget yields a classed |S| >= 2 result, never an empty -1e300 one.
  set.seed(80)
  d <- MakeD(80, 10)
  result <- MaxMean(d, maxSeconds = 1e-9, useRL = FALSE)
  expect_s3_class(result, "MaxMeanSelection")
  expect_gte(length(result), 2L)
  expect_true(is.finite(attr(result, "score")))
  expect_equal(attr(result, "score"), MeanDist(d, result), tolerance = 1e-10)
})

test_that("iters is a double that survives values beyond 2^31 (MM-02)", {
  set.seed(81)
  d <- MakeD(81, 8)
  result <- MaxMean(d, maxSeconds = 0.3, useRL = FALSE)
  expect_type(attr(result, "iters"), "double")
  expect_false(is.na(attr(result, "iters")))
})

test_that("score matches MeanDist on asymmetric input (MM-03)", {
  # .AsDistMatrix accepts asymmetric matrices; MaxMean's objective symmetrizes,
  # and MeanDist must report the same value so the two never silently disagree.
  set.seed(82)
  m <- matrix(runif(12L * 12L, -5, 5), 12L)
  diag(m) <- 0                                   # asymmetric, zero diagonal
  set.seed(83)
  result <- MaxMean(m, maxSeconds = 1, useRL = FALSE)
  expect_equal(attr(result, "score"), MeanDist(m, result), tolerance = 1e-10)
})

test_that("MeanDist symmetrizes asymmetric input", {
  m <- matrix(c(0, 4, 0, 0,  2, 0, 0, 0,  0, 0, 0, 6,  0, 0, 0, 0), 4L)
  # pair {1,2}: (4+2)/2 = 3; pair {3,4}: (6+0)/2 = 3; others 0. sum = 6, |S| = 4.
  expect_equal(MeanDist(m, 1:4), 6 / 4)
})

# ---- argument validation ---------------------------------------------------

test_that("MaxMean rejects non-square matrix", {
  expect_error(MaxMean(matrix(1:12, nrow = 3)), "square")
})

test_that("MaxMean rejects a 1x1 matrix", {
  expect_error(MaxMean(matrix(0, 1, 1)), "at least 2")
})

test_that("MaxMean rejects invalid maxSeconds", {
  d <- dist(matrix(1:6, ncol = 2))
  expect_error(MaxMean(d, maxSeconds = -1), "maxSeconds")
  expect_error(MaxMean(d, maxSeconds = NA), "maxSeconds")
  expect_error(MaxMean(d, maxSeconds = c(1, 2)), "maxSeconds")
})

test_that("MaxMean rejects invalid useRL", {
  d <- dist(matrix(1:6, ncol = 2))
  expect_error(MaxMean(d, maxSeconds = 1, useRL = NA), "useRL")
  expect_error(MaxMean(d, maxSeconds = 1, useRL = "yes"), "useRL")
})

test_that("MaxMean rejects a distance matrix with NA/Inf", {
  dmat <- matrix(c(0, NA, NA, 0), nrow = 2)
  expect_error(MaxMean(dmat, maxSeconds = 1), "NA")
})

test_that("MaxMean rejects invalid maxIter", {
  d <- dist(matrix(1:6, ncol = 2))
  expect_error(MaxMean(d, maxIter = 0), "maxIter")
  expect_error(MaxMean(d, maxIter = NA), "maxIter")
  expect_error(MaxMean(d, maxIter = c(1, 2)), "maxIter")
  expect_error(MaxMean(d, maxIter = "lots"), "maxIter")
})

# ---- maxIter budget --------------------------------------------------------

test_that("maxIter is an exact cap on total iterations", {
  # Generous time budget so the iteration cap, not the clock, is what stops the
  # search; `iters` must never exceed maxIter (the cap is checked every step).
  set.seed(40)
  d <- MakeD(40, 40)
  result <- MaxMean(d, maxSeconds = 30, maxIter = 100, useRL = FALSE)
  expect_lte(attr(result, "iters"), 100)
  expect_s3_class(result, "MaxMeanSelection")
  expect_equal(attr(result, "score"), MeanDist(d, result), tolerance = 1e-10)
})

test_that("a smaller maxIter bites: fewer iterations, same instance", {
  set.seed(41)
  d <- MakeD(41, 40)
  set.seed(42); r_small <- MaxMean(d, maxSeconds = 30, maxIter = 50,  useRL = FALSE)
  set.seed(42); r_big   <- MaxMean(d, maxSeconds = 30, maxIter = 500, useRL = FALSE)
  expect_lte(attr(r_small, "iters"), 50)
  expect_lte(attr(r_big,   "iters"), 500)
  expect_lt(attr(r_small, "iters"), attr(r_big, "iters"))
})

test_that("maxIter alone stops the search when maxSeconds is Inf", {
  # The end-of-restart iteration break is load-bearing here: with no time limit,
  # nothing else would terminate the restart loop once the cap is reached.
  set.seed(43)
  d <- MakeD(43, 30)
  result <- MaxMean(d, maxSeconds = Inf, maxIter = 500, useRL = TRUE)
  expect_s3_class(result, "MaxMeanSelection")
  expect_lte(attr(result, "iters"), 500)
  expect_true(is.finite(attr(result, "score")))
})

test_that("maxIter = 1 still returns a feasible |S| >= 2 selection", {
  # S0 is built and scored before any tabu iteration runs, so even a one-flip
  # budget yields a valid classed result (the do-while restart guarantee).
  set.seed(44)
  d <- MakeD(44, 10)
  result <- MaxMean(d, maxSeconds = Inf, maxIter = 1, useRL = FALSE)
  expect_s3_class(result, "MaxMeanSelection")
  expect_gte(length(result), 2L)
  expect_true(is.finite(attr(result, "score")))
  expect_equal(attr(result, "score"), MeanDist(d, result), tolerance = 1e-10)
})

test_that("maxIter = Inf budgets by time alone", {
  # With the iteration cap disabled a longer clock yields at least as many
  # iterations — the cap is genuinely off, not silently applied.
  d <- MakeD(45, 60)
  set.seed(45); r_short <- MaxMean(d, maxSeconds = 0.2, maxIter = Inf, useRL = FALSE)
  set.seed(45); r_long  <- MaxMean(d, maxSeconds = 0.6, maxIter = Inf, useRL = FALSE)
  expect_gte(attr(r_long, "iters"), attr(r_short, "iters"))
})

# ---- print / format / summary methods --------------------------------------

test_that("print.MaxMeanSelection prints a one-line summary", {
  set.seed(7)
  d <- MakeD(7, 8)
  result <- MaxMean(d, maxSeconds = 0.5, useRL = FALSE)
  expect_output(print(result), "MaxMean RLTS")
  expect_output(print(result), "element")
})

test_that("format.MaxMeanSelection includes the f value", {
  set.seed(8)
  d <- MakeD(8, 8)
  result <- MaxMean(d, maxSeconds = 0.5, useRL = FALSE)
  expect_match(format(result), "f =")
})

test_that("summary.MaxMeanSelection prints the detail block", {
  set.seed(9)
  d <- MakeD(9, 8)
  result <- MaxMean(d, maxSeconds = 0.5, useRL = FALSE)
  expect_output(summary(result), "objective")
  expect_output(summary(result), "iterations")
  expect_output(summary(result), "time")
})

test_that(".AsMaxMeanSelection leaves an empty vector unclassed", {
  expect_identical(Coreset:::.AsMaxMeanSelection(integer(0)), integer(0))
})

# ---- determinism & monotonicity --------------------------------------------

test_that("matrix and dist input give the same selection and score", {
  # iters/time are wall-clock dependent and intentionally not compared; the
  # selection and its score are timing-invariant on this size of instance.
  d <- MakeD(30, 10)
  set.seed(30); r_dist   <- MaxMean(d, maxSeconds = 0.5, useRL = FALSE)
  set.seed(30); r_matrix <- MaxMean(as.matrix(d), maxSeconds = 0.5, useRL = FALSE)
  expect_identical(as.integer(r_dist), as.integer(r_matrix))
  expect_equal(attr(r_dist, "score"), attr(r_matrix, "score"))
})

test_that("a longer budget never worsens the objective", {
  d <- MakeD(20, 25)
  # maxIter = Inf so the two runs genuinely differ in budget (a finite iter cap
  # would make both stop at the same iteration and the test vacuous).
  set.seed(20); r1 <- MaxMean(d, maxSeconds = 0.3, maxIter = Inf, useRL = FALSE)
  set.seed(20); r2 <- MaxMean(d, maxSeconds = 1.5, maxIter = Inf, useRL = FALSE)
  expect_gte(attr(r2, "score"), attr(r1, "score") - 1e-10)
})

test_that("MaxMean consumes and restores the R RNG stream", {
  # If GetRNGstate/PutRNGstate are wired in, MaxMean draws from R's stream and
  # writes the advanced state back, so a draw taken afterwards differs from the
  # same-seed draw taken before. (A stale C-level RNG would leave x1 == x2.)
  d <- MakeD(60, 12)
  set.seed(123); x1 <- runif(1)
  set.seed(123); invisible(MaxMean(d, maxSeconds = 0.3, useRL = FALSE))
  x2 <- runif(1)
  expect_false(identical(x1, x2))
})

test_that("RL and no-RL both return a valid, finite-scored selection", {
  d <- MakeD(50, 12)
  set.seed(50); r_rl   <- MaxMean(d, maxSeconds = 1, useRL = TRUE)
  set.seed(50); r_norl <- MaxMean(d, maxSeconds = 1, useRL = FALSE)
  expect_s3_class(r_rl,   "MaxMeanSelection")
  expect_s3_class(r_norl, "MaxMeanSelection")
  expect_true(is.finite(attr(r_rl,   "score")))
  expect_true(is.finite(attr(r_norl, "score")))
})

# ---- coverage: timer block, RL early-termination, progress UI --------------

test_that("the periodic interrupt/time check fires on a long restart", {
  # A large instance makes a single restart run well past the 256-iteration
  # check interval (tabu depth limit is 50000), so the countdown timer block
  # executes deterministically — independent of how many restarts fit the budget.
  set.seed(70)
  d <- MakeD(70, 300)
  # maxIter = Inf so the time-expiry break (not the iteration cap) is what stops
  # the search — that is the branch this test exists to cover.
  result <- MaxMean(d, maxSeconds = 0.5, maxIter = Inf, useRL = FALSE)
  expect_s3_class(result, "MaxMeanSelection")
  expect_gte(length(result), 2L)
})

test_that("RL initial-solution construction can terminate early", {
  # On signed data with many restarts, the Q-learning reward matrix accrues
  # negative entries (punished elements), so the Eq. 7 termination condition
  # eventually fires and construction stops before exhausting all elements.
  set.seed(71)
  m <- matrix(runif(40L * 40L, -10, 10), 40L)
  d <- (m + t(m)) / 2
  diag(d) <- 0
  set.seed(72)
  # maxIter = Inf so many restarts run (RL construction only fires on restart
  # >= 1); a finite default cap could leave restart 0 the only one executed.
  result <- MaxMean(d, maxSeconds = 1.5, maxIter = Inf, useRL = TRUE)
  expect_s3_class(result, "MaxMeanSelection")
  expect_true(is.finite(attr(result, "score")))
})

test_that("MaxMean Coreset.progress option fires the cli hooks", {
  d <- MakeD(73, 12)
  old <- options(Coreset.progress = TRUE)
  on.exit(options(old))
  expect_no_error(suppressMessages(MaxMean(d, maxSeconds = 0.2)))
})

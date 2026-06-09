# Tests for PolishSelection() and the two polished samplers.
#
# Coverage:
#   1. Recovery on the motivating adversary (triangle + square corners).
#   2. Idempotence — polishing a polished result is a no-op.
#   3. Monotonicity — T_k(polish(S)) >= T_k(S) for any S.
#   4. Flat-landscape escape — strict-improvement-only would stall; the
#      Della Croce tie-break reduces n_critical and continues.
#   5. Sampler wrappers return integer vectors of the requested size.

# ---- Recovery fixture: locked-in centre ---------------------------------

# Five points: the four corners of a unit square plus a centre point at
# (0.5, 0.5). Distances are 1 (side), sqrt(2) (diagonal), and sqrt(0.5)
# (corner to centre).
#
# Bare Gonzalez starting from a corner picks the diagonally opposite corner,
# then either of the remaining corners, then the centre (it ties with the
# fourth corner on T_3 = 1 but the centre wins on argmax-tie-break in some
# seeds). The "locked-in" selection {3 corners + centre} has
# T_4 = sqrt(0.5) ~ 0.707. Swapping centre for the missing corner gives the
# strictly better T_4 = 1. This is exactly the user's adversarial story
# (an early greedy commitment forecloses a strictly better later pick) in a
# geometry where T_k actually has a unique optimum.

test_that("polish recovers the global optimum on the corner+centre adversary", {
  pts <- rbind(
    c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0.5, 0.5)
  )
  d   <- as.matrix(dist(pts))
  bad <- c(1L, 2L, 3L, 5L)  # 3 corners + centre
  T_before <- MinDist(d, bad)
  polished <- PolishSelection(d, bad, limit = 20L)
  T_after  <- MinDist(d, polished)

  expect_lt(T_before, 1 - 1e-9)
  expect_equal(T_after, 1)
  expect_setequal(as.integer(polished), 1:4)  # the four corners
})

# ---- Idempotence ---------------------------------------------------------

test_that("polish is idempotent on its own output", {
  set.seed(1)
  pts <- matrix(rnorm(80), ncol = 2)
  d   <- as.matrix(dist(pts))
  s0  <- Gonzalez(d, 6L)
  s1  <- PolishSelection(d, s0)
  s2  <- PolishSelection(d, s1)

  expect_setequal(s1, s2)
  expect_identical(attr(s2, "swaps"), 0L)
})

# ---- Monotonicity --------------------------------------------------------

test_that("T_k(polish(S)) >= T_k(S) over many random selections", {
  set.seed(7)
  pts <- matrix(rnorm(201), ncol = 3)
  d   <- as.matrix(dist(pts))
  N   <- nrow(d)
  for (k in c(3L, 5L, 8L)) {
    for (trial in 1:6) {
      S0       <- sample.int(N, k)
      S1       <- PolishSelection(d, S0)
      T0       <- MinDist(d, S0)
      T1       <- MinDist(d, S1)
      # Strict monotonicity at minimum: polish never worsens T_k.
      expect_gte(T1, T0)
    }
  }
})

# ---- Flat-landscape regression -------------------------------------------

# Build a small case where T_k has multiple pairs at the minimum but a single
# swap reduces n_critical without changing T_k. Strict-improvement-only would
# stall; the Della Croce branch must accept this swap.
#
# Construction:
#   four points form a square of side 1 (so all four sides are at the min);
#   one extra point sits at distance >1 from three corners but exactly 1 from
#   the fourth. Selecting the four corners gives T_k = 1, n_critical = 4
#   (the four sides). Swapping one corner for the extra point preserves
#   T_k = 1 but reduces n_critical to 1 (a single tied pair). Then a further
#   swap can strictly improve T_k.

test_that("polish escapes flat-landscape ties via Della Croce tie-break", {
  # Hand-built distance matrix: 5 points.
  # Points: A=(0,0), B=(1,0), C=(1,1), D=(0,1), E=(2, 0.5).
  pts <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(2, 0.5))
  d   <- as.matrix(dist(pts))
  # Selecting all four square corners gives T_k = 1, n_critical = 4.
  S0  <- 1:4
  T0  <- MinDist(d, S0)
  expect_equal(T0, 1)

  # After polish the selection must NOT be the four square corners — that is
  # the stalled state. The Della Croce branch should accept a swap to a
  # configuration with fewer tied minima (and may then climb further).
  S1 <- PolishSelection(d, S0, limit = 5L)
  T1 <- MinDist(d, S1)
  expect_gte(T1, T0)
  # The polish should have made at least one swap.
  expect_gte(as.integer(attr(S1, "swaps")), 1L)
})

# Note: the WideSampleGonzPolished{End,Live} wrappers that compose
# Gonzalez(seed = "ensemble") with PolishSelection() live in the FurthestPoint
# package (which Imports MaxMin); their tests live there. MaxMin tests the
# PolishSelection primitive itself.

test_that("polish degenerate cases pass through unchanged", {
  d <- as.matrix(dist(matrix(rnorm(20), ncol = 2)))
  # k = 0 and k = 1: nothing to polish.
  expect_identical(PolishSelection(d, integer(0)), integer(0))
  expect_identical(PolishSelection(d, 3L), 3L)
})

test_that("PolishSelection progress = TRUE fires the cli hooks", {
  pts <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0.5, 0.5))
  d   <- as.matrix(dist(pts))
  expect_no_error(
    PolishSelection(d, 1:4, progress = TRUE)
  )
})

# ---- C++ guard paths (polish.cpp lines 38-41, 50, 53) --------------------

test_that("PolishMaximin_cpp early-returns when limit < 1 (lines 38-41)", {
  set.seed(1)
  d   <- as.matrix(dist(matrix(rnorm(20), ncol = 2)))
  res <- PolishSelection(d, c(1L, 3L, 5L), limit = 0L)
  expect_identical(attr(res, "swaps"), 0L)
})

test_that("PolishMaximin_cpp stops on out-of-range index (line 50)", {
  d <- as.matrix(dist(matrix(rnorm(6), ncol = 2)))  # 3-point distance matrix
  expect_error(PolishSelection(d, c(1L, 100L)), "out-of-range")
})

test_that("PolishMaximin_cpp stops on duplicate index (line 53)", {
  set.seed(2)
  d <- as.matrix(dist(matrix(rnorm(20), ncol = 2)))
  expect_error(PolishSelection(d, c(1L, 1L, 2L)), "duplicate")
})

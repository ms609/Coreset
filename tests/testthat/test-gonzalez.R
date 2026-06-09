# Tests for Gonzalez() seeding strategies and its distance-column oracle path.

make_data <- function(seed = 42, N = 60, dim = 4) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

single_seeds <- c("diameter", "anti_medoid", "medoid", "rowsum", "rownorm",
                  "peripheral", "random_furthest")

test_that("matrix and coordinate paths agree for every seed strategy", {
  dat <- make_data()
  for (s in single_seeds) {
    for (n in c(2L, 6L, 12L)) {
      # Seed the RNG identically so the random_furthest pivot matches on both
      # paths; deterministic strategies are unaffected by it.
      set.seed(1); mat <- Gonzalez(dat$d, n, seed = s)
      set.seed(1); pt  <- Gonzalez(n = n, points = dat$pts, seed = s)
      expect_identical(mat, pt, info = paste("seed", s, "n", n))
    }
  }
  # A centroid-free ensemble agrees across paths. Explicit pivots index points
  # directly and the distances are bit-identical on Euclidean data, so the
  # random-furthest starts coincide too.
  for (n in c(2L, 6L, 12L)) {
    mat <- Gonzalez(dat$d, n, seed = c("peripheral", "random_furthest"),
                    pivots = c(3L, 17L, 28L))
    pt  <- Gonzalez(n = n, points = dat$pts,
                    seed = c("peripheral", "random_furthest"),
                    pivots = c(3L, 17L, 28L))
    attributes(mat) <- NULL
    attributes(pt)  <- NULL
    expect_identical(mat, pt, info = paste("ensemble n", n))
  }
})

test_that("integer seed gives a bare pass from that index", {
  dat <- make_data()
  # integer seed selects the start index; any prior char strategy is irrelevant.
  a <- Gonzalez(dat$d, 8L, seed = 3L)
  b <- Gonzalez(dat$d, 8L, seed = 3L)
  expect_identical(a, b)
  expect_identical(a[[1]], 3L)
})

test_that("ensemble keeps the best anchor by T_k", {
  dat <- make_data()
  n   <- 8L
  anchors <- c("diameter", "anti_medoid", "rowsum", "rownorm")
  ens <- Gonzalez(dat$d, n, seed = anchors)
  ens_tk <- MinDist(dat$d, ens)
  for (s in anchors) {
    expect_gte(ens_tk + 1e-9, MinDist(dat$d, Gonzalez(dat$d, n, seed = s)))
  }
  expect_true(all(attr(ens, "winning_strategy") %in% anchors))
  expect_length(attr(ens, "strategy_results"), 4L)
  # Two-anchor ensemble works and has only those two entries.
  two <- Gonzalez(dat$d, n, seed = c("diameter", "rowsum"))
  expect_length(attr(two, "strategy_results"), 2L)
})

test_that("the default ensemble is the best-of-five O(N) seeds", {
  dat <- make_data()
  n   <- 8L
  # Matrix path: centroid is coordinate-only, so the default is peripheral plus
  # three random-furthest starts -> four strategies.
  mat <- Gonzalez(dat$d, n)
  expect_length(attr(mat, "strategy_results"), 4L)
  expect_identical(names(attr(mat, "strategy_results")),
                   c("peripheral", "random_furthest1", "random_furthest2",
                     "random_furthest3"))
  # Coordinate path adds centroid -> five strategies.
  pt <- Gonzalez(n = n, points = dat$pts)
  expect_length(attr(pt, "strategy_results"), 5L)
  expect_true("centroid" %in% names(attr(pt, "strategy_results")))
  # The default includes peripheral, so it is at least as good as it.
  expect_gte(MinDist(dat$d, mat) + 1e-9,
             MinDist(dat$d, Gonzalez(dat$d, n, seed = "peripheral")))
})

test_that("the default selection is reproducible under set.seed", {
  dat <- make_data()
  set.seed(1); a <- Gonzalez(dat$d, 8L)
  set.seed(1); b <- Gonzalez(dat$d, 8L)
  expect_identical(c(a), c(b))
})

test_that("pivots vector controls the random-furthest starts", {
  dat <- make_data()
  n   <- 8L
  # Unspecified draws three pivots (matrix default: peripheral + 3 random).
  expect_length(attr(Gonzalez(dat$d, n), "strategy_results"), 4L)
  # The vector's length sets the count (matrix default: peripheral + pivots).
  expect_length(attr(Gonzalez(dat$d, n, pivots = 1:6), "strategy_results"), 7L)
  # User-chosen pivots are honoured: each seeds at the point furthest from it.
  res <- Gonzalez(dat$d, n, pivots = c(2L, 30L))
  sr  <- attr(res, "strategy_results")
  expect_identical(sr[["random_furthest1"]]$s1,
                   as.integer(which.max(dat$d[, 2L])))
  expect_identical(sr[["random_furthest2"]]$s1,
                   as.integer(which.max(dat$d[, 30L])))
  # integer(0), NA, and NULL all disable the random starts equivalently: the
  # matrix default then reduces to peripheral alone.
  for (none in list(integer(0), NA, NULL)) {
    z <- Gonzalez(dat$d, n, pivots = none)
    expect_length(attr(z, "strategy_results"), 1L)
    expect_identical(attr(z, "winning_strategy"), "peripheral")
  }
  # Out-of-range pivots are rejected.
  expect_error(Gonzalez(dat$d, n, pivots = c(1L, 999L)), "pivots")
})

test_that("trivial cardinalities are handled", {
  dat <- make_data(N = 10)
  expect_identical(Gonzalez(dat$d, 0L), integer(0))
  expect_identical(Gonzalez(dat$d, 20L), seq_len(10L))
  expect_length(Gonzalez(dat$d, 1L, seed = "medoid"), 1L)
  # n == 1 under ensemble: all t_k are NA, first anchor wins.
  expect_length(Gonzalez(dat$d, 1L), 1L)
})

test_that("input validation", {
  dat <- make_data(N = 10)
  expect_error(Gonzalez(dat$d, -1L), "non-negative")
  expect_error(Gonzalez(dat$d, c(1L, 2L)), "single")
  expect_error(Gonzalez(dat$d, 3L, seed = "nope"), "arg")
  expect_error(Gonzalez("not a matrix", 3L), "dist|matrix")
})

test_that("explicitly naming centroid in a matrix ensemble warns and drops it", {
  dat <- make_data(N = 12)
  expect_warning(res <- Gonzalez(dat$d, 4L,
                                 seed = c("centroid", "peripheral"),
                                 pivots = integer(0)),
                 "coordinates")
  # Dropped, leaving peripheral alone.
  expect_identical(names(attr(res, "strategy_results")), "peripheral")
})

# ---- distance-column oracle path ----------------------------------------

test_that("the column-oracle path matches the matrix path given the same seed", {
  dat <- make_data()
  colFn <- function(i) dat$d[, i]
  for (n in c(2L, 5L, 15L)) {
    expect_identical(
      Gonzalez(colFn, n, N = nrow(dat$d), seed = 1L),
      Gonzalez(dat$d, n, seed = 1L)
    )
  }
})

test_that("column-oracle peripheral seed is deterministic and matrix-matched", {
  dat <- make_data()
  colFn <- function(i) dat$d[, i]
  # The default (ensemble) seed is unreachable from an oracle, so the path
  # falls back to the deterministic peripheral seed.
  s1 <- Gonzalez(colFn, 7L, N = nrow(dat$d))
  s2 <- Gonzalez(colFn, 7L, N = nrow(dat$d))
  expect_identical(s1, s2)
  # peripheral matrix seed should match Gonzalez(seed = "peripheral").
  expect_identical(s1, Gonzalez(dat$d, 7L, seed = "peripheral"))
})

test_that("column-oracle guards and contract", {
  dat <- make_data(N = 12)
  colFn <- function(i) dat$d[, i]
  expect_identical(Gonzalez(colFn, 0L, N = 12L), integer(0))
  expect_identical(Gonzalez(colFn, 20L, N = 12L), seq_len(12L))
  expect_length(Gonzalez(colFn, 1L, N = 12L, seed = 4L), 1L)
  # N is required on the oracle path: it cannot be inferred from the closure.
  expect_error(Gonzalez(colFn, 3L), "N")
  bad <- function(i) 1:3
  expect_error(Gonzalez(bad, 3L, N = 12L, seed = 1L), "length")
})

test_that("column-oracle warns on an unreachable named seed but not the default", {
  dat <- make_data(N = 12)
  colFn <- function(i) dat$d[, i]
  # A character/ensemble seed cannot be honoured from an oracle -> warn.
  expect_warning(Gonzalez(colFn, 4L, N = 12L, seed = "diameter"), "integer")
  expect_warning(Gonzalez(colFn, 4L, N = 12L, seed = c("diameter", "rowsum")),
                 "integer")
  # The default (unsupplied) seed and an integer seed are silent.
  expect_silent(Gonzalez(colFn, 4L, N = 12L))
  expect_silent(Gonzalez(colFn, 4L, N = 12L, seed = 1L))
})

# ---- .AsPointsMatrix validation (lines 40, 43, 46, 49-50) ------------------

test_that(".AsPointsMatrix coerces non-matrix, converts integer, and rejects bad input", {
  # Non-matrix coerced via as.matrix() (line 40): a plain vector becomes Nx1.
  r <- Gonzalez(n = 3L, points = 1:20, seed = 1L)
  expect_length(r, 3L)
  # Non-numeric matrix errors (line 43).
  expect_error(
    Gonzalez(n = 2L, points = matrix(c("a", "b", "c", "d"), 2L, 2L)),
    "numeric"
  )
  # Integer storage mode is silently coerced to double (line 46).
  r_int <- Gonzalez(n = 3L, points = matrix(1L:20L, ncol = 4L), seed = 1L)
  expect_length(r_int, 3L)
  # NA entries error (lines 49-50).
  expect_error(
    Gonzalez(n = 2L, points = matrix(c(1, 2, NA, 4), 2L, 2L)),
    "NA"
  )
})

# ---- .MaximinFromColumn progress (lines 351, 359) ---------------------------

test_that(".MaximinFromColumn progress = TRUE fires the cli hooks", {
  dat <- make_data()
  colFn <- function(i) dat$d[, i]
  expect_no_error(
    Gonzalez(colFn, 5L, N = nrow(dat$d), seed = 1L, progress = TRUE)
  )
})

# ---- .GonzalezColumn N/n/first validation (lines 311, 314, 323) ------------

test_that(".GonzalezColumn validates N < 1, n < 0, and first out of bounds", {
  dat <- make_data(N = 12)
  colFn <- function(i) dat$d[, i]
  # N < 1: line 311
  expect_error(Gonzalez(colFn, 3L, N = 0L),           "N")
  # n < 0: line 314
  expect_error(Gonzalez(colFn, -1L, N = 12L, seed = 1L), "n")
  # first out of bounds (> N): line 323
  expect_error(Gonzalez(colFn, 3L, N = 12L, seed = 15L), "first")
})

# ---- MaximinFrom_cpp stop on out-of-range first (maximin.cpp:16) ------------

test_that("MaximinFrom_cpp stops when seed index is out of range", {
  dat <- make_data(N = 10)
  # seed = 0L -> first = 0 < 1 -> Rcpp::stop in maximin.cpp line 16
  expect_error(Gonzalez(dat$d, 3L, seed = 0L), "first")
})

# ---- MaximinFromPoints_cpp stop on out-of-range first (maximin_points.cpp:57) --

test_that("MaximinFromPoints_cpp stops when seed index is out of range", {
  dat <- make_data(N = 10)
  # seed = 0L -> first = 0 < 1 -> Rcpp::stop in maximin_points.cpp line 57
  expect_error(Gonzalez(n = 3L, points = dat$pts, seed = 0L), "first")
})

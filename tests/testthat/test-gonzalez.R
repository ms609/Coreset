# Tests for Gonzalez() seeding strategies and its distance-column oracle path.

make_data <- function(seed = 42, N = 60, dim = 4) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

single_seeds <- c("diameter", "anti_medoid", "medoid", "rowsum", "rownorm",
                  "peripheral")

test_that("matrix and coordinate paths agree for every seed strategy", {
  dat <- make_data()
  for (s in single_seeds) {
    for (n in c(2L, 6L, 12L)) {
      mat <- Gonzalez(dat$d, n, seed = s)
      pt  <- Gonzalez(n = n, points = dat$pts, seed = s)
      expect_identical(mat, pt, info = paste("seed", s, "n", n))
    }
  }
  # Default (full ensemble) also agrees across paths.
  for (n in c(2L, 6L, 12L)) {
    mat <- Gonzalez(dat$d, n)
    pt  <- Gonzalez(n = n, points = dat$pts)
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
  ens <- Gonzalez(dat$d, n)   # default: full four-anchor ensemble
  ens_tk <- TkScore(dat$d, ens)
  for (s in c("diameter", "anti_medoid", "rowsum", "rownorm")) {
    expect_gte(ens_tk + 1e-9, TkScore(dat$d, Gonzalez(dat$d, n, seed = s)))
  }
  expect_true(all(attr(ens, "winning_strategy") %in%
                    c("diameter", "anti_medoid", "rowsum", "rownorm")))
  expect_length(attr(ens, "strategy_results"), 4L)
  # Two-anchor ensemble works and has only those two entries.
  two <- Gonzalez(dat$d, n, seed = c("diameter", "rowsum"))
  expect_length(attr(two, "strategy_results"), 2L)
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

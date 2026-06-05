# Tests for Gonzalez() seeding strategies and GonzalezColumn().

make_data <- function(seed = 42, N = 60, dim = 4) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

single_seeds <- c("diameter", "anti_medoid", "medoid", "rowsum", "rownorm",
                  "peripheral")

test_that("matrix and coordinate paths agree for every seed strategy", {
  dat <- make_data()
  for (s in c("ensemble", single_seeds)) {
    for (n in c(2L, 6L, 12L)) {
      mat <- Gonzalez(dat$d, n, seed = s)
      pt  <- Gonzalez(n = n, points = dat$pts, seed = s)
      attributes(mat) <- NULL
      attributes(pt)  <- NULL
      expect_identical(mat, pt, info = paste("seed", s, "n", n))
    }
  }
})

test_that("explicit `first` overrides `seed` and gives a bare pass", {
  dat <- make_data()
  # first index forced; seed argument must be ignored.
  a <- Gonzalez(dat$d, 8L, first = 3L, seed = "ensemble")
  b <- Gonzalez(dat$d, 8L, first = 3L, seed = "diameter")
  expect_identical(a, b)
  expect_identical(a[[1]], 3L)
})

test_that("ensemble keeps the best anchor by T_k", {
  dat <- make_data()
  n   <- 8L
  ens <- Gonzalez(dat$d, n, seed = "ensemble")
  ens_tk <- TkScore(dat$d, ens)
  for (s in c("diameter", "anti_medoid", "rowsum", "rownorm")) {
    expect_gte(ens_tk + 1e-9, TkScore(dat$d, Gonzalez(dat$d, n, seed = s)))
  }
  expect_true(attr(ens, "winning_strategy") %in%
                c("diameter", "anti_medoid", "rowsum", "rownorm"))
  expect_length(attr(ens, "strategy_results"), 4L)
})

test_that("trivial cardinalities are handled", {
  dat <- make_data(N = 10)
  expect_identical(Gonzalez(dat$d, 0L), integer(0))
  expect_identical(Gonzalez(dat$d, 20L), seq_len(10L))
  expect_length(Gonzalez(dat$d, 1L, seed = "medoid"), 1L)
  # n == 1 under ensemble: all t_k are NA, first anchor (diameter) wins.
  expect_length(Gonzalez(dat$d, 1L, seed = "ensemble"), 1L)
})

test_that("input validation", {
  dat <- make_data(N = 10)
  expect_error(Gonzalez(dat$d, -1L), "non-negative")
  expect_error(Gonzalez(dat$d, c(1L, 2L)), "single")
  expect_error(Gonzalez(dat$d, 3L, seed = "nope"))
  expect_error(Gonzalez("not a matrix", 3L), "dist|matrix")
})

# ---- GonzalezColumn -----------------------------------------------------

test_that("GonzalezColumn matches the matrix path given the same seed", {
  dat <- make_data()
  colFn <- function(i) dat$d[, i]
  for (n in c(2L, 5L, 15L)) {
    expect_identical(
      GonzalezColumn(colFn, nrow(dat$d), n, first = 1L),
      Gonzalez(dat$d, n, first = 1L)
    )
  }
})

test_that("GonzalezColumn peripheral seed is deterministic and matrix-matched", {
  dat <- make_data()
  colFn <- function(i) dat$d[, i]
  s1 <- GonzalezColumn(colFn, nrow(dat$d), 7L)
  s2 <- GonzalezColumn(colFn, nrow(dat$d), 7L)
  expect_identical(s1, s2)
  # peripheral matrix seed should match Gonzalez(seed = "peripheral").
  expect_identical(s1, Gonzalez(dat$d, 7L, seed = "peripheral"))
})

test_that("GonzalezColumn guards and contract", {
  dat <- make_data(N = 12)
  colFn <- function(i) dat$d[, i]
  expect_identical(GonzalezColumn(colFn, 12L, 0L), integer(0))
  expect_identical(GonzalezColumn(colFn, 12L, 20L), seq_len(12L))
  expect_length(GonzalezColumn(colFn, 12L, 1L, first = 4L), 1L)
  expect_error(GonzalezColumn("nope", 12L, 3L), "function")
  bad <- function(i) 1:3
  expect_error(GonzalezColumn(bad, 12L, 3L, first = 1L), "length")
})

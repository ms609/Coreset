# Tests for MaxMinSeed() and TkScore().

make_data <- function(seed = 7, N = 50, dim = 3) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

test_that("MaxMinSeed anchors match their definitions (matrix)", {
  dat <- make_data()
  d <- dat$d
  expect_identical(MaxMinSeed(d, method = "medoid"),
                   as.integer(which.min(rowSums(d))))
  expect_identical(MaxMinSeed(d, method = "rowsum"),
                   as.integer(which.max(rowSums(d))))
  expect_identical(MaxMinSeed(d, method = "rownorm"),
                   as.integer(which.max(rowSums(d ^ 2))))
  med <- which.min(rowSums(d))
  dd <- d[, med]; dd[med] <- -Inf
  expect_identical(MaxMinSeed(d, method = "anti_medoid"),
                   as.integer(which.max(dd)))
  d_off <- d; diag(d_off) <- -Inf
  expect_identical(MaxMinSeed(d, method = "diameter"),
                   as.integer(arrayInd(which.max(d_off), dim(d_off))[1L, 1L]))
})

test_that("MaxMinSeed coordinate path matches the matrix path", {
  dat <- make_data()
  for (m in c("peripheral", "diameter", "anti_medoid", "medoid",
              "rowsum", "rownorm")) {
    expect_identical(MaxMinSeed(points = dat$pts, method = m),
                     MaxMinSeed(dat$d, method = m), info = m)
  }
})

test_that("MaxMinSeed validates method", {
  dat <- make_data(N = 8)
  expect_error(MaxMinSeed(dat$d, method = "ensemble"))
  expect_error(MaxMinSeed(dat$d, method = "first"))
})

# ---- TkScore ------------------------------------------------------------

test_that("TkScore matrix and coordinate paths agree", {
  dat <- make_data()
  idx <- c(1L, 5L, 9L, 20L, 33L)
  expect_equal(TkScore(dat$d, idx),
               TkScore(idx = idx, points = dat$pts))
  # equals the brute-force minimum pairwise distance.
  sub <- dat$d[idx, idx]; diag(sub) <- Inf
  expect_equal(TkScore(dat$d, idx), min(sub))
})

test_that("TkScore returns NA for fewer than two points", {
  dat <- make_data(N = 6)
  expect_true(is.na(TkScore(dat$d, 3L)))
  expect_true(is.na(TkScore(idx = integer(0), points = dat$pts)))
})

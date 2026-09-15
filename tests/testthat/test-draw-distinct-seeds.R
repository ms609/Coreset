# Tests for DrawDistinctSeeds().

MakeData <- function(seed = 7, N = 60, dim = 3) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

test_that("each seed is the furthest element from some pivot", {
  dat <- MakeData()
  set.seed(1)
  seeds <- DrawDistinctSeeds(dat$d, nSeeds = 5)
  reachable <- unique(apply(dat$d, 2, which.max))
  expect_type(seeds, "integer")
  expect_true(all(seeds %in% reachable))
  expect_false(anyDuplicated(seeds) > 0)
  expect_false(is.unsorted(seeds))
  expect_lte(length(seeds), 5L)
})

test_that("matrix, dist, coordinate and column paths draw the same seeds", {
  dat <- MakeData()
  Column <- function(i) sqrt(colSums((t(dat$pts) - dat$pts[i, ]) ^ 2))
  ColumnNoSelf <- function(i) Column(i)[-i]
  set.seed(3); m <- DrawDistinctSeeds(dat$d, nSeeds = 4)
  set.seed(3); dd <- DrawDistinctSeeds(as.dist(dat$d), nSeeds = 4)
  set.seed(3); p <- DrawDistinctSeeds(points = dat$pts, nSeeds = 4)
  set.seed(3); f <- DrawDistinctSeeds(Column, N = nrow(dat$pts), nSeeds = 4)
  set.seed(3); g <- DrawDistinctSeeds(ColumnNoSelf, N = nrow(dat$pts), nSeeds = 4)
  expect_identical(dd, m)
  expect_identical(p, m)
  expect_identical(f, m)
  expect_identical(g, m)
})

test_that("seeds are FarFirst's random-furthest seeds", {
  dat <- MakeData()
  for (nSeeds in c(1L, 3L, 7L)) {
    set.seed(11)
    ff <- FarFirst(8L, dat$d, nSeeds = nSeeds)
    set.seed(11)
    seeds <- DrawDistinctSeeds(dat$d, nSeeds = nSeeds)
    scores <- vapply(seeds, function(s)
      attr(FarFirst(8L, dat$d, strategy = s), "score"), numeric(1))
    expect_identical(max(scores), attr(ff, "score"))

    set.seed(11)
    ffPts <- FarFirst(8L, points = dat$pts, nSeeds = nSeeds)
    set.seed(11)
    expect_identical(DrawDistinctSeeds(points = dat$pts, nSeeds = nSeeds), seeds)
    expect_equal(attr(ffPts, "score"), max(scores))
  }
})

test_that("a pool smaller than nSeeds returns every reachable seed", {
  # Points on a line: every pivot's furthest element is one of the two ends.
  pts <- matrix(c(1, 2, 3, 5, 8), ncol = 1)
  set.seed(1)
  expect_identical(DrawDistinctSeeds(points = pts, nSeeds = 4), c(1L, 5L))
})

test_that("maxDraws caps the pivots tried", {
  dat <- MakeData()
  set.seed(1)
  expect_length(DrawDistinctSeeds(dat$d, nSeeds = 10, maxDraws = 1), 1L)
  set.seed(1)
  pivot <- sample.int(nrow(dat$d), 1L)
  set.seed(1)
  expect_identical(DrawDistinctSeeds(dat$d, nSeeds = 10, maxDraws = 1),
                   as.integer(which.max(dat$d[, pivot])))
})

test_that("inputs are validated", {
  dat <- MakeData()
  Column <- function(i) dat$d[, i]
  expect_error(DrawDistinctSeeds(dat$d, nSeeds = 0), "nSeeds")
  expect_error(DrawDistinctSeeds(dat$d, nSeeds = NA), "nSeeds")
  expect_error(DrawDistinctSeeds(dat$d, nSeeds = 1:2), "nSeeds")
  expect_error(DrawDistinctSeeds(dat$d, maxDraws = 0), "maxDraws")
  expect_error(DrawDistinctSeeds(dat$d, maxDraws = c(1, 2)), "maxDraws")
  expect_error(DrawDistinctSeeds(Column), "`N`")
  expect_error(DrawDistinctSeeds(Column, N = NA), "`N`")
  expect_error(DrawDistinctSeeds(), "supply")
})

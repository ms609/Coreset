# Tests for the MaxMinSelection print + format S3 layer.

MakeData <- function(seed = 42, N = 60, dim = 4) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

test_that("FarFirst returns a self-describing MaxMinSelection that still indexes", {
  dat <- MakeData(N = 30)
  set.seed(1)
  sel <- FarFirst(6L, dat$d)
  expect_s3_class(sel, "MaxMinSelection")
  expect_type(sel, "integer")          # still an integer vector underneath
  expect_identical(attr(sel, "producer"), "FarFirst")
  # It indexes a matrix exactly like the bare vector it wraps.
  expect_identical(dim(dat$d[sel, sel]), c(6L, 6L))
})

test_that("the ensemble summary names the winning strategy and count", {
  dat <- MakeData(N = 30)
  set.seed(1)
  sel <- FarFirst(6L, dat$d)               # default: 8 random-furthest starts
  line <- format(sel)
  expect_match(line, "^6 elements \\([0-9 ]+\\) selected by ")
  expect_match(line, "farthest-first \\(best of [0-9]+ strategies")
  expect_match(line, "each at distance >= ")
})

test_that("a single pass omits the ensemble clause", {
  dat <- MakeData(N = 30)
  sel <- FarFirst(6L, dat$d, strategy = "diameter")
  line <- format(sel)
  expect_match(line, "selected by farthest-first, each at distance >= ")
  expect_false(grepl("best of", line))
})

test_that("a single element drops the distance clause (NA score)", {
  dat <- MakeData(N = 30)
  sel <- FarFirst(1L, dat$d, strategy = "medoid")
  line <- format(sel)
  expect_match(line, "^1 element \\([0-9]+\\) selected by farthest-first$")
  expect_false(grepl("distance", line))
})

test_that("DropAdd and Grasp report their own algorithm names", {
  dat <- MakeData(N = 30)
  da <- DropAdd(6L, dat$d)
  expect_s3_class(da, "MaxMinSelection")
  expect_match(format(da), "selected by DropAdd tabu search, each at distance >= ")
  set.seed(1)
  gr <- Grasp(d = dat$d, 6L, plateau = 20L, eliteSize = 4L)
  expect_s3_class(gr, "MaxMinSelection")
  expect_match(format(gr), "selected by GRASP with path-relinking")
})

test_that("the index list is truncated past the show limit", {
  dat <- MakeData(N = 60)
  sel <- FarFirst(25L, dat$d, strategy = 1L)   # 25 > the 20-index threshold
  line <- format(sel)
  expect_match(line, "\\.\\.\\. \\(\\+5 more\\)")
  # A short selection is shown in full, no ellipsis.
  expect_false(grepl("more", format(FarFirst(5L, dat$d, strategy = 1L))))
})

test_that("print returns its argument invisibly and emits the format line", {
  dat <- MakeData(N = 30)
  sel <- FarFirst(4L, dat$d, strategy = "peripheral")
  expect_output(ret <- withVisible(print(sel)), "selected by")
  expect_false(ret$visible)
  expect_identical(ret$value, sel)
})

test_that("an empty selection is left bare (no class)", {
  expect_identical(MaxMin:::.AsMaxMinSelection(integer(0), "FarFirst"),
                   integer(0))
  expect_identical(FarFirst(0L, MakeData(N = 10)$d), integer(0))
})

test_that(".FarFirstSelectedBy reports a tie among winning strategies", {
  fake <- structure(
    c(1L, 2L, 3L),
    score            = 1.5,
    winning_strategy = c("diameter", "rowsum"),
    strategy_results = list(diameter = 1, rowsum = 2, medoid = 3),
    producer         = "FarFirst",
    class            = "MaxMinSelection"
  )
  expect_match(format(fake),
               "best of 3 strategies, 2 tied: diameter, rowsum")
})

test_that(".MaxMinSelectedBy falls back for an unrecorded producer", {
  fake <- structure(c(1L, 2L), score = NA_real_, class = "MaxMinSelection")
  expect_match(format(fake), "selected by an unrecorded method")
})

test_that("ExactMaxMin prints proof status, indices and objective", {
  # Built directly so the test does not require the `highs` solver.
  proven <- .AsMaxMinSelection(
    structure(c(1L, 2L, 3L), score = 0.5, proven = TRUE,
              time_s = 0.1, solver = "highs", N = 10L, k = 3L),
    producer = "ExactMaxMin"
  )
  expect_s3_class(proven, "MaxMinSelection")
  expect_true(inherits(proven, "MaxMinSelection"))
  expect_match(format(proven),
               "^3 elements \\(1 2 3\\) selected by exact solver, proven optimal, each at distance >= 0.5$")

  incumbent <- proven
  attr(incumbent, "proven") <- FALSE
  expect_match(format(incumbent), "unproven incumbent, each at distance >= 0.5")

  expect_output(ret <- withVisible(print(proven)), "proven optimal")
  expect_false(ret$visible)
  expect_identical(ret$value, proven)
})

# ---- summary methods --------------------------------------------------------

test_that("summary of a FarFirst ensemble prints the per-strategy table", {
  dat <- MakeData(N = 30)
  set.seed(1)
  sel <- FarFirst(6L, dat$d)
  out <- capture.output(ret <- withVisible(summary(sel)))
  expect_false(ret$visible)
  expect_identical(ret$value, sel)
  expect_match(out[1], "selected by")                       # headline
  expect_true(any(grepl("strategies tried \\([0-9]+\\), best marked", out)))
  expect_true(any(grepl("strategy\\s+seed\\s+T_k", out)))
  expect_true(any(grepl("^  \\* random_furthest", out)))     # a winner is marked
  # Default is 8 starts; at least one row should appear.
  expect_gte(sum(grepl("^  [* ] random_furthest[0-9]", out)), 1L)
})

test_that("summary of a single FarFirst pass is just the headline", {
  dat <- MakeData(N = 30)
  out <- capture.output(summary(FarFirst(6L, dat$d, strategy = "diameter")))
  expect_length(out, 1L)
  expect_match(out, "farthest-first")
})

test_that("summary of DropAdd reports the secondary objective and effort", {
  dat <- MakeData(N = 30)
  out <- capture.output(summary(DropAdd(6L, dat$d)))
  expect_match(out[1], "DropAdd tabu search")
  expect_true(any(grepl("sum of pairwise distances:", out)))
  expect_true(any(grepl("iterations:", out)))
  expect_true(any(grepl("time:.* s", out)))
})

test_that("summary of Grasp reports iterations and path-relinking calls", {
  dat <- MakeData(N = 30)
  set.seed(1)
  out <- capture.output(summary(Grasp(d = dat$d, 6L, plateau = 20L, eliteSize = 4L)))
  expect_match(out[1], "GRASP")
  expect_true(any(grepl("refinement iterations:", out)))
  expect_true(any(grepl("path-relinking calls:", out)))
})

test_that("summary table tolerates NA T_k (all-NA ensemble)", {
  dat <- MakeData(N = 30)
  set.seed(1)
  out <- capture.output(summary(FarFirst(1L, dat$d)))   # m = 1 => every T_k NA
  expect_true(any(grepl("\\bNA\\b", out)))
})

test_that("summary of ExactMaxMin reports instance, objective and proof status", {
  proven <- .AsMaxMinSelection(
    structure(c(1L, 2L, 3L), score = 0.5, proven = TRUE,
              time_s = 0.1, solver = "highs", N = 10L, k = 3L),
    producer = "ExactMaxMin"
  )
  out <- capture.output(ret <- withVisible(summary(proven)))
  expect_false(ret$visible)
  expect_match(out[1], "proven optimal")
  expect_true(any(grepl("instance:\\s+n = 10, k = 3", out)))
  expect_true(any(grepl("objective:\\s+0.5 \\(proven optimal\\)", out)))
  expect_true(any(grepl("solver:\\s+highs", out)))

  incumbent <- proven; attr(incumbent, "proven") <- FALSE
  out2 <- capture.output(summary(incumbent))
  expect_true(any(grepl("objective:\\s+0.5 \\(lower bound \\(unproven\\)\\)", out2)))
})

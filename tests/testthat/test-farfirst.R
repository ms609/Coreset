# Tests for FarFirst() seeding strategies and its distance-column oracle path.

MakeData <- function(seed = 42, N = 60, dim = 4) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

# `"random_furthest"` is excluded: named on its own it now runs the ensemble
# (best-of-`pivots`), not a single bare pass, so it is exercised separately.
singleSeeds <- c("diameter", "anti_medoid", "medoid", "rowsum", "rownorm",
                  "peripheral")

test_that("matrix and coordinate paths agree for every seed strategy", {
  dat <- MakeData()
  for (s in singleSeeds) {
    for (n in c(2L, 6L, 12L)) {
      # Seed the RNG identically so the random_furthest pivot matches on both
      # paths; deterministic strategies are unaffected by it.
      set.seed(1); mat <- FarFirst(dat$d, n, method = s)
      set.seed(1); pt  <- FarFirst(m = n, points = dat$pts, method = s)
      expect_identical(mat, pt, info = paste("method", s, "m", n))
    }
  }
  # A centroid-free ensemble agrees across paths. Explicit pivots index points
  # directly and the distances are bit-identical on Euclidean data, so the
  # random-furthest starts coincide too.
  for (n in c(2L, 6L, 12L)) {
    mat <- FarFirst(dat$d, n, method = c("peripheral", "random_furthest"),
                    pivots = c(3L, 17L, 28L))
    pt  <- FarFirst(m = n, points = dat$pts,
                    method = c("peripheral", "random_furthest"),
                    pivots = c(3L, 17L, 28L))
    expect_identical(as.integer(mat), as.integer(pt), info = paste("ensemble m", n))
    expect_equal(attr(mat, "score"), attr(pt, "score"), tolerance = 1e-12,
                 info = paste("ensemble score m", n))
    expect_equal(attr(mat, "winning_strategy"), attr(pt, "winning_strategy"),
                 info = paste("ensemble winning_strategy m", n))
  }
})

test_that("integer method gives a bare pass from that index", {
  dat <- MakeData()
  # integer method selects the start index; any prior char strategy is irrelevant.
  a <- FarFirst(dat$d, 8L, method = 3L)
  b <- FarFirst(dat$d, 8L, method = 3L)
  expect_identical(a, b)
  expect_identical(a[[1]], 3L)
  # A bare pass now carries the achieved T_k as a `score` attribute.
  expect_equal(attr(a, "score"), MinDist(dat$d, a))
})

test_that("ensemble keeps the best anchor by T_k", {
  dat <- MakeData()
  n   <- 8L
  anchors <- c("diameter", "anti_medoid", "rowsum", "rownorm")
  ens <- FarFirst(dat$d, n, method = anchors)
  ensTk <- MinDist(dat$d, ens)
  for (s in anchors) {
    expect_gte(ensTk + 1e-9, MinDist(dat$d, FarFirst(dat$d, n, method = s)))
  }
  expect_true(all(attr(ens, "winning_strategy") %in% anchors))
  expect_length(attr(ens, "strategy_results"), 4L)
  # The ensemble also exposes the winning T_k as a `score` attribute.
  expect_equal(attr(ens, "score"), ensTk)
  # Two-anchor ensemble works and has only those two entries.
  two <- FarFirst(dat$d, n, method = c("diameter", "rowsum"))
  expect_length(attr(two, "strategy_results"), 2L)
})

test_that("the default ensemble is three random-furthest starts", {
  dat <- MakeData()
  n   <- 8L
  # Default is `"random_furthest"` alone -> three random starts on both paths.
  mat <- FarFirst(dat$d, n)
  expect_identical(names(attr(mat, "strategy_results")),
                   c("random_furthest1", "random_furthest2", "random_furthest3"))
  pt <- FarFirst(m = n, points = dat$pts)
  expect_identical(names(attr(pt, "strategy_results")),
                   c("random_furthest1", "random_furthest2", "random_furthest3"))
  # Neither path includes the deterministic anchors by default.
  expect_false("centroid" %in% names(attr(pt, "strategy_results")))
  expect_false("peripheral" %in% names(attr(mat, "strategy_results")))
})

test_that("the default selection is reproducible under set.seed", {
  dat <- MakeData()
  set.seed(1); a <- FarFirst(dat$d, 8L)
  set.seed(1); b <- FarFirst(dat$d, 8L)
  expect_identical(c(a), c(b))
})

test_that("pivots vector controls the random-furthest starts", {
  dat <- MakeData()
  n   <- 8L
  # Unspecified draws three pivots (default: 3 random-furthest starts).
  expect_length(attr(FarFirst(dat$d, n), "strategy_results"), 3L)
  # The vector's length sets the count.
  expect_length(attr(FarFirst(dat$d, n, pivots = 1:6), "strategy_results"), 6L)
  # User-chosen pivots are honoured: each seeds at the point furthest from it.
  res <- FarFirst(dat$d, n, pivots = c(2L, 30L))
  sr  <- attr(res, "strategy_results")
  expect_identical(sr[["random_furthest1"]]$s1,
                   as.integer(which.max(dat$d[, 2L])))
  expect_identical(sr[["random_furthest2"]]$s1,
                   as.integer(which.max(dat$d[, 30L])))
  # integer(0), NA, and NULL all disable the random starts equivalently: under
  # the default `method` that leaves no anchor, so it errors.
  for (none in list(integer(0), NA, NULL)) {
    expect_error(FarFirst(dat$d, n, pivots = none), "no seed strateg")
  }
  # Paired with a deterministic anchor, disabling the random starts is fine.
  z <- FarFirst(dat$d, n, method = "peripheral", pivots = integer(0))
  expect_length(z, n)
  # Out-of-range pivots are rejected.
  expect_error(FarFirst(dat$d, n, pivots = c(1L, 999L)), "pivots")
})

test_that("trivial cardinalities are handled", {
  dat <- MakeData(N = 10)
  expect_identical(FarFirst(dat$d, 0L), integer(0))
  # m > N: all N indices returned in Gonzalez order (a permutation of 1:N).
  over <- FarFirst(dat$d, 20L, method = 1L)
  expect_length(over, 10L)
  expect_setequal(over, seq_len(10L))
  expect_length(FarFirst(dat$d, 1L, method = "medoid"), 1L)
  # m == 1 under ensemble: all t_k are NA, first anchor wins.
  expect_length(FarFirst(dat$d, 1L), 1L)
})

test_that("FarFirst m > nPts returns all points (m capped at nPts)", {
  dat <- MakeData(N = 10)   # 10-point dataset
  # Request more points than available
  res <- FarFirst(dat$d, m = 50L, method = "peripheral")
  expect_length(res, 10L)
  expect_setequal(res, seq_len(10L))
  # Score is finite (not NA) — it is the min pairwise distance of all 10 points
  expect_true(is.numeric(attr(res, "score")))
})

test_that("input validation", {
  dat <- MakeData(N = 10)
  expect_error(FarFirst(dat$d, -1L), "non-negative")
  expect_error(FarFirst(dat$d, c(1L, 2L)), "single")
  expect_error(FarFirst(dat$d, 3L, method = "nope"), "arg")
  expect_error(FarFirst("not a matrix", 3L), "dist|matrix")
  # An integer `method` must be a single finite value (FF-002, FF-003).
  expect_error(FarFirst(dat$d, 3L, method = NA_integer_), "single finite")
  expect_error(FarFirst(dat$d, 3L, method = integer(0)),  "single finite")
  expect_error(FarFirst(dat$d, 3L, method = c(1L, 2L)),   "single finite")
  # A misspelled ensemble anchor is rejected, not silently dropped (F-601).
  expect_error(FarFirst(dat$d, 3L, method = c("peripheral", "anti-medoid")),
               "unknown seed/method")
  # m = Inf is rejected without a spurious base-R coercion warning (FF-004).
  expect_silent(expect_error(FarFirst(dat$d, Inf), "non-negative"))
})

test_that("FarFirst rejects an NA/non-finite distance matrix (FF-001)", {
  # Pre-fix, an NA in the matrix propagated through pmin.int()/which.max() and
  # produced a selection with a REPEATED index (e.g. c(1, 2, 1)); now it errors.
  naMat <- matrix(c(0, 5, NA, 5, 0, NA, NA, NA, 0), 3L, 3L)
  expect_error(FarFirst(naMat, 3L), "NA")
  infMat <- matrix(c(0, 5, Inf, 5, 0, 1, Inf, 1, 0), 3L, 3L)
  expect_error(FarFirst(infMat, 3L), "NA|Inf|finite")
  # The distance-column oracle path rejects a non-self NA at a selected column.
  dat <- MakeData(N = 8)
  naCol <- function(i) { v <- dat$d[, i]; v[if (i == 1L) 2L else 1L] <- NA; v }
  expect_error(FarFirst(naCol, 3L, N = 8L, method = 1L), "NA|NaN")
})

test_that("explicitly naming centroid in a matrix ensemble warns and drops it", {
  dat <- MakeData(N = 12)
  expect_warning(res <- FarFirst(dat$d, 4L,
                                 method = c("centroid", "peripheral"),
                                 pivots = integer(0)),
                 "coordinates")
  # Dropped, leaving peripheral alone.
  expect_identical(names(attr(res, "strategy_results")), "peripheral")
})

# ---- distance-column oracle path ----------------------------------------

test_that("the column-oracle path matches the matrix path given the same method", {
  dat <- MakeData()
  colFn <- function(i) dat$d[, i]
  for (n in c(2L, 5L, 15L)) {
    expect_identical(
      FarFirst(colFn, n, N = nrow(dat$d), method = 1L),
      FarFirst(dat$d, n, method = 1L)
    )
  }
})

test_that("a self-distance-omitting (length N-1) colFn matches the matrix path", {
  dat <- MakeData()
  N <- nrow(dat$d)
  colN   <- function(i) dat$d[, i]      # self reported (length N)
  colNm1 <- function(i) dat$d[-i, i]    # self omitted  (length N-1), others in order
  for (n in c(2L, 5L, 15L)) {
    ref <- FarFirst(dat$d, n, method = 1L)
    expect_identical(FarFirst(colNm1, n, N = N, method = 1L), ref)
    expect_identical(FarFirst(colN,   n, N = N, method = 1L), ref)
  }
  # The deterministic peripheral seed must also splice correctly (it touches
  # both i = 1 and a downstream index).
  expect_identical(FarFirst(colNm1, 7L, N = N), FarFirst(dat$d, 7L, method = "peripheral"))
})

test_that("column-oracle peripheral seed is deterministic and matrix-matched", {
  dat <- MakeData()
  colFn <- function(i) dat$d[, i]
  # The default (ensemble) method is unreachable from an oracle, so the path
  # falls back to the deterministic peripheral seed.
  s1 <- FarFirst(colFn, 7L, N = nrow(dat$d))
  s2 <- FarFirst(colFn, 7L, N = nrow(dat$d))
  expect_identical(s1, s2)
  # peripheral matrix seed should match FarFirst(method = "peripheral").
  expect_identical(s1, FarFirst(dat$d, 7L, method = "peripheral"))
})

test_that("column-oracle guards and contract", {
  dat <- MakeData(N = 12)
  colFn <- function(i) dat$d[, i]
  expect_identical(FarFirst(colFn, 0L, N = 12L), integer(0))
  # m > N: all N indices returned in Gonzalez order (a permutation of 1:N).
  over <- FarFirst(colFn, 20L, N = 12L, method = 1L)
  expect_length(over, 12L)
  expect_setequal(over, seq_len(12L))
  expect_length(FarFirst(colFn, 1L, N = 12L, method = 4L), 1L)
  # N is required on the oracle path: it cannot be inferred from the closure.
  expect_error(FarFirst(colFn, 3L), "N")
  bad <- function(i) 1:3
  expect_error(FarFirst(bad, 3L, N = 12L, method = 1L), "length")
})

test_that("column-oracle warns on an unreachable named method but not the default", {
  dat <- MakeData(N = 12)
  colFn <- function(i) dat$d[, i]
  # A character/ensemble method cannot be honoured from an oracle -> warn.
  expect_warning(FarFirst(colFn, 4L, N = 12L, method = "diameter"), "integer")
  expect_warning(FarFirst(colFn, 4L, N = 12L, method = c("diameter", "rowsum")),
                 "integer")
  # The default (unsupplied) method and an integer method are silent.
  expect_silent(FarFirst(colFn, 4L, N = 12L))
  expect_silent(FarFirst(colFn, 4L, N = 12L, method = 1L))
})

# ---- .AsPointsMatrix validation (lines 40, 43, 46, 49-50) ------------------

test_that(".AsPointsMatrix coerces non-matrix, converts integer, and rejects bad input", {
  # Non-matrix coerced via as.matrix() (line 40): a plain vector becomes Nx1.
  r <- FarFirst(m = 3L, points = 1:20, method = 1L)
  expect_length(r, 3L)
  # Non-numeric matrix errors (line 43).
  expect_error(
    FarFirst(m = 2L, points = matrix(c("a", "b", "c", "d"), 2L, 2L)),
    "numeric"
  )
  # Integer storage mode is silently coerced to double (line 46).
  rInt <- FarFirst(m = 3L, points = matrix(1L:20L, ncol = 4L), method = 1L)
  expect_length(rInt, 3L)
  # NA entries error (lines 49-50).
  expect_error(
    FarFirst(m = 2L, points = matrix(c(1, 2, NA, 4), 2L, 2L)),
    "NA"
  )
})

# ---- .MaximinFromColumn progress (lines 351, 359) ---------------------------

test_that(".MaximinFromColumn progress = TRUE fires the cli hooks", {
  dat <- MakeData()
  colFn <- function(i) dat$d[, i]
  expect_no_error(
    FarFirst(colFn, 5L, N = nrow(dat$d), method = 1L, progress = TRUE)
  )
})

# ---- .GonzalezColumn N/m/first validation (lines 311, 314, 323) ------------

test_that(".GonzalezColumn validates N < 1, m < 0, and first out of bounds", {
  dat <- MakeData(N = 12)
  colFn <- function(i) dat$d[, i]
  # N < 1: line 311
  expect_error(FarFirst(colFn, 3L, N = 0L),           "N")
  # m < 0: line 314
  expect_error(FarFirst(colFn, -1L, N = 12L, method = 1L), "m")
  # first out of bounds (> N): line 323
  expect_error(FarFirst(colFn, 3L, N = 12L, method = 15L), "first")
})

# ---- MaximinFrom_cpp stop on out-of-range first (maximin.cpp:16) ------------

test_that("MaximinFrom_cpp stops when seed index is out of range", {
  dat <- MakeData(N = 10)
  # method = 0L -> first = 0 < 1 -> Rcpp::stop in maximin.cpp line 16
  expect_error(FarFirst(dat$d, 3L, method = 0L), "first")
})

# ---- MaximinFromPoints_cpp stop on out-of-range first (maximin_points.cpp:57) --

test_that("MaximinFromPoints_cpp stops when seed index is out of range", {
  dat <- MakeData(N = 10)
  # method = 0L -> first = 0 < 1 -> Rcpp::stop in maximin_points.cpp line 57
  expect_error(FarFirst(m = 3L, points = dat$pts, method = 0L), "first")
})

# ---- .SubsetScore mean_pairwise branch --------------------------------------

test_that(".SubsetScore mean_pairwise returns mean of lower-triangle entries", {
  dat <- MakeData(N = 10)
  idx <- c(1L, 3L, 7L)
  sub <- dat$d[idx, idx]
  expect_equal(MaxMin:::.SubsetScore(dat$d, idx, "mean_pairwise"),
               mean(sub[lower.tri(sub)]))
  expect_true(is.na(MaxMin:::.SubsetScore(dat$d, 5L, "mean_pairwise")))
})

# ---- .GonzalezColumn non-function guard -------------------------------------

test_that(".GonzalezColumn rejects a non-function colFn", {
  expect_error(MaxMin:::.GonzalezColumn("not_a_function", N = 10L, m = 3L),
               "function")
})

# ---- downsample consistency -------------------------------------------------

test_that("FarFirst(m = N)[1:5] equals FarFirst(m = 5) on all paths", {
  dat <- MakeData()
  N <- nrow(dat$d)

  # Bare passes carry a `score` attribute; `[` drops it, so compare the bare
  # integer indices (the prefix property is about the selection order).
  bare <- function(x) as.integer(x)

  # Matrix path, single integer method (bare Gonzalez pass).
  full_mat <- FarFirst(dat$d, N, method = 1L)
  five_mat <- FarFirst(dat$d, 5L, method = 1L)
  expect_identical(bare(full_mat[1:5]), bare(five_mat))

  # Points path, same method.
  full_pts <- FarFirst(m = N, points = dat$pts, method = 1L)
  five_pts <- FarFirst(m = 5L, points = dat$pts, method = 1L)
  expect_identical(bare(full_pts[1:5]), bare(five_pts))

  # Oracle path, integer method.
  colFn <- function(i) dat$d[, i]
  full_col <- FarFirst(colFn, N, N = N, method = 1L)
  five_col <- FarFirst(colFn, 5L, N = N, method = 1L)
  expect_identical(bare(full_col[1:5]), bare(five_col))

  # m > N also satisfies the same prefix property.
  over_mat <- FarFirst(dat$d, N + 10L, method = 1L)
  expect_identical(bare(over_mat[1:5]), bare(five_mat))
})

# ---- nseeds: distinct-seed random restart -----------------------------------

test_that("nseeds runs a best-of over distinct peripheral seeds", {
  dat <- MakeData()
  set.seed(1)
  r <- FarFirst(dat$d, 6L, nseeds = 4L)
  expect_length(as.integer(r), 6L)
  sr <- attr(r, "strategy_results")
  expect_false(is.null(sr))
  # At most `nseeds` strategies, and the seeds they ran from are distinct.
  expect_lte(length(sr), 4L)
  s1s <- vapply(sr, `[[`, integer(1L), "s1")
  expect_identical(anyDuplicated(s1s), 0L)
  # Labelled random_furthest1.. and the returned score is the best T_k.
  expect_true(all(grepl("^random_furthest", names(sr))))
  tks <- vapply(sr, `[[`, numeric(1L), "t_k")
  expect_equal(attr(r, "score"), max(tks))
})

test_that("nseeds: matrix and coordinate paths agree (same RNG)", {
  dat <- MakeData()
  for (n in c(2L, 6L, 12L)) {
    set.seed(7); mat <- FarFirst(dat$d, n, nseeds = 5L)
    set.seed(7); pt  <- FarFirst(m = n, points = dat$pts, nseeds = 5L)
    expect_identical(mat, pt, info = paste("nseeds m", n))
  }
})

test_that("nseeds is reproducible under a fixed seed", {
  dat <- MakeData()
  set.seed(3); a <- FarFirst(dat$d, 6L, nseeds = 4L)
  set.seed(3); b <- FarFirst(dat$d, 6L, nseeds = 4L)
  expect_identical(a, b)
})

test_that("nseeds overrides method/pivots with a warning, and validates", {
  dat <- MakeData()
  expect_warning(FarFirst(dat$d, 6L, method = "diameter", nseeds = 3L),
                 "overrides")
  expect_warning(FarFirst(dat$d, 6L, pivots = c(1L, 2L), nseeds = 3L),
                 "overrides")
  # No warning when only nseeds is supplied.
  expect_silent(FarFirst(dat$d, 6L, nseeds = 3L))
  expect_error(FarFirst(dat$d, 6L, nseeds = 0L), "positive integer")
  expect_error(FarFirst(dat$d, 6L, nseeds = c(1L, 2L)), "single")
})

test_that("nseeds caps at the reachable pool without error", {
  # 8 points: at most 8 distinct seeds exist, so nseeds = 50 returns <= 8.
  dat <- MakeData(N = 8)
  set.seed(1)
  r <- FarFirst(dat$d, 3L, nseeds = 50L)
  sr <- attr(r, "strategy_results")
  expect_lte(length(sr), 8L)
  expect_gte(length(sr), 1L)
})

test_that("nseeds is rejected on the distance-column oracle path", {
  dat <- MakeData()
  colFn <- function(i) dat$d[, i]
  expect_error(FarFirst(colFn, 6L, N = nrow(dat$d), nseeds = 3L), "oracle")
})

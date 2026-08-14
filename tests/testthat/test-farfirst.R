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
      set.seed(1); mat <- FarFirst(n, dat$d, strategy = s)
      set.seed(1); pt  <- FarFirst(k = n, points = dat$pts, strategy = s)
      expect_identical(mat, pt, info = paste("strategy", s, "k", n))
    }
  }
  # A anti_centroid-free ensemble agrees across paths. With the same RNG seed the
  # random-furthest seeds coincide on both paths (distances are bit-identical on
  # Euclidean data).
  for (n in c(2L, 6L, 12L)) {
    set.seed(5)
    mat <- FarFirst(n, dat$d, strategy = c("peripheral", "random_furthest"))
    set.seed(5)
    pt  <- FarFirst(k = n, points = dat$pts,
                    strategy = c("peripheral", "random_furthest"))
    expect_identical(as.integer(mat), as.integer(pt), info = paste("ensemble k", n))
    expect_equal(attr(mat, "score"), attr(pt, "score"), tolerance = 1e-12,
                 info = paste("ensemble score k", n))
    expect_equal(attr(mat, "winning_strategy"), attr(pt, "winning_strategy"),
                 info = paste("ensemble winning_strategy k", n))
  }
})

test_that("integer strategy gives a bare pass from that index", {
  dat <- MakeData()
  # integer strategy selects the start index; any prior char strategy is irrelevant.
  a <- FarFirst(8L, dat$d, strategy = 3L)
  b <- FarFirst(8L, dat$d, strategy = 3L)
  expect_identical(a, b)
  expect_identical(a[[1]], 3L)
  # A bare pass now carries the achieved T_k as a `score` attribute.
  expect_equal(attr(a, "score"), MinDist(d = dat$d, a))
})

test_that("ensemble keeps the best anchor by T_k", {
  dat <- MakeData()
  n   <- 8L
  anchors <- c("diameter", "anti_medoid", "rowsum", "rownorm")
  ens <- FarFirst(n, dat$d, strategy = anchors)
  ensTk <- MinDist(d = dat$d, ens)
  for (s in anchors) {
    expect_gte(ensTk + 1e-9, MinDist(d = dat$d, FarFirst(n, dat$d, strategy = s)))
  }
  expect_true(all(attr(ens, "winning_strategy") %in% anchors))
  expect_length(attr(ens, "strategy_results"), 4L)
  # The ensemble also exposes the winning T_k as a `score` attribute.
  expect_equal(attr(ens, "score"), ensTk)
  # Two-anchor ensemble works and has only those two entries.
  two <- FarFirst(n, dat$d, strategy = c("diameter", "rowsum"))
  expect_length(attr(two, "strategy_results"), 2L)
})

test_that("the default ensemble is three random-furthest starts", {
  dat <- MakeData()
  n   <- 8L
  # Default is `"random_furthest"` alone -> three random starts on both paths.
  mat <- FarFirst(n, dat$d)
  expect_identical(names(attr(mat, "strategy_results")),
                   paste0("random_furthest", 1:3))
  pt <- FarFirst(k = n, points = dat$pts)
  expect_identical(names(attr(pt, "strategy_results")),
                   paste0("random_furthest", 1:3))
  # Neither path includes the deterministic anchors by default.
  expect_false("anti_centroid" %in% names(attr(pt, "strategy_results")))
  expect_false("peripheral" %in% names(attr(mat, "strategy_results")))
})

test_that("the default selection is reproducible under set.seed", {
  dat <- MakeData()
  set.seed(1); a <- FarFirst(8L, dat$d)
  set.seed(1); b <- FarFirst(8L, dat$d)
  expect_identical(c(a), c(b))
})

test_that("nSeeds controls the number of random-furthest starts", {
  dat <- MakeData()
  n   <- 8L
  # Default is 3 random-furthest starts.
  expect_length(attr(FarFirst(n, dat$d), "strategy_results"), 3L)
  # nSeeds sets the count.
  expect_length(attr(FarFirst(n, dat$d, nSeeds = 4L), "strategy_results"), 4L)
  expect_length(attr(FarFirst(n, dat$d, nSeeds = 1L), "strategy_results"), 1L)
})

test_that("trivial cardinalities are handled", {
  dat <- MakeData(N = 10)
  expect_identical(FarFirst(0L, dat$d), integer(0))
  # k > N: all N indices returned in Gonzalez order (a permutation of 1:N).
  over <- FarFirst(20L, dat$d, strategy = 1L)
  expect_length(over, 10L)
  expect_setequal(over, seq_len(10L))
  expect_length(FarFirst(1L, dat$d, strategy = "medoid"), 1L)
  # k == 1 under ensemble: all t_k are NA, first anchor wins.
  expect_length(FarFirst(1L, dat$d), 1L)
})

test_that("FarFirst k > nPts returns all points (k capped at nPts)", {
  dat <- MakeData(N = 10)   # 10-point dataset
  # Request more points than available
  res <- FarFirst(k = 50L, dat$d, strategy = "peripheral")
  expect_length(res, 10L)
  expect_setequal(res, seq_len(10L))
  # Score is finite (not NA) â€” it is the min pairwise distance of all 10 points
  expect_true(is.numeric(attr(res, "score")))
})

test_that("input validation", {
  dat <- MakeData(N = 10)
  expect_error(FarFirst(-1L, dat$d), "non-negative")
  expect_error(FarFirst(c(1L, 2L), dat$d), "single")
  expect_error(FarFirst(3L, dat$d, strategy = "nope"), "arg")
  expect_error(FarFirst(3L, "not a matrix"), "dist|matrix")
  # An integer `strategy` must be a single finite value (FF-002, FF-003).
  expect_error(FarFirst(3L, dat$d, strategy = NA_integer_), "single finite")
  expect_error(FarFirst(3L, dat$d, strategy = integer(0)),  "single finite")
  expect_error(FarFirst(3L, dat$d, strategy = c(1L, 2L)),   "single finite")
  # A misspelled ensemble anchor is rejected, not silently dropped (F-601).
  expect_error(FarFirst(3L, dat$d, strategy = c("peripheral", "anti-medoid")),
               "unknown strateg")
  # k = Inf is rejected without a spurious base-R coercion warning (FF-004).
  expect_silent(expect_error(FarFirst(Inf, dat$d), "non-negative"))
})

test_that("FarFirst rejects an NA/non-finite distance matrix (FF-001)", {
  # Pre-fix, an NA in the matrix propagated through pmin.int()/which.max() and
  # produced a selection with a REPEATED index (e.g. c(1, 2, 1)); now it errors.
  naMat <- matrix(c(0, 5, NA, 5, 0, NA, NA, NA, 0), 3L, 3L)
  expect_error(FarFirst(3L, naMat), "NA")
  infMat <- matrix(c(0, 5, Inf, 5, 0, 1, Inf, 1, 0), 3L, 3L)
  expect_error(FarFirst(3L, infMat), "NA|Inf|finite")
  # The distance-column oracle path rejects a non-self NA at a selected column.
  dat <- MakeData(N = 8)
  naCol <- function(i) { v <- dat$d[, i]; v[if (i == 1L) 2L else 1L] <- NA; v }
  expect_error(FarFirst(3L, naCol, N = 8L, strategy = 1L), "NA|NaN")
})

test_that("N is warned on matrix/coordinate paths but not on the oracle path", {
  dat <- MakeData(N = 12)
  # Matrix path: N is ignored -> warn.
  expect_warning(FarFirst(4L, dat$d, N = 12L, strategy = 1L), "ignored")
  # Coordinate path: N is ignored -> warn.
  expect_warning(FarFirst(4L, points = dat$pts, N = 12L, strategy = 1L), "ignored")
  # Oracle path: N is required and used -> no warning.
  colFn <- function(i) dat$d[, i]
  expect_no_warning(FarFirst(4L, colFn, N = nrow(dat$d), strategy = 1L))
})

test_that("explicitly naming anti_centroid in a matrix ensemble warns and drops it", {
  dat <- MakeData(N = 12)
  expect_warning(res <- FarFirst(4L, dat$d,
                                 strategy = c("anti_centroid", "peripheral")),
                 "coordinates")
  # Dropped, leaving peripheral alone.
  expect_identical(names(attr(res, "strategy_results")), "peripheral")
})

# ---- distance-column oracle path ----------------------------------------

test_that("the column-oracle path matches the matrix path given the same strategy", {
  dat <- MakeData()
  colFn <- function(i) dat$d[, i]
  for (n in c(2L, 5L, 15L)) {
    expect_identical(
      FarFirst(n, colFn, N = nrow(dat$d), strategy = 1L),
      FarFirst(n, dat$d, strategy = 1L)
    )
  }
})

test_that("a self-distance-omitting (length N-1) colFn matches the matrix path", {
  dat <- MakeData()
  N <- nrow(dat$d)
  colN   <- function(i) dat$d[, i]      # self reported (length N)
  colNm1 <- function(i) dat$d[-i, i]    # self omitted  (length N-1), others in order
  for (n in c(2L, 5L, 15L)) {
    ref <- FarFirst(n, dat$d, strategy = 1L)
    expect_identical(FarFirst(n, colNm1, N = N, strategy = 1L), ref)
    expect_identical(FarFirst(  n, colN, N = N, strategy = 1L), ref)
  }
  # The deterministic peripheral seed must also splice correctly (it touches
  # both i = 1 and a downstream index).
  expect_identical(FarFirst(7L, colNm1, N = N), FarFirst(7L, dat$d, strategy = "peripheral"))
})

test_that("column-oracle peripheral seed is deterministic and matrix-matched", {
  dat <- MakeData()
  colFn <- function(i) dat$d[, i]
  # The default (ensemble) strategy is unreachable from an oracle, so the path
  # falls back to the deterministic peripheral seed.
  s1 <- FarFirst(7L, colFn, N = nrow(dat$d))
  s2 <- FarFirst(7L, colFn, N = nrow(dat$d))
  expect_identical(s1, s2)
  # peripheral matrix seed should match FarFirst(strategy = "peripheral").
  expect_identical(s1, FarFirst(7L, dat$d, strategy = "peripheral"))
})

test_that("column-oracle guards and contract", {
  dat <- MakeData(N = 12)
  colFn <- function(i) dat$d[, i]
  expect_identical(FarFirst(0L, colFn, N = 12L), integer(0))
  # k > N: all N indices returned in Gonzalez order (a permutation of 1:N).
  over <- FarFirst(20L, colFn, N = 12L, strategy = 1L)
  expect_length(over, 12L)
  expect_setequal(over, seq_len(12L))
  expect_length(FarFirst(1L, colFn, N = 12L, strategy = 4L), 1L)
  # N is required on the oracle path: it cannot be inferred from the closure.
  expect_error(FarFirst(3L, colFn), "N")
  bad <- function(i) 1:3
  expect_error(FarFirst(3L, bad, N = 12L, strategy = 1L), "length")
})

test_that("column-oracle warns on an unreachable named strategy but not the default", {
  dat <- MakeData(N = 12)
  colFn <- function(i) dat$d[, i]
  # A character/ensemble strategy cannot be honoured from an oracle -> warn.
  expect_warning(FarFirst(4L, colFn, N = 12L, strategy = "diameter"), "integer")
  expect_warning(FarFirst(4L, colFn, N = 12L, strategy = c("diameter", "rowsum")),
                 "integer")
  # The default (unsupplied) strategy and an integer strategy are silent.
  expect_silent(FarFirst(4L, colFn, N = 12L))
  expect_silent(FarFirst(4L, colFn, N = 12L, strategy = 1L))
})

# ---- .AsPointsMatrix validation (lines 40, 43, 46, 49-50) ------------------

test_that(".AsPointsMatrix coerces non-matrix, converts integer, and rejects bad input", {
  # Non-matrix coerced via as.matrix() (line 40): a plain vector becomes Nx1.
  r <- FarFirst(k = 3L, points = 1:20, strategy = 1L)
  expect_length(r, 3L)
  # Non-numeric matrix errors (line 43).
  expect_error(
    FarFirst(k = 2L, points = matrix(c("a", "b", "c", "d"), 2L, 2L)),
    "numeric"
  )
  # Integer storage mode is silently coerced to double (line 46).
  rInt <- FarFirst(k = 3L, points = matrix(1L:20L, ncol = 4L), strategy = 1L)
  expect_length(rInt, 3L)
  # NA entries error (lines 49-50).
  expect_error(
    FarFirst(k = 2L, points = matrix(c(1, 2, NA, 4), 2L, 2L)),
    "NA"
  )
})

# ---- .MaximinFromColumn progress (lines 351, 359) ---------------------------

test_that(".MaximinFromColumn progress bar fires without error", {
  dat <- MakeData()
  colFn <- function(i) dat$d[, i]
  old <- options(MaxMin.progress = TRUE)
  on.exit(options(old))
  expect_no_error(FarFirst(5L, colFn, N = nrow(dat$d), strategy = 1L))
})

# ---- .GonzalezColumn N/k/first validation (lines 311, 314, 323) ------------

test_that(".GonzalezColumn validates N < 1, k < 0, and first out of bounds", {
  dat <- MakeData(N = 12)
  colFn <- function(i) dat$d[, i]
  # N < 1: line 311
  expect_error(FarFirst(3L, colFn, N = 0L),           "N")
  # k < 0: line 314
  expect_error(FarFirst(-1L, colFn, N = 12L, strategy = 1L), "k")
  # first out of bounds (> N): line 323
  expect_error(FarFirst(3L, colFn, N = 12L, strategy = 15L), "first")
})

# ---- MaximinFrom_cpp stop on out-of-range first (maximin.cpp:16) ------------

test_that("MaximinFrom_cpp stops when seed index is out of range", {
  dat <- MakeData(N = 10)
  # strategy = 0L -> first = 0 < 1 -> Rcpp::stop in maximin.cpp line 16
  expect_error(FarFirst(3L, dat$d, strategy = 0L), "first")
})

# ---- MaximinFromPoints_cpp stop on out-of-range first (maximin_points.cpp:57) --

test_that("MaximinFromPoints_cpp stops when seed index is out of range", {
  dat <- MakeData(N = 10)
  # strategy = 0L -> first = 0 < 1 -> Rcpp::stop in maximin_points.cpp line 57
  expect_error(FarFirst(k = 3L, points = dat$pts, strategy = 0L), "first")
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
  expect_error(MaxMin:::.GonzalezColumn("not_a_function", N = 10L, k = 3L),
               "function")
})

# ---- downsample consistency -------------------------------------------------

test_that("FarFirst(k = N)[1:5] equals FarFirst(k = 5) on all paths", {
  dat <- MakeData()
  N <- nrow(dat$d)

  # Bare passes carry a `score` attribute; `[` drops it, so compare the bare
  # integer indices (the prefix property is about the selection order).
  bare <- function(x) as.integer(x)

  # Matrix path, single integer method (bare Gonzalez pass).
  full_mat <- FarFirst(N, dat$d, strategy = 1L)
  five_mat <- FarFirst(5L, dat$d, strategy = 1L)
  expect_identical(bare(full_mat[1:5]), bare(five_mat))

  # Points path, same method.
  full_pts <- FarFirst(k = N, points = dat$pts, strategy = 1L)
  five_pts <- FarFirst(k = 5L, points = dat$pts, strategy = 1L)
  expect_identical(bare(full_pts[1:5]), bare(five_pts))

  # Oracle path, integer method.
  colFn <- function(i) dat$d[, i]
  full_col <- FarFirst(N, colFn, N = N, strategy = 1L)
  five_col <- FarFirst(5L, colFn, N = N, strategy = 1L)
  expect_identical(bare(full_col[1:5]), bare(five_col))

  # k > N also satisfies the same prefix property.
  over_mat <- FarFirst(N + 10L, dat$d, strategy = 1L)
  expect_identical(bare(over_mat[1:5]), bare(five_mat))
})

# ---- nSeeds: distinct-seed random restart -----------------------------------

test_that("nSeeds runs a best-of over distinct peripheral seeds", {
  dat <- MakeData()
  set.seed(1)
  r <- FarFirst(6L, dat$d, nSeeds = 4L)
  expect_length(as.integer(r), 6L)
  sr <- attr(r, "strategy_results")
  expect_false(is.null(sr))
  # At most `nSeeds` strategies, and the seeds they ran from are distinct.
  expect_lte(length(sr), 4L)
  s1s <- vapply(sr, `[[`, integer(1L), "s1")
  expect_identical(anyDuplicated(s1s), 0L)
  # Labelled random_furthest1.. and the returned score is the best T_k.
  expect_true(all(grepl("^random_furthest", names(sr))))
  tks <- vapply(sr, `[[`, numeric(1L), "t_k")
  expect_equal(attr(r, "score"), max(tks))
})

test_that("nSeeds: matrix and coordinate paths agree (same RNG)", {
  dat <- MakeData()
  for (n in c(2L, 6L, 12L)) {
    set.seed(7); mat <- FarFirst(n, dat$d, nSeeds = 5L)
    set.seed(7); pt  <- FarFirst(k = n, points = dat$pts, nSeeds = 5L)

    # Indices, seeds and winners must agree exactly; the matrix and coordinate
    # kernels accumulate T_k via different summation orders, so a per-strategy
    # t_k can differ in its last bit even when every selection is identical --
    # compare those with tolerance rather than via expect_identical().
    info <- paste("nSeeds k", n)
    expect_identical(as.integer(mat), as.integer(pt), info = info)
    expect_equal(attr(mat, "score"), attr(pt, "score"), tolerance = 1e-9,
                 info = paste(info, "score"))
    expect_identical(attr(mat, "winning_strategy"), attr(pt, "winning_strategy"),
                      info = paste(info, "winning_strategy"))
    srMat <- attr(mat, "strategy_results")
    srPt  <- attr(pt, "strategy_results")
    expect_identical(names(srMat), names(srPt), info = paste(info, "names"))
    for (nm in names(srMat)) {
      expect_identical(srMat[[nm]]$s1, srPt[[nm]]$s1, info = paste(info, nm, "s1"))
      expect_identical(srMat[[nm]]$idx, srPt[[nm]]$idx, info = paste(info, nm, "idx"))
      expect_equal(srMat[[nm]]$t_k, srPt[[nm]]$t_k, tolerance = 1e-9,
                   info = paste(info, nm, "t_k"))
    }
  }
})

test_that("nSeeds is reproducible under a fixed seed", {
  dat <- MakeData()
  set.seed(3); a <- FarFirst(6L, dat$d, nSeeds = 4L)
  set.seed(3); b <- FarFirst(6L, dat$d, nSeeds = 4L)
  expect_identical(a, b)
})

test_that("nSeeds validates its argument", {
  dat <- MakeData()
  expect_error(FarFirst(6L, dat$d, nSeeds = 0L), "positive integer")
  expect_error(FarFirst(6L, dat$d, nSeeds = c(1L, 2L)), "single")
})

test_that("nSeeds caps at the reachable pool without error", {
  # 8 points: at most 8 distinct seeds exist, so nSeeds = 50 returns <= 8.
  dat <- MakeData(N = 8)
  set.seed(1)
  r <- FarFirst(3L, dat$d, nSeeds = 50L)
  sr <- attr(r, "strategy_results")
  expect_lte(length(sr), 8L)
  expect_gte(length(sr), 1L)
})

test_that("nSeeds is silently ignored on the distance-column oracle path", {
  dat <- MakeData()
  colFn <- function(i) dat$d[, i]
  expect_no_error(FarFirst(6L, colFn, N = nrow(dat$d), nSeeds = 3L))
})


test_that("FarFirst results are invariant to mc.cores (small anchors)", {
  # The O(N^2) anchor primitives parallelise at any non-trivial size; results
  # must be bit-identical at every thread count. CRAN caps tests at 2 cores.
  set.seed(21)
  pts <- matrix(rnorm(200L * 3L), ncol = 3L)
  d <- as.matrix(dist(pts))
  old <- options(mc.cores = NULL)
  on.exit(options(old), add = TRUE)
  anchors <- c("diameter", "anti_medoid", "medoid", "rowsum", "rownorm")
  options(mc.cores = 1L)
  a1 <- FarFirst(12L, d, strategy = anchors)
  p1 <- FarFirst(12L, points = pts, strategy = anchors)
  options(mc.cores = 2L)
  a2 <- FarFirst(12L, d, strategy = anchors)
  p2 <- FarFirst(12L, points = pts, strategy = anchors)
  expect_identical(a1, a2)
  expect_identical(p1, p2)
})

test_that("FarFirst results are invariant to mc.cores (parallel greedy pass)", {
  # N above the kernels' thread-engagement threshold (~32k), so the fused
  # pmin/argmax sweep really runs chunked; the chunk merge must reproduce the
  # serial first-maximum rule exactly.
  set.seed(22)
  pts <- matrix(rnorm(40000L * 2L), ncol = 2L)
  old <- options(mc.cores = NULL)
  on.exit(options(old), add = TRUE)
  options(mc.cores = 1L)
  s1 <- FarFirst(50L, points = pts, strategy = 5L)
  options(mc.cores = 2L)
  s2 <- FarFirst(50L, points = pts, strategy = 5L)
  expect_identical(s1, s2)
})

test_that("AllFinite_cpp matches the anyNA/is.finite guard semantics", {
  # .AsDistMatrix's finite check runs through this single-pass scan; it must
  # flag exactly what `anyNA(x) || any(!is.finite(x))` flags, for double and
  # integer storage, serially and through the parallel OR-reduction.
  af <- MaxMin:::AllFinite_cpp
  expect_true(af(c(0, 1.5, -2), 1L))
  expect_false(af(c(0, NA_real_), 1L))
  expect_false(af(c(0, NaN), 1L))
  expect_false(af(c(0, Inf), 1L))
  expect_false(af(c(0, -Inf), 1L))
  expect_true(af(c(1L, 3L), 1L))
  expect_false(af(c(1L, NA_integer_), 1L))
  expect_error(af("a", 1L), "storage")
  # The parallel reduction engages past 2^20 elements (CRAN 2-core cap).
  big <- runif(1200000L)
  expect_true(af(big, 2L))
  big[[999983L]] <- NaN
  expect_false(af(big, 2L))
})

test_that("points pass is cross-path identical at every dimension block shape", {
  # The fused update sweeps up to four dimensions per block, with distinct
  # code paths for a single-block dim (<= 4) and each final-block width
  # (dim mod 4). The matrix path never touches that code, so cross-path
  # agreement certifies every instantiation.
  for (dim in c(1:5, 8:11)) {
    set.seed(dim)
    pts <- matrix(rnorm(120L * dim), ncol = dim)
    d <- as.matrix(dist(pts))
    pm <- FarFirst(40L, d, strategy = 5L)
    pp <- FarFirst(40L, points = pts, strategy = 5L)
    expect_identical(as.integer(pm), as.integer(pp))
    expect_identical(attr(pm, "score"), attr(pp, "score"))
  }
})

test_that("multi-seed kernels reproduce the single-seed kernels exactly", {
  # The ensemble path runs its restarts through MaximinMulti*_cpp; each column
  # must equal the corresponding single-seed pass, indices and T_k alike, on
  # both the restart-parallel arm (nSeeds > 1, threads > 1) and the serial one.
  dat <- MakeData(seed = 7, N = 150L, dim = 6L)
  seeds <- c(3L, 41L, 99L, 128L)
  for (k in c(1L, 2L, 25L)) {
    refM <- lapply(seeds, function(s) MaxMin:::.MaximinFrom(dat$d, k, s))
    refP <- lapply(seeds, function(s) MaxMin:::.MaximinFromPoints(dat$pts, k, s))
    for (nThread in c(1L, 2L)) {
      gotM <- MaxMin:::MaximinMultiFrom_cpp(dat$d, k, seeds, nThread)
      gotP <- MaxMin:::MaximinMultiFromPoints_cpp(dat$pts, k, seeds, nThread)
      for (j in seq_along(seeds)) {
        expect_identical(gotM[["idx"]][, j], as.integer(refM[[j]]))
        expect_identical(gotM[["t_k"]][[j]], attr(refM[[j]], "t_k"))
        expect_identical(gotP[["idx"]][, j], as.integer(refP[[j]]))
        expect_identical(gotP[["t_k"]][[j]], attr(refP[[j]], "t_k"))
      }
    }
  }
})

test_that("multi-seed kernels reject out-of-range arguments", {
  dat <- MakeData(N = 20L)
  expect_error(MaxMin:::MaximinMultiFrom_cpp(dat$d, 0L, 1L, 1L), "'n' must be")
  expect_error(MaxMin:::MaximinMultiFrom_cpp(dat$d, 21L, 1L, 1L), "'n' must be")
  expect_error(MaxMin:::MaximinMultiFrom_cpp(dat$d, 5L, integer(0), 1L),
               "at least one seed")
  expect_error(MaxMin:::MaximinMultiFrom_cpp(dat$d, 5L, c(1L, 21L), 1L),
               "'firsts' must lie")
  expect_error(MaxMin:::MaximinMultiFrom_cpp(dat$d, 5L, 0L, 1L),
               "'firsts' must lie")
  expect_error(MaxMin:::MaximinMultiFromPoints_cpp(dat$pts, 0L, 1L, 1L),
               "'n' must be")
  expect_error(MaxMin:::MaximinMultiFromPoints_cpp(dat$pts, 5L, integer(0), 1L),
               "at least one seed")
  expect_error(MaxMin:::MaximinMultiFromPoints_cpp(dat$pts, 5L, c(1L, 99L), 1L),
               "'firsts' must lie")
})

test_that("restarts fall back to per-pass threading above the size threshold", {
  # Past ~32k points a single pass can occupy every thread itself, so the
  # kernel keeps the seed loop serial and threads the sweep instead. Only the
  # coordinate path reaches that size in practice. CRAN caps tests at 2 cores.
  skip_on_cran()
  set.seed(23)
  pts <- matrix(rnorm(33000L * 2L), ncol = 2L)
  seeds <- c(11L, 20001L)
  ref <- lapply(seeds, function(s) MaxMin:::.MaximinFromPoints(pts, 20L, s))
  got <- MaxMin:::MaximinMultiFromPoints_cpp(pts, 20L, seeds, 2L)
  for (j in seq_along(seeds)) {
    expect_identical(got[["idx"]][, j], as.integer(ref[[j]]))
    expect_identical(got[["t_k"]][[j]], attr(ref[[j]], "t_k"))
  }
})

test_that("default ensemble is invariant to mc.cores (restart-parallel arm)", {
  # The restarts run one per thread below the size threshold; the winner and
  # every strategy_results record must be identical at any thread count.
  old <- options(mc.cores = NULL)
  on.exit(options(old), add = TRUE)
  set.seed(24)
  pts <- matrix(rnorm(400L * 3L), ncol = 3L)
  d <- as.matrix(dist(pts))
  options(mc.cores = 1L)
  set.seed(5); m1 <- FarFirst(20L, d, nSeeds = 6L)
  set.seed(5); p1 <- FarFirst(20L, points = pts, nSeeds = 6L)
  options(mc.cores = 2L)
  set.seed(5); m2 <- FarFirst(20L, d, nSeeds = 6L)
  set.seed(5); p2 <- FarFirst(20L, points = pts, nSeeds = 6L)
  expect_identical(m1, m2)
  expect_identical(p1, p2)
})

test_that("an ensemble whose anchors share a seed keeps one record per label", {
  # Distinct seeds are solved once and mapped back to labels, so anchors that
  # collide on a seed still get their own strategy_results entry, and a label
  # that ties the best score is still reported as a winner.
  set.seed(25)
  pts <- matrix(rnorm(80L * 2L), ncol = 2L)
  d <- as.matrix(dist(pts))
  anchors <- c("diameter", "anti_medoid", "medoid", "rowsum", "rownorm",
               "peripheral")
  seeds <- vapply(anchors, function(a) MaxMin:::.PickPoint(d, a), integer(1L))
  expect_true(anyDuplicated(seeds) > 0L)     # the collision this test needs
  r <- FarFirst(8L, d, strategy = anchors)
  sr <- attr(r, "strategy_results")
  expect_identical(names(sr), anchors)
  expect_identical(vapply(sr, `[[`, integer(1L), "s1"), seeds)
  # Colliding anchors must report the same record, not merely the same score.
  dup <- anchors[duplicated(seeds) | duplicated(seeds, fromLast = TRUE)]
  expect_identical(sr[[dup[[1L]]]]$idx, sr[[dup[[2L]]]]$idx)
  expect_true(all(dup %in% attr(r, "winning_strategy")) ||
                !any(dup %in% attr(r, "winning_strategy")))
})

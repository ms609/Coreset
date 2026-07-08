# Tests for the composable-coreset `maxCandidates` thinning on DropAdd()/Grasp()
# and its shared helpers (.ResolveCap, .FarFirstThin).

MakeData <- function(seed = 42, N = 20, dim = 4) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

# ---- .ResolveCap: every branch ---------------------------------------------

test_that(".ResolveCap validates and normalises the cap", {
  R <- MaxMin:::.ResolveCap
  # Non-scalar / NA -> error.
  expect_error(R(c(1, 2), 10L, 3L), "single number")
  expect_error(R(NA_real_, 10L, 3L), "single number")
  # 0 / Inf disable thinning -> NA.
  expect_identical(R(0, 10L, 3L), NA_integer_)
  expect_identical(R(Inf, 10L, 3L), NA_integer_)
  # Negative / non-integer -> error.
  expect_error(R(-1, 10L, 3L), "positive integer")
  expect_error(R(3.5, 10L, 3L), "positive integer")
  # Cap at or above n is not binding -> NA (also exercises the overflow-safe
  # double comparison with a huge value).
  expect_identical(R(10L, 10L, 3L), NA_integer_)
  expect_identical(R(1e12, 10L, 3L), NA_integer_)
  # Positive cap below k -> error.
  expect_error(R(2L, 10L, 5L), ">= k")
  # Binding cap -> the integer coreset size.
  expect_identical(R(6L, 10L, 3L), 6L)
})

# ---- .FarFirstThin: the index round-trip (critical) ------------------------

test_that(".FarFirstThin maps coreset-local indices back to original space", {
  dat <- MakeData()
  m <- 8L
  core <- as.integer(FarFirst(m, d = dat$d, strategy = "peripheral"))
  # The peripheral coreset is NOT the first m rows, so a verbatim (no-remap)
  # return would be wrong -- this is exactly what the round-trip must fix.
  expect_false(identical(core, seq_len(m)))

  # A deterministic stub solver that "selects" coreset-local positions 1 and 3.
  stub <- function(d, points) {
    structure(c(1L, 3L), score = 42, secondary = 7,
              class = "MaxMinSelection", producer = "stub")
  }
  expect_warning(
    res <- MaxMin:::.FarFirstThin(2L, m, d = dat$d, points = NULL,
                                  RunOnSubset = stub, label = "stub"),
    "thinned 20 candidates to 8"
  )
  # Returned indices are original-space (members of the coreset), sorted, and
  # NOT the coreset-local positions.
  expect_setequal(as.integer(res), core[c(1L, 3L)])
  expect_identical(as.integer(res), sort(core[c(1L, 3L)]))
  expect_true(all(as.integer(res) %in% core))
  expect_false(identical(as.integer(res), c(1L, 3L)))
  # Class and the (original-space-valid) attributes are carried across.
  expect_s3_class(res, "MaxMinSelection")
  expect_identical(attr(res, "score"), 42)
  expect_identical(attr(res, "secondary"), 7)
  expect_identical(attr(res, "producer"), "stub")
})

test_that(".FarFirstThin draws no session RNG (deterministic coreset)", {
  dat <- MakeData()
  stub <- function(d, points) {
    structure(c(1L, 2L), score = 0, class = "MaxMinSelection", producer = "stub")
  }
  set.seed(123)
  before <- .Random.seed
  suppressWarnings(
    MaxMin:::.FarFirstThin(2L, 6L, d = dat$d, points = NULL,
                           RunOnSubset = stub, label = "stub")
  )
  # The peripheral seed + the no-RNG stub must leave the stream untouched.
  expect_identical(.Random.seed, before)
})

# ---- DropAdd integration ----------------------------------------------------

test_that("DropAdd thins when maxCandidates binds (matrix path)", {
  dat <- MakeData()
  core <- as.integer(FarFirst(8L, d = dat$d, strategy = "peripheral"))
  expect_warning(r <- DropAdd(3L, d = dat$d, maxCandidates = 8L), "thinned")
  expect_length(r, 3L)
  expect_true(all(as.integer(r) %in% core))             # original-space, in coreset
  expect_equal(attr(r, "score"), MinDist(d = dat$d, r)) # score valid in original space
  expect_s3_class(r, "MaxMinSelection")
  expect_identical(attr(r, "producer"), "DropAdd")
})

test_that("DropAdd thins when maxCandidates binds (points path)", {
  dat <- MakeData()
  expect_warning(r <- DropAdd(3L, points = dat$pts, maxCandidates = 8L), "thinned")
  expect_length(r, 3L)
  expect_equal(attr(r, "score"), MinDist(points = dat$pts, idx = r))
})

test_that("DropAdd maxCandidates is a no-op when off or non-binding", {
  dat <- MakeData()
  expect_no_warning(off  <- DropAdd(3L, d = dat$d, maxCandidates = 0L))
  expect_no_warning(inf  <- DropAdd(3L, d = dat$d, maxCandidates = Inf))
  expect_no_warning(over <- DropAdd(3L, d = dat$d, maxCandidates = 100L)) # >= n = 20
  # Default cap (46340) is far above n, so the default call never thins.
  expect_no_warning(def  <- DropAdd(3L, d = dat$d))
  expect_identical(as.integer(off), as.integer(over))
  expect_identical(as.integer(off), as.integer(def))
})

test_that("DropAdd errors when maxCandidates < k", {
  dat <- MakeData()
  expect_error(DropAdd(5L, d = dat$d, maxCandidates = 3L), ">= k")
})

test_that("DropAdd thinned result is deterministic and original-space", {
  dat <- MakeData()
  a <- suppressWarnings(DropAdd(3L, d = dat$d, maxCandidates = 8L))
  b <- suppressWarnings(DropAdd(3L, d = dat$d, maxCandidates = 8L))
  # Selection + score are deterministic (the volatile `time_s` attribute aside).
  expect_identical(as.integer(a), as.integer(b))
  expect_equal(attr(a, "score"), attr(b, "score"))
  # The selection is a subset of the original-space coreset.
  core <- as.integer(FarFirst(8L, d = dat$d, strategy = "peripheral"))
  expect_true(all(as.integer(a) %in% core))
})

test_that("DropAdd points and matrix paths share the same coreset", {
  dat <- MakeData()
  core_d <- as.integer(FarFirst(8L, d = dat$d, strategy = "peripheral"))
  core_p <- as.integer(FarFirst(8L, points = dat$pts, strategy = "peripheral"))
  expect_identical(core_d, core_p)                      # coreset agrees bit-for-bit
  rd <- suppressWarnings(DropAdd(3L, d = dat$d, maxCandidates = 8L))
  rp <- suppressWarnings(DropAdd(3L, points = dat$pts, maxCandidates = 8L))
  # The subproblem differs only in representation; scores agree within the
  # documented DropAdd points-vs-matrix construction-seed tolerance (5%).
  expect_gte(attr(rp, "score"), attr(rd, "score") * 0.95)
  expect_gte(attr(rd, "score"), attr(rp, "score") * 0.95)
})

# ---- Grasp integration ------------------------------------------------------

test_that("Grasp thins when maxCandidates binds", {
  dat <- MakeData()
  core <- as.integer(FarFirst(8L, d = dat$d, strategy = "peripheral"))
  set.seed(5)
  expect_warning(r <- Grasp(3L, dat$d, plateau = 20L, maxCandidates = 8L), "thinned")
  expect_length(r, 3L)
  expect_true(all(as.integer(r) %in% core))
  expect_equal(attr(r, "score"), MinDist(d = dat$d, r))
  expect_identical(attr(r, "producer"), "Grasp")
})

test_that("Grasp maxCandidates off/non-binding and < k", {
  dat <- MakeData()
  expect_no_warning(Grasp(3L, dat$d, plateau = 20L, maxCandidates = 0L))
  expect_no_warning(Grasp(3L, dat$d, plateau = 20L))     # default 2000 >> n
  expect_error(Grasp(5L, dat$d, maxCandidates = 3L), ">= k")
})

test_that("Grasp thinned result is reproducible under a fixed seed", {
  dat <- MakeData()
  set.seed(5); a <- suppressWarnings(Grasp(3L, dat$d, plateau = 20L, maxCandidates = 8L))
  set.seed(5); b <- suppressWarnings(Grasp(3L, dat$d, plateau = 20L, maxCandidates = 8L))
  expect_identical(as.integer(a), as.integer(b))
  expect_equal(attr(a, "score"), attr(b, "score"))
})

# ---- degenerate distances (zeros / ties) -----------------------------------

test_that("thinning tolerates degenerate (tied/zero) distances", {
  set.seed(1)
  pts <- matrix(rnorm(30), ncol = 3)
  pts <- rbind(pts, pts)                 # 20 points: every point is duplicated
  d <- as.matrix(dist(pts))              # carries exact zeros and ties
  # Coreset construction and the m x m subproblem must not error on zeros/ties.
  expect_no_error(suppressWarnings(DropAdd(3L, d = d, maxCandidates = 8L)))
  expect_no_error(suppressWarnings(DropAdd(3L, points = pts, maxCandidates = 8L)))
  expect_no_error(suppressWarnings(Grasp(3L, d, plateau = 20L, maxCandidates = 8L)))
})

# ---- ExactMaxMin: internal warm-starts are pinned to no thinning ------------

test_that("ExactMaxMin does not emit a thinning warning on a small instance", {
  skip_if_not_installed("highs")
  dat <- MakeData(N = 12)
  # n = 12 is well below both caps; the maxCandidates = 0L pins on the internal
  # Grasp/DropAdd warm-starts keep it that way regardless of their defaults.
  expect_no_warning(ExactMaxMin(3L, dat$d, maxSeconds = 30))
})

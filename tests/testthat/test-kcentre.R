# Discrete k-centre: covering-radius score, CDSh heuristic, exact covering IP.

# Brute-force optimum covering radius over all k-subsets (tiny n only).
brute_kcentre <- function(d, k) {
  n <- nrow(d)
  combs <- utils::combn(n, k)
  radii <- apply(combs, 2L, function(idx) {
    max(apply(d[, idx, drop = FALSE], 1L, min))
  })
  min(radii)
}

# Covering radius of a centre set, the slow-but-obvious reference.
ref_radius <- function(d, idx) {
  max(apply(d[, idx, drop = FALSE], 1L, min))
}

# ----- KCentreRadius --------------------------------------------------------

test_that("KCentreRadius matches the brute-force covering radius", {
  set.seed(11)
  pts <- matrix(rnorm(40), ncol = 2)
  d <- as.matrix(stats::dist(pts))
  for (idx in list(3L, c(1L, 5L), c(2L, 8L, 15L, 19L))) {
    expect_equal(KCentreRadius(d, idx), ref_radius(d, idx))
  }
})

test_that("KCentreRadius matrix and coordinate paths agree on Euclidean data", {
  set.seed(12)
  pts <- matrix(rnorm(90), ncol = 3)
  d <- stats::dist(pts)
  idx <- c(1L, 7L, 13L, 22L)
  expect_equal(KCentreRadius(d, idx), KCentreRadius(points = pts, idx = idx))
})

test_that("KCentreRadius accepts a single centre and a dist object", {
  set.seed(13)
  pts <- matrix(rnorm(30), ncol = 2)
  d <- stats::dist(pts)
  dm <- as.matrix(d)
  expect_equal(KCentreRadius(d, 4L), max(dm[, 4L]))
})

test_that("KCentreRadius rejects NA, duplicate and empty index sets", {
  d <- as.matrix(stats::dist(matrix(rnorm(20), ncol = 2)))
  expect_error(KCentreRadius(d, c(1L, NA_integer_)), "NA")
  expect_error(KCentreRadius(d, c(2L, 2L)), "duplicate")
  expect_error(KCentreRadius(d, integer(0)), "at least one")
})

# ----- KCentre (CDSh) -------------------------------------------------------

test_that("KCentre returns a valid centre set whose radius attribute is exact", {
  set.seed(21)
  d <- as.matrix(stats::dist(matrix(rnorm(120), ncol = 2)))
  res <- KCentre(d, 5L)
  expect_s3_class(res, "KCentreSelection")
  expect_lte(length(res), 5L)
  expect_false(anyDuplicated(as.integer(res)) > 0L)
  # The attached radius is the true covering radius of the returned centres.
  expect_equal(attr(res, "radius"), KCentreRadius(d, as.integer(res)))
})

test_that("CDSh covers at least as tightly as Gonzalez (the SOTA-beat)", {
  # Across several clustered instances CDSh never does worse than a Gonzalez
  # (farthest-first) pass, and is strictly better on average.
  cdsh <- numeric(0)
  gonz <- numeric(0)
  for (s in 1:6) {
    set.seed(100 + s)
    centres <- matrix(rnorm(8), ncol = 2) * 6
    pts <- centres[rep(1:4, each = 30), ] + matrix(rnorm(240), ncol = 2)
    d <- as.matrix(stats::dist(pts))
    cdsh <- c(cdsh, KCentreRadius(d, as.integer(KCentre(d, 6L))))
    gonz <- c(gonz, KCentreRadius(d, as.integer(FarFirst(d, 6L, method = "peripheral"))))
  }
  expect_true(all(cdsh <= gonz + 1e-9))   # never worse
  expect_lt(mean(cdsh), mean(gonz))        # strictly better on average
})

test_that("KCentre handles the k = n and k = 1 boundaries", {
  d <- as.matrix(stats::dist(matrix(rnorm(20), ncol = 2)))
  n <- nrow(d)
  all_centres <- KCentre(d, n)
  expect_equal(sort(as.integer(all_centres)), seq_len(n))
  expect_equal(attr(all_centres, "radius"), 0)

  one <- KCentre(d, 1L)
  expect_lte(length(one), 1L)
  expect_equal(attr(one, "radius"), KCentreRadius(d, as.integer(one)))
})

test_that("KCentre rejects bad k and non-finite distances", {
  d <- as.matrix(stats::dist(matrix(rnorm(20), ncol = 2)))
  expect_error(KCentre(d, nrow(d) + 1L), "1 <= k")
  expect_error(KCentre(d, 0L), "positive integer")
  d[1L, 2L] <- Inf
  expect_error(KCentre(d, 3L), "NA/NaN/Inf")
})

# ----- ExactKCentre ---------------------------------------------------------

test_that("ExactKCentre proves the brute-force optimum on small instances", {
  skip_if_not_installed("highs")
  skip_if_not_installed("Matrix")
  for (seed in c(31L, 32L, 33L)) {
    set.seed(seed)
    pts <- matrix(rnorm(20), ncol = 2)          # n = 10
    d <- as.matrix(stats::dist(pts))
    for (k in c(1L, 2L, 3L)) {
      res <- ExactKCentre(d, k, progress = FALSE)
      expect_s3_class(res, "KCentreExact")
      expect_true(res$proven)
      expect_lte(res$n_centres, k)
      expect_equal(res$radius, brute_kcentre(d, k))
      # Returned centres genuinely achieve the proven radius.
      expect_equal(KCentreRadius(d, res$indices), res$radius)
    }
  }
})

test_that("ExactKCentre is no worse than CDSh, which is no worse than optimal", {
  skip_if_not_installed("highs")
  skip_if_not_installed("Matrix")
  set.seed(41)
  d <- as.matrix(stats::dist(matrix(rnorm(24), ncol = 2)))   # n = 12
  k <- 3L
  exact <- ExactKCentre(d, k, progress = FALSE)
  heur <- KCentreRadius(d, as.integer(KCentre(d, k)))
  expect_lte(exact$radius, heur + 1e-9)       # exact optimum <= heuristic
  expect_equal(exact$radius, brute_kcentre(d, k))
})

test_that("ExactKCentre handles k = n trivially", {
  skip_if_not_installed("highs")
  skip_if_not_installed("Matrix")
  d <- as.matrix(stats::dist(matrix(rnorm(16), ncol = 2)))
  res <- ExactKCentre(d, nrow(d), progress = FALSE)
  expect_true(res$proven)
  expect_equal(res$radius, 0)
})

# ----- US-spelling aliases --------------------------------------------------

test_that("US-spelling aliases are identical bindings to the UK functions", {
  expect_identical(KCenter, KCentre)
  expect_identical(KCenterRadius, KCentreRadius)
  expect_identical(ExactKCenter, ExactKCentre)
})

test_that("US-spelling aliases give the same results as the UK functions", {
  set.seed(51)
  d <- as.matrix(stats::dist(matrix(rnorm(80), ncol = 2)))
  expect_equal(as.integer(KCenter(d, 4L)), as.integer(KCentre(d, 4L)))
  idx <- c(1L, 5L, 9L)
  expect_equal(KCenterRadius(d, idx), KCentreRadius(d, idx))

  skip_if_not_installed("highs")
  skip_if_not_installed("Matrix")
  d10 <- d[1:10, 1:10]
  expect_equal(ExactKCenter(d10, 3L, progress = FALSE)$radius,
               ExactKCentre(d10, 3L, progress = FALSE)$radius)
})

# ----- red-team regressions (Round 8: KC-001..007) --------------------------

test_that("KCentre is never worse than Gonzalez on small random instances (KC-001)", {
  # The binary search's feasibility predicate is non-monotone and could land
  # ABOVE the Gonzalez 2-approximation; the Gonzalez floor now forbids that.
  for (seed in 1:25) {
    set.seed(seed)
    n <- sample(8:50, 1L)
    d <- as.matrix(stats::dist(matrix(rnorm(n * 2L), ncol = 2L)))
    k <- sample(2:min(8L, n - 1L), 1L)
    cdsh <- KCentreRadius(d, as.integer(KCentre(d, k)))
    gonz <- KCentreRadius(d, as.integer(FarFirst(d, k, method = "peripheral")))
    expect_lte(cdsh, gonz + 1e-9)
  }
})

test_that("KCentre never exceeds twice the optimum (restored 2-approx; KC-001)", {
  skip_if_not_installed("highs")
  skip_if_not_installed("Matrix")
  for (seed in 1:6) {
    set.seed(200L + seed)
    d <- as.matrix(stats::dist(matrix(rnorm(40L), ncol = 2L)))   # n = 20 (exhaustive)
    opt <- ExactKCentre(d, 4L, progress = FALSE)$radius
    heur <- KCentreRadius(d, as.integer(KCentre(d, 4L)))
    if (opt > 0) expect_lte(heur, 2 * opt + 1e-9)
  }
})

test_that("KCentre and ExactKCentre reject asymmetric distance matrices (KC-002)", {
  d <- as.matrix(stats::dist(matrix(rnorm(20L), ncol = 2L)))
  d[1L, 2L] <- d[1L, 2L] + 1                       # break symmetry (still finite)
  expect_error(KCentre(d, 3L), "symmetric")
  skip_if_not_installed("highs")
  skip_if_not_installed("Matrix")
  expect_error(ExactKCentre(d, 3L, progress = FALSE), "symmetric")
})

test_that("KCentreRadius rejects out-of-range indices on both paths (KC-003)", {
  set.seed(71)
  pts <- matrix(rnorm(30L), ncol = 2L)
  d <- as.matrix(stats::dist(pts))
  expect_error(KCentreRadius(d, 0L), "1, nrow")
  expect_error(KCentreRadius(d, nrow(d) + 1L), "1, nrow")
  expect_error(KCentreRadius(points = pts, idx = 0L), "1, nrow")
})

test_that("ExactKCentre returns a consistent unproven incumbent under a tiny budget (KC-004, KC-005)", {
  skip_if_not_installed("highs")
  skip_if_not_installed("Matrix")
  set.seed(61)
  d <- as.matrix(stats::dist(matrix(rnorm(60L), ncol = 2L)))     # n = 30
  res <- ExactKCentre(d, 4L, maxSeconds = 1e-6, progress = FALSE)
  expect_false(res$proven)
  expect_lte(res$n_centres, 4L)
  # The reported radius is the true covering radius of the returned centres,
  # proven or not (KC-004).
  expect_equal(res$radius, KCentreRadius(d, res$indices))
})

test_that("KCentre collapses duplicate centres with a consistent radius (KC-007)", {
  d <- matrix(1, 8L, 8L); diag(d) <- 0            # all pairwise distances equal
  res <- KCentre(d, 3L)
  expect_lte(length(res), 3L)
  expect_false(anyDuplicated(as.integer(res)) > 0L)
  expect_equal(attr(res, "radius"), KCentreRadius(d, as.integer(res)))
})

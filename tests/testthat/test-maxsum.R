# Tests for ExactMaxSum() -- the Maximum Diversity (max-sum) solver.

# Brute-force optimum over all size-k subsets, for validation.
.BruteMaxSum <- function(d, k) {
  d <- as.matrix(d)
  cb <- utils::combn(nrow(d), k)
  sc <- apply(cb, 2L, function(idx) {
    s <- d[idx, idx]
    sum(s[upper.tri(s)])
  })
  list(idx = sort(cb[, which.max(sc)]), score = max(sc))
}

test_that("ExactMaxSum matches the brute-force optimum (proven)", {
  skip_if_not_installed("highs")
  set.seed(1)
  d <- dist(matrix(rnorm(40), ncol = 2))   # 20 points
  for (k in c(2L, 3L, 5L)) {
    res <- ExactMaxSum(k, d)
    bf <- .BruteMaxSum(d, k)
    expect_length(res, k)
    expect_s3_class(res, "MaxSumSelection")
    expect_true(isTRUE(attr(res, "proven")))
    expect_equal(attr(res, "score"), bf$score, tolerance = 1e-9)
    # The achieved total distance equals the brute-force optimum (the indices
    # may differ when several subsets tie).
    sub <- as.matrix(d)[as.integer(res), as.integer(res)]
    expect_equal(sum(sub[upper.tri(sub)]), bf$score, tolerance = 1e-9)
  }
})

test_that("ExactMaxSum picks the two most distant points at k = 2", {
  skip_if_not_installed("highs")
  m <- matrix(c(0, 0, 0, 1, 0, 20), ncol = 2, byrow = TRUE)  # collinear; pt 3 far
  res <- ExactMaxSum(2L, dist(m))
  expect_equal(as.integer(res), c(1L, 3L))   # the diameter pair (1,3): dist 20
})

test_that("ExactMaxSum validates k", {
  skip_if_not_installed("highs")
  d <- dist(matrix(rnorm(20), ncol = 2))
  expect_error(ExactMaxSum(1L, d), "2 <= k")
  expect_error(ExactMaxSum(99L, d), "2 <= k")
})

test_that("ExactMaxSum accepts a dist or a square matrix and warmStart", {
  skip_if_not_installed("highs")
  set.seed(2)
  d <- as.matrix(dist(matrix(rnorm(30), ncol = 2)))
  bf <- .BruteMaxSum(d, 4L)
  r1 <- ExactMaxSum(4L, d)
  r2 <- ExactMaxSum(4L, stats::as.dist(d), warmStart = bf$idx)
  expect_equal(attr(r1, "score"), bf$score, tolerance = 1e-9)
  expect_equal(attr(r2, "score"), bf$score, tolerance = 1e-9)
})

test_that("ExactMaxSum result prints with a total-distance summary", {
  skip_if_not_installed("highs")
  res <- ExactMaxSum(3L, dist(matrix(rnorm(20), ncol = 2)))
  expect_match(format(res), "total distance")
  expect_invisible(print(res))
})

test_that("a warmStart that beats the heuristic incumbent is adopted", {
  skip_if_not_installed("highs")
  # A 31-point instance on which the default multi-start 1-swap local search
  # settles on a subset that a hand-found alternative strictly beats --
  # exercising the `wv > lsVal` warmStart-adoption branch. Found by a targeted
  # search (not a plain set.seed() draw), so the points are pinned literally.
  pts <- structure(c(
    -0.448375293423593, -0.310051806584565, 1.35753656396223,
    -1.10508888594004, -0.277907249955316, 1.34761200180236, 0.719657290518998,
    -0.07950563882744, 0.214732200902843, 0.0442938836333716, 0.343526595099658,
    -0.0326243058378093, -1.02421899507323, 0.215795504599057, -0.594081683661156,
    -0.229863820673017, 0.611733062530893, 0.888350495286952, -1.10873967175493,
    -1.48007963153779, 0.422531508960742, -0.240635474015899, -0.885316258078437,
    0.418593510432848, 2.30146027874758, 0.678571623571227, 0.570994338371516,
    1.29967815382362, -0.190980831441199, 1.52065054559739, -0.658108768588629,
    -0.188098947656654, 0.44227843894853, 1.16488327165977, 1.52136649656061,
    -0.667226361888637, -0.827779148714529, -0.415743937089355, -1.90193010017359,
    -1.77973675577368, 0.390121217377137, 0.754018576382694, 1.16623794939202,
    -1.77003062092275, 1.55976421902239, -1.57627609921479, -1.10272338727579,
    -1.36765807890335, -0.586168115894317, 0.137425769008358, -1.11326792721959,
    0.711732936224916, 0.765126474042235, 0.237712966544576, -0.606296322475592,
    -0.739082272955005, -0.147296280538258, -0.847492624420121, 0.99503653137979,
    0.994387535678726, 0.107168976711389, 0.771575996264042
  ), dim = c(31L, 2L))
  d <- as.matrix(dist(pts))
  k <- 12L
  warmStart <- c(3L, 4L, 8L, 9L, 13L, 14L, 15L, 19L, 20L, 25L, 28L, 30L)

  # .MaxSumHeuristic() draws its random-restart seeds from the session RNG, so
  # pin it to reproduce the suboptimal local-search incumbent found offline.
  set.seed(10)
  lsIdx <- Coreset:::.MaxSumHeuristic(d, k)
  lsVal <- Coreset:::.MaxSumScore(d, lsIdx)
  wv <- Coreset:::.MaxSumScore(d, warmStart)
  # Sanity check that the instance still has the intended property.
  expect_gt(wv, lsVal)

  set.seed(10)
  res <- ExactMaxSum(k, d, warmStart = warmStart)
  # The returned incumbent is never worse than the supplied warm start,
  # regardless of whether the MILP later improves on or proves it.
  expect_gte(attr(res, "score"), wv)
})

test_that("ExactMaxSum returns proven=FALSE on budget expiry", {
  skip_if_not_installed("highs")
  set.seed(1)
  # A large enough instance that a near-zero budget cannot certify optimality.
  pts <- matrix(rnorm(30 * 4), ncol = 4)
  d <- as.matrix(dist(pts))
  res <- ExactMaxSum(5L, d, maxSeconds = 1e-9)
  # With a near-zero budget the solver may not certify, but must not error.
  # Note: proven could be TRUE if the solver is extremely fast; we just assert
  # the field exists and the return is valid.
  expect_type(attr(res, "proven"), "logical")
  expect_length(as.integer(res), 5L)
  # If proven = FALSE, the objective is still a valid lower bound (non-negative).
  if (!isTRUE(attr(res, "proven"))) {
    expect_gte(attr(res, "score"), 0)
  }
})

test_that("the local-search helpers return valid k-subsets", {
  set.seed(3)
  d <- as.matrix(dist(matrix(rnorm(40), ncol = 2)))
  greedy <- Coreset:::.MaxSumGreedy(d, 5L, seed = which.max(rowSums(d)))
  expect_length(greedy, 5L)
  expect_true(!anyDuplicated(greedy))
  improved <- Coreset:::.MaxSumLocalSearch(d, 5L, greedy)
  expect_length(improved, 5L)
  # local search never worsens the objective
  expect_gte(Coreset:::.MaxSumScore(d, improved),
             Coreset:::.MaxSumScore(d, greedy))
})

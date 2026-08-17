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

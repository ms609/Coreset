# Tests for MaxEntropy() -- maximum-entropy (maxdet) subset selection.

test_that("MaxEntropy returns a well-formed selection", {
  set.seed(1)
  d <- dist(matrix(rnorm(40), ncol = 2))
  sel <- MaxEntropy(4L, d)

  expect_s3_class(sel, "MaxEntropySelection")
  expect_length(sel, 4L)
  expect_type(as.integer(sel), "integer")
  expect_identical(as.integer(sel), sort(as.integer(sel)))   # returned sorted
  expect_true(all(sel >= 1L & sel <= 20L))
  expect_false(anyDuplicated(sel) > 0L)
  expect_false(is.null(attr(sel, "logDet")))
  expect_equal(attr(sel, "score"), attr(sel, "logDet"))
  expect_true(attr(sel, "negMass") >= 0)
})

test_that("MaxEntropy accepts a dist or a square matrix identically", {
  set.seed(2)
  d <- dist(matrix(rnorm(30), ncol = 3))
  expect_identical(as.integer(MaxEntropy(3L, d)),
                   as.integer(MaxEntropy(3L, as.matrix(d))))
})

test_that("MaxEntropy is deterministic (no random seed)", {
  set.seed(3)
  d <- dist(matrix(rnorm(60), ncol = 2))
  a <- MaxEntropy(5L, d)
  b <- MaxEntropy(5L, d)
  expect_identical(as.integer(a), as.integer(b))
})

test_that("the greedy C++ selector matches the R pivoted-Cholesky reference", {
  # Independent R port of the greedy (max residual conditional variance).
  greedyRef <- function(kp, k, seed) {
    n <- nrow(kp); dd <- diag(kp); L <- matrix(0, n, k)
    avail <- rep(TRUE, n); perm <- integer(k)
    for (t in seq_len(k)) {
      j <- if (t == 1L) seed else which.max(ifelse(avail, dd, -Inf))
      perm[t] <- j; avail[j] <- FALSE
      ljt <- sqrt(max(dd[j], 0))
      L[j, t] <- ljt
      rows <- which(avail)
      if (length(rows) && ljt > 0) {
        prev <- if (t > 1L) L[rows, seq_len(t - 1L), drop = FALSE] %*%
          L[j, seq_len(t - 1L)] else 0
        L[rows, t] <- (kp[rows, j] - prev) / ljt
        dd[rows] <- pmax(dd[rows] - L[rows, t]^2, 0)
      }
    }
    perm
  }
  set.seed(4)
  d <- as.matrix(dist(matrix(rnorm(50), ncol = 2)))
  kp <- .MaxEntropyRepair(.MaxEntropyKernel(d), "clip")
  seed <- which.min(rowSums(kp))
  for (k in 2:5) {
    expect_identical(MaxEntropyGreedy_cpp(kp, k, as.integer(seed)),
                     as.integer(greedyRef(kp, k, seed)))
  }
})

# Independent R Cholesky log-determinant (mirrors src/maxentropy.cpp SubLogDet):
# -Inf for a non-positive-definite block.
cholLogDet <- function(kp, idx) {
  sub <- (kp[idx, idx, drop = FALSE] + t(kp[idx, idx, drop = FALSE])) / 2
  r <- tryCatch(chol(sub), error = function(e) NULL)
  if (is.null(r)) -Inf else 2 * sum(log(diag(r)))
}

test_that("the exact C++ enumeration matches a brute-force log-det argmax", {
  set.seed(5)
  d <- as.matrix(dist(matrix(rnorm(24), ncol = 2)))   # n = 12
  kp <- .MaxEntropyRepair(.MaxEntropyKernel(d), "clip")
  bruteForce <- function(kp, k) {
    combos <- utils::combn(nrow(kp), k)
    vals <- apply(combos, 2L, function(idx) cholLogDet(kp, idx))
    sort(combos[, which.max(vals)])
  }
  for (k in 2:4) {
    expect_identical(as.integer(MaxEntropyExact_cpp(kp, k)),
                     as.integer(bruteForce(kp, k)))
  }
})

test_that("selection and reported log det use one consistent metric (ME-003)", {
  set.seed(8)
  d <- as.matrix(dist(matrix(rnorm(30), ncol = 2)))
  kp <- .MaxEntropyRepair(.MaxEntropyKernel(d), "clip")
  for (k in 2:5) {
    sel <- MaxEntropy(k, d)
    expect_equal(attr(sel, "logDet"), cholLogDet(kp, as.integer(sel)))
    expect_equal(MaxEntropyLogDet_cpp(kp, as.integer(sel)),
                 cholLogDet(kp, as.integer(sel)))
  }
})

test_that("out-of-range seed does not crash the greedy core (ME-001)", {
  # A seed beyond [1, n] must be ignored (fall back to the diagonal argmax),
  # not indexed out of bounds -- which previously segfaulted the R session.
  K <- diag(5)
  expect_length(MaxEntropyGreedy_cpp(K, 2L, 99L), 2L)
  expect_length(MaxEntropyGreedy_cpp(K, 2L, 0L), 2L)
  expect_length(MaxEntropyGreedy_cpp(K, 2L, -3L), 2L)
})

test_that("k beyond the distinct-point count warns and reports -Inf (ME-002)", {
  base <- rbind(c(0, 0), c(10, 0), c(0, 10))      # 3 distinct points
  d <- dist(rbind(base, base))                     # n = 6, only 3 distinct
  expect_warning(sel <- MaxEntropy(4L, d), "distinct")
  expect_true(is.infinite(attr(sel, "logDet")))
  expect_length(sel, 4L)
  # k within the distinct count: no warning, finite log det.
  expect_no_warning(sel3 <- MaxEntropy(3L, d))
  expect_true(is.finite(attr(sel3, "logDet")))
})

test_that("non-finite distances are rejected (ME-004)", {
  d <- as.matrix(dist(matrix(c(0, 0, 10, 0, 0, 10, 10, 10), ncol = 2, byrow = TRUE)))
  dNA <- d; dNA[1, 2] <- dNA[2, 1] <- NA
  expect_error(MaxEntropy(2L, dNA), "finite")
  dInf <- d; dInf[1, 2] <- dInf[2, 1] <- Inf
  expect_error(MaxEntropy(2L, dInf), "finite")
})

test_that("auto mode uses exact for small choose(n, k), greedy otherwise", {
  set.seed(6)
  d <- dist(matrix(rnorm(30), ncol = 2))              # n = 15
  expect_true(attr(MaxEntropy(3L, d), "exact"))       # choose(15, 3) = 455
  expect_false(attr(MaxEntropy(3L, d, maxCombos = 100), "exact"))
  # Forcing exact beyond the ceiling is an error, not a silent fallback.
  expect_error(MaxEntropy(3L, d, exact = TRUE, maxCombos = 100), "maxCombos")
})

test_that("MaxEntropy is density-blind: a duplicated point is never co-selected", {
  # Five well-separated points; duplicate point 1 several times. The maxdet
  # optimum must never take two zero-distance copies (det -> 0), so each
  # selection contains at most one member of the duplicate cluster.
  base <- rbind(c(0, 0), c(10, 0), c(0, 10), c(10, 10), c(5, 5))
  pts <- rbind(base, base[rep(1L, 4L), ])             # 4 extra copies of point 1
  d <- dist(pts)
  dupes <- c(1L, 6L, 7L, 8L, 9L)
  for (k in 2:4) {
    sel <- MaxEntropy(k, d, exact = FALSE)            # greedy
    expect_lte(sum(sel %in% dupes), 1L)
    sel2 <- MaxEntropy(k, d, exact = TRUE)            # exact
    expect_lte(sum(sel2 %in% dupes), 1L)
  }
})

test_that("MaxEntropy honours k edge cases and rejects bad k", {
  set.seed(7)
  d <- dist(matrix(rnorm(20), ncol = 2))              # n = 10
  expect_length(MaxEntropy(1L, d), 1L)
  expect_identical(as.integer(MaxEntropy(10L, d)), 1:10)   # k == n: all points
  expect_error(MaxEntropy(0L, d), "1 <= k")
  expect_error(MaxEntropy(11L, d), "1 <= k")
})

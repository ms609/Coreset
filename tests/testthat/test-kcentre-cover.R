# Tests for CoverDecide_cpp, the exhaustive covering decision behind
# ExactKCentre(). The kernel is exhaustive, so agreement with an oracle is
# exact -- "infeasible" asserts that NO cover of size <= k exists, which is
# stronger than any solver status.

Decide <- function(d, r, k, maxSeconds = 30) {
  Coreset:::CoverDecide_cpp(d, r, k, maxSeconds)
}

# Brute force: the true minimum number of centres covering everything within r.
BruteMinCover <- function(d, r) {
  n <- nrow(d)
  cov <- d <= r
  for (size in seq_len(n)) {
    combos <- utils::combn(n, size)
    for (ci in seq_len(ncol(combos))) {
      if (all(apply(cov[combos[, ci], , drop = FALSE], 2L, any))) {
        return(size)
      }
    }
  }
  n + 1L                                                            # nocov
}

ValidCover <- function(d, sel, r) {
  length(sel) >= 1L && !is.unsorted(sel, strictly = TRUE) &&
    all(sel >= 1L) && all(sel <= nrow(d)) &&
    KCentreRadius(d = d, sel) <= r
}

test_that("CoverDecide agrees exactly with brute force over random probes", {
  probes <- 0L
  for (seed in 1:12) {
    set.seed(seed)
    n <- 9L
    d <- as.matrix(stats::dist(matrix(stats::rnorm(2L * n), ncol = 2L)))
    radii <- sort(unique(d[upper.tri(d)]))
    for (r in radii[seq(1L, length(radii), length.out = 6L)]) {
      truth <- BruteMinCover(d, r)
      for (k in 1:4) {
        got <- Decide(d, r, k, 30)
        expect_identical(got$status,
                         if (truth <= k) "feasible" else "infeasible")
        if (identical(got$status, "feasible")) {
          expect_lte(length(got$witness), k)
          expect_true(ValidCover(d, got$witness, r))
        }
        probes <- probes + 1L
      }
    }
  }
  expect_gt(probes, 200L)
})

test_that("CoverDecide agrees with a set-cover IP past brute force's reach", {
  skip_if_not_installed("highs")
  skip_if_not_installed("Matrix")
  set.seed(77)
  n <- 30L
  d <- as.matrix(stats::dist(matrix(stats::rnorm(2L * n), ncol = 2L)))
  radii <- sort(unique(d[upper.tri(d)]))
  IpMin <- function(r) {
    cover <- which(d <= r, arr.ind = TRUE)
    A <- Matrix::sparseMatrix(i = cover[, 2L], j = cover[, 1L], x = 1,
                              dims = c(n, n))
    res <- highs::highs_solve(
      L = rep.int(1, n), lower = rep.int(0, n), upper = rep.int(1, n), A = A,
      lhs = rep.int(1, n), rhs = rep.int(Inf, n), types = rep.int("I", n),
      maximum = FALSE, control = list(threads = 1L, time_limit = 60))
    length(which(res$primal_solution > 0.5))
  }
  for (r in radii[seq(2L, length(radii), length.out = 8L)]) {
    truth <- IpMin(r)
    for (k in c(1L, truth - 1L, truth, truth + 1L)) {
      if (k < 1L) next
      got <- Decide(d, r, k, 60)
      expect_identical(got$status,
                       if (truth <= k) "feasible" else "infeasible")
    }
  }
})

test_that("unit propagation forces centres and can refute on its own", {
  # Three far-apart pairs: each point's only cover within r is its own pair,
  # so three centres are forced and k < 3 is refuted without any search.
  d <- as.matrix(stats::dist(matrix(c(0, 0, 0.1, 0, 100, 0, 100.1, 0,
                                      200, 0, 200.1, 0), ncol = 2L,
                                    byrow = TRUE)))
  r <- 0.5
  expect_identical(Decide(d, r, 3L)$status, "feasible")
  expect_identical(Decide(d, r, 2L)$status, "infeasible")
  expect_identical(Decide(d, r, 1L)$status, "infeasible")
  expect_length(Decide(d, r, 3L)$witness, 3L)
  # Refutation by propagation alone visits no search node.
  expect_identical(Decide(d, r, 2L)$nodes, 0)
})

test_that("dominance keeps exactly one of a set of identical points", {
  # Four coincident points plus one far away: the duplicates have identical
  # neighbourhoods, so the tie-break must retain one, not none and not all.
  pts <- rbind(matrix(0, nrow = 4L, ncol = 2L), c(10, 0))
  d <- as.matrix(stats::dist(pts))
  expect_identical(Decide(d, 0.5, 2L)$status, "feasible")
  expect_identical(Decide(d, 0.5, 1L)$status, "infeasible")
  got <- Decide(d, 0.5, 2L)
  expect_true(ValidCover(d, got$witness, 0.5))
  # A radius covering everything needs one centre, and any point serves.
  expect_identical(Decide(d, 10, 1L)$status, "feasible")
})

test_that("independent components are budgeted separately", {
  # Two identical 5-point rings 1000 apart: each needs 2 centres at this
  # radius, so the instance needs 4 and refutes 3.
  ring <- function(shift) {
    a <- seq(0, 2 * pi, length.out = 6L)[-6L]
    cbind(cos(a) + shift, sin(a))
  }
  d <- as.matrix(stats::dist(rbind(ring(0), ring(1000))))
  r <- 1.2
  one <- BruteMinCover(as.matrix(stats::dist(ring(0))), r)
  expect_identical(Decide(d, r, 2L * one)$status, "feasible")
  expect_identical(Decide(d, r, 2L * one - 1L)$status, "infeasible")
  expect_true(ValidCover(d, Decide(d, r, 2L * one)$witness, r))
})

# Two probes that genuinely need the search, with their radii pinned by index
# so the tests do not pay for a bisection. Both instances come from the
# Mersenne-Twister stream, so the ordering of the distinct distances -- and
# hence the index -- is the same everywhere.
Radii <- function(d) sort(unique(d[upper.tri(d)]))

test_that("CoverDecide stops inside the search when the budget expires", {
  # 200 points in five dimensions: k = 14 cannot cover at this radius, and the
  # exhaustive proof runs to ~10^8 nodes, so the deadline lands inside the
  # search rather than inside the reduction.
  set.seed(2)
  d <- as.matrix(stats::dist(matrix(stats::runif(1000L), ncol = 5L)))
  r <- Radii(d)[[2079L]]                  # one below the k = 14 optimum
  got <- Decide(d, r, 14L, 0.2)           # ~250x short of what the proof needs
  expect_identical(got$status, "inconclusive")
  expect_gt(got$nodes, 0)                 # it did search before giving up
})

test_that("CoverDecide reports inconclusive before it starts", {
  set.seed(9)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(400L), ncol = 2L)))
  # No budget at all: the deadline is past before the first reduction pass.
  got <- Decide(d, Radii(d)[[50L]], 12L, -1)
  expect_identical(got$status, "inconclusive")
  expect_identical(got$nodes, 0)
})

test_that("a refutation the reduction cannot make runs and is exhaustive", {
  # Below the optimum by one radius: the reduction leaves a real instance and
  # the search has to prove the refutation itself.
  set.seed(1)
  d <- as.matrix(stats::dist(matrix(stats::runif(500L), ncol = 2L)))
  radii <- Radii(d)
  deep <- Decide(d, radii[[3568L]], 10L, 60)
  expect_identical(deep$status, "infeasible")
  expect_gt(deep$nodes, 1000)
  # The next radius up is the optimum, and its witness covers.
  at <- Decide(d, radii[[3569L]], 10L, 60)
  expect_identical(at$status, "feasible")
  expect_true(ValidCover(d, at$witness, radii[[3569L]]))
})

test_that("CoverDecide is deterministic", {
  set.seed(5)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(60L), ncol = 2L)))
  radii <- sort(unique(d[upper.tri(d)]))
  r <- radii[[round(length(radii) / 3)]]
  a <- Decide(d, r, 4L)
  b <- Decide(d, r, 4L)
  expect_identical(a$status, b$status)
  expect_identical(a$witness, b$witness)
  expect_identical(a$nodes, b$nodes)
})

test_that("degenerate radii are decided without search", {
  set.seed(3)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(30L), ncol = 2L)))
  n <- nrow(d)
  # Every point covered by any single centre.
  expect_identical(Decide(d, max(d), 1L)$status, "feasible")
  # Below every off-diagonal distance each point covers only itself, so n
  # centres are needed and n - 1 is refuted.
  tiny <- min(d[d > 0]) / 2
  expect_identical(Decide(d, tiny, n)$status, "feasible")
  expect_identical(Decide(d, tiny, n - 1L)$status, "infeasible")
})

test_that("ExactKCentre needs neither highs nor Matrix", {
  set.seed(11)
  d <- as.matrix(stats::dist(matrix(stats::rnorm(30L), ncol = 2L)))
  res <- ExactKCentre(3L, d)
  expect_true(attr(res, "proven"))
  expect_identical(attr(res, "solver"), "cover search")
  expect_equal(attr(res, "radius"), KCentreRadius(d = d, as.integer(res)))
})

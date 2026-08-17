# Tests for ThresholdDecide_cpp and the .MaxISVerdict probe it backs.
#
# Oracle: brute-force clique enumeration on tiny graphs. The kernel's contract
# is exact: "feasible" iff the complement graph H holds a k-clique, and a
# returned witness is one. A second oracle at sizes past brute force solves
# the equivalent packing IP with `highs`.

# Does the symmetric adjacency matrix contain a clique of size k?
.HasKClique <- function(hAdj, k) {
  n <- nrow(hAdj)
  cmb <- utils::combn(n, k)
  for (j in seq_len(ncol(cmb))) {
    s <- cmb[, j]
    if (all(hAdj[s, s][upper.tri(diag(k))])) {
      return(TRUE)
    }
  }
  FALSE
}

# H edges at one threshold, the probe input .MaxISVerdict takes.
.ProbeArgs <- function(d, lambda) {
  ut <- which(upper.tri(d))
  rc <- arrayInd(ut, dim(d))
  far <- d[ut] >= lambda
  list(hi = rc[far, 1L], hj = rc[far, 2L])
}

# Distance matrix realising a given H edge list: H-pairs far, others near.
.HDist <- function(n, hi, hj, far = 10, near = 1) {
  d <- matrix(near, n, n)
  diag(d) <- 0
  for (e in seq_along(hi)) {
    d[hi[e], hj[e]] <- far
    d[hj[e], hi[e]] <- far
  }
  d
}

.Decide <- function(hi, hj, n, k, maxSeconds = 60) {
  Coreset:::ThresholdDecide_cpp(as.integer(hi), as.integer(hj), n, k, maxSeconds)
}

# Is `cl` a clique of hAdj?
.IsClique <- function(hAdj, cl) {
  if (length(cl) < 2L) {
    return(TRUE)
  }
  sub <- hAdj[cl, cl, drop = FALSE]
  all(sub[upper.tri(sub)])
}

test_that("the triangle scans match the R idioms they replace", {
  set.seed(7)
  for (n in c(2L, 3L, 40L)) {
    d <- as.matrix(stats::dist(matrix(stats::rnorm(n * 2L), ncol = 2L)))
    ut <- which(upper.tri(d))
    rc <- arrayInd(ut, dim(d))
    ud <- d[ut]
    for (lambda in c(-Inf, stats::quantile(ud, c(0, 0.4, 1), names = FALSE),
                     Inf)) {
      keep <- ud >= lambda
      expect_identical(Coreset:::TriangleAtLeast_cpp(d, lambda), ud[keep])
      h <- Coreset:::EdgesAtLeast_cpp(d, lambda)
      expect_identical(h$hi, rc[keep, 1L])
      expect_identical(h$hj, rc[keep, 2L])
    }
  }
})

test_that("ThresholdDecide_cpp agrees with brute force on random probes", {
  for (seed in 1:20) {
    set.seed(seed)
    n <- sample(8:12, 1)
    pts <- matrix(stats::rnorm(n * 3L), ncol = 3L)
    d <- as.matrix(stats::dist(pts))
    ud <- d[upper.tri(d)]
    for (k in 2:5) {
      for (p in c(0.3, 0.6, 0.9)) {
        lambda <- stats::quantile(ud, p, names = FALSE)
        pa <- .ProbeArgs(d, lambda)
        hi <- pa$hi
        hj <- pa$hj
        hAdj <- matrix(FALSE, n, n)
        hAdj[cbind(hi, hj)] <- TRUE
        hAdj[cbind(hj, hi)] <- TRUE
        expected <- .HasKClique(hAdj, k)
        got <- .Decide(hi, hj, n, k)
        info <- sprintf("seed=%d n=%d k=%d p=%.1f", seed, n, k, p)

        expect_identical(got$status,
                         if (expected) "feasible" else "infeasible",
                         info = info)
        if (expected) {
          cl <- got$witness
          expect_length(cl, k)
          expect_false(is.unsorted(cl, strictly = TRUE), info = info)
          expect_true(.IsClique(hAdj, cl), info = info)
        } else {
          expect_identical(got$witness, integer(0), info = info)
        }
      }
    }
  }
})

test_that("ThresholdDecide_cpp handles targeted structures", {
  # Path: peel cascades to an empty core.
  expect_identical(.Decide(1:5, 2:6, 6L, 4L)$status, "infeasible")

  # Complete tripartite K[2,2,2]: the core survives the peel, but chi = 3, so
  # the colour bound refutes k = 4 at the root.
  hi <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4)
  hj <- c(3, 4, 5, 6, 3, 4, 5, 6, 5, 6, 5, 6)
  expect_identical(.Decide(hi, hj, 6L, 4L)$status, "infeasible")
  # ... and hosts a triangle at k = 3.
  got <- .Decide(hi, hj, 6L, 3L)
  expect_identical(got$status, "feasible")
  expect_length(got$witness, 3L)

  # K5 at k = 4: the search stops at k rather than growing the clique.
  cmb <- utils::combn(5L, 2L)
  expect_length(.Decide(cmb[1, ], cmb[2, ], 5L, 4L)$witness, 4L)

  # Single edge at k = 2 is itself the witness.
  expect_identical(.Decide(1L, 2L, 2L, 2L)$witness, 1:2)

  # No H edges at all: nothing can survive.
  expect_identical(.Decide(integer(0), integer(0), 4L, 2L)$status, "infeasible")

  # C5 at k = 3: chi = 3 cannot refute it at the root, so the refutation is
  # the search's own -- omega = 2.
  expect_identical(.Decide(c(1, 2, 3, 4, 5), c(2, 3, 4, 5, 1), 5L, 3L)$status,
                   "infeasible")

  # C5 bridged to a triangle: one component, whose triangle a greedy descent
  # from the highest-degree vertex misses.
  got <- .Decide(c(1, 2, 3, 4, 5, 6, 7, 6, 5), c(2, 3, 4, 5, 1, 7, 8, 8, 6),
                 8L, 3L)
  expect_identical(got$status, "feasible")
  expect_identical(got$witness, c(6L, 7L, 8L))

  # Two disjoint C5s: both components must be refuted in turn.
  expect_identical(.Decide(c(1:5, 6:10), c(2:5, 1, 7:10, 6), 10L, 3L)$status,
                   "infeasible")

  # A C5 component refuted first, then a bridged-triangle component found:
  # the search continues past a refuted component.
  got <- .Decide(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 11, 10),
                 c(2, 3, 4, 5, 1, 7, 8, 9, 10, 6, 12, 13, 13, 11), 13L, 3L)
  expect_identical(got$status, "feasible")
  expect_identical(got$witness, c(11L, 12L, 13L))

  # Determinism: identical calls give identical output.
  set.seed(99)
  n <- 12L
  m <- matrix(stats::runif(n * n), n, n)
  hAdj <- (m + t(m)) > 1
  hi <- row(hAdj)[upper.tri(hAdj) & hAdj]
  hj <- col(hAdj)[upper.tri(hAdj) & hAdj]
  expect_identical(.Decide(hi, hj, n, 3L), .Decide(hi, hj, n, 3L))
})

test_that("ThresholdDecide_cpp reports inconclusive when the budget expires", {
  # The deadline is only checked every 1024 search nodes, so triggering it
  # needs a single H-component whose exhaustive search takes >= 1024 Expand()
  # calls. This instance reliably does (verified well above that floor) while
  # still resolving in well under a second at a real budget.
  set.seed(1350)
  n <- 300L
  m <- matrix(stats::runif(n * n), n, n)
  hAdj <- (m + t(m)) > 1
  diag(hAdj) <- FALSE
  hi <- row(hAdj)[upper.tri(hAdj) & hAdj]
  hj <- col(hAdj)[upper.tri(hAdj) & hAdj]

  # A generous budget still proves infeasibility outright (`-O0 --coverage`
  # builds run slower, so this stays well clear of the real deadline).
  full <- .Decide(hi, hj, n, 16L, maxSeconds = 60)
  expect_identical(full$status, "infeasible")
  expect_identical(full$witness, integer(0))

  # A budget that expires before the first deadline check can complete the
  # search leaves the probe genuinely undecided.
  expired <- .Decide(hi, hj, n, 16L, maxSeconds = 0)
  expect_identical(expired$status, "inconclusive")
  expect_identical(expired$witness, integer(0))
})

test_that("ThresholdDecide_cpp agrees with the packing IP past brute force", {
  testthat::skip_if_not_installed("highs")
  # alpha(G) >= k is the same question, asked algebraically: one row per
  # G-edge, x binary, maximise sum x.
  IPFeasible <- function(hAdj, k) {
    n <- nrow(hAdj)
    gi <- row(hAdj)[upper.tri(hAdj) & !hAdj]
    gj <- col(hAdj)[upper.tri(hAdj) & !hAdj]
    nEdge <- length(gi)
    if (nEdge == 0L) {
      return(TRUE)
    }
    res <- highs::highs_solve(
      L = rep.int(1, n), lower = rep.int(0, n), upper = rep.int(1, n),
      A = Matrix::sparseMatrix(i = rep.int(seq_len(nEdge), 2L),
                               j = c(gi, gj), x = 1, dims = c(nEdge, n)),
      lhs = rep.int(-Inf, nEdge), rhs = rep.int(1, nEdge),
      types = rep.int("I", n), maximum = TRUE,
      control = list(threads = 1L, time_limit = 60))
    res$objective_value >= k - 1e-6
  }
  for (seed in 1:6) {
    set.seed(seed)
    n <- 30L
    m <- matrix(stats::runif(n * n), n, n)
    hAdj <- (m + t(m)) > 1.1                    # ~ 40% density
    diag(hAdj) <- FALSE
    hi <- row(hAdj)[upper.tri(hAdj) & hAdj]
    hj <- col(hAdj)[upper.tri(hAdj) & hAdj]
    for (k in c(3L, 5L, 7L)) {
      got <- .Decide(hi, hj, n, k)
      info <- sprintf("seed=%d k=%d", seed, k)
      expect_identical(got$status,
                       if (IPFeasible(hAdj, k)) "feasible" else "infeasible",
                       info = info)
      if (length(got$witness)) {
        expect_true(.IsClique(hAdj, got$witness), info = info)
      }
    }
  }
})

test_that(".MaxISVerdict decides probes and validates its witness", {
  # Three far-apart clusters of two: H is K[2,2,2].
  hi <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4)
  hj <- c(3, 4, 5, 6, 3, 4, 5, 6, 5, 6, 5, 6)
  d <- .HDist(6L, hi, hj)
  pa <- .ProbeArgs(d, 5)

  # k = 4: chi = 3 refutes at the root -- the certifying-probe fast path.
  v <- Coreset:::.MaxISVerdict(d, 6L, pa$hi, pa$hj, 5, 4L, 60)
  expect_identical(v$verdict, "infeasible")
  expect_identical(v$witness, integer(0))

  # k = 3: a witness, whose pairwise distances all clear lambda.
  v <- Coreset:::.MaxISVerdict(d, 6L, pa$hi, pa$hj, 5, 3L, 60)
  expect_identical(v$verdict, "feasible")
  expect_length(v$witness, 3L)
  sub <- d[v$witness, v$witness]
  expect_true(all(sub[upper.tri(sub)] >= 5))

  # lambda below every distance: empty G, all vertices independent.
  paLow <- .ProbeArgs(d, 0.5)
  v <- Coreset:::.MaxISVerdict(d, 6L, paLow$hi, paLow$hj, 0.5, 3L, 60)
  expect_identical(v$verdict, "feasible")
  expect_identical(v$witness, 1:6)

  # A budget that cannot buy any search is inconclusive, never a verdict.
  v <- Coreset:::.MaxISVerdict(d, 6L, pa$hi, pa$hj, 5, 3L, 0)
  expect_identical(v$verdict, "inconclusive")
})

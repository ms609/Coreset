# Tests for ThresholdReduce_cpp and the reduced .MaxISVerdict probe.
#
# Oracle: brute-force clique enumeration on tiny graphs. The kernel reduces
# the complement graph H of a threshold probe; its contract is
# verdict-preservation (never lose a k-clique), not alpha-preservation, so
# the properties checked are: returned cliques are real, refutations agree
# with brute force, colourings are proper, and the emitted model covers
# every G-edge among its variables.

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

# Upper-triangle probe inputs for .MaxISVerdict from a distance matrix.
.ProbeArgs <- function(d, lambda) {
  ut <- which(upper.tri(d))
  rc <- arrayInd(ut, dim(d))
  list(ui = rc[, 1L], uj = rc[, 2L], gEdge = d[ut] < lambda)
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

test_that("ThresholdReduce_cpp agrees with brute force on random probes", {
  Reduce1 <- function(hi, hj, n, k) {
    MaxMin:::ThresholdReduce_cpp(as.integer(hi), as.integer(hj), n, k)
  }
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
        hi <- pa$ui[!pa$gEdge]
        hj <- pa$uj[!pa$gEdge]
        hAdj <- matrix(FALSE, n, n)
        hAdj[cbind(hi, hj)] <- TRUE
        hAdj[cbind(hj, hi)] <- TRUE
        expected <- .HasKClique(hAdj, k)
        red <- Reduce1(hi, hj, n, k)
        info <- sprintf("seed=%d n=%d k=%d p=%.1f", seed, n, k, p)

        if (length(red$clique)) {
          # A returned clique is a real k-clique of H.
          cl <- red$clique
          expect_length(cl, k)
          expect_false(is.unsorted(cl, strictly = TRUE))
          expect_true(all(hAdj[cbind(rep(cl, each = k), rep(cl, k))] |
                            rep(cl, each = k) == rep(cl, k)), info = info)
          expect_true(expected, info = info)
        } else if (length(red$comps) == 0L) {
          # A combinatorial refutation must agree with brute force.
          expect_false(expected, info = info)
        } else {
          # Verdict-preservation: a k-clique, if any, survives in some
          # component; the model covers every G-edge among its variables.
          inComp <- vapply(red$comps, function(cp) {
            .HasKClique(hAdj[cp$vars, cp$vars, drop = FALSE], k)
          }, logical(1))
          expect_equal(any(inComp), expected, info = info)
          for (cp in red$comps) {
            vars <- cp$vars
            expect_false(is.unsorted(vars, strictly = TRUE))
            nv <- length(vars)
            # Proper colouring of H[vars]; class-mates are G-adjacent.
            for (a in seq_len(nv - 1L)) {
              for (b in seq(a + 1L, nv)) {
                u <- vars[a]
                v <- vars[b]
                covered <- if (hAdj[u, v]) {
                  cp$cls[a] != cp$cls[b]                # proper colouring
                } else if (cp$cls[a] == cp$cls[b]) {
                  TRUE                                  # class row covers it
                } else {
                  any(cp$ei == a & cp$ej == b)          # pairwise row needed
                }
                expect_true(covered, info = sprintf("%s pair %d-%d", info, u, v))
              }
            }
            # Emitted rows are exactly cross-colour G-edges.
            expect_true(all(!hAdj[cbind(vars[cp$ei], vars[cp$ej])]))
            expect_true(all(cp$cls[cp$ei] != cp$cls[cp$ej]))
          }
        }
      }
    }
  }
})

test_that("ThresholdReduce_cpp handles targeted structures", {
  Reduce1 <- function(hi, hj, n, k) {
    MaxMin:::ThresholdReduce_cpp(as.integer(hi), as.integer(hj), n, k)
  }
  # Path: peel cascades to an empty core.
  red <- Reduce1(1:5, 2:6, 6L, 4L)
  expect_length(red$clique, 0L)
  expect_length(red$comps, 0L)

  # Complete tripartite K[2,2,2]: core survives peel but chi = 3 < k = 4.
  hi <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4)
  hj <- c(3, 4, 5, 6, 3, 4, 5, 6, 5, 6, 5, 6)
  red <- Reduce1(hi, hj, 6L, 4L)
  expect_length(red$clique, 0L)
  expect_length(red$comps, 0L)
  # ... and hosts a triangle at k = 3, found greedily.
  red <- Reduce1(hi, hj, 6L, 3L)
  expect_length(red$clique, 3L)

  # K5 at k = 4: the greedy clique stops at k.
  cmb <- utils::combn(5L, 2L)
  red <- Reduce1(cmb[1, ], cmb[2, ], 5L, 4L)
  expect_length(red$clique, 4L)

  # Single edge at k = 2 is itself the witness.
  red <- Reduce1(1L, 2L, 2L, 2L)
  expect_identical(red$clique, 1:2)

  # No H edges at all: nothing can survive.
  red <- Reduce1(integer(0), integer(0), 4L, 2L)
  expect_length(red$clique, 0L)
  expect_length(red$comps, 0L)

  # Determinism: identical calls give identical output.
  set.seed(99)
  n <- 12L
  m <- matrix(stats::runif(n * n), n, n)
  hAdj <- (m + t(m)) > 1
  hi <- row(hAdj)[upper.tri(hAdj) & hAdj]
  hj <- col(hAdj)[upper.tri(hAdj) & hAdj]
  expect_identical(Reduce1(hi, hj, n, 3L), Reduce1(hi, hj, n, 3L))
})

test_that(".MaxISVerdict decides reduction-only probes without the solver", {
  # Three far-apart clusters of two: H is K[2,2,2].
  hi <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4)
  hj <- c(3, 4, 5, 6, 3, 4, 5, 6, 5, 6, 5, 6)
  d <- .HDist(6L, hi, hj)
  pa <- .ProbeArgs(d, 5)

  # k = 4: chi = 3 refutes with no IP -- the certifying-probe fast path.
  v <- MaxMin:::.MaxISVerdict(d, 6L, pa$ui, pa$uj, pa$gEdge, 5, 4L, 60)
  expect_identical(v$verdict, "infeasible")
  expect_identical(v$witness, integer(0))

  # k = 3: the greedy clique is the witness, again with no IP.
  v <- MaxMin:::.MaxISVerdict(d, 6L, pa$ui, pa$uj, pa$gEdge, 5, 3L, 60)
  expect_identical(v$verdict, "feasible")
  expect_length(v$witness, 3L)
  sub <- d[v$witness, v$witness]
  expect_true(all(sub[upper.tri(sub)] >= 5))

  # lambda below every distance: empty G, all vertices independent.
  paLow <- .ProbeArgs(d, 0.5)
  v <- MaxMin:::.MaxISVerdict(d, 6L, paLow$ui, paLow$uj, paLow$gEdge, 0.5, 3L, 60)
  expect_identical(v$verdict, "feasible")
  expect_identical(v$witness, 1:6)
})

test_that(".MaxISVerdict IP branch decides probes the reduction cannot", {
  testthat::skip_if_not_installed("highs")
  # C5: chi = 3 cannot refute k = 3, but omega = 2 -- the IP must refute.
  d <- .HDist(5L, c(1, 2, 3, 4, 5), c(2, 3, 4, 5, 1))
  pa <- .ProbeArgs(d, 5)
  v <- MaxMin:::.MaxISVerdict(d, 5L, pa$ui, pa$uj, pa$gEdge, 5, 3L, 60)
  expect_identical(v$verdict, "infeasible")

  # C5 bridged to a triangle: one component, greedy clique misses the
  # triangle, the IP finds it.
  hi <- c(1, 2, 3, 4, 5, 6, 7, 6, 5)
  hj <- c(2, 3, 4, 5, 1, 7, 8, 8, 6)
  d <- .HDist(8L, hi, hj)
  pa <- .ProbeArgs(d, 5)
  v <- MaxMin:::.MaxISVerdict(d, 8L, pa$ui, pa$uj, pa$gEdge, 5, 3L, 60)
  expect_identical(v$verdict, "feasible")
  expect_identical(sort(v$witness), c(6L, 7L, 8L))

  # Two disjoint C5s: both components must be refuted in turn.
  d <- .HDist(10L, c(1:5, 6:10), c(2:5, 1, 7:10, 6))
  pa <- .ProbeArgs(d, 5)
  v <- MaxMin:::.MaxISVerdict(d, 10L, pa$ui, pa$uj, pa$gEdge, 5, 3L, 60)
  expect_identical(v$verdict, "infeasible")

  # A C5 component refuted first, then a bridged-triangle component found
  # feasible: exercises the continue-past-a-refuted-component path.
  hi <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 11, 10)
  hj <- c(2, 3, 4, 5, 1, 7, 8, 9, 10, 6, 12, 13, 13, 11)
  d <- .HDist(13L, hi, hj)
  pa <- .ProbeArgs(d, 5)
  v <- MaxMin:::.MaxISVerdict(d, 13L, pa$ui, pa$uj, pa$gEdge, 5, 3L, 60)
  expect_identical(v$verdict, "feasible")
  expect_identical(sort(v$witness), c(11L, 12L, 13L))
})

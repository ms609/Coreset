# Tests for DropAdd (Porumbel, Hao & Glover 2011).
#
# The validation oracle is Geo 100 1 with m=10: published best-known MaxMin
# objective is 89.37 and the instance is small enough that the algorithm
# should reach it within a 10-second budget on any modern machine.

# Helper: brute-force MaxMin objective for a selection.
.MaxminBrute <- function(S, dmat) {
  if (length(S) < 2L) return(NA_real_)
  sub <- dmat[S, S]
  diag(sub) <- Inf
  min(sub)
}

# Locate the bundled Geo loader. Skip tests that need it if it isn't reachable
# (e.g. installed-package CMD check without dev/ in tree).
.TryGeoLoader <- function() {
  pkgRoot <- tryCatch(rprojroot::find_root(rprojroot::is_r_package),
                       error = function(e) NULL)
  if (is.null(pkgRoot)) return(NULL)
  loader <- file.path(pkgRoot, "dev", "benchmarks", "mdplib", "loader.R")
  if (!file.exists(loader)) return(NULL)
  env <- new.env(parent = baseenv())
  source(loader, local = env, chdir = TRUE)
  # mdplib_root() relies on sys.frame()$ofile which becomes NULL after the
  # source() call returns; injecting an explicit absolute path here avoids
  # the lookup entirely.
  env$.MDPLIB_INSTANCE_DIR <- file.path(
    pkgRoot, "dev", "benchmarks", "mdplib", "extracted", "instances"
  )
  env
}

# ---------------------------------------------------------------------------
# 1. Smoke test
# ---------------------------------------------------------------------------
test_that("DropAdd smoke: 20 pts in 5-D, k=5", {
  set.seed(7)
  pts <- matrix(rnorm(20 * 5), ncol = 5)
  d <- dist(pts)
  res <- DropAdd(d, k = 5L, maxSeconds = 1)
  expect_length(res, 5L)
  expect_equal(length(unique(res)), 5L)
  expect_true(all(res %in% seq_len(20L)))
  expect_gt(attr(res, "score"), 0)
  expect_gte(attr(res, "iters"), 5L)
  expect_true(attr(res, "time_s") >= 0)
})

# ---------------------------------------------------------------------------
# 2. Construction-only sanity
# ---------------------------------------------------------------------------
test_that("DropAdd construction-only result matches its MaxMin score", {
  set.seed(11)
  pts <- matrix(rnorm(30 * 4), ncol = 4)
  dmat <- as.matrix(dist(pts))
  res <- DropAdd(dmat, k = 6L, maxIter = 0L)
  expect_length(res, 6L)
  expect_equal(attr(res, "iters"), 0L)
  # Objective stored equals the actual MaxMin over the returned indices.
  expect_equal(attr(res, "score"), .MaxminBrute(res, dmat))
})

# ---------------------------------------------------------------------------
# 3. Improvement over construction on a real benchmark instance
# ---------------------------------------------------------------------------
test_that("DropAdd never worsens the constructive solution", {
  geoEnv <- .TryGeoLoader()
  skip_if(is.null(geoEnv), "Geo loader not available")
  path <- file.path(geoEnv$.MDPLIB_INSTANCE_DIR, "Geo", "Geo 100 1.txt")
  skip_if_not(file.exists(path), "Geo 100 1 instance not present")
  geo <- geoEnv$read_mdplib_geo(path)
  dmat <- geoEnv$mdplib_geo_dist(geo)

  cons <- DropAdd(dmat, k = 10L, maxIter = 0L)
  full <- DropAdd(dmat, k = 10L, plateau = 2000L)
  expect_gte(attr(full, "score"), attr(cons, "score") - 1e-9)
})

# ---------------------------------------------------------------------------
# 4. Validation oracle: Geo 100 1, m=10
# ---------------------------------------------------------------------------
# The proven optimum on this instance is 89.37 (Della Croce, Grosso & Locatelli
# 2009 prove all n=100 MMDP optima; cross-method best-known in Marti et al.
# 2008/2022). Porumbel's drop-add search excludes the just-dropped point x# from
# the add candidates for one iteration (Add X(k) = Z - X(k); paper p.281), which
# forces turnover so the deterministic FIFO+greedy escapes local optima and
# reaches the optimum (89.3701) well within a 10 s budget (Della Croce et al.
# report n=100 optima are reachable in 2-3 s). Test 5 directly asserts the x#
# exclusion (adds != drops), the invariant a prior frozen-at-86.99 version broke.
test_that("DropAdd reaches the Geo 100 1 m=10 optimum (89.37)", {
  geoEnv <- .TryGeoLoader()
  skip_if(is.null(geoEnv), "Geo loader not available")
  path <- file.path(geoEnv$.MDPLIB_INSTANCE_DIR, "Geo", "Geo 100 1.txt")
  skip_if_not(file.exists(path), "Geo 100 1 instance not present")
  geo <- geoEnv$read_mdplib_geo(path)
  dmat <- geoEnv$mdplib_geo_dist(geo)

  # Time-budgeted oracle: run the full 10 s, as before, to confirm the search
  # reaches the proven optimum (this is a correctness oracle, not a frozen
  # result, so a wall-clock budget is appropriate here).
  res <- DropAdd(dmat, k = 10L, maxSeconds = 10, plateau = 100000000L)
  bestKnown <- 89.37
  expect_gte(attr(res, "score"), 0.999 * bestKnown)
  # Returned indices truly achieve the reported objective.
  expect_equal(attr(res, "score"), .MaxminBrute(res, dmat),
               tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# 5. FIFO invariant: across the first k main-loop iterations, every initially
# selected point is dropped exactly once.
# ---------------------------------------------------------------------------
test_that("DropAdd FIFO drops each initial point once in first k iterations", {
  set.seed(2026)
  pts <- matrix(rnorm(10 * 3), ncol = 3)
  dmat <- as.matrix(dist(pts))
  k <- 4L

  # Capture the constructive S via the internal helper.
  cons <- MaxMin:::.DropAddConstruct(dmat, k)
  sInit <- cons$S
  expect_length(sInit, k)
  expect_equal(length(unique(sInit)), k)

  # Drive the PRODUCTION C++ loop for exactly 2*k iterations and capture the
  # drop/add sequences via the internal-only .DropAddTrace() scaffolding (the
  # public DropAdd() deliberately exposes no trace hook).
  tr <- MaxMin:::.DropAddTrace(dmat, k = k, maxSeconds = 60, maxIter = 2L * k)
  expect_equal(tr$iters, 2L * k)
  expect_length(tr$drops, 2L * k)
  # First k drops are exactly the constructive members, in some order
  # (FIFO order is determined by iter_stamp, which was 1..k in construction
  # so drop order is the construction order).
  expect_setequal(tr$drops[seq_len(k)], sInit)
  expect_equal(tr$drops[seq_len(k)], sInit)  # actually in order

  # Across 2*k iterations, every dropped point was previously in S.
  # No index can be dropped twice in the first k iterations.
  expect_equal(anyDuplicated(tr$drops[seq_len(k)]), 0L)

  # Porumbel p.281: the just-dropped point is excluded from the add candidates
  # for that iteration, so the added point is never the one just dropped. This
  # is the diversification invariant; without it the search freezes.
  expect_true(all(tr$adds != tr$drops))
})

# ---------------------------------------------------------------------------
# 6. Time budget honoured on a 200-point Ran-format instance
# ---------------------------------------------------------------------------
test_that("DropAdd respects maxSeconds within reasonable slack", {
  set.seed(99)
  pts <- matrix(runif(200 * 5), ncol = 5)
  dmat <- as.matrix(dist(pts))
  t0 <- Sys.time()
  # Disable stagnation so the wall-clock ceiling is the binding criterion.
  res <- DropAdd(dmat, k = 20L, maxSeconds = 0.05, plateau = 100000000L)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_lte(attr(res, "time_s"), 1.5)
  expect_lte(elapsed, 2.0)
  expect_length(res, 20L)
})

# ---------------------------------------------------------------------------
# 8. seed parameter, input validation, and progress output
# ---------------------------------------------------------------------------
test_that("DropAdd is deterministic and validates inputs", {
  set.seed(11)
  dmat <- as.matrix(dist(matrix(rnorm(20 * 3), ncol = 3)))

  # The search is RNG-free, so repeated calls are identical. (The old `seed`
  # argument was a documented no-op and has been removed from the API.)
  set.seed(1)
  r1 <- DropAdd(dmat, k = 4L, maxIter = 0L)
  set.seed(999)
  r2 <- DropAdd(dmat, k = 4L, maxIter = 0L)
  attr(r1, "time_s") <- attr(r2, "time_s") <- NULL
  expect_identical(r1, r2, label = "DropAdd is RNG-independent: different seeds give identical result")

  # k validation
  expect_error(DropAdd(dmat, k = 1L),   "2 <= k")
  expect_error(DropAdd(dmat, k = 25L),  "2 <= k")
  # plateau validation
  expect_error(DropAdd(dmat, k = 4L, plateau = 0L),  "plateau")
  # maxIter validation
  expect_error(DropAdd(dmat, k = 4L, maxIter = -1L),  "maxIter")
  # maxSeconds validation
  expect_error(DropAdd(dmat, k = 4L, maxSeconds = 0),   "maxSeconds")
  expect_error(DropAdd(dmat, k = 4L, maxSeconds = NA_real_), "maxSeconds")
})

test_that("DropAdd progress = TRUE fires the cli hooks", {
  dmat <- as.matrix(dist(matrix(rnorm(15 * 2), ncol = 2)))
  expect_no_error(suppressMessages(
    DropAdd(dmat, k = 3L, maxIter = 2L, progress = TRUE)
  ))
})

test_that("DropAdd rejects an NA distance matrix (FF-001 / T7-08)", {
  # The matrix path coerces via .AsDistMatrix, which now guards against NA so
  # the search cannot silently accept a corrupt matrix.
  naMat <- matrix(c(0, 5, NA, 5, 0, NA, NA, NA, 0), 3L, 3L)
  expect_error(DropAdd(naMat, k = 3L), "NA")
})

# ---------------------------------------------------------------------------
# 9. Degenerate and stopping-criterion corner cases (C++ path).
# Tie-breaking / caseA / k=2 record branches are exercised by tests 10a-10c.
# ---------------------------------------------------------------------------
test_that("DropAdd performs no iterations when k == n", {
  dmat <- as.matrix(dist(matrix(rnorm(5 * 2), ncol = 2)))
  res  <- DropAdd(dmat, k = 5L)
  # All 5 points selected, no drop-add move exists.
  expect_length(res, 5L)
  expect_equal(attr(res, "iters"), 0L)
})

test_that("DropAdd time budget halts when both other criteria are disabled", {
  set.seed(1)
  dmat <- as.matrix(dist(matrix(rnorm(30 * 3), ncol = 3)))
  # .Machine$integer.max disables both other stopping criteria; only the
  # wall-clock budget can end the loop.
  res <- expect_returns_within(
    DropAdd(dmat, k = 5L,
            maxIter = .Machine$integer.max, plateau = .Machine$integer.max,
            maxSeconds = 0.001),
    limit = 5)
  expect_gte(attr(res, "iters"), 1L)    # at least one iteration ran
  expect_lte(attr(res, "time_s"), 0.5)
})

# ---------------------------------------------------------------------------
# 10. C++ path coverage: tie-breaking and minDist-count branches
#
# Tests 10a-10c exercise specific branches in dropadd.cpp that require exact
# distance ties.
# ---------------------------------------------------------------------------

test_that("DropAdd C++ construction covers sumDist tie-break (lines 82-85)", {
  # 5-point geometry: A=(0,0), B=(3,0), C=(0,2), P=(1,1), Q=(2,1).
  # Construction: B (max row-sum) -> C -> A.
  # Then P and Q both have minDist=sqrt(2) to {B,C,A}, but
  # sumDist[Q]=sqrt(2)+2*sqrt(5) > sumDist[P]=2*sqrt(2)+sqrt(5),
  # so the tie-break updates xNew to Q (lines 83-85 fire).
  # When A is added, d(P,A)=sqrt(2)==minDist[P] -> line 103 fires too.
  pts <- rbind(c(0,0), c(3,0), c(0,2), c(1,1), c(2,1))
  d   <- as.matrix(dist(pts))
  res <- DropAdd(d, k = 4L, maxIter = 0L)
  expect_length(res, 4L)
  expect_equal(attr(res, "iters"), 0L)
})

test_that("DropAdd C++ k=2 covers DROP else-branch (lines 214-215)", {
  # Rhombus: (0,0),(1,1),(2,0),(1,-1).  All edges=sqrt(2), diagonals=2.
  # Construction: {(0,0),(2,0)} (endpoints of the longer diagonal).
  # Drop (0,0): the remaining selected (2,0) had (0,0) as its only peer,
  # so minDistCount[(2,0)] drops to 0.  In the recompute, the self-mask
  # skips the only surviving S member, leaving mns=Inf -> else branch fires.
  ptsRh <- rbind(c(0,0), c(1,1), c(2,0), c(1,-1))
  dRh   <- as.matrix(dist(ptsRh))
  res <- DropAdd(dRh, k = 2L, maxIter = 4L)
  expect_length(res, 2L)
})

test_that("DropAdd C++ main-loop ADD covers tie-break (255-258) and equality (277)", {
  # 7-point geometry ensures:
  #   * Lines 255-258: after dropping ANC, X=(1,0) and Y=(3,0) tie on
  #     minDist=1 but sumDist[Y]>sumDist[X], so Y wins via the update.
  #   * Line 277: when Y is added, d(Y,W)=0.5=minDist[W] -> count increments.
  #
  # Points: P(0,0), Q(4,0), R(0,4), ANC(10,10), X(1,0), Y(3,0), W(3.5,0).
  # Construction: ANC -> P -> Q -> R  (ANC has max row-sum).
  pts7 <- rbind(c(0,0), c(4,0), c(0,4), c(10,10), c(1,0), c(3,0), c(3.5,0))
  d7   <- as.matrix(dist(pts7))
  res <- DropAdd(d7, k = 4L, maxIter = 1L)
  expect_length(res, 4L)
  expect_equal(attr(res, "iters"), 1L)
})

# ---------------------------------------------------------------------------
# 7. Frozen golden-value regression for the C++ search trajectory.
#
# Replaces the former "C++ vs hand-maintained R reference" bit-identity test.
# That comparison was inherently FP-fragile (R accumulates dist/rowSums in
# long double, the C++ port in double), so it flickered run-to-run on random
# coordinate instances. Here we freeze the C++ output against a fixed *integer*
# symmetric distance matrix: integers <= 2^53 are exact in double and the C++
# path uses no long double, so every recorded value is reproducible and
# portable across platforms. Any change to the search trajectory -- indices,
# objective, secondary, or run length -- trips this test. Correctness (as
# opposed to mere stability) is anchored separately by the Geo 100 1 optimum
# oracle (test 4).
#
# To regenerate after an intentional algorithm change, re-run the block that
# produced these values (same seed/matrix) and paste the new expectations.
# ---------------------------------------------------------------------------

# Fixed integer symmetric distance matrix (24 points), built deterministically.
.DropAddGoldenMatrix <- function() {
  set.seed(2026)
  n <- 24L
  M <- matrix(0, n, n)
  M[upper.tri(M)] <- sample(5:995, n * (n - 1L) / 2L, replace = TRUE)
  M <- M + t(M)
  storage.mode(M) <- "double"
  M
}

test_that("DropAdd C++ trajectory matches frozen golden values (maxIter)", {
  M <- .DropAddGoldenMatrix()
  # idx | score | secondary | iters, frozen from the production C++ path.
  golden <- list(
    "0"   = list(idx = c(1,7,8,9,13,18,21,23),    score = 304, sec = 19133, it = 0L),
    "1"   = list(idx = c(1,7,8,9,13,18,21,23),    score = 304, sec = 19133, it = 1L),
    "5"   = list(idx = c(1,7,8,9,13,18,21,23),    score = 304, sec = 19133, it = 5L),
    "50"  = list(idx = c(6,8,15,16,19,20,21,23),  score = 330, sec = 16975, it = 50L),
    "200" = list(idx = c(6,8,15,16,19,20,21,23),  score = 330, sec = 16975, it = 200L)
  )
  for (mi in names(golden)) {
    g   <- golden[[mi]]
    res <- DropAdd(M, k = 8L, maxIter = as.integer(mi))
    expect_identical(as.integer(res), as.integer(g$idx), info = paste("maxIter", mi))
    expect_equal(attr(res, "score"),     g$score, tolerance = 1e-12, info = mi)
    expect_equal(attr(res, "secondary"), g$sec,   tolerance = 1e-12, info = mi)
    expect_identical(attr(res, "iters"), g$it,             info = mi)
  }
})

test_that("DropAdd C++ trajectory matches frozen golden values (plateau)", {
  M <- .DropAddGoldenMatrix()
  # The run length under plateau is itself an output, so frozen `iters` confirms
  # the deterministic stagnation criterion fires exactly where expected.
  golden <- list(
    "5"   = list(idx = c(1,7,8,9,13,18,21,23),    score = 304, sec = 19133, it = 5L),
    "25"  = list(idx = c(6,8,15,16,19,20,21,23),  score = 330, sec = 16975, it = 44L),
    "100" = list(idx = c(6,8,15,16,19,20,21,23),  score = 330, sec = 16975, it = 119L)
  )
  for (pl in names(golden)) {
    g   <- golden[[pl]]
    res <- DropAdd(M, k = 8L, plateau = as.integer(pl))
    expect_identical(as.integer(res), as.integer(g$idx), info = paste("plateau", pl))
    expect_equal(attr(res, "score"),     g$score, tolerance = 1e-12, info = pl)
    expect_equal(attr(res, "secondary"), g$sec,   tolerance = 1e-12, info = pl)
    expect_identical(attr(res, "iters"), g$it,             info = pl)
  }
})

test_that("DropAdd secondary attribute equals upper-triangle distance sum", {
  set.seed(2026)
  pts <- matrix(rnorm(20 * 3), ncol = 3)
  dmat <- as.matrix(dist(pts))
  res <- DropAdd(dmat, k = 5L, maxIter = 20L)
  # Brute-force the secondary: upper-triangle sum of distances within selection.
  sub <- dmat[res, res]
  expected_secondary <- sum(sub[upper.tri(sub)])
  expect_equal(attr(res, "secondary"), expected_secondary, tolerance = 1e-10)
})

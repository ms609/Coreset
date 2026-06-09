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
test_that("DropAdd smoke: 20 pts in 5-D, m=5", {
  set.seed(7)
  pts <- matrix(rnorm(20 * 5), ncol = 5)
  d <- dist(pts)
  res <- DropAdd(d, m = 5L, timeBudgetS = 1)
  expect_length(res, 5L)
  expect_equal(length(unique(res)), 5L)
  expect_true(all(res %in% seq_len(20L)))
  expect_gt(attr(res, "objective"), 0)
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
  res <- DropAdd(dmat, m = 6L, maxIter = 0L)
  expect_length(res, 6L)
  expect_equal(attr(res, "iters"), 0L)
  # Objective stored equals the actual MaxMin over the returned indices.
  expect_equal(attr(res, "objective"), .MaxminBrute(res, dmat))
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

  cons <- DropAdd(dmat, m = 10L, maxIter = 0L)
  full <- DropAdd(dmat, m = 10L, plateau = 2000L)
  expect_gte(attr(full, "objective"), attr(cons, "objective") - 1e-9)
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
  res <- DropAdd(dmat, m = 10L, timeBudgetS = 10, plateau = 100000000L)
  bestKnown <- 89.37
  expect_gte(attr(res, "objective"), 0.999 * bestKnown)
  # Returned indices truly achieve the reported objective.
  expect_equal(attr(res, "objective"), .MaxminBrute(res, dmat),
               tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# 4b. Streamlined-record self-consistency under .verify
# ---------------------------------------------------------------------------
test_that("DropAdd .verify=TRUE passes silently on a real instance", {
  geoEnv <- .TryGeoLoader()
  skip_if(is.null(geoEnv), "Geo loader not available")
  path <- file.path(geoEnv$.MDPLIB_INSTANCE_DIR, "Geo", "Geo 100 1.txt")
  skip_if_not(file.exists(path), "Geo 100 1 instance not present")
  geo <- geoEnv$read_mdplib_geo(path)
  dmat <- geoEnv$mdplib_geo_dist(geo)

  # Cap iters to keep verification cheap; the inner check is O(n*m) per iter.
  expect_silent(
    DropAdd(dmat, m = 10L, timeBudgetS = 30, maxIter = 100L,
              .verify = TRUE)
  )
})

# ---------------------------------------------------------------------------
# 5. FIFO invariant: across the first m main-loop iterations, every initially
# selected point is dropped exactly once.
# ---------------------------------------------------------------------------
test_that("DropAdd FIFO drops each initial point once in first m iterations", {
  set.seed(2026)
  pts <- matrix(rnorm(10 * 3), ncol = 3)
  dmat <- as.matrix(dist(pts))
  m <- 4L

  # Capture the constructive S via the internal helper.
  cons <- MaxMin:::.DropAddConstruct(dmat, m)
  sInit <- cons$S
  expect_length(sInit, m)
  expect_equal(length(unique(sInit)), m)

  # Drive the PRODUCTION loop via the .trace hook for exactly 2*m iterations.
  traceEnv <- new.env()
  res <- DropAdd(dmat, m = m, timeBudgetS = 60, maxIter = 2L * m,
                   .trace = traceEnv)
  expect_equal(attr(res, "iters"), 2L * m)
  expect_length(traceEnv$drops, 2L * m)
  # First m drops are exactly the constructive members, in some order
  # (FIFO order is determined by iter_stamp, which was 1..m in construction
  # so drop order is the construction order).
  expect_setequal(traceEnv$drops[seq_len(m)], sInit)
  expect_equal(traceEnv$drops[seq_len(m)], sInit)  # actually in order

  # Across 2*m iterations, every dropped point was previously in S.
  # No index can be dropped twice in the first m iterations.
  expect_equal(anyDuplicated(traceEnv$drops[seq_len(m)]), 0L)

  # Porumbel p.281: the just-dropped point is excluded from the add candidates
  # for that iteration, so the added point is never the one just dropped. This
  # is the diversification invariant; without it the search freezes.
  expect_true(all(traceEnv$adds != traceEnv$drops))
})

# ---------------------------------------------------------------------------
# 6. Time budget honoured on a 200-point Ran-format instance
# ---------------------------------------------------------------------------
test_that("DropAdd respects timeBudgetS within reasonable slack", {
  set.seed(99)
  pts <- matrix(runif(200 * 5), ncol = 5)
  dmat <- as.matrix(dist(pts))
  t0 <- Sys.time()
  # Disable stagnation so the wall-clock ceiling is the binding criterion.
  res <- DropAdd(dmat, m = 20L, timeBudgetS = 1, plateau = 100000000L)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_lte(attr(res, "time_s"), 1.5)
  expect_lte(elapsed, 2.0)
  expect_length(res, 20L)
})

# ---------------------------------------------------------------------------
# 8. seed parameter, input validation, and progress output
# ---------------------------------------------------------------------------
test_that("DropAdd seed parameter and input validation", {
  set.seed(11)
  dmat <- as.matrix(dist(matrix(rnorm(20 * 3), ncol = 3)))

  # seed path
  r1 <- DropAdd(dmat, m = 4L, maxIter = 0L, seed = 1L)
  r2 <- DropAdd(dmat, m = 4L, maxIter = 0L, seed = 1L)
  expect_identical(r1, r2)

  # m validation
  expect_error(DropAdd(dmat, m = 1L),   "2 <= m")
  expect_error(DropAdd(dmat, m = 25L),  "2 <= m")
  # plateau validation
  expect_error(DropAdd(dmat, m = 4L, plateau = 0L),  "plateau")
  # maxIter validation
  expect_error(DropAdd(dmat, m = 4L, maxIter = -1L),  "maxIter")
  # timeBudgetS validation
  expect_error(DropAdd(dmat, m = 4L, timeBudgetS = 0),   "timeBudgetS")
  expect_error(DropAdd(dmat, m = 4L, timeBudgetS = NA_real_), "timeBudgetS")
})

test_that("DropAdd progress = TRUE fires the cli hooks", {
  dmat <- as.matrix(dist(matrix(rnorm(15 * 2), ncol = 2)))
  expect_no_error(DropAdd(dmat, m = 3L, maxIter = 2L, progress = TRUE))
})

# ---------------------------------------------------------------------------
# 9. R-path corner cases: tie-breaking, caseA updates, m=2 else branch
# ---------------------------------------------------------------------------
test_that("DropAdd construction tie-breaking and caseA update (.verify)", {
  # Unit square: all sides=1, diagonals=sqrt(2).
  # Construction for m=3: seed=P1, add diagonal P3 (caseA fires for P2 and P4),
  # then P2 and P4 are tied (lines 51-53 of .DropAddConstruct).
  ptsSq <- rbind(c(0,0), c(1,0), c(1,1), c(0,1))
  dmat   <- as.matrix(dist(ptsSq))
  # .verify=TRUE routes through R path with brute-force record checks
  res <- DropAdd(dmat, m = 3L, maxIter = 4L, .verify = TRUE)
  expect_length(res, 3L)
})

test_that("DropAdd m=2 else-branch and ADD tie-breaking (.verify)", {
  # Rhombus: all edges=sqrt(2), diagonals=2.
  # m=2: after drop of seed, two candidates are equidistant from remaining
  # selected point (lines 360-366 else branch + lines 391-393 ADD tie).
  ptsRh <- rbind(c(0,0), c(1,1), c(2,0), c(1,-1))
  dmat   <- as.matrix(dist(ptsRh))
  res <- DropAdd(dmat, m = 2L, maxIter = 4L, .verify = TRUE)
  expect_length(res, 2L)
})

test_that("DropAdd caseA ADD update fires in main loop (.verify)", {
  # Unit square m=3: first main-loop ADD of P4 has d[P4,P1]=1 = minDist[P1]=1
  # (caseA, line 409).
  ptsSq <- rbind(c(0,0), c(1,0), c(1,1), c(0,1))
  dmat   <- as.matrix(dist(ptsSq))
  expect_no_error(
    DropAdd(dmat, m = 3L, maxIter = 3L, .verify = TRUE)
  )
})

test_that("DropAdd effectiveMax = 0 when m == n (.verify)", {
  dmat <- as.matrix(dist(matrix(rnorm(5 * 2), ncol = 2)))
  res  <- DropAdd(dmat, m = 5L, .verify = TRUE)
  # All 5 points selected, no loop iterations.
  expect_length(res, 5L)
  expect_equal(attr(res, "iters"), 0L)
})

test_that("DropAdd .trace + .verify together cover init and loop writes", {
  set.seed(2026)
  dmat <- as.matrix(dist(matrix(rnorm(12 * 2), ncol = 2)))
  te   <- new.env()
  res  <- DropAdd(dmat, m = 4L, maxIter = 5L, .verify = TRUE, .trace = te)
  # .trace should have been initialised on R path (lines 295-296) and
  # written during the loop (lines 438-439).
  expect_equal(length(te$drops), attr(res, "iters"))
  expect_equal(length(te$adds),  attr(res, "iters"))
})

# ---------------------------------------------------------------------------
# 7. R reference loop and C++ port are bit-identical on a fixed iter budget.
#
# .verify = TRUE routes through the original R reference loop (with brute-
# force record assertions); the default routes through DropAdd_cpp. The
# two paths must produce identical indices, objective, secondary, and iters
# for any iter budget.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 10. C++ path coverage: tie-breaking and minDist-count branches
#
# Tests 10a-10c exercise specific branches in dropadd.cpp that require exact
# distance ties; they use the default (C++) path (no .verify).
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
  res <- DropAdd(d, m = 4L, maxIter = 0L)
  expect_length(res, 4L)
  expect_equal(attr(res, "iters"), 0L)
})

test_that("DropAdd C++ m=2 covers DROP else-branch (lines 214-215)", {
  # Rhombus: (0,0),(1,1),(2,0),(1,-1).  All edges=sqrt(2), diagonals=2.
  # Construction: {(0,0),(2,0)} (endpoints of the longer diagonal).
  # Drop (0,0): the remaining selected (2,0) had (0,0) as its only peer,
  # so minDistCount[(2,0)] drops to 0.  In the recompute, the self-mask
  # skips the only surviving S member, leaving mns=Inf -> else branch fires.
  ptsRh <- rbind(c(0,0), c(1,1), c(2,0), c(1,-1))
  dRh   <- as.matrix(dist(ptsRh))
  res <- DropAdd(dRh, m = 2L, maxIter = 4L)
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
  res <- DropAdd(d7, m = 4L, maxIter = 1L)
  expect_length(res, 4L)
  expect_equal(attr(res, "iters"), 1L)
})

# ---------------------------------------------------------------------------
# 7. R reference loop and C++ port are bit-identical.
#
# .verify = TRUE routes through the original R reference loop (with brute-
# force record assertions); the default routes through DropAdd_cpp. The
# two paths must produce identical indices, objective, secondary, and iters
# for any iter budget.
# ---------------------------------------------------------------------------

test_that("DropAdd R reference loop and C++ port are bit-identical", {
  set.seed(2026)
  pts <- matrix(rnorm(60 * 4), ncol = 4)
  dmat <- as.matrix(dist(pts))

  for (maxIter in c(0L, 1L, 5L, 50L, 200L)) {
    outR <- DropAdd(dmat, m = 8L, maxIter = maxIter, .verify = TRUE)
    outC <- DropAdd(dmat, m = 8L, maxIter = maxIter)
    expect_identical(outR,                        outC)
    expect_identical(attr(outR, "objective"),  attr(outC, "objective"))
    expect_identical(attr(outR, "secondary"),  attr(outC, "secondary"))
    expect_identical(attr(outR, "iters"),      attr(outC, "iters"))
  }

  # Same parity under the stagnation stopping rule (plateau binds, no
  # iteration cap): the run-length is itself an output, so identical iters
  # confirms both paths take the deterministic criterion identically.
  for (mni in c(5L, 25L, 100L)) {
    outR <- DropAdd(dmat, m = 8L, plateau = mni, .verify = TRUE)
    outC <- DropAdd(dmat, m = 8L, plateau = mni)
    expect_identical(outR,                        outC)
    expect_identical(attr(outR, "objective"),  attr(outC, "objective"))
    expect_identical(attr(outR, "secondary"),  attr(outC, "secondary"))
    expect_identical(attr(outR, "iters"),      attr(outC, "iters"))
  }
})

# Tests for DropAddTS (Porumbel, Hao & Glover 2011).
#
# The validation oracle is Geo 100 1 with m=10: published best-known MaxMin
# objective is 89.37 and the instance is small enough that the algorithm
# should reach it within a 10-second budget on any modern machine.

# Helper: brute-force MaxMin objective for a selection.
.maxmin_brute <- function(S, dmat) {
  if (length(S) < 2L) return(NA_real_)
  sub <- dmat[S, S]
  diag(sub) <- Inf
  min(sub)
}

# Locate the bundled Geo loader. Skip tests that need it if it isn't reachable
# (e.g. installed-package CMD check without dev/ in tree).
.try_geo_loader <- function() {
  pkg_root <- tryCatch(rprojroot::find_root(rprojroot::is_r_package),
                       error = function(e) NULL)
  if (is.null(pkg_root)) return(NULL)
  loader <- file.path(pkg_root, "dev", "benchmarks", "mdplib", "loader.R")
  if (!file.exists(loader)) return(NULL)
  env <- new.env(parent = baseenv())
  source(loader, local = env, chdir = TRUE)
  # mdplib_root() relies on sys.frame()$ofile which becomes NULL after the
  # source() call returns; injecting an explicit absolute path here avoids
  # the lookup entirely.
  env$.MDPLIB_INSTANCE_DIR <- file.path(
    pkg_root, "dev", "benchmarks", "mdplib", "extracted", "instances"
  )
  env
}

# ---------------------------------------------------------------------------
# 1. Smoke test
# ---------------------------------------------------------------------------
test_that("DropAddTS smoke: 20 pts in 5-D, m=5", {
  set.seed(7)
  pts <- matrix(rnorm(20 * 5), ncol = 5)
  d <- dist(pts)
  res <- DropAddTS(d, m = 5L, time_budget_s = 1)
  expect_length(res$indices, 5L)
  expect_equal(length(unique(res$indices)), 5L)
  expect_true(all(res$indices %in% seq_len(20L)))
  expect_gt(res$objective, 0)
  expect_gte(res$iters, 5L)
  expect_true(res$time_s >= 0)
})

# ---------------------------------------------------------------------------
# 2. Construction-only sanity
# ---------------------------------------------------------------------------
test_that("DropAddTS construction-only result matches its MaxMin score", {
  set.seed(11)
  pts <- matrix(rnorm(30 * 4), ncol = 4)
  dmat <- as.matrix(dist(pts))
  res <- DropAddTS(dmat, m = 6L, max_iter = 0L)
  expect_length(res$indices, 6L)
  expect_equal(res$iters, 0L)
  # Objective stored equals the actual MaxMin over the returned indices.
  expect_equal(res$objective, .maxmin_brute(res$indices, dmat))
})

# ---------------------------------------------------------------------------
# 3. Improvement over construction on a real benchmark instance
# ---------------------------------------------------------------------------
test_that("DropAddTS never worsens the constructive solution", {
  geo_env <- .try_geo_loader()
  skip_if(is.null(geo_env), "Geo loader not available")
  path <- file.path(geo_env$.MDPLIB_INSTANCE_DIR, "Geo", "Geo 100 1.txt")
  skip_if_not(file.exists(path), "Geo 100 1 instance not present")
  geo <- geo_env$read_mdplib_geo(path)
  dmat <- geo_env$mdplib_geo_dist(geo)

  cons <- DropAddTS(dmat, m = 10L, max_iter = 0L)
  full <- DropAddTS(dmat, m = 10L, max_no_improve = 2000L)
  expect_gte(full$objective, cons$objective - 1e-9)
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
test_that("DropAddTS reaches the Geo 100 1 m=10 optimum (89.37)", {
  geo_env <- .try_geo_loader()
  skip_if(is.null(geo_env), "Geo loader not available")
  path <- file.path(geo_env$.MDPLIB_INSTANCE_DIR, "Geo", "Geo 100 1.txt")
  skip_if_not(file.exists(path), "Geo 100 1 instance not present")
  geo <- geo_env$read_mdplib_geo(path)
  dmat <- geo_env$mdplib_geo_dist(geo)

  # Time-budgeted oracle: run the full 10 s, as before, to confirm the search
  # reaches the proven optimum (this is a correctness oracle, not a frozen
  # result, so a wall-clock budget is appropriate here).
  res <- DropAddTS(dmat, m = 10L, time_budget_s = 10, max_no_improve = 100000000L)
  best_known <- 89.37
  expect_gte(res$objective, 0.999 * best_known)
  # Returned indices truly achieve the reported objective.
  expect_equal(res$objective, .maxmin_brute(res$indices, dmat),
               tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# 4b. Streamlined-record self-consistency under .verify
# ---------------------------------------------------------------------------
test_that("DropAddTS .verify=TRUE passes silently on a real instance", {
  geo_env <- .try_geo_loader()
  skip_if(is.null(geo_env), "Geo loader not available")
  path <- file.path(geo_env$.MDPLIB_INSTANCE_DIR, "Geo", "Geo 100 1.txt")
  skip_if_not(file.exists(path), "Geo 100 1 instance not present")
  geo <- geo_env$read_mdplib_geo(path)
  dmat <- geo_env$mdplib_geo_dist(geo)

  # Cap iters to keep verification cheap; the inner check is O(n*m) per iter.
  expect_silent(
    DropAddTS(dmat, m = 10L, time_budget_s = 30, max_iter = 100L,
              .verify = TRUE)
  )
})

# ---------------------------------------------------------------------------
# 5. FIFO invariant: across the first m main-loop iterations, every initially
# selected point is dropped exactly once.
# ---------------------------------------------------------------------------
test_that("DropAddTS FIFO drops each initial point once in first m iterations", {
  set.seed(2026)
  pts <- matrix(rnorm(10 * 3), ncol = 3)
  dmat <- as.matrix(dist(pts))
  m <- 4L

  # Capture the constructive S via the internal helper.
  cons <- MaxMin:::.DropAddConstruct(dmat, m)
  S_init <- cons$S
  expect_length(S_init, m)
  expect_equal(length(unique(S_init)), m)

  # Drive the PRODUCTION loop via the .trace hook for exactly 2*m iterations.
  trace_env <- new.env()
  res <- DropAddTS(dmat, m = m, time_budget_s = 60, max_iter = 2L * m,
                   .trace = trace_env)
  expect_equal(res$iters, 2L * m)
  expect_length(trace_env$drops, 2L * m)
  # First m drops are exactly the constructive members, in some order
  # (FIFO order is determined by iter_stamp, which was 1..m in construction
  # so drop order is the construction order).
  expect_setequal(trace_env$drops[seq_len(m)], S_init)
  expect_equal(trace_env$drops[seq_len(m)], S_init)  # actually in order

  # Across 2*m iterations, every dropped point was previously in S.
  # No index can be dropped twice in the first m iterations.
  expect_equal(anyDuplicated(trace_env$drops[seq_len(m)]), 0L)

  # Porumbel p.281: the just-dropped point is excluded from the add candidates
  # for that iteration, so the added point is never the one just dropped. This
  # is the diversification invariant; without it the search freezes.
  expect_true(all(trace_env$adds != trace_env$drops))
})

# ---------------------------------------------------------------------------
# 6. Time budget honoured on a 200-point Ran-format instance
# ---------------------------------------------------------------------------
test_that("DropAddTS respects time_budget_s within reasonable slack", {
  set.seed(99)
  pts <- matrix(runif(200 * 5), ncol = 5)
  dmat <- as.matrix(dist(pts))
  t0 <- Sys.time()
  # Disable stagnation so the wall-clock ceiling is the binding criterion.
  res <- DropAddTS(dmat, m = 20L, time_budget_s = 1, max_no_improve = 100000000L)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_lte(res$time_s, 1.5)
  expect_lte(elapsed, 2.0)
  expect_length(res$indices, 20L)
})

# ---------------------------------------------------------------------------
# 8. seed parameter, input validation, and progress output
# ---------------------------------------------------------------------------
test_that("DropAddTS seed parameter and input validation", {
  set.seed(11)
  dmat <- as.matrix(dist(matrix(rnorm(20 * 3), ncol = 3)))

  # seed path
  r1 <- DropAddTS(dmat, m = 4L, max_iter = 0L, seed = 1L)
  r2 <- DropAddTS(dmat, m = 4L, max_iter = 0L, seed = 1L)
  expect_identical(r1$indices, r2$indices)

  # m validation
  expect_error(DropAddTS(dmat, m = 1L),   "2 <= m")
  expect_error(DropAddTS(dmat, m = 25L),  "2 <= m")
  # max_no_improve validation
  expect_error(DropAddTS(dmat, m = 4L, max_no_improve = 0L),  "max_no_improve")
  # max_iter validation
  expect_error(DropAddTS(dmat, m = 4L, max_iter = -1L),  "max_iter")
  # time_budget_s validation
  expect_error(DropAddTS(dmat, m = 4L, time_budget_s = 0),   "time_budget_s")
  expect_error(DropAddTS(dmat, m = 4L, time_budget_s = NA_real_), "time_budget_s")
})

test_that("DropAddTS progress = TRUE fires the cli hooks", {
  dmat <- as.matrix(dist(matrix(rnorm(15 * 2), ncol = 2)))
  expect_no_error(DropAddTS(dmat, m = 3L, max_iter = 2L, progress = TRUE))
})

# ---------------------------------------------------------------------------
# 9. R-path corner cases: tie-breaking, caseA updates, m=2 else branch
# ---------------------------------------------------------------------------
test_that("DropAddTS construction tie-breaking and caseA update (.verify)", {
  # Unit square: all sides=1, diagonals=sqrt(2).
  # Construction for m=3: seed=P1, add diagonal P3 (caseA fires for P2 and P4),
  # then P2 and P4 are tied (lines 51-53 of .DropAddConstruct).
  pts_sq <- rbind(c(0,0), c(1,0), c(1,1), c(0,1))
  dmat   <- as.matrix(dist(pts_sq))
  # .verify=TRUE routes through R path with brute-force record checks
  res <- DropAddTS(dmat, m = 3L, max_iter = 4L, .verify = TRUE)
  expect_length(res$indices, 3L)
})

test_that("DropAddTS m=2 else-branch and ADD tie-breaking (.verify)", {
  # Rhombus: all edges=sqrt(2), diagonals=2.
  # m=2: after drop of seed, two candidates are equidistant from remaining
  # selected point (lines 360-366 else branch + lines 391-393 ADD tie).
  pts_rh <- rbind(c(0,0), c(1,1), c(2,0), c(1,-1))
  dmat   <- as.matrix(dist(pts_rh))
  res <- DropAddTS(dmat, m = 2L, max_iter = 4L, .verify = TRUE)
  expect_length(res$indices, 2L)
})

test_that("DropAddTS caseA ADD update fires in main loop (.verify)", {
  # Unit square m=3: first main-loop ADD of P4 has d[P4,P1]=1 = min_dist[P1]=1
  # (caseA, line 409).
  pts_sq <- rbind(c(0,0), c(1,0), c(1,1), c(0,1))
  dmat   <- as.matrix(dist(pts_sq))
  expect_no_error(
    DropAddTS(dmat, m = 3L, max_iter = 3L, .verify = TRUE)
  )
})

test_that("DropAddTS effective_max = 0 when m == n (.verify)", {
  dmat <- as.matrix(dist(matrix(rnorm(5 * 2), ncol = 2)))
  res  <- DropAddTS(dmat, m = 5L, .verify = TRUE)
  # All 5 points selected, no loop iterations.
  expect_length(res$indices, 5L)
  expect_equal(res$iters, 0L)
})

test_that("DropAddTS .trace + .verify together cover init and loop writes", {
  set.seed(2026)
  dmat <- as.matrix(dist(matrix(rnorm(12 * 2), ncol = 2)))
  te   <- new.env()
  res  <- DropAddTS(dmat, m = 4L, max_iter = 5L, .verify = TRUE, .trace = te)
  # .trace should have been initialised on R path (lines 295-296) and
  # written during the loop (lines 438-439).
  expect_equal(length(te$drops), res$iters)
  expect_equal(length(te$adds),  res$iters)
})

# ---------------------------------------------------------------------------
# 7. R reference loop and C++ port are bit-identical on a fixed iter budget.
#
# .verify = TRUE routes through the original R reference loop (with brute-
# force record assertions); the default routes through DropAddTS_cpp. The
# two paths must produce identical indices, objective, secondary, and iters
# for any iter budget.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 10. C++ path coverage: tie-breaking and min_dist-count branches
#
# Tests 10a-10c exercise specific branches in dropadd.cpp that require exact
# distance ties; they use the default (C++) path (no .verify).
# ---------------------------------------------------------------------------

test_that("DropAddTS C++ construction covers sum_dist tie-break (lines 82-85)", {
  # 5-point geometry: A=(0,0), B=(3,0), C=(0,2), P=(1,1), Q=(2,1).
  # Construction: B (max row-sum) -> C -> A.
  # Then P and Q both have min_dist=sqrt(2) to {B,C,A}, but
  # sum_dist[Q]=sqrt(2)+2*sqrt(5) > sum_dist[P]=2*sqrt(2)+sqrt(5),
  # so the tie-break updates x_new to Q (lines 83-85 fire).
  # When A is added, d(P,A)=sqrt(2)==min_dist[P] -> line 103 fires too.
  pts <- rbind(c(0,0), c(3,0), c(0,2), c(1,1), c(2,1))
  d   <- as.matrix(dist(pts))
  res <- DropAddTS(d, m = 4L, max_iter = 0L)
  expect_length(res$indices, 4L)
  expect_equal(res$iters, 0L)
})

test_that("DropAddTS C++ m=2 covers DROP else-branch (lines 214-215)", {
  # Rhombus: (0,0),(1,1),(2,0),(1,-1).  All edges=sqrt(2), diagonals=2.
  # Construction: {(0,0),(2,0)} (endpoints of the longer diagonal).
  # Drop (0,0): the remaining selected (2,0) had (0,0) as its only peer,
  # so min_dist_count[(2,0)] drops to 0.  In the recompute, the self-mask
  # skips the only surviving S member, leaving mns=Inf -> else branch fires.
  pts_rh <- rbind(c(0,0), c(1,1), c(2,0), c(1,-1))
  d_rh   <- as.matrix(dist(pts_rh))
  res <- DropAddTS(d_rh, m = 2L, max_iter = 4L)
  expect_length(res$indices, 2L)
})

test_that("DropAddTS C++ main-loop ADD covers tie-break (255-258) and equality (277)", {
  # 7-point geometry ensures:
  #   * Lines 255-258: after dropping ANC, X=(1,0) and Y=(3,0) tie on
  #     min_dist=1 but sum_dist[Y]>sum_dist[X], so Y wins via the update.
  #   * Line 277: when Y is added, d(Y,W)=0.5=min_dist[W] -> count increments.
  #
  # Points: P(0,0), Q(4,0), R(0,4), ANC(10,10), X(1,0), Y(3,0), W(3.5,0).
  # Construction: ANC -> P -> Q -> R  (ANC has max row-sum).
  pts7 <- rbind(c(0,0), c(4,0), c(0,4), c(10,10), c(1,0), c(3,0), c(3.5,0))
  d7   <- as.matrix(dist(pts7))
  res <- DropAddTS(d7, m = 4L, max_iter = 1L)
  expect_length(res$indices, 4L)
  expect_equal(res$iters, 1L)
})

# ---------------------------------------------------------------------------
# 7. R reference loop and C++ port are bit-identical.
#
# .verify = TRUE routes through the original R reference loop (with brute-
# force record assertions); the default routes through DropAddTS_cpp. The
# two paths must produce identical indices, objective, secondary, and iters
# for any iter budget.
# ---------------------------------------------------------------------------

test_that("DropAddTS R reference loop and C++ port are bit-identical", {
  set.seed(2026)
  pts <- matrix(rnorm(60 * 4), ncol = 4)
  dmat <- as.matrix(dist(pts))

  for (max_iter in c(0L, 1L, 5L, 50L, 200L)) {
    out_R <- DropAddTS(dmat, m = 8L, max_iter = max_iter, .verify = TRUE)
    out_C <- DropAddTS(dmat, m = 8L, max_iter = max_iter)
    expect_identical(out_R$indices,   out_C$indices)
    expect_identical(out_R$objective, out_C$objective)
    expect_identical(out_R$secondary, out_C$secondary)
    expect_identical(out_R$iters,     out_C$iters)
  }

  # Same parity under the stagnation stopping rule (max_no_improve binds, no
  # iteration cap): the run-length is itself an output, so identical iters
  # confirms both paths take the deterministic criterion identically.
  for (mni in c(5L, 25L, 100L)) {
    out_R <- DropAddTS(dmat, m = 8L, max_no_improve = mni, .verify = TRUE)
    out_C <- DropAddTS(dmat, m = 8L, max_no_improve = mni)
    expect_identical(out_R$indices,   out_C$indices)
    expect_identical(out_R$objective, out_C$objective)
    expect_identical(out_R$secondary, out_C$secondary)
    expect_identical(out_R$iters,     out_C$iters)
  }
})

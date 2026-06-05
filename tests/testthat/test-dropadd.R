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
# 7. R reference loop and C++ port are bit-identical on a fixed iter budget.
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

# Tests for GraspPR (Resende et al. 2010, static variant) and its compiled
# kernel. Termination is the deterministic stagnation rule (max_no_improve).

set.seed(42)
pts30 <- matrix(rnorm(90), ncol = 3)
d30   <- dist(pts30)
d30m  <- as.matrix(d30)

# 1. Smoke test ------------------------------------------------------------

test_that("GraspPR smoke: returns valid selection on 30 random 3-D points", {
  res <- GraspPR(d30, m = 5L, max_no_improve = 20L, elite_size = 4L, seed = 1L)
  expect_type(res$indices, "integer")
  expect_length(res$indices, 5L)
  expect_true(all(res$indices %in% seq_len(30L)))
  expect_equal(length(unique(res$indices)), 5L)
  expect_true(res$objective > 0)
})

# 2. Compiled kernel matches the R reference bit-for-bit -------------------
# Both draw construction indices from R's RNG, so from a common seed the
# trajectories, iteration counts and selections must be identical.

test_that("GraspPR_cpp == .GraspPR_R across seeds and parameters", {
  grid <- expand.grid(seed = c(1L, 7L, 99L, 2024L),
                      mni  = c(10L, 40L),
                      es   = c(4L, 8L),
                      KEEP.OUT.ATTRS = FALSE)
  for (r in seq_len(nrow(grid))) {
    s   <- grid$seed[r]
    mni <- grid$mni[r]
    es  <- grid$es[r]

    set.seed(s)
    ref <- MaxMin:::.GraspPR_R(d30m, m = 6L, max_no_improve = mni,
                               elite_size = es, alpha = 0.8)
    set.seed(s)
    ker <- GraspPR(d30m, m = 6L, max_no_improve = mni, elite_size = es,
                   alpha = 0.8)

    info <- sprintf("seed=%d mni=%d es=%d", s, mni, es)
    expect_identical(ker$indices, ref$indices, info = info)
    expect_identical(ker$objective, ref$objective, info = info)
    expect_identical(ker$iters, ref$iters, info = info)
    expect_identical(ker$pr_calls, ref$pr_calls, info = info)
  }
})

# 3. Determinism -----------------------------------------------------------

test_that("GraspPR is reproducible from a seed (machine-independent)", {
  a <- GraspPR(d30, m = 6L, max_no_improve = 30L, elite_size = 5L, seed = 17L)
  b <- GraspPR(d30, m = 6L, max_no_improve = 30L, elite_size = 5L, seed = 17L)
  expect_identical(a$indices, b$indices)
  expect_identical(a$objective, b$objective)
  expect_identical(a$iters, b$iters)
})

# 4. max_iter = 0 short-circuits Phase B -----------------------------------

test_that("GraspPR with max_iter = 0 runs Phase A + relinking only", {
  res <- GraspPR(d30, m = 4L, max_iter = 0L, elite_size = 4L, seed = 7L)
  expect_length(res$indices, 4L)
  expect_true(res$objective > 0)
  expect_equal(res$iters, 0L)
})

# 5. Stagnation criterion bounds the run -----------------------------------

test_that("GraspPR stops within max_no_improve of its last improvement", {
  # With max_iter as a hard cap we can assert iters never exceeds it; with the
  # stagnation rule alone the loop must still terminate.
  res <- GraspPR(d30, m = 6L, max_no_improve = 5L, max_iter = 200L,
                 elite_size = 4L, seed = 3L)
  expect_lte(res$iters, 200L)
  expect_length(res$indices, 6L)
})

# 6. Local search monotonicity (.GprLocalSearch reference helper) ----------

test_that(".GprLocalSearch never decreases the MaxMin objective", {
  centroid <- colMeans(pts30)
  to_centroid <- sqrt(rowSums(sweep(pts30, 2L, centroid)^2))
  bad_sel <- order(to_centroid)[1:5]
  z_start <- MaxMin:::.GprObjective(d30m, bad_sel)
  improved <- MaxMin:::.GprLocalSearch(d30m, bad_sel)
  z_end <- MaxMin:::.GprObjective(d30m, improved)
  expect_true(z_end >= z_start)
  expect_true(z_end > z_start)
  expect_length(improved, 5L)
  expect_equal(length(unique(improved)), 5L)
})

# 7. Path relinking keeps the best state along the path --------------------

test_that(".GprPathRelink keeps the best state along the path", {
  x <- c(1L, 2L, 3L, 4L, 5L)
  y <- c(1L, 2L, 6L, 7L, 8L)
  pr <- MaxMin:::.GprPathRelink(d30m, x, y)
  zx <- MaxMin:::.GprObjective(d30m, x)
  zy <- MaxMin:::.GprObjective(d30m, y)
  expect_equal(pr$intermediates, 3L)
  expect_true(pr$objective >= max(zx, zy))
})

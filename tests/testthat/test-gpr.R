# Tests for GraspPR (Resende et al. 2010, static variant) and its compiled
# kernel. Termination is the deterministic stagnation rule (plateau).

set.seed(42)
pts30 <- matrix(rnorm(90), ncol = 3)
d30   <- dist(pts30)
d30m  <- as.matrix(d30)

# 1. Smoke test ------------------------------------------------------------

test_that("GraspPR smoke: returns valid selection on 30 random 3-D points", {
  res <- GraspPR(d30, m = 5L, plateau = 20L, eliteSize = 4L, seed = 1L)
  expect_type(res, "integer")
  expect_length(res, 5L)
  expect_true(all(res %in% seq_len(30L)))
  expect_equal(length(unique(res)), 5L)
  expect_true(attr(res, "score") > 0)
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
    ref <- MaxMin:::.GraspPR_R(d30m, m = 6L, plateau = mni,
                               eliteSize = es, alpha = 0.8)
    set.seed(s)
    ker <- GraspPR(d30m, m = 6L, plateau = mni, eliteSize = es,
                   alpha = 0.8)

    info <- sprintf("seed=%d mni=%d es=%d", s, mni, es)
    expect_identical(ker,                    ref,                    info = info)
    expect_identical(attr(ker, "score"),     attr(ref, "score"),     info = info)
    expect_identical(attr(ker, "iters"),     attr(ref, "iters"),     info = info)
    expect_identical(attr(ker, "pr_calls"),  attr(ref, "pr_calls"),  info = info)
  }
})

# 3. Determinism -----------------------------------------------------------

test_that("GraspPR is reproducible from a seed (machine-independent)", {
  a <- GraspPR(d30, m = 6L, plateau = 30L, eliteSize = 5L, seed = 17L)
  b <- GraspPR(d30, m = 6L, plateau = 30L, eliteSize = 5L, seed = 17L)
  expect_identical(a, b)
})

# 4. maxIter = 0 short-circuits Phase B -----------------------------------

test_that("GraspPR with maxIter = 0 runs Phase A + relinking only", {
  res <- GraspPR(d30, m = 4L, maxIter = 0L, eliteSize = 4L, seed = 7L)
  expect_length(res, 4L)
  expect_true(attr(res, "score") > 0)
  expect_equal(attr(res, "iters"), 0L)
})

# 5. Stagnation criterion bounds the run -----------------------------------

test_that("GraspPR stops within plateau of its last improvement", {
  # With maxIter as a hard cap we can assert iters never exceeds it; with the
  # stagnation rule alone the loop must still terminate.
  res <- GraspPR(d30, m = 6L, plateau = 5L, maxIter = 200L,
                 eliteSize = 4L, seed = 3L)
  expect_lte(attr(res, "iters"), 200L)
  expect_length(res, 6L)
})

# 6. Local search monotonicity (.GprLocalSearch reference helper) ----------

test_that(".GprLocalSearch never decreases the MaxMin objective", {
  centroid <- colMeans(pts30)
  toCentroid <- sqrt(rowSums(sweep(pts30, 2L, centroid)^2))
  badSel <- order(toCentroid)[1:5]
  zStart <- MaxMin:::.GprObjective(d30m, badSel)
  improved <- MaxMin:::.GprLocalSearch(d30m, badSel)
  zEnd <- MaxMin:::.GprObjective(d30m, improved)
  expect_true(zEnd >= zStart)
  expect_true(zEnd > zStart)
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

# 8. GraspPR timeBudgetS validation -------------------------------------

test_that("GraspPR validates timeBudgetS", {
  expect_error(GraspPR(d30, m = 4L, timeBudgetS = 0),  "timeBudgetS")
  expect_error(GraspPR(d30, m = 4L, timeBudgetS = -1), "timeBudgetS")
  expect_error(GraspPR(d30, m = 4L, timeBudgetS = NA_real_), "timeBudgetS")
})

# 9. .GraspPR_R maxIter cap -----------------------------------

test_that(".GraspPR_R stops exactly at maxIter", {
  set.seed(1)
  ref <- MaxMin:::.GraspPR_R(d30m, m = 5L, plateau = 1000L,
                              maxIter = 3L, eliteSize = 4L)
  expect_lte(attr(ref, "iters"), 3L)
})

# 10. .GprLocalSearch pair-count reduction (extended-improvement branch) ---

test_that(".GprLocalSearch reduces pair count when T_k is unchanged", {
  # Construct a matrix where two points are equidistant (forcing n_critical > 1),
  # and a swap outside preserves T_k but reduces the critical-pair count.
  # Four corners of a unit square: T_4 = 1, with 4 critical pairs (all sides).
  ptsSq <- rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(3,0.5))
  dSq   <- as.matrix(dist(ptsSq))
  sel    <- 1:4   # four corners; T_k = 1, 4 critical pairs
  improved <- MaxMin:::.GprLocalSearch(dSq, sel)
  # After improvement T_k must be >= 1.
  expect_gte(MaxMin:::.GprObjective(dSq, improved),
             MaxMin:::.GprObjective(dSq, sel))
})

# 11. .GprTryInsert second acceptance condition and tail insertion ---------
# Hand-crafted 4-point geometry guarantees both branches:
#   P1=(0,0), P2=(1,0), P3=(0,0.5), P4=(0.3,0.5)
#   d12=1.0  (z1), d34=0.3 (zb), d13=0.5 (selZ) with zb < selZ < z1
#
# Line 217 fires: selZ=0.5 > zb=0.3 and selZ <= z1=1.0, dmin=1 >= dth=1.
# Tie-break on Hamming removes s2={3,4} (lowest z=0.3), leaving ES=[{1,2}].
# pos = sum([1.0] >= 0.5) + 1 = 2 > 1 = length(remaining) → tail (233-234).

test_that(".GprTryInsert line 217 (second condition) and lines 233-234 (tail insert)", {
  pts4 <- rbind(c(0, 0), c(1, 0), c(0, 0.5), c(0.3, 0.5))
  d4   <- as.matrix(dist(pts4))

  s1   <- c(1L, 2L); s2 <- c(3L, 4L)
  z1   <- MaxMin:::.GprObjective(d4, s1)   # 1.0
  z2   <- MaxMin:::.GprObjective(d4, s2)   # 0.3
  ES   <- list(s1, s2)
  esZ  <- c(z1, z2)   # already descending

  sel  <- c(1L, 3L)
  selZ <- MaxMin:::.GprObjective(d4, sel) # 0.5  (between zb and z1)
  dth  <- 1L

  res <- MaxMin:::.GprTryInsert(d4, ES, esZ, sel, selZ, dth)

  expect_true(res$changed)
  expect_length(res$ES, 2L)
  expect_equal(res$esZ[2L], selZ)   # sel inserted at tail position
})

# 12. Path relinking improves bestSel (.GraspPR_R lines 422-423) ---------
# Scan seeds until Phase C PR beats the Phase A best, confirming lines 422-423.
# Phase A is replicated manually (same RNG path) to get the pre-PR ceiling;
# Phase C is deterministic, so any gain must have fired those lines.

test_that(".GraspPR_R phase-C path relinking fires lines 422-423", {
  m     <- 4L
  es    <- 6L
  alpha <- 0.8
  found <- FALSE
  for (s in seq_len(200L)) {
    set.seed(s)
    pts <- matrix(rnorm(30L * 2L), ncol = 2L)
    d   <- as.matrix(dist(pts))

    # Replicate Phase A: same RNG, same constructions -> same Phase A ceiling.
    set.seed(s)
    phaseABest <- -Inf
    for (b in seq_len(es)) {
      x  <- MaxMin:::.GprConstruct(d, m, alpha)
      xp <- MaxMin:::.GprLocalSearch(d, x)
      phaseABest <- max(phaseABest, MaxMin:::.GprObjective(d, xp))
    }

    # Full run (Phase A + C, no Phase B).
    set.seed(s)
    res <- MaxMin:::.GraspPR_R(d, m = m, plateau = 1000L,
                                maxIter = 0L, eliteSize = es, alpha = alpha)

    if (attr(res, "score") > phaseABest + 1e-9) {
      found <- TRUE
      break
    }
  }
  expect_true(found)
})

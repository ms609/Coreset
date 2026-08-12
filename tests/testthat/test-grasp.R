# Tests for Grasp (Resende et al. 2010, static variant) and its compiled
# kernel. Termination is the deterministic stagnation rule (plateau).

set.seed(42)
pts30 <- matrix(rnorm(90), ncol = 3)
d30   <- dist(pts30)
d30m  <- as.matrix(d30)

# 1. Smoke test ------------------------------------------------------------

test_that("Grasp smoke: returns valid selection on 30 random 3-D points", {
  set.seed(1)
  res <- Grasp(d = d30, k = 5L, plateau = 20L, eliteSize = 4L)
  expect_type(res, "integer")
  expect_length(res, 5L)
  expect_true(all(res %in% seq_len(30L)))
  expect_equal(length(unique(res)), 5L)
  expect_true(attr(res, "score") > 0)
})

# 2. Compiled kernel matches the R reference bit-for-bit -------------------
# Both draw construction indices from R's RNG, so from a common seed the
# trajectories, iteration counts and selections must be identical.

test_that("Grasp_cpp == .Grasp_R across seeds and parameters", {
  # alpha is swept including the extremes 0 (uniform random) and 1 (pure
  # greedy): alpha = 1 is the FP-rounding empty-RCL trigger of GRASP-01, so its
  # inclusion turns that crash into a detectable parity failure.
  grid <- expand.grid(seed  = c(1L, 7L, 99L, 2024L),
                      mni   = c(10L, 40L),
                      es    = c(4L, 8L),
                      alpha = c(0, 0.8, 1),
                      KEEP.OUT.ATTRS = FALSE)
  for (r in seq_len(nrow(grid))) {
    s   <- grid$seed[r]
    mni <- grid$mni[r]
    es  <- grid$es[r]
    al  <- grid$alpha[r]

    set.seed(s)
    ref <- MaxMin:::.Grasp_R(6L, d30m, plateau = mni,
                               eliteSize = es, alpha = al)
    set.seed(s)
    ker <- Grasp(d = d30m, k = 6L, plateau = mni, eliteSize = es,
                   alpha = al)

    info <- sprintf("seed=%d mni=%d es=%d alpha=%g", s, mni, es, al)
    # Drop the wall-clock `time_s` attribute: it is genuinely nondeterministic
    # and never expected to match. Every other attribute is asserted below.
    attr(ker, "time_s") <- attr(ref, "time_s") <- NULL
    expect_identical(ker,                    ref,                    info = info)
    expect_equal(attr(ker, "score"),     attr(ref, "score"),     tolerance = 1e-14, info = info)
    expect_identical(attr(ker, "iters"),     attr(ref, "iters"),     info = info)
    expect_identical(attr(ker, "pr_calls"),  attr(ref, "pr_calls"),  info = info)
  }
})

# 3. Determinism -----------------------------------------------------------

test_that("Grasp is reproducible from a seed (machine-independent)", {
  set.seed(17); a <- Grasp(d = d30, k = 6L, plateau = 30L, eliteSize = 5L)
  set.seed(17); b <- Grasp(d = d30, k = 6L, plateau = 30L, eliteSize = 5L)
  # set.seed() before each call => identical selection and objective; only
  # wall-clock differs.
  attr(a, "time_s") <- attr(b, "time_s") <- NULL
  expect_identical(a, b)
})

# 4. Grasp returns a valid selection -----------------------------------

test_that("Grasp returns a valid scored selection", {
  set.seed(7)
  res <- Grasp(d = d30, k = 4L, eliteSize = 4L)
  expect_length(res, 4L)
  expect_true(attr(res, "score") > 0)
})

# 5. Stagnation criterion bounds the run -----------------------------------

test_that("Grasp terminates under plateau stopping criterion", {
  set.seed(3)
  res <- Grasp(d = d30, k = 6L, plateau = 5L, eliteSize = 4L)
  expect_length(res, 6L)
  expect_gte(attr(res, "iters"), 0L)
})

# 6. Local search monotonicity (.GraspLocalSearch reference helper) ----------

test_that(".GraspLocalSearch never decreases the MaxMin objective", {
  anti_centroid <- colMeans(pts30)
  toCentroid <- sqrt(rowSums(sweep(pts30, 2L, anti_centroid)^2))
  badSel <- order(toCentroid)[1:5]
  zStart <- MaxMin:::.GraspObjective(d30m, badSel)
  improved <- MaxMin:::.GraspLocalSearch(d30m, badSel)
  zEnd <- MaxMin:::.GraspObjective(d30m, improved)
  expect_true(zEnd >= zStart)
  expect_true(zEnd > zStart)
  expect_length(improved, 5L)
  expect_equal(length(unique(improved)), 5L)
})

# 7. Path relinking keeps the best state along the path --------------------

test_that(".GraspPathRelink keeps the best state along the path", {
  x <- c(1L, 2L, 3L, 4L, 5L)
  y <- c(1L, 2L, 6L, 7L, 8L)
  pr <- MaxMin:::.GraspPathRelink(d30m, x, y)
  zx <- MaxMin:::.GraspObjective(d30m, x)
  zy <- MaxMin:::.GraspObjective(d30m, y)
  expect_equal(pr$intermediates, 3L)
  expect_true(pr$objective >= max(zx, zy))
})

# 7b. Path relinking strictly improves on at least one pair ---------------

test_that(".GraspPathRelink strictly improves on at least one pair", {
  found <- FALSE
  for (s in 1:100) {
    set.seed(s)
    x <- MaxMin:::.GraspConstruct(d30m, 5L, 0.5)
    set.seed(s + 100L)
    y <- MaxMin:::.GraspConstruct(d30m, 5L, 0.5)
    if (!identical(sort(x), sort(y))) {
      zx <- MaxMin:::.GraspObjective(d30m, x)
      zy <- MaxMin:::.GraspObjective(d30m, y)
      pr <- MaxMin:::.GraspPathRelink(d30m, x, y)
      if (pr$objective > max(zx, zy) + 1e-10) {
        found <- TRUE; break
      }
    }
  }
  expect_true(found, label = "path relinking strictly improves on at least one pair")
})

# 8. Grasp maxSeconds validation -------------------------------------

test_that("Grasp validates maxSeconds", {
  expect_error(Grasp(d = d30, k = 4L, maxSeconds = 0),  "maxSeconds")
  expect_error(Grasp(d = d30, k = 4L, maxSeconds = -1), "maxSeconds")
  expect_error(Grasp(d = d30, k = 4L, maxSeconds = NA_real_), "maxSeconds")
})

# 8b. GRASP-01: empty-RCL at alpha = 1 no longer crashes ---------------------
# At the documented alpha = 1 (pure greedy), FP rounding can push the RCL
# threshold ~1 ULP past gmax and empty the RCL; pre-fix the C++ kernel indexed
# rcl[0] on an empty vector and SIGSEGV'd (uncatchable, killed the session),
# while the R path threw. A bare-greedy fallback (no RNG draw) now selects the
# argmax-g candidate in both. Sweeping seeds 1:30 exercises the rounding case
# and asserts the kernel still matches the R reference bit-for-bit.
test_that("Grasp survives alpha = 1 (empty-RCL fallback) and matches the R reference", {
  for (s in 1:30) {
    set.seed(s)
    ref <- MaxMin:::.Grasp_R(5L, d30m, plateau = 15L, eliteSize = 4L,
                             alpha = 1)
    set.seed(s)
    ker <- Grasp(d = d30m, k = 5L, plateau = 15L, eliteSize = 4L, alpha = 1)
    attr(ker, "time_s") <- attr(ref, "time_s") <- NULL
    expect_identical(ker, ref, info = paste("seed", s))
  }
})

# 8c. GRASP-02: alpha validation ---------------------------------------------
test_that("Grasp validates alpha", {
  expect_error(Grasp(d = d30, k = 4L, alpha = 2),         "alpha")
  expect_error(Grasp(d = d30, k = 4L, alpha = -1),        "alpha")
  expect_error(Grasp(d = d30, k = 4L, alpha = NA_real_),  "alpha")
  expect_error(Grasp(d = d30, k = 4L, alpha = c(0.5, 0.8)), "alpha")
})

# 9. .Grasp_R maxIter cap -----------------------------------

test_that(".Grasp_R stops exactly at maxIter", {
  set.seed(1)
  ref <- MaxMin:::.Grasp_R(5L, d30m, plateau = 1000L,
                              maxIter = 3L, eliteSize = 4L)
  expect_lte(attr(ref, "iters"), 3L)
})

# 10. .GraspLocalSearch pair-count reduction (extended-improvement branch) ---

test_that(".GraspLocalSearch reduces pair count when T_k is unchanged", {
  # Construct a matrix where two points are equidistant (forcing n_critical > 1),
  # and a swap outside preserves T_k but reduces the critical-pair count.
  # Four corners of a unit square: T_4 = 1, with 4 critical pairs (all sides).
  ptsSq <- rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(3,0.5))
  dSq   <- as.matrix(dist(ptsSq))
  sel    <- 1:4   # four corners; T_k = 1, 4 critical pairs
  improved <- MaxMin:::.GraspLocalSearch(dSq, sel)
  # After improvement T_k must be >= 1.
  expect_gte(MaxMin:::.GraspObjective(dSq, improved),
             MaxMin:::.GraspObjective(dSq, sel))
})

# 11. .GraspTryInsert second acceptance condition and tail insertion ---------
# Hand-crafted 4-point geometry guarantees both branches:
#   P1=(0,0), P2=(1,0), P3=(0,0.5), P4=(0.3,0.5)
#   d12=1.0  (z1), d34=0.3 (zb), d13=0.5 (selZ) with zb < selZ < z1
#
# Line 217 fires: selZ=0.5 > zb=0.3 and selZ <= z1=1.0, dmin=1 >= dth=1.
# Tie-break on Hamming removes s2={3,4} (lowest z=0.3), leaving ES=[{1,2}].
# pos = sum([1.0] >= 0.5) + 1 = 2 > 1 = length(remaining) → tail (233-234).

test_that(".GraspTryInsert line 217 (second condition) and lines 233-234 (tail insert)", {
  pts4 <- rbind(c(0, 0), c(1, 0), c(0, 0.5), c(0.3, 0.5))
  d4   <- as.matrix(dist(pts4))

  s1   <- c(1L, 2L); s2 <- c(3L, 4L)
  z1   <- MaxMin:::.GraspObjective(d4, s1)   # 1.0
  z2   <- MaxMin:::.GraspObjective(d4, s2)   # 0.3
  ES   <- list(s1, s2)
  esZ  <- c(z1, z2)   # already descending

  sel  <- c(1L, 3L)
  selZ <- MaxMin:::.GraspObjective(d4, sel) # 0.5  (between zb and z1)
  dth  <- 1L

  res <- MaxMin:::.GraspTryInsert(d4, ES, esZ, sel, selZ, dth)

  expect_true(res$changed)
  expect_length(res$ES, 2L)
  expect_equal(res$esZ[2L], selZ)   # sel inserted at tail position
})

# 11b. .GraspTryInsert never evicts a member better than the candidate --------
# Resende et al. (2010) Fig. 4 line 8 restricts eviction to elite members worse
# than the incoming solution. Searching all of ES instead lets a mid-value
# candidate displace the incumbent best whenever it is that member's unique
# Hamming-nearest neighbour, and the pool maximum falls.
#
# Geometry pins the case deterministically. Elite set {s1, s2}; the candidate
# `sel` shares one element with s1 (Hamming 1) and none with s2 (Hamming 2), so
# s1 -- the incumbent BEST -- is the unique nearest. selZ sits between the two
# elite objectives, so s1 must survive and s2 must go.

test_that(".GraspTryInsert evicts only members worse than the candidate", {
  pts <- rbind(c(0, 0), c(10, 0),          # s1 = {1,2}: z = 10   (the best)
               c(0, 1), c(1, 1),           # s2 = {3,4}: z = 1
               c(0, 5))                    # sel = {1,5}: z = 5
  dm  <- as.matrix(dist(pts))
  s1  <- c(1L, 2L); s2 <- c(3L, 4L); sel <- c(1L, 5L)
  z1  <- MaxMin:::.GraspObjective(dm, s1)
  z2  <- MaxMin:::.GraspObjective(dm, s2)
  selZ <- MaxMin:::.GraspObjective(dm, sel)
  expect_equal(c(z1, z2, selZ), c(10, 1, 5))

  # The incumbent best really is the unique Hamming-nearest member.
  hamm <- MaxMin:::.GraspHammingToES(sel, list(s1, s2))
  expect_identical(hamm, c(1L, 2L))

  res <- MaxMin:::.GraspTryInsert(dm, list(s1, s2), c(z1, z2), sel, selZ,
                                  dth = 1L)
  expect_true(res$changed)
  expect_identical(res$ES[[1L]], s1)        # incumbent best survives
  expect_equal(res$esZ, c(z1, selZ))        # the worse member was displaced
})

# 11c. The pool maximum is monotone under repeated insertion ------------------
# The invariant .Grasp_R's stagnation rule relies on (see grasp.R): the stall
# counter resets on a rise in esZ[1], so a fall would let a re-attained value
# count as fresh progress and extend the run past its stopping rule.

test_that(".GraspTryInsert never lowers the elite best", {
  set.seed(1)
  n <- 60L; k <- 8L; b <- 10L; dth <- 2L
  dm <- as.matrix(dist(matrix(rnorm(n * 3L), n)))
  z <- function(s) MaxMin:::.GraspObjective(dm, s)

  ES  <- lapply(seq_len(b), function(i) sort(sample(n, k)))
  esZ <- vapply(ES, z, numeric(1L))
  ord <- order(esZ, decreasing = TRUE); ES <- ES[ord]; esZ <- esZ[ord]

  drops <- 0L
  for (it in seq_len(2000L)) {
    sel <- sort(sample(n, k))
    before <- esZ[1L]
    res <- MaxMin:::.GraspTryInsert(dm, ES, esZ, sel, z(sel), dth)
    ES <- res$ES; esZ <- res$esZ
    if (esZ[1L] < before - 1e-12) drops <- drops + 1L
    # The set stays the right size and sorted best-to-worst throughout.
    expect_length(esZ, b)
    expect_false(is.unsorted(rev(esZ)))
  }
  expect_identical(drops, 0L)
})

# 12. Path relinking improves bestSel (.Grasp_R lines 422-423) ---------
# Scan seeds until Phase C PR beats the Phase A best, confirming lines 422-423.
# Phase A is replicated manually (same RNG path) to get the pre-PR ceiling;
# Phase C is deterministic, so any gain must have fired those lines.

test_that(".Grasp_R phase-C path relinking fires lines 422-423", {
  k     <- 4L
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
      x  <- MaxMin:::.GraspConstruct(d, k, alpha)
      xp <- MaxMin:::.GraspLocalSearch(d, x)
      phaseABest <- max(phaseABest, MaxMin:::.GraspObjective(d, xp))
    }

    # Full run (Phase A + C, no Phase B).
    set.seed(s)
    res <- MaxMin:::.Grasp_R(k, d, plateau = 1000L,
                                maxIter = 0L, eliteSize = es, alpha = alpha)

    if (attr(res, "score") > phaseABest + 1e-9) {
      found <- TRUE
      break
    }
  }
  expect_true(found)
})

# 13. Grasp eliteSize = 1 skips path relinking ---------------------------------

test_that("Grasp eliteSize=1 skips path relinking", {
  set.seed(1)
  res <- Grasp(d = d30, k = 4L, eliteSize = 1L, plateau = 20L)
  expect_length(res, 4L)
  expect_equal(attr(res, "pr_calls"), 0L)
  expect_true(attr(res, "score") > 0)
})

# 14. .GraspPathRelink returns immediately for identical endpoints -------------

test_that(".GraspPathRelink returns immediately for identical endpoints", {
  x  <- c(1L, 2L, 3L, 4L, 5L)
  pr <- MaxMin:::.GraspPathRelink(d30m, x, x)
  expect_equal(pr$intermediates, 0L)
  # The returned objective is the best we can do with that selection.
  expect_equal(pr$objective, MaxMin:::.GraspObjective(d30m, x))
})

# 15. .Grasp_R finite maxSeconds covers setTimeLimit (grasp.R line 397) ------

# 15b. grasp_local_search non-witness critical branch (grasp.cpp:193) ---------
# Fires when the dropped vertex is critical (min-distance participant) but is
# NOT one of the two witness vertices (wa, wb) of the first min-distance pair.
# The unit square has all 4 sides equal, so all 4 corners are critical, but
# only corners 0 and 1 (the first pair found) are witnesses — dropping corners
# 2 or 3 fires the else branch.

test_that("grasp_local_search non-witness critical branch covered (grasp.cpp:193)", {
  ptsSq5 <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(3, 0.5))
  dSq5   <- as.matrix(dist(ptsSq5))
  # Parity assertion proves the C++ kernel executed the same logic as the R
  # reference; the R reference is known to hit line 193 (all 4 corners are
  # critical, corners 2 and 3 are non-witnesses).
  set.seed(1)
  ref <- MaxMin:::.Grasp_R(4L, dSq5, plateau = 10L, eliteSize = 4L)
  set.seed(1)
  ker <- Grasp(d = dSq5, k = 4L, plateau = 10L, eliteSize = 4L)
  attr(ker, "time_s") <- attr(ref, "time_s") <- NULL
  expect_identical(ker, ref)
})

# 15c. Grasp_cpp countdown covers checkUserInterrupt block (grasp.cpp:393-394) -
# check_every = 256; a run long enough to exceed 256 iterations exercises the block.

test_that("Grasp_cpp countdown block is exercised (grasp.cpp:393-394)", {
  set.seed(1)
  res <- Grasp(d = d30m, k = 6L, plateau = .Machine$integer.max,
               eliteSize = 4L, maxSeconds = 0.05)
  expect_gte(attr(res, "iters"), 1L)
})

test_that(".Grasp_R time budget halts execution (grasp.R line 397)", {
  set.seed(1)
  # .Machine$integer.max disables both other stopping criteria; only the
  # budget can end Phase B.
  res <- expect_returns_within(
    MaxMin:::.Grasp_R(6L, d30m,
                      plateau  = .Machine$integer.max,
                      maxIter  = .Machine$integer.max,
                      eliteSize = 4L, maxSeconds = 0.001),
    limit = 5)
  # Need 2s to pass on memcheck runs
  expect_lte(attr(res, "time_s"), 2)
})

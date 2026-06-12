# Tests for MaxMinSeed() and MinDist().

MakeData <- function(seed = 7, N = 50, dim = 3) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

test_that("MaxMinSeed anchors match their definitions (matrix)", {
  dat <- MakeData()
  d <- dat$d
  expect_identical(MaxMinSeed(d, strategy = "medoid"),
                   as.integer(which.min(rowSums(d))))
  expect_identical(MaxMinSeed(d, strategy = "rowsum"),
                   as.integer(which.max(rowSums(d))))
  expect_identical(MaxMinSeed(d, strategy = "rownorm"),
                   as.integer(which.max(rowSums(d ^ 2))))
  med <- which.min(rowSums(d))
  dd <- d[, med]; dd[med] <- -Inf
  expect_identical(MaxMinSeed(d, strategy = "anti_medoid"),
                   as.integer(which.max(dd)))
  dOff <- d; diag(dOff) <- -Inf
  expect_identical(MaxMinSeed(d, strategy = "diameter"),
                   as.integer(arrayInd(which.max(dOff), dim(dOff))[1L, 1L]))
})

test_that("MaxMinSeed coordinate path matches the matrix path", {
  dat <- MakeData()
  for (k in c("peripheral", "random_furthest", "diameter", "anti_medoid",
              "medoid", "rowsum", "rownorm")) {
    # Seed identically so the random_furthest pivot matches across paths.
    set.seed(1); pt  <- MaxMinSeed(points = dat$pts, strategy = k)
    set.seed(1); mat <- MaxMinSeed(dat$d, strategy = k)
    expect_identical(pt, mat, info = k)
  }
})

test_that("centroid seed is the farthest point from the coordinate mean", {
  dat <- MakeData()
  mu  <- colMeans(dat$pts)
  d2  <- rowSums((dat$pts - rep(mu, each = nrow(dat$pts))) ^ 2)
  expect_identical(MaxMinSeed(points = dat$pts, strategy = "centroid"),
                   as.integer(which.max(d2)))
  # It has no distance-matrix form.
  expect_error(MaxMinSeed(dat$d, strategy = "centroid"), "coordinates")
  expect_error(FarFirst(5L, dat$d, strategy = "centroid"), "coordinates")
})

test_that("random_furthest is reproducible under set.seed", {
  dat <- MakeData()
  set.seed(1); a <- MaxMinSeed(points = dat$pts, strategy = "random_furthest")
  set.seed(1); b <- MaxMinSeed(points = dat$pts, strategy = "random_furthest")
  expect_identical(a, b)
})

test_that("MaxMinSeed validates method", {
  dat <- MakeData(N = 8)
  expect_error(MaxMinSeed(dat$d, strategy = "ensemble"))
  expect_error(MaxMinSeed(dat$d, strategy = "first"))
})

# ---- MinDist ------------------------------------------------------------

test_that("MinDist matrix and coordinate paths agree", {
  dat <- MakeData()
  idx <- c(1L, 5L, 9L, 20L, 33L)
  expect_equal(MinDist(dat$d, idx),
               MinDist(idx = idx, points = dat$pts))
  # equals the brute-force minimum pairwise distance.
  sub <- dat$d[idx, idx]; diag(sub) <- Inf
  expect_equal(MinDist(dat$d, idx), min(sub))
})

test_that("MinDist returns NA for fewer than two points", {
  dat <- MakeData(N = 6)
  expect_true(is.na(MinDist(dat$d, 3L)))
  expect_true(is.na(MinDist(idx = integer(0), points = dat$pts)))
})

test_that("MinDist rejects NA and duplicate idx (F-604/F-605)", {
  dat <- MakeData(N = 8)
  # NA was silently NA on the matrix path but an error on the points path.
  expect_error(MinDist(dat$d, c(1L, NA, 3L)), "NA")
  expect_error(MinDist(idx = c(1L, NA), points = dat$pts), "NA")
  # Duplicate indices made any selection score 0 (self-distance survives).
  expect_error(MinDist(dat$d, c(1L, 2L, 2L)), "duplicate")
  expect_error(MinDist(idx = c(1L, 1L), points = dat$pts), "duplicate")
})

# ---- "first" seed (.MaxMinSeed line 16; .MaxMinSeedPoints line 55) ---------

test_that("Gonzalez strategy='first' uses index 1 as anchor (both paths)", {
  dat <- MakeData()
  rMat <- FarFirst(4L, dat$d, strategy = "first")
  expect_identical(rMat[1L], 1L)
  rPts <- FarFirst(k = 4L, points = dat$pts, strategy = "first")
  expect_identical(rPts[1L], 1L)
})

# ---- Degenerate diameter (.MaxMinSeed line 31; .MaxMinSeedPoints line 68;
#       .GonzEnsemble AnchorSeed line 206; .GonzEnsembleFromPoints line 328) -

test_that("diameter anchor returns 1 on zero-distance data", {
  # All points at the origin -> all pairwise distances are 0 -> dMax <= 0.
  ptsDegen <- matrix(0, nrow = 5L, ncol = 2L)
  dDegen   <- as.matrix(dist(ptsDegen))
  expect_identical(MaxMinSeed(dDegen, strategy = "diameter"), 1L)
  expect_identical(MaxMinSeed(points = ptsDegen, strategy = "diameter"), 1L)
  # Ensemble paths: AnchorSeed("diameter") hits the degenerate branch.
  ensM  <- FarFirst(2L, dDegen, strategy = c("diameter", "rowsum"))
  ensPt <- FarFirst(k = 2L, points = ptsDegen, strategy = c("diameter", "rowsum"))
  expect_length(ensM,  2L)
  expect_length(ensPt, 2L)
})

# ---- Non-first anchor wins (.GonzEnsemble lines 241-242;
#       .GonzEnsembleFromPoints lines 362-363) --------------------------------

test_that(".GonzEnsemble non-first anchor wins when it is the better one", {
  # Scan seeds until we find an instance where two strategies yield different
  # T_k; construct a 2-anchor ensemble with the LOSER first and WINNER second so
  # the update branch (bestI <- i, bestTk <- tk) fires.
  found <- FALSE
  for (s in seq_len(50L)) {
    set.seed(s)
    pts <- matrix(rnorm(80L * 3L), ncol = 3L)
    d   <- as.matrix(dist(pts))
    n   <- 6L
    tkD  <- MinDist(d, FarFirst(n, d, strategy = "diameter"))
    tkAm <- MinDist(d, FarFirst(n, d, strategy = "anti_medoid"))
    if (isTRUE(all.equal(tkD, tkAm))) next
    loser  <- if (tkAm > tkD) "diameter"    else "anti_medoid"
    winner <- if (tkAm > tkD) "anti_medoid" else "diameter"
    bestTk <- max(tkD, tkAm)
    ensM  <- FarFirst(n, d, strategy = c(loser, winner))
    ensPt <- FarFirst(k = n, points = pts, strategy = c(loser, winner))
    expect_gte(MinDist(d, ensM),  bestTk - 1e-9)
    expect_gte(MinDist(d, ensPt), bestTk - 1e-9)
    found <- TRUE
    break
  }
  expect_true(found)
})

# ---- .MaxMinSeed / .MaxMinSeedPoints unknown-method fallthrough -------------

test_that(".MaxMinSeed and .MaxMinSeedPoints stop on an unknown method", {
  dat <- MakeData()
  expect_error(MaxMin:::.MaxMinSeed(dat$d, "nope"), "Unknown")
  expect_error(MaxMin:::.MaxMinSeedPoints(dat$pts, "nope"), "Unknown")
})

# ---- .GonzEnsemble internal guards (bypassed by FarFirst() validation) -----

test_that(".GonzEnsemble validates k and anchors and handles trivial k", {
  dat <- MakeData()
  d   <- dat$d
  n   <- nrow(d)
  expect_error(MaxMin:::.GonzEnsemble(d, -1L, "peripheral"), "non-negative")
  expect_error(MaxMin:::.GonzEnsemble(d, 3L, character(0)), "at least one")
  # k >= nPts: all points returned with score = NA_real_ (FF-005)
  trivial <- MaxMin:::.GonzEnsemble(d, n, "peripheral")
  expect_identical(as.integer(trivial), seq_len(n))
  expect_true(is.na(attr(trivial, "score")))
  expect_identical(MaxMin:::.GonzEnsemble(d, 0L, "peripheral"), integer(0))
})

# ---- .GonzEnsembleFromPoints internal guards --------------------------------

test_that(".GonzEnsembleFromPoints validates k and anchors and handles trivial k", {
  dat  <- MakeData()
  pts  <- dat$pts
  nPts <- nrow(pts)
  expect_error(MaxMin:::.GonzEnsembleFromPoints(pts, -1L, "peripheral"), "non-negative")
  expect_error(MaxMin:::.GonzEnsembleFromPoints(pts, 3L, character(0)), "at least one")
  # k >= nPts: all points returned with score = NA_real_ (FF-005)
  trivial <- MaxMin:::.GonzEnsembleFromPoints(pts, nPts, "peripheral")
  expect_identical(as.integer(trivial), seq_len(nPts))
  expect_true(is.na(attr(trivial, "score")))
  expect_identical(MaxMin:::.GonzEnsembleFromPoints(pts, 0L, "peripheral"), integer(0))
})

# ---- .ExpandAnchors empty-specs guard (no anchors produce no specs) ----------

test_that(".ExpandAnchors stops when no specs are produced", {
  # Calling the internal directly with empty rfSeeds and no deterministic anchors
  # triggers the guard; nseeds = 0 can't be reached via FarFirst (validated >= 1).
  expect_error(
    MaxMin:::.ExpandAnchors("random_furthest", integer(0), function(n) 1L),
    "no seed strateg"
  )
})

# ---- k = 1 ensemble on the points path (.GonzEnsembleFromPoints line 360,
#       368) -----------------------------------------------------------------

test_that(".GonzEnsembleFromPoints k=1 propagates NA t_k through all strategies", {
  dat <- MakeData()
  # Each strategy returns one point -> t_k = NA -> is.na(tk) next fires for
  # every non-first strategy (line 360); is.na(bestTk) selects bestI (line 368).
  res <- FarFirst(k = 1L, points = dat$pts)
  expect_length(res, 1L)
  strat <- attr(res, "strategy_results")
  expect_true(all(is.na(vapply(strat, `[[`, numeric(1L), "t_k"))))
})

# ---- ensemble AnchorSeed branches not covered by single-strategy FarFirst ----
# Single-string strategy calls go through .MaxMinSeedPoints directly; the
# AnchorSeed() closure inside .GonzEnsemble/.GonzEnsembleFromPoints is only
# reached when those names appear in a multi-anchor ensemble call.

test_that(".GonzEnsemble medoid branch covered via multi-anchor ensemble", {
  dat <- MakeData()
  # medoid in a 2-anchor ensemble -> AnchorSeed("medoid") = GetMedoid() fires
  res <- FarFirst(5L, dat$d, strategy = c("medoid", "peripheral"))
  expect_length(res, 5L)
  expect_true(all(attr(res, "winning_strategy") %in% c("medoid", "peripheral")))
})

test_that(".GonzEnsembleFromPoints medoid/centroid/rownorm branches covered via ensemble", {
  dat <- MakeData()
  # Three anchors only reachable through the points-path ensemble AnchorSeed
  res <- FarFirst(k = 5L, points = dat$pts,
                  strategy = c("medoid", "centroid", "rownorm"))
  expect_length(res, 5L)
  expect_true(all(attr(res, "winning_strategy") %in%
                  c("medoid", "centroid", "rownorm")))
})

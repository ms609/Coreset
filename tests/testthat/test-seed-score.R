# Tests for MaxMinSeed() and MinDist().

MakeData <- function(seed = 7, N = 50, dim = 3) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

test_that("MaxMinSeed anchors match their definitions (matrix)", {
  dat <- MakeData()
  d <- dat$d
  expect_identical(MaxMinSeed(d, method = "medoid"),
                   as.integer(which.min(rowSums(d))))
  expect_identical(MaxMinSeed(d, method = "rowsum"),
                   as.integer(which.max(rowSums(d))))
  expect_identical(MaxMinSeed(d, method = "rownorm"),
                   as.integer(which.max(rowSums(d ^ 2))))
  med <- which.min(rowSums(d))
  dd <- d[, med]; dd[med] <- -Inf
  expect_identical(MaxMinSeed(d, method = "anti_medoid"),
                   as.integer(which.max(dd)))
  dOff <- d; diag(dOff) <- -Inf
  expect_identical(MaxMinSeed(d, method = "diameter"),
                   as.integer(arrayInd(which.max(dOff), dim(dOff))[1L, 1L]))
})

test_that("MaxMinSeed coordinate path matches the matrix path", {
  dat <- MakeData()
  for (m in c("peripheral", "random_furthest", "diameter", "anti_medoid",
              "medoid", "rowsum", "rownorm")) {
    # Seed identically so the random_furthest pivot matches across paths.
    set.seed(1); pt  <- MaxMinSeed(points = dat$pts, method = m)
    set.seed(1); mat <- MaxMinSeed(dat$d, method = m)
    expect_identical(pt, mat, info = m)
  }
})

test_that("centroid seed is the farthest point from the coordinate mean", {
  dat <- MakeData()
  mu  <- colMeans(dat$pts)
  d2  <- rowSums((dat$pts - rep(mu, each = nrow(dat$pts))) ^ 2)
  expect_identical(MaxMinSeed(points = dat$pts, method = "centroid"),
                   as.integer(which.max(d2)))
  # It has no distance-matrix form.
  expect_error(MaxMinSeed(dat$d, method = "centroid"), "coordinates")
  expect_error(Gonzalez(dat$d, 5L, seed = "centroid"), "coordinates")
})

test_that("random_furthest is reproducible under set.seed", {
  dat <- MakeData()
  set.seed(1); a <- MaxMinSeed(points = dat$pts, method = "random_furthest")
  set.seed(1); b <- MaxMinSeed(points = dat$pts, method = "random_furthest")
  expect_identical(a, b)
})

test_that("MaxMinSeed validates method", {
  dat <- MakeData(N = 8)
  expect_error(MaxMinSeed(dat$d, method = "ensemble"))
  expect_error(MaxMinSeed(dat$d, method = "first"))
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

# ---- "first" seed (.MaxMinSeed line 16; .MaxMinSeedPoints line 55) ---------

test_that("Gonzalez seed='first' uses index 1 as anchor (both paths)", {
  dat <- MakeData()
  rMat <- Gonzalez(dat$d, 4L, seed = "first")
  expect_identical(rMat[1L], 1L)
  rPts <- Gonzalez(n = 4L, points = dat$pts, seed = "first")
  expect_identical(rPts[1L], 1L)
})

# ---- Degenerate diameter (.MaxMinSeed line 31; .MaxMinSeedPoints line 68;
#       .GonzEnsemble AnchorSeed line 206; .GonzEnsembleFromPoints line 328) -

test_that("diameter anchor returns 1 on zero-distance data", {
  # All points at the origin -> all pairwise distances are 0 -> dMax <= 0.
  ptsDegen <- matrix(0, nrow = 5L, ncol = 2L)
  dDegen   <- as.matrix(dist(ptsDegen))
  expect_identical(MaxMinSeed(dDegen, method = "diameter"), 1L)
  expect_identical(MaxMinSeed(points = ptsDegen, method = "diameter"), 1L)
  # Ensemble paths: AnchorSeed("diameter") hits the degenerate branch.
  ensM  <- Gonzalez(dDegen, 2L, seed = c("diameter", "rowsum"))
  ensPt <- Gonzalez(n = 2L, points = ptsDegen, seed = c("diameter", "rowsum"))
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
    tkD  <- MinDist(d, Gonzalez(d, n, seed = "diameter"))
    tkAm <- MinDist(d, Gonzalez(d, n, seed = "anti_medoid"))
    if (isTRUE(all.equal(tkD, tkAm))) next
    loser  <- if (tkAm > tkD) "diameter"    else "anti_medoid"
    winner <- if (tkAm > tkD) "anti_medoid" else "diameter"
    bestTk <- max(tkD, tkAm)
    ensM  <- Gonzalez(d, n, seed = c(loser, winner))
    ensPt <- Gonzalez(n = n, points = pts, seed = c(loser, winner))
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

# ---- .GonzEnsemble internal guards (bypassed by Gonzalez() validation) -----

test_that(".GonzEnsemble validates n and anchors and handles trivial n", {
  dat <- MakeData()
  d   <- dat$d
  n   <- nrow(d)
  expect_error(MaxMin:::.GonzEnsemble(d, -1L, "peripheral"), "non-negative")
  expect_error(MaxMin:::.GonzEnsemble(d, 3L, character(0)), "at least one")
  expect_identical(MaxMin:::.GonzEnsemble(d, n, "peripheral"), seq_len(n))
  expect_identical(MaxMin:::.GonzEnsemble(d, 0L, "peripheral"), integer(0))
})

# ---- .GonzEnsembleFromPoints internal guards --------------------------------

test_that(".GonzEnsembleFromPoints validates n and anchors and handles trivial n", {
  dat  <- MakeData()
  pts  <- dat$pts
  nPts <- nrow(pts)
  expect_error(MaxMin:::.GonzEnsembleFromPoints(pts, -1L, "peripheral"), "non-negative")
  expect_error(MaxMin:::.GonzEnsembleFromPoints(pts, 3L, character(0)), "at least one")
  expect_identical(MaxMin:::.GonzEnsembleFromPoints(pts, nPts, "peripheral"),
                   seq_len(nPts))
  expect_identical(MaxMin:::.GonzEnsembleFromPoints(pts, 0L, "peripheral"), integer(0))
})

# ---- .ExpandAnchors empty-specs guard (only random_furthest + no pivots) ----

test_that(".GonzEnsemble stops when random_furthest only and pivots is empty", {
  dat <- MakeData()
  expect_error(
    MaxMin:::.GonzEnsemble(dat$d, 4L, "random_furthest", integer(0)),
    "no seed strateg"
  )
  expect_error(
    MaxMin:::.GonzEnsembleFromPoints(dat$pts, 4L, "random_furthest", integer(0)),
    "no seed strateg"
  )
})

# ---- n = 1 ensemble on the points path (.GonzEnsembleFromPoints line 360,
#       368) -----------------------------------------------------------------

test_that(".GonzEnsembleFromPoints n=1 propagates NA t_k through all strategies", {
  dat <- MakeData()
  # Each strategy returns one point -> t_k = NA -> is.na(tk) next fires for
  # every non-first strategy (line 360); is.na(bestTk) selects bestI (line 368).
  res <- Gonzalez(n = 1L, points = dat$pts)
  expect_length(res, 1L)
  strat <- attr(res, "strategy_results")
  expect_true(all(is.na(vapply(strat, `[[`, numeric(1L), "t_k"))))
})

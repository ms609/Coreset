# Tests for MaxMinSeed() and TkScore().

make_data <- function(seed = 7, N = 50, dim = 3) {
  set.seed(seed)
  pts <- matrix(rnorm(N * dim), ncol = dim)
  list(pts = pts, d = as.matrix(dist(pts)))
}

test_that("MaxMinSeed anchors match their definitions (matrix)", {
  dat <- make_data()
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
  d_off <- d; diag(d_off) <- -Inf
  expect_identical(MaxMinSeed(d, method = "diameter"),
                   as.integer(arrayInd(which.max(d_off), dim(d_off))[1L, 1L]))
})

test_that("MaxMinSeed coordinate path matches the matrix path", {
  dat <- make_data()
  for (m in c("peripheral", "random_furthest", "diameter", "anti_medoid",
              "medoid", "rowsum", "rownorm")) {
    expect_identical(MaxMinSeed(points = dat$pts, method = m),
                     MaxMinSeed(dat$d, method = m), info = m)
  }
})

test_that("centroid seed is the farthest point from the coordinate mean", {
  dat <- make_data()
  mu  <- colMeans(dat$pts)
  d2  <- rowSums((dat$pts - rep(mu, each = nrow(dat$pts))) ^ 2)
  expect_identical(MaxMinSeed(points = dat$pts, method = "centroid"),
                   as.integer(which.max(d2)))
  # It has no distance-matrix form.
  expect_error(MaxMinSeed(dat$d, method = "centroid"), "coordinates")
  expect_error(Gonzalez(dat$d, 5L, seed = "centroid"), "coordinates")
})

test_that("random_furthest is reproducible and RNG-isolated", {
  dat <- make_data()
  set.seed(1);   a <- MaxMinSeed(points = dat$pts, method = "random_furthest")
  set.seed(404); b <- MaxMinSeed(points = dat$pts, method = "random_furthest")
  expect_identical(a, b)                      # ambient-independent
  set.seed(3); u1 <- runif(1)
  set.seed(3); invisible(MaxMinSeed(dat$d, method = "random_furthest"))
  u2 <- runif(1)
  expect_identical(u1, u2)                     # caller's stream untouched
})

test_that("MaxMinSeed validates method", {
  dat <- make_data(N = 8)
  expect_error(MaxMinSeed(dat$d, method = "ensemble"))
  expect_error(MaxMinSeed(dat$d, method = "first"))
})

# ---- TkScore ------------------------------------------------------------

test_that("TkScore matrix and coordinate paths agree", {
  dat <- make_data()
  idx <- c(1L, 5L, 9L, 20L, 33L)
  expect_equal(TkScore(dat$d, idx),
               TkScore(idx = idx, points = dat$pts))
  # equals the brute-force minimum pairwise distance.
  sub <- dat$d[idx, idx]; diag(sub) <- Inf
  expect_equal(TkScore(dat$d, idx), min(sub))
})

test_that("TkScore returns NA for fewer than two points", {
  dat <- make_data(N = 6)
  expect_true(is.na(TkScore(dat$d, 3L)))
  expect_true(is.na(TkScore(idx = integer(0), points = dat$pts)))
})

# ---- "first" seed (.MaxMinSeed line 16; .MaxMinSeedPoints line 55) ---------

test_that("Gonzalez seed='first' uses index 1 as anchor (both paths)", {
  dat <- make_data()
  r_mat <- Gonzalez(dat$d, 4L, seed = "first")
  expect_identical(r_mat[1L], 1L)
  r_pts <- Gonzalez(n = 4L, points = dat$pts, seed = "first")
  expect_identical(r_pts[1L], 1L)
})

# ---- Degenerate diameter (.MaxMinSeed line 31; .MaxMinSeedPoints line 68;
#       .GonzEnsemble anchor_seed line 206; .GonzEnsembleFromPoints line 328) -

test_that("diameter anchor returns 1 on zero-distance data", {
  # All points at the origin -> all pairwise distances are 0 -> d_max <= 0.
  pts_degen <- matrix(0, nrow = 5L, ncol = 2L)
  d_degen   <- as.matrix(dist(pts_degen))
  expect_identical(MaxMinSeed(d_degen, method = "diameter"), 1L)
  expect_identical(MaxMinSeed(points = pts_degen, method = "diameter"), 1L)
  # Ensemble paths: anchor_seed("diameter") hits the degenerate branch.
  ens_m  <- Gonzalez(d_degen, 2L, seed = c("diameter", "rowsum"))
  ens_pt <- Gonzalez(n = 2L, points = pts_degen, seed = c("diameter", "rowsum"))
  expect_length(ens_m,  2L)
  expect_length(ens_pt, 2L)
})

# ---- Non-first anchor wins (.GonzEnsemble lines 241-242;
#       .GonzEnsembleFromPoints lines 362-363) --------------------------------

test_that(".GonzEnsemble non-first anchor wins when it is the better one", {
  # Scan seeds until we find an instance where two strategies yield different
  # T_k; construct a 2-anchor ensemble with the LOSER first and WINNER second so
  # the update branch (best_i <- i, best_tk <- tk) fires.
  found <- FALSE
  for (s in seq_len(50L)) {
    set.seed(s)
    pts <- matrix(rnorm(80L * 3L), ncol = 3L)
    d   <- as.matrix(dist(pts))
    n   <- 6L
    tk_d  <- TkScore(d, Gonzalez(d, n, seed = "diameter"))
    tk_am <- TkScore(d, Gonzalez(d, n, seed = "anti_medoid"))
    if (isTRUE(all.equal(tk_d, tk_am))) next
    loser  <- if (tk_am > tk_d) "diameter"    else "anti_medoid"
    winner <- if (tk_am > tk_d) "anti_medoid" else "diameter"
    best_tk <- max(tk_d, tk_am)
    ens_m  <- Gonzalez(d, n, seed = c(loser, winner))
    ens_pt <- Gonzalez(n = n, points = pts, seed = c(loser, winner))
    expect_gte(TkScore(d, ens_m),  best_tk - 1e-9)
    expect_gte(TkScore(d, ens_pt), best_tk - 1e-9)
    found <- TRUE
    break
  }
  expect_true(found)
})

# ---- n = 1 ensemble on the points path (.GonzEnsembleFromPoints line 360,
#       368) -----------------------------------------------------------------

test_that(".GonzEnsembleFromPoints n=1 propagates NA t_k through all strategies", {
  dat <- make_data()
  # Each strategy returns one point -> t_k = NA -> is.na(tk) next fires for
  # every non-first strategy (line 360); is.na(best_tk) selects best_i (line 368).
  res <- Gonzalez(n = 1L, points = dat$pts)
  expect_length(res, 1L)
  strat <- attr(res, "strategy_results")
  expect_true(all(is.na(vapply(strat, `[[`, numeric(1L), "t_k"))))
})

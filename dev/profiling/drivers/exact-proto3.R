# Prototype v3: heuristic-seeded search. DropAdd gives a proven-feasible lower
# bound (every realized distance <= it is feasible); gallop UP from there to the
# first infeasible threshold, then binary-search the (small) bracket. Collapses
# ~16 IP solves to O(log gap) -- ~1-3 when the heuristic is optimal.
# Layered on sparse-A + cheap-edges (v2). Plus an objective_bound cutoff probe.
set.seed(5813)
suppressMessages({ library(MaxMin); library(Matrix); library(highs) })
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")

# feasibility oracle: does G(lambda) admit an IS of size >= m? returns witness.
# cutoff: stop highs as soon as it finds obj >= m (we don't need the true max).
make_oracle <- function(d, n, ui, uj, ud, m, cutoff = TRUE) {
  function(lambda, remaining) {
    if (remaining <= 0) return(list(v = "inconclusive", w = integer(0)))
    e <- ud < lambda; nE <- sum(e)
    if (nE == 0L) return(list(v = "feasible", w = seq_len(n)))
    A <- sparseMatrix(i = rep(seq_len(nE), 2L), j = c(ui[e], uj[e]), x = 1, dims = c(nE, n))
    ctrl <- list(threads = 1L, time_limit = remaining)
    if (cutoff) ctrl$objective_bound <- m            # stop once max IS reaches m
    res <- highs_solve(L = rep.int(1,n), lower = rep.int(0,n), upper = rep.int(1,n),
      A = A, lhs = rep.int(-Inf,nE), rhs = rep.int(1,nE), types = rep.int("I",n),
      maximum = TRUE, control = ctrl)
    sel <- which(res$primal_solution > 0.5)
    ok <- if (length(sel) < 2L) TRUE else { sub <- d[sel,sel]; !any(sub[upper.tri(sub)] < lambda) }
    if (ok && length(sel) >= m) return(list(v = "feasible", w = sel))
    if (identical(res$status_message, "Optimal") && ok && length(sel) < m)
      return(list(v = "infeasible", w = integer(0)))
    list(v = "inconclusive", w = integer(0))
  }
}

Exact_v3 <- function(d, m, timeBudgetS = 600, cutoff = TRUE, seedMethod = "dropadd") {
  t0 <- proc.time()[[3L]]; Elapsed <- function() proc.time()[[3L]] - t0
  d <- as.matrix(d); n <- nrow(d); m <- as.integer(m)
  ut <- which(upper.tri(d)); rc <- arrayInd(ut, dim(d))
  ui <- rc[,1L]; uj <- rc[,2L]; ud <- d[ut]
  cand <- sort(unique(ud)); nCand <- length(cand)
  feas <- make_oracle(d, n, ui, uj, ud, m, cutoff)

  # ---- heuristic seed: a provably-achievable m-subset -> lower bound ----
  scv <- function(idx) { s <- d[idx, idx]; diag(s) <- Inf; min(s) }
  cands_seed <- list()
  if (seedMethod %in% c("grasp", "best"))
    cands_seed$grasp <- sort(as.integer(MaxMin::Grasp(d, m, plateau = 50L, seed = 1L)))
  if (seedMethod %in% c("dropadd", "best"))
    cands_seed$dropadd <- sort(as.integer(MaxMin::DropAdd(d = d, m = m, plateau = 512L)))
  vals <- vapply(cands_seed, scv, numeric(1))
  Sh <- cands_seed[[which.max(vals)]]; LB <- max(vals)
  i0 <- findInterval(LB, cand)                       # cand[i0] == LB (realized)
  best <- i0; bestW <- Sh; nProbe <- 0L; inconcl <- FALSE

  # ---- gallop up from i0 to first infeasible index ----
  step <- 1L; loF <- i0; hiX <- NA_integer_; probe <- i0 + 1L
  while (probe <= nCand) {
    rem <- timeBudgetS - Elapsed(); if (rem <= 0) { inconcl <- TRUE; break }
    v <- feas(cand[probe], rem); nProbe <- nProbe + 1L
    if (v$v == "feasible") { loF <- probe; best <- probe; bestW <- v$w
      step <- step * 2L; probe <- probe + step
    } else if (v$v == "infeasible") { hiX <- probe; break } else { inconcl <- TRUE; break }
  }
  if (is.na(hiX) && !inconcl) hiX <- nCand + 1L       # all feasible to the top
  # ---- binary search the bracket (loF, hiX) ----
  if (!inconcl) {
    lo <- loF + 1L; hi <- min(hiX - 1L, nCand)
    while (lo <= hi) {
      rem <- timeBudgetS - Elapsed(); if (rem <= 0) { inconcl <- TRUE; break }
      mid <- (lo + hi) %/% 2L
      v <- feas(cand[mid], rem); nProbe <- nProbe + 1L
      if (v$v == "feasible") { best <- mid; bestW <- v$w; lo <- mid + 1L }
      else if (v$v == "infeasible") { hi <- mid - 1L } else { inconcl <- TRUE; break }
    }
  }
  idx <- sort(bestW[seq_len(m)]); s2 <- d[idx, idx]; diag(s2) <- Inf
  list(indices = idx, objective = min(s2), proven = !inconcl,
       time_s = Elapsed(), n = n, m = as.integer(m), nProbe = nProbe, LB = LB)
}

bench <- function(case, k = 10L) {
  pts <- as.matrix(cases[[case]][["points"]]); storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts)); n <- nrow(d)
  t <- proc.time()[[3L]]; r0 <- MaxMin::ExactMaxMin(d, k, timeBudgetS = 600, progress = FALSE); t0 <- proc.time()[[3L]] - t
  t <- proc.time()[[3L]]; rg <- Exact_v3(d, k, cutoff = FALSE, seedMethod = "grasp"); tg <- proc.time()[[3L]] - t
  t <- proc.time()[[3L]]; rb <- Exact_v3(d, k, cutoff = FALSE, seedMethod = "best");  tb <- proc.time()[[3L]] - t
  ok <- isTRUE(all.equal(r0$objective, rg$objective)) && isTRUE(all.equal(r0$objective, rb$objective)) &&
        r0$proven && rg$proven && rb$proven
  cat(sprintf("%-16s n=%4d | installed %.2fs (16 probes) | v3-grasp %.2fs (%d probes, %.1fx) | v3-best %.2fs (%d probes, %.1fx) | LB=%.4f opt=%.4f | match=%s\n",
    case, n, t0, tg, rg$nProbe, t0/tg, tb, rb$nProbe, t0/tb, rg$LB, r0$objective, ok))
}
for (cs in c("tc22_penguins", "tc11_ionosphere")) bench(cs)

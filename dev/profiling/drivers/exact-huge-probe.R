# Is exact even feasible for the huge datasets? Probe spam (n=4601), the
# smallest of the four. Strong GRASP warm start (does it land the optimum?),
# then time the SINGLE infeasibility probe just above the warm-start LB -- the
# one solve GRASP cannot remove. If that one probe closes fast, exact is viable;
# if it can't close in a few minutes, no number of GRASP restarts helps.
suppressMessages({ library(MaxMin); library(Matrix); library(highs) })
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
case <- Sys.getenv("FP_HUGE_CASE", "tc21_spam")
k <- 10L; probe_cap <- as.numeric(Sys.getenv("FP_PROBE_CAP", "200"))

pts <- as.matrix(cases[[case]][["points"]]); storage.mode(pts) <- "double"; N <- nrow(pts)
t <- proc.time()[[3L]]; d <- as.matrix(stats::dist(pts)); cat(sprintf("%s n=%d: dist matrix %.1fs (%.2f GB)\n", case, N, proc.time()[[3L]]-t, N*N*8/1e9))

ut <- which(upper.tri(d)); rc <- arrayInd(ut, dim(d)); ui <- rc[,1L]; uj <- rc[,2L]; ud <- d[ut]
cand <- sort(unique(ud)); nCand <- length(cand)
cat(sprintf("  %d distinct distances\n", nCand))

scv <- function(idx){ s <- d[idx,idx]; diag(s) <- Inf; min(s) }
set.seed(1)
t <- proc.time()[[3L]]
gr <- vapply(1:20, function(s) scv(Grasp(d, k, plateau = 50L)), numeric(1))
da <- scv(DropAdd(d = d, m = k, plateau = 512L))
LB <- max(gr, da); cat(sprintf("  warm start (20 Grasp + DropAdd): LB=%.6f  %.1fs  (best Grasp=%.6f, DropAdd=%.6f)\n", LB, proc.time()[[3L]]-t, max(gr), da))

i0 <- findInterval(LB, cand)
cat(sprintf("  LB index i0=%d / %d  (=%.4f%% up the distance range)\n", i0, nCand, 100*i0/nCand))

probe <- function(idx, cap) {
  lambda <- cand[idx]; e <- ud < lambda; nE <- sum(e)
  A <- sparseMatrix(i = rep(seq_len(nE),2L), j = c(ui[e], uj[e]), x = 1, dims = c(nE, N))
  t <- proc.time()[[3L]]
  res <- highs_solve(L = rep.int(1,N), lower = rep.int(0,N), upper = rep.int(1,N),
    A = A, lhs = rep.int(-Inf,nE), rhs = rep.int(1,nE), types = rep.int("I",N),
    maximum = TRUE, control = list(threads = 1L, time_limit = cap))
  el <- proc.time()[[3L]] - t
  alpha <- sum(res$primal_solution > 0.5)
  list(nE = nE, alpha = alpha, status = res$status_message, el = el)
}

cat(sprintf("\n  probing cand[i0+1] (just above LB) with %.0fs cap ...\n", probe_cap))
r <- probe(i0 + 1L, probe_cap)
cat(sprintf("  nEdge=%d (%.1f%% of pairs)  alpha=%d  status=%s  %.1fs\n",
            r$nE, 100*r$nE/nCand, r$alpha, r$status, r$el))
verdict <- if (r$status == "Optimal" && r$alpha < k) "INFEASIBLE (LB is the optimum -> 1 probe proves it)" else
           if (r$alpha >= k) "FEASIBLE (LB suboptimal -> gallop would continue, more probes)" else
           "INCONCLUSIVE (one probe exceeded the cap -> exact impractical at this n)"
cat(sprintf("  => %s\n", verdict))

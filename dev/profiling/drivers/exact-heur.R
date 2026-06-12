# How tight are the heuristics vs the exact optimum at k=10? Decides whether a
# heuristic-seeded exact search can prove optimality in ~1 probe.
suppressMessages({ library(MaxMin) })
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
sc <- function(d, idx) { s <- d[idx, idx]; diag(s) <- Inf; min(s) }
for (cs in c("tc22_penguins", "tc11_ionosphere")) {
  pts <- as.matrix(cases[[cs]][["points"]]); storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts)); k <- 10L
  opt <- MaxMin::ExactMaxMin(d, k, maxSeconds = 600)$objective
  da512  <- sc(d, MaxMin::DropAdd(d = d, m = k, plateau = 512L))
  da5000 <- sc(d, MaxMin::DropAdd(d = d, m = k, plateau = 5000L))
  gp <- sc(d, MaxMin::Grasp(d, k, plateau = 50L, seed = 1L))
  # best of a few seeds / methods
  best <- max(da512, da5000, gp)
  cat(sprintf("%-16s opt=%.5f | dropadd512=%.5f dropadd5000=%.5f grasp=%.5f | best=%.5f gap=%.3f%% %s\n",
    cs, opt, da512, da5000, gp, best, 100*(opt-best)/opt,
    if (abs(best-opt) < 1e-9) "<-- HEURISTIC OPTIMAL" else ""))
}

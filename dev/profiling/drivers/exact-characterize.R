# Does the heuristic-seed speedup GENERALISE? Characterise, across every case
# with exact ground truth (the 12 small manuscript cases, n<=240) + the two
# smallest targets, at the k values the job will use (2,4,6,10):
#   (1) how often the best heuristic == the exact optimum  (=> 1-probe win),
#   (2) that the new seeded search returns the EXACT optimum every time, even
#       when the heuristic misses (correctness independent of seed quality),
#   (3) the probe count distribution (worst case vs binary search's ~16).
suppressMessages({ library(MaxMin); library(Matrix); library(highs) })
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
source("dev/profiling/drivers/exact-proto3.R", local = TRUE)   # defines Exact_v3, make_oracle

scv <- function(d, idx) { s <- d[idx, idx]; diag(s) <- Inf; min(s) }
# strongest cheap heuristic seed: Grasp (several seeds) + DropAdd, take the best
best_heur <- function(d, m) {
  vs <- c(
    vapply(1:3, function(s) scv(d, MaxMin::Grasp(d, m, plateau = 50L, seed = s)), numeric(1)),
    scv(d, MaxMin::DropAdd(d = d, m = m, plateau = 512L))
  )
  max(vs)
}

small <- c("tc20_zoo","tc7_ring","tc9_iris","tc8_highdim_gaussians","tc3_multiscale",
           "tc1_uniform","tc2_two_unequal","tc6_density_gradient","tc16_sonar",
           "tc5_clusters_outliers","tc10_glass","tc4_hierarchical",
           "tc22_penguins","tc11_ionosphere")
ks <- c(2L, 4L, 6L, 10L)

rows <- list()
for (cs in small) {
  pts <- as.matrix(cases[[cs]][["points"]]); storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts)); n <- nrow(d)
  for (k in ks) {
    opt <- MaxMin::ExactMaxMin(d, k, timeBudgetS = 600, progress = FALSE)
    bh  <- best_heur(d, k)
    rv  <- Exact_v3(d, k, cutoff = FALSE, seedMethod = "best")   # Grasp+DropAdd seed
    rows[[length(rows)+1L]] <- data.frame(
      case = cs, n = n, k = k,
      exact = round(opt$objective, 5),
      heur  = round(bh, 5),
      heur_opt = abs(bh - opt$objective) < 1e-9,
      new_correct = isTRUE(all.equal(rv$objective, opt$objective)) && rv$proven,
      probes = rv$nProbe, stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, rows)
cat("\n================ per-(case,k) ================\n")
print(R, row.names = FALSE)
cat("\n================ summary ================\n")
cat(sprintf("heuristic == exact optimum : %d / %d  (%.0f%%)\n",
    sum(R$heur_opt), nrow(R), 100*mean(R$heur_opt)))
cat(sprintf("new search returns EXACT   : %d / %d  (correctness)\n",
    sum(R$new_correct), nrow(R)))
cat(sprintf("probe count: min %d  median %d  max %d   (binary search baseline ~16-17)\n",
    min(R$probes), as.integer(median(R$probes)), max(R$probes)))
cat(sprintf("worst case probes when heuristic MISSES: %d\n",
    if (any(!R$heur_opt)) max(R$probes[!R$heur_opt]) else 0L))

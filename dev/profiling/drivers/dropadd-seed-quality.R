# Issue #1: does DropAdd()'s TS initializer choice matter?
#
# Current default: matrix kernel seeds at the max-row-sum point ("rowsum",
# argmax_i sum_j d_ij). Question: does any of the 7 single-anchor seed
# strategies used by the Gonzalez ensemble ("peripheral", "random_furthest",
# "diameter", "anti_medoid", "medoid", "rowsum", "rownorm") give DropAdd's
# *final* (post tabu-search) objective a meaningfully better starting point?
# Ground truth: ExactMaxMin() proven optima over the reference case grid.
#
# Usage: R --vanilla --no-echo < dev/profiling/drivers/dropadd-seed-quality.R
lib <- Sys.getenv("FP_LIB", "")
if (nzchar(lib)) library(Coreset, lib.loc = lib) else pkgload::load_all(".", quiet = TRUE)
options(Coreset.progress = FALSE)

load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")

cs <- c("tc1_uniform", "tc2_two_unequal", "tc3_multiscale",
        "tc4_hierarchical", "tc5_clusters_outliers", "tc6_density_gradient",
        "tc9_iris", "tc10_glass", "tc16_sonar", "tc22_penguins")
ks <- c(2L, 5L, 10L, 20L)

strategies <- c("peripheral", "random_furthest", "diameter",
                "anti_medoid", "medoid", "rowsum", "rownorm")

Scv <- function(d, idx) {
  s <- d[idx, idx]
  diag(s) <- Inf
  min(s)
}

rows <- list()
for (caseName in cs) {
  pts <- as.matrix(cases[[caseName]][["points"]])
  storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts))
  n <- nrow(d)
  for (k in ks) {
    if (k >= n) next
    set.seed(42L)
    t0 <- proc.time()[[3L]]
    opt <- ExactMaxMin(k, d, maxSeconds = 120)
    optWall <- proc.time()[[3L]] - t0
    optScore <- attr(opt, "score")
    proven <- isTRUE(attr(opt, "proven"))

    # Default DropAdd (unseeded -- current rowsum construction default).
    set.seed(1L)
    defIdx <- DropAdd(k, d)
    defScore <- Scv(d, defIdx)

    stratScores <- setNames(numeric(length(strategies)), strategies)
    stratTk     <- setNames(numeric(length(strategies)), strategies)
    for (strat in strategies) {
      seedIdx <- Coreset:::.PickPoint(d, strat)
      stratTk[[strat]] <- attr(FarFirst(k, d, strategy = as.numeric(seedIdx)), "score")
      set.seed(1L)
      r <- DropAdd(k, d, seed = seedIdx)
      stratScores[[strat]] <- Scv(d, r)
    }
    # "Ensemble-initializer" arm: the anchor the 7-way FarFirst ensemble would
    # itself pick (largest construction T_k), scored by its DropAdd outcome.
    ensemblePickScore <- stratScores[[names(which.max(stratTk))]]
    bestOf7Score      <- max(stratScores)

    rows[[length(rows) + 1L]] <- c(
      list(case = caseName, n = n, k = k, optScore = optScore,
           proven = proven, optWall = optWall, defaultScore = defScore,
           ensemblePickScore = ensemblePickScore, bestOf7Score = bestOf7Score),
      as.list(stratScores)
    )
    cat(sprintf(
      "%-22s n=%4d k=%2d opt=%.6f(proven=%s,%.1fs) default=%.6f  ensemble_pick=%.6f  best_of_7=%.6f\n",
      caseName, n, k, optScore, proven, optWall, defScore, ensemblePickScore, bestOf7Score))
  }
}

df <- do.call(rbind.data.frame, rows)
saveRDS(df, "dev/profiling/dropadd-seed-quality.rds")

cat("\n=== Gap to proven optimum (optScore - score), lower is better ===\n")
gapCols <- c("defaultScore", strategies, "ensemblePickScore", "bestOf7Score")
gaps <- sapply(gapCols, function(cn) df$optScore - df[[cn]])
gapSummary <- data.frame(
  strategy = gapCols,
  mean_gap = colMeans(gaps),
  max_gap  = apply(gaps, 2, max),
  n_optimal = colSums(abs(gaps) < 1e-9)
)
print(gapSummary[order(gapSummary$mean_gap), ], row.names = FALSE)

cat("\n=== Cases where any alt strategy beat default (default not optimal) ===\n")
worseThanBest <- df$defaultScore < apply(df[strategies], 1, max) - 1e-9
print(df[worseThanBest, c("case", "n", "k", "optScore", "defaultScore")])

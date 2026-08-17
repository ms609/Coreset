# A/B validation: MaxMean() vs published best-known f_best on the MDPI Type-I
# n=500 benchmark instances (Marti/Duarte; best-known from Lai & Hao 2016,
# reproduced in Nijimbere et al. 2020 Table 1).
#
# RESULT (2026-06-17, useRL = TRUE, 30 s/instance, set.seed(1)): MaxMean reaches
# the published best-known value on all three instances tested (gap < 5e-7, i.e.
# floating-point noise):
#     MDPI1_500  81.277044  MATCH      MDPI2_500  78.610216  MATCH
#     MDPI3_500  76.300787  MATCH
# This validates the RL-initialisation fidelity fixes (red-team MM-04/05/06).
#
# Instances are NOT committed (the full set is ~64 MB). To reproduce:
#   1. Download https://grafo.etsii.urjc.es/optsicom/edp/edp_instances.zip
#   2. Extract typeI/MDPI{1,2,3}_500.txt into dev/_mdpi/
#   3. Rscript dev/ab-mdpi.R           (or pass instance names as args)
#      AB_BUDGET=100 matches the paper's n=500 time limit.
#
# Instance format: lines "i j d_ij" (1-based, upper triangle, signed distances).
.libPaths(c(".agent-mm", .libPaths()))
suppressMessages(library("Coreset"))

best_known <- c(MDPI1_500 = 81.277044, MDPI2_500 = 78.610216, MDPI3_500 = 76.300787)

ReadInstance <- function(path) {
  tab <- scan(path, what = list(i = integer(), j = integer(), d = double()),
              quiet = TRUE)
  n <- max(tab$i, tab$j)
  d <- matrix(0, n, n)
  idx_ij <- cbind(tab$i, tab$j)
  idx_ji <- cbind(tab$j, tab$i)
  d[idx_ij] <- tab$d
  d[idx_ji] <- tab$d
  d
}

budget <- as.numeric(Sys.getenv("AB_BUDGET", unset = "30"))
args <- commandArgs(trailingOnly = TRUE)
insts <- if (length(args)) args else names(best_known)

cat(sprintf("budget = %g s/instance\n", budget))
cat(sprintf("%-12s %12s %12s %10s %8s %6s\n",
            "instance", "best_known", "MaxMean_f", "gap", "iters", "|S|"))
for (nm in insts) {
  d <- ReadInstance(file.path("dev/_mdpi", paste0(nm, ".txt")))
  set.seed(1)
  # maxIter = Inf: budget purely by time (the n=500 paper protocol). Without
  # this the default maxIter = 1000 would stop far short of convergence.
  res <- MaxMean(d, maxSeconds = budget, maxIter = Inf, useRL = TRUE)
  f  <- attr(res, "score")
  bk <- best_known[[nm]]
  cat(sprintf("%-12s %12.6f %12.6f %+10.2e %8.3g %6d   %s\n",
              nm, bk, f, f - bk, attr(res, "iters"), length(res),
              if (f >= bk - 1e-4) "MATCH" else if (f >= bk - 0.05) "~near" else "BELOW"))
}

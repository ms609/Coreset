# Quality-vs-time frontier harness for Grasp (round 4).
#
# Forked from grasp-battery.R (which asserts bit-identity and remains the
# verifier for identity-preserving changes). This driver measures what the
# battery cannot: the two comparisons that matter once a change is allowed to
# alter the search trajectory (dev/profiling/ROUND4_PLAN.md section 5):
#
#   (a) equal plateau -- wall time and objective at comparable search effort.
#       For an identity-preserving change the objectives must agree exactly,
#       which doubles as an identity spot-check at canonical n = 2000 shapes
#       the battery's n <= 500 grid never reaches.
#   (b) equal budget -- objective attained within a fixed maxSeconds
#       wall-clock budget: the quality-per-unit-time frontier. Run in two
#       sub-modes: plateau disabled ("open", pure iteration-throughput
#       quality) and canonical plateau 512 ("canon", the manuscript shape).
#       When the budget binds, Grasp_cpp skips phase C entirely, so a faster
#       kernel may finish phase B inside the budget and reach path relinking
#       while the slower build gets none; iters and pr_calls are recorded so
#       that mechanism stays visible in the comparison.
#
# GRASP is randomised, so every cell runs several seeds, and the comparison
# reports the mean AND the worst per-seed delta -- a lever that raises the
# mean while occasionally collapsing is a regression for a manuscript claim.
# Budget-gated runs are wall-clock-gated and therefore not exactly
# reproducible run to run; the seed spread is the signal, not any single run.
#
# Usage: Rscript grasp-frontier.R <out.rds>            # capture
#        Rscript grasp-frontier.R <new.rds> <base.rds> # capture + compare
#
# Wall-clock on this box gives *relative* ratios for choosing between builds,
# not quotable figures (those belong on Hamilton).
library(Coreset)
args <- commandArgs(trailingOnly = TRUE)
outFile <- args[[1L]]

PLATEAU_OFF <- .Machine$integer.max
CANON <- 512L

# Instances: canonical synthetic coreset shape (n = 2000, dim = 10), plus one
# real case from FurthestPoint thinned to 2000 by farthest-first, mirroring
# the canonical pipeline's maxCandidates = 2000 coreset. Data is fixed per
# instance (independent of the run seed) so the seed spread is pure algorithm
# randomness on one instance.
MakeSynthetic <- function(n, dim, dataSeed) {
  set.seed(dataSeed)
  as.matrix(dist(matrix(rnorm(n * dim), n)))
}
instances <- list(syn2000 = MakeSynthetic(2000L, 10L, 42L))

casesPaths <- c("../furthest-point/data/cases.rda",
                "../../../furthest-point/data/cases.rda",
                "C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
casesPath <- Filter(file.exists, casesPaths)
if (length(casesPath)) {
  casesEnv <- new.env()
  load(casesPath[[1L]], envir = casesEnv)
  pts <- casesEnv$cases$tc14_satellite$points
  set.seed(42)
  core <- FarFirst(2000L, points = pts)
  instances$sat2000 <- as.matrix(dist(pts[core, , drop = FALSE]))
} else {
  message("cases.rda not found; running synthetic instance only")
}

RunCell <- function(d, k, plateau, budget, seed) {
  set.seed(seed)
  t0 <- proc.time()[[3L]]
  o <- Coreset:::Grasp_cpp(d, k, plateau, .Machine$integer.max, 10L, 0.8,
                          budget)
  list(time = proc.time()[[3L]] - t0, obj = as.numeric(o$objective),
       iters = as.numeric(o$iters), pr = as.numeric(o$pr_calls))
}

Cells <- function(inst, mode, k, plateau, budget, seeds) {
  expand.grid(inst = inst, mode = mode, k = k, plateau = plateau,
              budget = budget, seed = seeds,
              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

grid <- rbind(
  # (a) equal plateau: identity + throughput, both ladder ends.
  Cells("syn2000", "plateau", 100L, c(8L, 64L, 256L), Inf, 1:3),
  Cells("syn2000", "plateau",  10L, c(64L, 256L),     Inf, 1:3),
  # (b) equal budget, plateau off: pure quality-per-second.
  Cells("syn2000", "open",  100L, PLATEAU_OFF, c(1, 3, 8),    1:5),
  Cells("syn2000", "open",   10L, PLATEAU_OFF, c(0.25, 1),    1:5),
  # (b) equal budget at the canonical plateau: the manuscript shape.
  Cells("syn2000", "canon", 100L, CANON,       c(1, 3, 8),    1:5),
  Cells("syn2000", "canon",  10L, CANON,       c(0.25, 1),    1:5)
)
if (!is.null(instances$sat2000)) {
  grid <- rbind(
    grid,
    Cells("sat2000", "plateau", 100L, 64L,   Inf,     1:3),
    Cells("sat2000", "canon",   100L, CANON, c(1, 3), 1:5)
  )
}

res <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  res[[i]] <- RunCell(instances[[g$inst]], g$k, g$plateau, g$budget, g$seed)
}
saveRDS(list(grid = grid, res = res), outFile)
cat("captured", nrow(grid), "cells ->", outFile, "\n")

if (length(args) >= 2L) {
  base <- readRDS(args[[2L]])
  stopifnot(identical(dim(base$grid), dim(grid)))
  Pull <- function(rr, f) vapply(rr, `[[`, numeric(1L), f)
  newObj <- Pull(res, "obj");        baseObj <- Pull(base$res, "obj")
  newT   <- Pull(res, "time");       baseT   <- Pull(base$res, "time")
  newIt  <- Pull(res, "iters");      baseIt  <- Pull(base$res, "iters")
  newPr  <- Pull(res, "pr");         basePr  <- Pull(base$res, "pr")

  cat("\n=== FRONTIER: (a) equal plateau ===\n")
  pa <- grid$mode == "plateau"
  cellKey <- paste(grid$inst, "k", grid$k, "plateau", grid$plateau)
  for (key in unique(cellKey[pa])) {
    sub <- pa & cellKey == key
    ident <- all(newObj[sub] == baseObj[sub])
    cat(sprintf("%-32s speed %5.2fx  objectives %s\n", key,
                sum(baseT[sub]) / max(sum(newT[sub]), 1e-9),
                if (ident) "IDENTICAL" else "*** DIFFER ***"))
  }

  cat("\n=== FRONTIER: (b) equal budget (dT_k = new - base; + is better) ===\n")
  pb <- grid$mode %in% c("open", "canon")
  cellKey <- paste(grid$inst, grid$mode, "k", grid$k, "budget", grid$budget)
  for (key in unique(cellKey[pb])) {
    sub <- pb & cellKey == key
    d <- newObj[sub] - baseObj[sub]
    rel <- d / baseObj[sub]
    cat(sprintf(paste0("%-36s mean %+.3e (%+6.2f%%)  worst %+.3e  ",
                       "win/tie/loss %d/%d/%d  iters %.0f->%.0f  pr %.0f->%.0f\n"),
                key, mean(d), 100 * mean(rel), min(d),
                sum(d > 0), sum(d == 0), sum(d < 0),
                mean(baseIt[sub]), mean(newIt[sub]),
                mean(basePr[sub]), mean(newPr[sub])))
  }
}

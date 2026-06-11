# dev/profiling/drivers/exact.R
# Area 4 (user-requested): ExactMaxMin node-packing solver.
# Pure-R orchestration around highs MILP solves. Question: where does the wall
# time go -- R-side construction (edge enumeration + DENSE constraint matrix)
# vs the highs IP solves -- and is the dense `A` the scaling wall?
#
# Representative workload: penguins (n=342), the smallest of the six datasets
# the manuscript currently cannot solve exactly. k=10 (the manuscript regime).
#
# Override lib via FP_LIB env (timestamped profiling install); default installed.
set.seed(5813)
lib <- Sys.getenv("FP_LIB", "")
if (nzchar(lib)) library(MaxMin, lib.loc = lib) else library(MaxMin)

load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
case <- Sys.getenv("FP_CASE", "tc22_penguins")
k    <- as.integer(Sys.getenv("FP_K", "10"))
pts  <- as.matrix(cases[[case]][["points"]]); storage.mode(pts) <- "double"
d    <- as.matrix(stats::dist(pts))
cat(sprintf("case=%s n=%d k=%d  (%d distinct distances)\n",
            case, nrow(d), k, length(unique(d[upper.tri(d)]))))

gc(reset = TRUE)
t0 <- proc.time()[[3L]]
r  <- MaxMin::ExactMaxMin(d, m = k, maxSeconds = 600, progress = FALSE)
el <- proc.time()[[3L]] - t0
g  <- gc()
cat(sprintf("Elapsed: %.2f s   obj=%.6f  proven=%s\n", el, r$objective, r$proven))
cat(sprintf("R-side peak: Ncells %.1f Mb / Vcells %.1f Mb\n", g[1, 6], g[2, 6]))

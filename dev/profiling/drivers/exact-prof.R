# Rprof harness for ExactMaxMin (area 4). Line + memory profiling on the
# penguins n=342 k=10 workload; summarise self/total by function and by line.
set.seed(5813)
library(MaxMin)
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
pts <- as.matrix(cases[["tc22_penguins"]][["points"]]); storage.mode(pts) <- "double"
d   <- as.matrix(stats::dist(pts))

src <- file.path(find.package("MaxMin"), "R")   # for line numbers if available
out <- tempfile(fileext = ".out")
Rprof(out, line.profiling = TRUE, memory.profiling = TRUE, interval = 0.005)
r <- MaxMin::ExactMaxMin(d, m = 10L, timeBudgetS = 600, progress = FALSE)
Rprof(NULL)
cat(sprintf("obj=%.6f proven=%s\n\n", r$objective, r$proven))

s <- summaryRprof(out, lines = "show", memory = "both")
cat("===== BY SELF TIME (function) =====\n")
bs <- summaryRprof(out, memory = "both")$by.self
print(head(bs, 18))
cat("\n===== BY LINE (self) =====\n")
print(head(s$by.self, 20))

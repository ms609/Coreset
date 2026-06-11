# Per-probe cost instrumentation for the binary search, + highs option listing.
# Where is the solver time spent: a few hard (infeasible) probes near the
# optimum, or spread across many?
suppressMessages({ library(MaxMin); library(Matrix); library(highs) })
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")

cat("===== highs solver options (early-termination relevant) =====\n")
op <- tryCatch(highs_available_solver_options(), error = function(e) NULL)
if (!is.null(op)) {
  nm <- names(op); if (is.data.frame(op)) {
    keep <- grepl("gap|bound|target|limit|sol", op$option, ignore.case = TRUE)
    print(op[keep, intersect(c("option","type","default","range"), nm)])
  } else print(op)
}

probe_log <- function(case, k = 10L) {
  pts <- as.matrix(cases[[case]][["points"]]); storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts)); n <- nrow(d)
  ut <- which(upper.tri(d)); rc <- arrayInd(ut, dim(d))
  ui <- rc[,1L]; uj <- rc[,2L]; ud <- d[ut]
  cand <- sort(unique(ud)); nCand <- length(cand)
  lo <- 2L; hi <- nCand; rows <- list()
  while (lo <= hi) {
    mid <- (lo + hi) %/% 2L; lambda <- cand[mid]
    e <- ud < lambda; nE <- sum(e)
    A <- sparseMatrix(i = rep(seq_len(nE), 2L), j = c(ui[e], uj[e]), x = 1, dims = c(nE, n))
    t <- proc.time()[[3L]]
    res <- highs_solve(L = rep.int(1,n), lower = rep.int(0,n), upper = rep.int(1,n),
      A = A, lhs = rep.int(-Inf,nE), rhs = rep.int(1,nE), types = rep.int("I",n),
      maximum = TRUE, control = list(threads = 1L, time_limit = 600))
    st <- proc.time()[[3L]] - t
    alpha <- sum(res$primal_solution > 0.5)
    feas <- alpha >= k
    rows[[length(rows)+1L]] <- data.frame(mid=mid, lambda=round(lambda,4), nEdge=nE,
      alpha=alpha, feasible=feas, solve_s=round(st,3))
    if (feas) lo <- mid + 1L else hi <- mid - 1L
  }
  R <- do.call(rbind, rows)
  cat(sprintf("\n===== %s  n=%d k=%d  (%d probes, total solve %.2fs) =====\n",
              case, n, k, nrow(R), sum(R$solve_s)))
  print(R, row.names = FALSE)
  cat(sprintf("  feasible probes: %d (%.2fs) | infeasible probes: %d (%.2fs)\n",
      sum(R$feasible), sum(R$solve_s[R$feasible]),
      sum(!R$feasible), sum(R$solve_s[!R$feasible])))
}
for (cs in c("tc22_penguins", "tc11_ionosphere")) probe_log(cs)

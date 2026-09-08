# Does the Kuo--Glover--Dhir contribution cap pay for being exact? COUNTS ONLY
# -- this script makes no wall-clock claim; timing for this package belongs on
# Hamilton.
#
# Two arms of the same model: `U(k)` caps a_i by the k largest distances from
# i, `U(k-1)` by the k - 1 largest, which is what a_i actually sums over. For
# each arm it records the branch-and-bound node count, the root relaxation
# value, and the certified objective, so the tighter arm can be checked for
# reaching the same optimum rather than a cheaper wrong one.
#
# Usage: Rscript maxsum-ucap.R [out.rds]
suppressMessages(library("Coreset"))
stopifnot(requireNamespace("highs", quietly = TRUE),
          requireNamespace("Matrix", quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)
outFile <- if (length(args)) args[[1L]] else "maxsum-ucap.rds"

Dist <- function(n, dim, seed) {
  set.seed(seed)
  as.matrix(dist(matrix(rnorm(n * dim), ncol = dim)))
}
Score <- function(d, idx) {
  s <- d[idx, idx, drop = FALSE]
  sum(s[upper.tri(s)])
}

# The model ExactMaxSum() builds, with the cap width and the integrality
# switch exposed.
Kgd <- function(d, k, terms, relax, timeLimit) {
  n <- nrow(d)
  U <- vapply(seq_len(n), function(i) {
    sum(sort(d[i, -i], decreasing = TRUE)[seq_len(min(terms, n - 1L))])
  }, double(1))
  ri <- rep(1L, n); ci <- seq_len(n); vi <- rep(1, n)
  rA <- 1L + seq_len(n)
  ri <- c(ri, rA, rA); ci <- c(ci, n + seq_len(n), seq_len(n))
  vi <- c(vi, rep(1, n), -U)
  rB <- 1L + n + seq_len(n)
  ri <- c(ri, rB, rep(rB, each = n))
  ci <- c(ci, n + seq_len(n), rep(seq_len(n), times = n))
  vi <- c(vi, rep(1, n), -as.vector(d))
  A <- Matrix::sparseMatrix(i = ri, j = ci, x = vi,
                            dims = c(1L + 2L * n, 2L * n))
  res <- highs::highs_solve(
    L = c(rep(0, n), rep(0.5, n)), lower = rep(0, 2L * n),
    upper = c(rep(1, n), U), A = A,
    lhs = c(k, rep(-Inf, 2L * n)), rhs = c(k, rep(0, 2L * n)),
    types = if (relax) rep("C", 2L * n) else c(rep("I", n), rep("C", n)),
    maximum = TRUE,
    control = list(threads = 1L, time_limit = timeLimit))
  sel <- sort(which(res$primal_solution[seq_len(n)] > 0.5))
  nodes <- res$info[["mip_node_count"]]
  list(status = res$status_message,
       value = if (length(sel) == k) Score(d, sel) else NA_real_,
       bound = res$objective_value,
       nodes = if (is.null(nodes)) NA_real_ else as.numeric(nodes))
}

cases <- expand.grid(k = c(5L, 10L), n = c(30L, 40L, 50L, 60L))
out <- NULL
for (i in seq_len(nrow(cases))) {
  n <- cases$n[[i]]; k <- cases$k[[i]]
  d <- Dist(n, 5L, seed = 1L)
  loose <- Kgd(d, k, k, FALSE, 120)
  tight <- Kgd(d, k, k - 1L, FALSE, 120)
  out <- rbind(out, data.frame(
    n = n, k = k,
    looseStatus = loose$status, looseNodes = loose$nodes,
    looseRoot = Kgd(d, k, k, TRUE, 120)$bound, looseValue = loose$value,
    tightStatus = tight$status, tightNodes = tight$nodes,
    tightRoot = Kgd(d, k, k - 1L, TRUE, 120)$bound, tightValue = tight$value))
  cat(sprintf("n=%3d k=%2d | %-18s %9.0f nodes | %-18s %9.0f nodes | same=%s\n",
              n, k, loose$status, loose$nodes, tight$status, tight$nodes,
              isTRUE(all.equal(loose$value, tight$value))))
  saveRDS(out, outFile)
}
cat("Written to ", outFile, "\n", sep = "")

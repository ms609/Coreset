# verify.R <libpath> <tag>
# Per-lib correctness + timing battery for the maximin kernel changes
# (squared-space pass + free T_k). Writes dev/profiling/verify-<tag>.rds.
#
#   Correctness recorded per case:
#     - matrix-path selection and points-path selection (for cross-path identity
#       within this build, and old-vs-new identity across builds)
#     - for ensembles: each strategy's reported t_k vs a MinDist() reference
#   Timing: end-to-end FarFirst full_ms at canonical (N, dim, n).
args <- commandArgs(trailingOnly = TRUE)
libpath <- normalizePath(args[[1]])
tag     <- args[[2]]
suppressMessages({ library(MaxMin, lib.loc = libpath); library(bench) })

bench1 <- function(thunk, iters = 15L) {
  m <- bench::mark(thunk(), iterations = iters, check = FALSE, filter_gc = FALSE)
  as.numeric(m$median) * 1000
}

# ---- correctness battery ----------------------------------------------------
cases <- list()
add <- function(...) cases[[length(cases) + 1L]] <<- list(...)

fixedPivots <- c(1L, 7L, 13L)
detAnchors  <- c("peripheral", "diameter", "medoid", "anti_medoid",
                 "rowsum", "rownorm")

corr <- list()
for (dseed in 1:6) {
  for (N in c(150L, 400L)) {
    for (dim in c(2L, 3L, 8L)) {
      set.seed(dseed)
      pts <- matrix(rnorm(N * dim), ncol = dim)
      d   <- as.matrix(dist(pts))
      ns  <- unique(c(2L, 5L, as.integer(round(N * 0.1)),
                      as.integer(round(N * 0.5)), N - 1L))
      strategies <- c(
        lapply(detAnchors, function(s) list(kind = "single", seed = s)),
        list(list(kind = "int", seed = 1L)),
        list(list(kind = "ens", seed = c("diameter", "rowsum", "anti_medoid"))),
        list(list(kind = "ens", seed = "random_furthest"))
      )
      for (n in ns) {
        for (st in strategies) {
          key <- paste(dseed, N, dim, n, st$kind,
                       paste(st$seed, collapse = "+"), sep = "|")
          matSel <- FarFirst(d, n, method = st$seed, pivots = fixedPivots)
          ptSel  <- FarFirst(m = n, points = pts, method = st$seed,
                             pivots = fixedPivots)
          # t_k fidelity: compare every reported strategy t_k to MinDist().
          tkErr <- NA_real_
          sr <- attr(matSel, "strategy_results")
          if (!is.null(sr)) {
            errs <- vapply(sr, function(r) {
              ref <- if (length(r$idx) >= 2L) MinDist(d, r$idx) else NA_real_
              if (is.na(r$t_k) && is.na(ref)) 0 else abs(r$t_k - ref)
            }, numeric(1L))
            tkErr <- max(errs, na.rm = TRUE)
          }
          corr[[key]] <- list(
            mat = as.integer(matSel), pts = as.integer(ptSel),
            cross_ident = identical(as.integer(matSel), as.integer(ptSel)),
            tkErr = tkErr
          )
        }
      }
    }
  }
}

# ---- timing -----------------------------------------------------------------
set.seed(5813)
N <- 6000L
timing <- list()
for (dim in c(2L, 10L)) {
  pts <- matrix(rnorm(N * dim), ncol = dim)
  d   <- as.matrix(dist(pts))
  piv <- c(11L, 222L, 3333L)
  for (n in c(60L, 300L, 1200L, 3000L)) {
    fp <- function() FarFirst(m = n, points = pts, method = "random_furthest",
                              pivots = piv)
    fm <- function() FarFirst(d, n, method = "random_furthest", pivots = piv)
    timing[[paste(dim, n, sep = "|")]] <- list(
      dim = dim, n = n, ratio = n / N,
      points_ms = round(bench1(fp), 2),
      matrix_ms = round(bench1(fm), 2)
    )
  }
}

saveRDS(list(tag = tag, corr = corr, timing = timing),
        sprintf("dev/profiling/verify-%s.rds", tag))
cat(sprintf("[%s] %d correctness cases, %d timing cases written\n",
            tag, length(corr), length(timing)))

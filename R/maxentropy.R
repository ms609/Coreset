# maxentropy.R
#
# Maximum-entropy (maxdet) subset selection: choose the k-subset maximising
# log det K_S for a similarity kernel K built from a distance matrix -- the
# Shewry & Wynn (1987) maximum-entropy sampling criterion, equivalently the
# maximum-a-posteriori mode of a determinantal point process (Kulesza & Taskar
# 2012). The combinatorial core (greedy pivoted Cholesky and exact enumeration)
# is in src/maxentropy.cpp; the kernel construction and its positive-
# semidefinite repair are orchestrated here, on top of the LAPACK partial
# eigendecomposition in src/symeigen.cpp, so the repair stays a transparent
# modelling step.
#
# A redundant tree adds zero volume (det -> 0) and is never co-selected, so the
# objective is exactly density-blind. Because the exact argmax is NP-hard the
# greedy selector is the workhorse; exact enumeration is used only where
# choose(n, k) is small enough to certify the optimum.

# ----- kernel + PSD repair ---------------------------------------------------

# RBF / Gaussian similarity kernel from a symmetric distance matrix. `sigma`
# defaults to the median of the POSITIVE distances -- robust when many pairs
# are exact duplicates, where the median over all pairs would collapse toward
# zero. The median and the kernel are computed in C++ over one triangle and
# equal `median(d[d > 0])` and `exp(-(d ^ 2) / (2 * sigma ^ 2))` bit for bit.
.MaxEntropyKernel <- function(d, sigma = NULL) {
  if (is.null(sigma)) {
    sigma <- MedianPositiveUpper_cpp(d)
    if (is.na(sigma)) sigma <- 1
  }
  k <- RbfKernel_cpp(d, sigma)
  dimnames(k) <- dimnames(d)
  attr(k, "sigma") <- sigma
  # Return:
  k
}

# Positive-semidefinite repair of the kernel, and the negative-eigenvalue mass:
#
#   clip     -- eigen-clip negatives to zero (the nearest PSD matrix)
#   shift    -- diagonal loading by |lambda_min| (preserves eigenvectors)
#   truncate -- retain the top positive dimensions holding `keep` of the
#               positive eigen-mass (a lossy low-rank embedding)
#
# negMass is the fraction of total |eigenvalue| carried by the negatives -- the
# magnitude of the repair, reported so a caller sees when the Euclidean
# approximation is doing real work.
#
# Cost: this is MaxEntropy()'s O(n^3) step. A kernel that is positive-definite
# is certified by one Cholesky factorisation (n^3 / 3) and returned as is.
# For clip the factorisation is of ks + delta I, delta = min(4 n eps ||ks||,
# tol), so a kernel that is positive-definite to within round-off -- genuine
# Euclidean distances in low dimension, whose computed spectrum has many
# round-off negatives -- is certified too: any negative eigenvalue is then
# within the resolution of a dense eigen-solver, the clip would change ks
# only by round-off, and negMass is 0 by its definition. shift gets no such
# margin (delta = 0): its repair of a round-off-indefinite kernel is a ridge
# of `tol`, not a round-off change, so it keeps the plain positive-definite
# test. Otherwise SymEigenPartial_cpp() tridiagonalises once
# (4/3 n^3), takes every eigenvalue from the tridiagonal form (O(n^2)), and
# computes only the eigenvectors the repair needs: the smaller side of zero for
# clip, none for shift, and for truncate the smaller of the kept and dropped
# sets. The repaired kernel is a rank-p update of ks (n^2 p) -- subtracting the
# dropped components or rebuilding from the kept ones.
#
# `symmetric = TRUE` promises that `k` is exactly symmetric, as
# .MaxEntropyKernel() returns it, and skips the averaging (two n^2 passes);
# `k` is then used as is, attributes included.
.MaxEntropyPrepare <- function(k, method = c("clip", "shift", "truncate"),
                               keep = 0.99, tol = 1e-9, symmetric = FALSE) {
  method <- match.arg(method)
  if (symmetric) {
    ks <- k
  } else {
    ks <- (k + t(k)) / 2
    attributes(ks) <- attributes(ks)["dim"]      # a bare matrix on every path
  }
  if (method != "truncate" &&
      CholCertificate_cpp(ks, if (method == "clip") tol else 0)) {
    # Return: numerically positive-definite; nothing to repair.
    return(list(kp = ks, negMass = 0))
  }
  eig <- SymEigenPartial_cpp(ks, switch(method, clip = 1L, shift = 0L,
                                        truncate = 2L), keep)
  lam <- eig[["values"]]                                   # ascending
  neg <- lam[lam < -tol]
  negMass <- if (length(lam)) sum(abs(neg)) / sum(abs(lam)) else 0
  kp <- switch(method,
    shift = if (length(lam) && lam[[1]] < 0) {
      ks + (-lam[[1]] + tol) * diag(nrow(ks))
    } else {
      ks
    },
    # clip: subtract the negative eigen-component (kp = ks - V L- V^T, a rank-m
    # update) or, when the negatives are the majority, rebuild from the
    # positive one (kp = V L+ V^T). truncate: likewise from whichever of the
    # kept and dropped components is the smaller set.
    clip = ,
    truncate = if (eig[["side"]] == 0L) {
      ks
    } else if (eig[["side"]] < 0L) {
      RankUpdate_cpp(ks, eig[["vectors"]], -eig[["pvalues"]])
    } else {
      RankUpdate_cpp(NULL, eig[["vectors"]], eig[["pvalues"]])
    }
  )
  # Return:
  list(kp = kp, negMass = negMass)
}

# Back-compat thin wrappers (used by the tests and any external callers).
.MaxEntropyNegMass <- function(k, tol = 1e-9) {
  # Return:
  .MaxEntropyPrepare(k, "clip", tol = tol)$negMass
}
.MaxEntropyRepair <- function(k, method = c("clip", "shift", "truncate"),
                              keep = 0.99, tol = 1e-9) {
  # Return:
  .MaxEntropyPrepare(k, match.arg(method), keep, tol)$kp
}

# ----- exported solver ------------------------------------------------------

#' Maximum-entropy (maxdet) subset selection
#'
#' `MaxEntropy()` selects the `k` points that maximise the log-determinant of
#' their kernel block, \eqn{\log\det K_S}. This is equivalent to finding the set
#' of `k` points that span the largest volume, which corresponds to
#' the maximum-entropy sampling criterion \insertCite{Shewry1987}{Coreset} and
#' the maximum-_a-posteriori_ mode of a determinantal point process
#' \insertCite{Kulesza2012}{Coreset}.
#'
#' A radial-basis kernel \eqn{K_{ij} = \exp(-d_{ij}^2 / 2\sigma^2)} is built from
#' the supplied distances and repaired to a positive-semidefinite matrix.
#' The exact argmax is NP-hard \insertCite{Kulesza2012}{Coreset}.
#' A greedy approximation is built by pivoted Cholesky, adding at each step
#' the point of largest residual conditional variance.
#' Ties are broken by selecting the more peripheral point.
#'
#' Large instances run faster when \R is linked to an optimised BLAS, such as
#' OpenBLAS, Intel MKL or Apple Accelerate.
#'
#' @param k Integer specifying target selection size, \eqn{1 \le k \le n}.
#' @param d `dist` object or square numeric distance matrix over the \eqn{n}
#'   points.
#' @param sigma Optional numeric specifying kernel bandwidth;
#' defaults to the median positive distance.
#' @param repair Character selecting a positive semi-definite repair method:
#' `"clip"` (nearest), `"shift"` (diagonal loading) or `"truncate"` (low-rank
#' embedding).
#' @param exact Logical: `TRUE` uses explicit enumeration, failing with an error
#' if `maxCombos` is exceeded; `FALSE` uses the greedy approximation.
#' `NA` uses exact enumeration when `choose(n, k) <= maxCombos`,
#' greedy otherwise.
#' @param maxCombos Integer specifying ceiling on `choose(n, k)` for exact
#' enumeration.
#' @return `MaxEntropy()` returns an integer vector of length `k` (sorted
#'   ascending) with class `"MaxEntropySelection"`, carrying attributes:
#'   \describe{
#'     \item{score}{The retained \eqn{\log\det K_S} of the selection.
#'       `-Inf` is returned for a degenerate selection where `k` exceeds the
#'       number of distinct points.}
#'     \item{negMass}{Fraction of spectral mass removed by the positive
#'       semi-definite repair.}
#'     \item{sigma, repair, exact}{The bandwidth, repair, and whether the
#'       optimum was certified by enumeration.}
#'     \item{seed, N, k}{The peripheral seed index, instance size, target size.}
#'   }
#' @references \insertAllCited{}
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(40), ncol = 2)
#' MaxEntropy(4L, dist(pts))
#' @export
MaxEntropy <- function(k, d, sigma = NULL,
                       repair = c("clip", "shift", "truncate"),
                       exact = NA, maxCombos = 3e5L) {
  repair <- match.arg(repair)
  d <- .ExactAsMatrix(d)
  if (!all(is.finite(d))) {
    stop("`d` must contain only finite values (no NA, NaN, or Inf)")
  }
  n <- nrow(d)
  k <- as.integer(k)
  if (is.na(k) || k < 1L || k > n) {
    stop("`k` must satisfy 1 <= k <= nrow(d)")
  }
  nDistinct <- DistinctRows_cpp(d)                 # == sum(!duplicated(d))
  if (k > nDistinct) {
    warning(sprintf(paste0("`k` (%d) exceeds the number of distinct points (%d); ",
                           "the selection must repeat near-identical points and ",
                           "its log-determinant is -Inf (degenerate)."),
                    k, nDistinct))
  }

  kern <- .MaxEntropyKernel(d, sigma)
  sigmaUsed <- attr(kern, "sigma")
  prep <- .MaxEntropyPrepare(kern, repair, symmetric = TRUE)
  negMass <- prep$negMass
  kp <- prep$kp

  cnk <- choose(n, k)
  useExact <- if (is.na(exact)) is.finite(cnk) && cnk <= maxCombos else isTRUE(exact)
  if (useExact && (!is.finite(cnk) || cnk > maxCombos)) {
    stop(sprintf(paste0("Exact enumeration needs choose(n, k) = %.3g <= ",
                        "maxCombos = %.3g; raise `maxCombos` or use the greedy ",
                        "(`exact = FALSE`)."), cnk, maxCombos))
  }

  seed <- which.min(rowSums(kp))
  if (useExact) {
    idx <- MaxEntropyExact_cpp(kp, k)
  } else {
    idx <- sort(MaxEntropyGreedy_cpp(kp, k, as.integer(seed)))
  }

  # Report the Cholesky log-determinant.
  # When k exceeds the distinct-point count the selection must repeat a point,
  # so its Gram matrix is singular and the true log-determinant is -Inf.
  logDet <- if (k > nDistinct) -Inf else MaxEntropyLogDet_cpp(kp, idx)

  # Return:
  structure(idx,
            score = logDet, negMass = negMass, sigma = sigmaUsed,
            repair = repair, exact = useExact, seed = as.integer(seed),
            N = n, k = k,
            producer = "MaxEntropy", class = "MaxEntropySelection")
}

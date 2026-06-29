# maxentropy.R
#
# Maximum-entropy (maxdet) subset selection: choose the k-subset maximising
# log det K_S for a similarity kernel K built from a distance matrix -- the
# Shewry & Wynn (1987) maximum-entropy sampling criterion, equivalently the
# maximum-a-posteriori mode of a determinantal point process (Kulesza & Taskar
# 2012). The combinatorial core (greedy pivoted Cholesky and exact enumeration)
# is in src/maxentropy.cpp; the kernel construction and its positive-
# semidefinite repair are kept here in R, reusing the same LAPACK eigendecom-
# position `eigen()` provides so the repair is a transparent modelling step.
#
# A redundant tree adds zero volume (det -> 0) and is never co-selected, so the
# objective is exactly density-blind. Because the exact argmax is NP-hard the
# greedy selector is the workhorse; exact enumeration is used only where
# choose(n, k) is small enough to certify the optimum.

# ----- kernel + PSD repair (in R; LAPACK via eigen()) -----------------------

# RBF / Gaussian similarity kernel from a distance matrix. `sigma` defaults to
# the median of the POSITIVE distances -- robust when many pairs are exact
# duplicates, where the median over all pairs would collapse toward zero.
.MaxEntropyKernel <- function(d, sigma = NULL) {
  if (is.null(sigma)) {
    pos <- d[d > 0]
    sigma <- if (length(pos)) stats::median(pos) else 1
  }
  k <- exp(-(d ^ 2) / (2 * sigma ^ 2))
  attr(k, "sigma") <- sigma
  # Return:
  k
}

# One eigendecomposition serving both outputs the wrapper needs from the kernel:
# the positive-semidefinite repair and the negative-eigenvalue mass. eigen() is
# the wrapper's dominant cost at large n, so it is computed ONCE here (rather than
# once for the repair and again, values-only, for negMass), and the default clip
# repair is reconstructed with a symmetric rank-k product (tcrossprod / dsyrk),
# ~3.6x faster than the general V (Lambda V^T) matmul and numerically identical.
#
#   clip     -- eigen-clip negatives to zero (the nearest PSD matrix)
#   shift    -- diagonal loading by |lambda_min| (preserves eigenvectors)
#   truncate -- retain the top positive dimensions holding `keep` of the
#               positive eigen-mass (a lossy low-rank embedding)
# negMass is the fraction of total |eigenvalue| carried by the negatives -- the
# magnitude of the repair, reported so a caller sees when the Euclidean
# approximation is doing real work.
.MaxEntropyPrepare <- function(k, method = c("clip", "shift", "truncate"),
                               keep = 0.99, tol = 1e-9) {
  method <- match.arg(method)
  ks <- (k + t(k)) / 2
  e <- eigen(ks, symmetric = TRUE)
  lam <- e$values
  neg <- lam[lam < -tol]
  negMass <- if (length(lam)) sum(abs(neg)) / sum(abs(lam)) else 0
  if (method == "shift") {
    lmin <- min(lam)
    kp <- if (lmin < 0) ks + (-lmin + tol) * diag(nrow(ks)) else ks
  } else {
    lam[lam < 0] <- 0
    if (method == "truncate") {
      pos <- lam[lam > 0]
      if (length(pos)) {
        cum <- cumsum(sort(pos, decreasing = TRUE)) / sum(pos)
        thresh <- sort(pos, decreasing = TRUE)[which(cum >= keep)[1]]
        lam[lam < thresh] <- 0
      }
      kp <- e$vectors %*% (lam * t(e$vectors))
    } else {
      # clip: nearest PSD = (V sqrt(Lambda))(V sqrt(Lambda))^T by a symmetric
      # rank-k update; scale eigenvector column m by sqrt(lambda_m) (column-major,
      # so `each = n` aligns each scalar with its column).
      kp <- tcrossprod(e$vectors * rep(sqrt(lam), each = nrow(ks)))
    }
  }
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
#' `MaxEntropy()` selects the `k`-subset of points that maximises the
#' log-determinant of its kernel block, \eqn{\log\det K_S} -- the spanned volume
#' of the selection, the maximum-entropy sampling criterion
#' \insertCite{Shewry1987}{MaxMin} and the maximum-a-posteriori mode of a
#' determinantal point process \insertCite{Kulesza2012}{MaxMin}. A redundant
#' point lies in the span of those already chosen, adds zero volume, and is
#' never taken, so the selection is exactly density-blind.
#'
#' A radial-basis kernel \eqn{K_{ij} = \exp(-d_{ij}^2 / 2\sigma^2)} is built from
#' the supplied distances (`sigma` defaulting to the median positive distance)
#' and repaired to a positive-semidefinite matrix, because a general distance is
#' not of negative type and \eqn{\log\det} requires it. The exact argmax is
#' NP-hard \insertCite{Kulesza2012}{MaxMin}, so the selection is built greedily
#' by pivoted Cholesky -- at each step adding the point of largest residual
#' conditional variance -- with exact enumeration substituted where
#' \eqn{\binom{n}{k}} does not exceed `maxCombos`. The greedy first pivot is tied
#' on a unit-diagonal kernel and is broken deterministically by the most
#' peripheral point (least total similarity), so no random seed is used.
#'
#' @param k Integer target selection size, \eqn{1 \le k \le n}.
#' @param d A `dist` object or square numeric distance matrix over the `n`
#'   points.
#' @param sigma Optional kernel bandwidth; defaults to the median positive
#'   distance.
#' @param repair PSD repair for the kernel: `"clip"` (nearest PSD; the default),
#'   `"shift"` (diagonal loading) or `"truncate"` (low-rank embedding).
#' @param exact Logical, or `NA` (the default) to choose automatically: use the
#'   exact enumeration when `choose(n, k) <= maxCombos`, otherwise the greedy.
#'   `TRUE` forces enumeration (error if it exceeds `maxCombos`); `FALSE` forces
#'   the greedy.
#' @param maxCombos Numeric ceiling on `choose(n, k)` for exact enumeration.
#' @return `MaxEntropy()` returns an integer vector of length `k` (sorted
#'   ascending) with class `"MaxEntropySelection"`, carrying attributes:
#'   \describe{
#'     \item{logDet, score}{The retained \eqn{\log\det K_S} of the selection;
#'       `-Inf` for a degenerate selection (one forced to repeat near-identical
#'       points because `k` exceeds the number of distinct points, which also
#'       warns).}
#'     \item{negMass}{Fraction of spectral mass removed by the PSD repair.}
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
                       exact = NA, maxCombos = 3e5) {
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
  nDistinct <- sum(!duplicated(d))
  if (k > nDistinct) {
    warning(sprintf(paste0("`k` (%d) exceeds the number of distinct points (%d); ",
                           "the selection must repeat near-identical points and ",
                           "its log-determinant is -Inf (degenerate)."),
                    k, nDistinct))
  }

  kern <- .MaxEntropyKernel(d, sigma)
  sigmaUsed <- attr(kern, "sigma")
  prep <- .MaxEntropyPrepare(kern, repair)
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

  # Report the Cholesky log-determinant the exact selector maximises -- the same
  # quantity for both paths. When k exceeds the distinct-point count the selection
  # must repeat a point (pigeonhole), so its Gram matrix is singular and the true
  # log-determinant is -Inf; report that honestly rather than the large finite
  # value the clip-repaired kernel's softened duplicates would otherwise yield.
  logDet <- if (k > nDistinct) -Inf else MaxEntropyLogDet_cpp(kp, idx)

  # Return:
  structure(idx,
            logDet = logDet, score = logDet, negMass = negMass,
            sigma = sigmaUsed, repair = repair, exact = useExact,
            seed = as.integer(seed), N = n, k = k,
            producer = "MaxEntropy", class = "MaxEntropySelection")
}

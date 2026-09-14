# Tests for the MaxEntropy() kernel construction and its PSD repair: the C++
# kernel, the LAPACK partial eigendecomposition, and .MaxEntropyPrepare()
# against a full-eigen() reference.

# Deterministic fixtures. d^1.2 breaks the negative-type property of Euclidean
# distances, so the RBF kernel becomes indefinite: `dNeg` has 10 of 20 negative
# eigenvalues (the negative side is returned), `dPos` 31 of 40 (the positive
# side is the smaller one). `dPsd` is Euclidean and numerically positive-
# definite, so it takes the Cholesky fast path.
.Fixtures <- function() {
  set.seed(13)
  dPsd <- as.matrix(dist(matrix(rnorm(60), ncol = 3)))
  dNeg <- dPsd ^ 1.2
  dPos <- as.matrix(dist(matrix(rnorm(80), ncol = 2))) ^ 1.2
  list(psd = dPsd, neg = dNeg, pos = dPos)
}
.Sym <- function(k) {
  ks <- (k + t(k)) / 2
  attributes(ks) <- attributes(ks)["dim"]
  ks
}

# The kernel and repair as MaxEntropy() 1.0.0 computed them, in R.
.RefKernel <- function(d, sigma = NULL) {
  if (is.null(sigma)) {
    pos <- d[d > 0]
    sigma <- if (length(pos)) stats::median(pos) else 1
  }
  k <- exp(-(d ^ 2) / (2 * sigma ^ 2))
  attr(k, "sigma") <- sigma
  k
}
.RefPrepare <- function(k, method, keep = 0.99, tol = 1e-9) {
  ks <- .Sym(k)
  e <- eigen(ks, symmetric = TRUE)
  lam <- e$values
  neg <- lam[lam < -tol]
  negMass <- sum(abs(neg)) / sum(abs(lam))
  if (method == "shift") {
    lmin <- min(lam)
    kp <- if (lmin < 0) ks + (-lmin + tol) * diag(nrow(ks)) else ks
  } else {
    lam[lam < 0] <- 0
    if (method == "truncate") {
      pos <- lam[lam > 0]
      cum <- cumsum(sort(pos, decreasing = TRUE)) / sum(pos)
      thresh <- sort(pos, decreasing = TRUE)[which(cum >= keep)[1]]
      lam[lam < thresh] <- 0
    }
    kp <- e$vectors %*% (lam * t(e$vectors))
  }
  list(kp = kp, negMass = negMass)
}

test_that(".MaxEntropyKernel is bit-identical to its R formula", {
  for (s in 1:6) {
    set.seed(s)
    n <- 4L + s                                   # even and odd pair counts
    d <- as.matrix(dist(matrix(rnorm(n * 2), ncol = 2)))
    if (s %% 2 == 0) d[1, 2] <- d[2, 1] <- 0      # a zero pair is not "positive"
    expect_identical(.MaxEntropyKernel(d), .RefKernel(d))
    expect_identical(.MaxEntropyKernel(d, sigma = 0.7), .RefKernel(d, sigma = 0.7))
  }
  # No positive distance at all: sigma falls back to 1.
  z <- matrix(0, 3, 3)
  expect_identical(.MaxEntropyKernel(z), .RefKernel(z))
  expect_identical(attr(.MaxEntropyKernel(z), "sigma"), 1)
  expect_identical(.MaxEntropyKernel(matrix(0, 1, 1)), .RefKernel(matrix(0, 1, 1)))
  expect_true(is.na(MedianPositiveUpper_cpp(matrix(0, 0, 0))))
  expect_error(MedianPositiveUpper_cpp(matrix(1, 2, 3)), "square")
  expect_error(RbfKernel_cpp(matrix(1, 2, 3), 1), "square")
})

test_that("SymEigenPartial_cpp reproduces eigen() on each side of the spectrum", {
  fx <- .Fixtures()
  ks <- .Sym(.MaxEntropyKernel(fx$neg))
  n <- nrow(ks)
  vals <- rev(eigen(ks, symmetric = TRUE, only.values = TRUE)$values)
  m <- sum(vals < 0)
  expect_identical(m, 10L)

  v0 <- SymEigenPartial_cpp(ks, 0L, 0.99)          # values only
  expect_equal(v0$values, vals)
  expect_identical(v0$side, 0L)
  expect_identical(dim(v0$vectors), c(n, 0L))

  v1 <- SymEigenPartial_cpp(ks, 1L, 0.99)          # negative side (the minority)
  expect_identical(v1$side, -1L)
  expect_identical(dim(v1$vectors), c(n, m))
  expect_equal(v1$pvalues, vals[seq_len(m)])
  expect_equal(ks %*% v1$vectors, v1$vectors %*% diag(v1$pvalues, m))
  expect_equal(crossprod(v1$vectors), diag(m))

  # Majority-negative spectrum: the positive side is returned instead.
  ksPos <- .Sym(.MaxEntropyKernel(fx$pos))
  valsPos <- rev(eigen(ksPos, symmetric = TRUE, only.values = TRUE)$values)
  mPos <- sum(valsPos < 0)
  expect_gt(mPos, nrow(ksPos) / 2)
  v1p <- SymEigenPartial_cpp(ksPos, 1L, 0.99)
  expect_identical(v1p$side, 1L)
  expect_identical(ncol(v1p$vectors), nrow(ksPos) - mPos)
  expect_equal(v1p$pvalues, valsPos[valsPos >= 0])
  expect_equal(ksPos %*% v1p$vectors, v1p$vectors %*% diag(v1p$pvalues))

  # Positive-definite kernel: no negative eigenvalue, no vectors.
  ksPsd <- .Sym(.MaxEntropyKernel(fx$psd))
  vPsd <- SymEigenPartial_cpp(ksPsd, 1L, 0.99)
  expect_identical(vPsd$side, 0L)
  expect_identical(ncol(vPsd$vectors), 0L)
  expect_gt(min(vPsd$values), 0)

  # Truncate: top r eigenvalues holding `keep` of the positive mass, by the
  # same rule the R reference applies.
  pos <- vals[vals > 0]
  for (keep in c(0.5, 0.9, 0.99)) {
    r <- which(cumsum(sort(pos, decreasing = TRUE)) / sum(pos) >= keep)[1]
    v2 <- SymEigenPartial_cpp(ks, 2L, keep)
    expect_identical(v2$side, 1L)
    expect_identical(ncol(v2$vectors), r)
    expect_equal(v2$pvalues, vals[n - r + seq_len(r)])
    expect_equal(ks %*% v2$vectors, v2$vectors %*% diag(v2$pvalues, r))
  }
  expect_identical(ncol(SymEigenPartial_cpp(ks, 2L, 2)$vectors), length(pos))
  expect_identical(ncol(SymEigenPartial_cpp(ks, 2L, 0)$vectors), 1L)
  expect_identical(ncol(SymEigenPartial_cpp(-diag(3), 2L, 0.99)$vectors), 0L)

  # Edge cases.
  e0 <- SymEigenPartial_cpp(matrix(0, 0, 0), 1L, 0.99)
  expect_length(e0$values, 0L)
  # 1 x 1 negative: the empty positive side is the smaller one, so the clip
  # repair rebuilds from nothing (kp = 0).
  e1 <- SymEigenPartial_cpp(matrix(-2, 1, 1), 1L, 0.99)
  expect_equal(e1$values, -2)
  expect_identical(e1$side, 1L)
  expect_identical(dim(e1$vectors), c(1L, 0L))
  expect_identical(.MaxEntropyPrepare(matrix(-2, 1, 1), "clip")$kp, matrix(0, 1, 1))
  e2 <- SymEigenPartial_cpp(matrix(3, 1, 1), 1L, 0.99)
  expect_identical(e2$side, 0L)
  expect_equal(e2$values, 3)
  e3 <- SymEigenPartial_cpp(matrix(3, 1, 1), 2L, 0.99)
  expect_equal(e3$pvalues, 3)
  expect_equal(abs(e3$vectors), matrix(1, 1, 1))
  expect_error(SymEigenPartial_cpp(matrix(0, 2, 3), 0L, 0.99), "square")
  expect_error(SymEigenPartial_cpp(ks, 3L, 0.99), "mode")
})

test_that("RankUpdate_cpp is the explicit rank-p product", {
  set.seed(12)
  V <- matrix(rnorm(24), 8, 3)
  lam <- c(0.5, 2, 1e-3)
  base <- crossprod(matrix(rnorm(64), 8))
  expect_equal(RankUpdate_cpp(NULL, V, lam), V %*% diag(lam) %*% t(V))
  r <- RankUpdate_cpp(base, V, lam)
  expect_equal(r, base + V %*% diag(lam) %*% t(V))
  expect_identical(r, t(r))
  # A negative weight is roundoff on the wrong side of the cut: treated as 0.
  expect_equal(RankUpdate_cpp(NULL, V, c(0.5, -1, 1e-3)),
               V[, c(1, 3)] %*% diag(c(0.5, 1e-3)) %*% t(V[, c(1, 3)]))
  V0 <- V[, 0, drop = FALSE]
  expect_identical(RankUpdate_cpp(base, V0, numeric(0)), base)
  expect_identical(RankUpdate_cpp(NULL, V0, numeric(0)), matrix(0, 8, 8))
  expect_error(RankUpdate_cpp(NULL, V, c(1, 2)), "one entry per column")
  expect_error(RankUpdate_cpp(matrix(0, 3, 3), V, lam), "nrow")
})

test_that(".MaxEntropyPrepare matches the full-eigendecomposition reference", {
  fx <- .Fixtures()
  for (d in fx) {
    kern <- .MaxEntropyKernel(d)
    for (method in c("clip", "shift", "truncate")) {
      ref <- .RefPrepare(kern, method)
      new <- .MaxEntropyPrepare(kern, method)
      expect_equal(new$kp, ref$kp, tolerance = 1e-12)
      expect_equal(new$negMass, ref$negMass)
      expect_identical(new$kp, t(new$kp))
    }
    expect_equal(.MaxEntropyPrepare(kern, "truncate", keep = 0.6)$kp,
                 .RefPrepare(kern, "truncate", keep = 0.6)$kp, tolerance = 1e-12)
  }
  # The clip repair leaves a PSD kernel, whichever side it was built from.
  for (d in fx[c("neg", "pos")]) {
    kp <- .MaxEntropyPrepare(.MaxEntropyKernel(d), "clip")$kp
    expect_gte(min(eigen(kp, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
  }
  # Positive-definite kernel: returned as symmetrised, negMass exactly 0.
  kern <- .MaxEntropyKernel(fx$psd)
  for (method in c("clip", "shift")) {
    prep <- .MaxEntropyPrepare(kern, method)
    expect_identical(prep$kp, .Sym(kern))
    expect_identical(prep$negMass, 0)
  }
  # Positive-SEMIdefinite kernel (a zero pivot fails the Cholesky certificate,
  # yet no eigenvalue is negative): nothing to clip or shift.
  psd <- diag(c(1, 0))
  for (method in c("clip", "shift")) {
    prep <- .MaxEntropyPrepare(psd, method)
    expect_identical(prep$kp, psd)
    expect_identical(prep$negMass, 0)
  }
  expect_identical(.MaxEntropyPrepare(psd, "truncate")$kp, psd)
  # Back-compat wrappers.
  kern <- .MaxEntropyKernel(fx$neg)
  expect_identical(.MaxEntropyNegMass(kern), .MaxEntropyPrepare(kern, "clip")$negMass)
  expect_identical(.MaxEntropyRepair(kern, "shift"), .MaxEntropyPrepare(kern, "shift")$kp)
  expect_identical(.MaxEntropyRepair(kern), .MaxEntropyPrepare(kern, "clip")$kp)
})

test_that("MaxEntropy selects identically on PSD and indefinite kernels", {
  # End-to-end: the selection is what the 1.0.0 (full-eigen) pipeline picks,
  # for k within the numerical rank, on every repair.
  fx <- .Fixtures()
  for (d in fx) {
    n <- nrow(d)
    for (repair in c("clip", "shift", "truncate")) {
      kern <- .RefKernel(d)
      kp <- .RefPrepare(kern, repair)$kp
      seed <- which.min(rowSums(kp))
      for (k in c(2L, 4L)) {
        ref <- if (choose(n, k) <= 3e5) MaxEntropyExact_cpp(kp, k) else
          sort(MaxEntropyGreedy_cpp(kp, k, as.integer(seed)))
        sel <- MaxEntropy(k, d, repair = repair)
        expect_identical(as.integer(sel), as.integer(ref))
        expect_identical(attr(sel, "repair"), repair)
      }
    }
  }
  expect_gt(attr(MaxEntropy(3L, fx$neg), "negMass"), 0)
  expect_identical(attr(MaxEntropy(3L, fx$psd), "negMass"), 0)
})

test_that("the column-sweep greedy is bit-identical to the per-row form", {
  # The 1.0.0 greedy accumulated each row's dot product against the pivot row
  # in its own scalar loop; the column sweep keeps that summation order.
  greedyRows <- function(kp, k, seed) {
    n <- nrow(kp); dd <- diag(kp); L <- matrix(0, n, k)
    avail <- rep(TRUE, n); perm <- integer(k)
    for (t in seq_len(k)) {
      j <- if (t == 1L) seed else which.max(ifelse(avail, dd, -Inf))
      perm[t] <- j; avail[j] <- FALSE
      ljt <- sqrt(max(dd[j], 0)); L[j, t] <- ljt
      if (ljt > 0) for (row in which(avail)) {
        prev <- 0
        for (s in seq_len(t - 1L)) prev <- prev + L[row, s] * L[j, s]
        val <- (kp[row, j] - prev) / ljt
        L[row, t] <- val; dd[row] <- max(dd[row] - val * val, 0)
      }
    }
    perm
  }
  # k stays within the numerical rank of each repaired kernel (the residual
  # variance at every pick is >= 1e-3, three orders above round-off): past it
  # the picks are noise, and a compiler that contracts a * b + c into an FMA
  # (arm64) would legitimately break bit-identity with the interpreter there.
  fx <- .Fixtures()
  kMax <- c(psd = 10L, neg = 10L, pos = 9L)
  for (nm in names(kMax)) {
    kp <- .MaxEntropyPrepare(.MaxEntropyKernel(fx[[nm]]), "clip")$kp
    seed <- which.min(rowSums(kp))
    for (k in c(1L, 6L, kMax[[nm]])) {
      expect_identical(MaxEntropyGreedy_cpp(kp, k, as.integer(seed)),
                       as.integer(greedyRows(kp, k, seed)))
    }
  }
  expect_identical(MaxEntropyGreedy_cpp(diag(4), 6L, 2L), c(2L, 1L, 3L, 4L))  # k > n
  expect_length(MaxEntropyGreedy_cpp(diag(4), 0L, 1L), 0L)
})

test_that("DistinctRows_cpp counts rows as duplicated() does", {
  base <- rbind(c(0, 0), c(10, 0), c(0, 10))
  d <- as.matrix(dist(rbind(base, base, base[1, ])))          # 3 distinct of 7
  expect_identical(DistinctRows_cpp(d), sum(!duplicated(d)))
  expect_identical(DistinctRows_cpp(d), 3L)
  set.seed(2)
  m <- matrix(rnorm(40), 8)
  expect_identical(DistinctRows_cpp(m), 8L)
  m2 <- m[c(1:8, 3, 5, 5), ]
  expect_identical(DistinctRows_cpp(m2), sum(!duplicated(m2)))
  z <- matrix(c(0, -0, 0, -0), 2)                             # -0 == 0
  expect_identical(DistinctRows_cpp(z), sum(!duplicated(z)))
  expect_identical(DistinctRows_cpp(z), 1L)
  expect_identical(DistinctRows_cpp(matrix(0, 0, 0)), 0L)
  expect_identical(DistinctRows_cpp(matrix(0, 1, 1)), 1L)
  expect_identical(DistinctRows_cpp(matrix(c(1, 1, 2), 3)), 2L)   # single column
})

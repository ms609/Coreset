// symeigen.cpp
//
// LAPACK-backed partial symmetric eigendecomposition and rank-p reconstruction
// for the maximum-entropy kernel repair (R/maxentropy.R). Every repair needs
// the FULL spectrum (negMass, lambda_min, the truncation threshold) but only a
// SUBSET of the eigenvectors: those of the negative eigenvalues (clip -- or of
// the positive ones when they are the minority), or the top r (truncate).
// The matrix is tridiagonalised once (dsytrd, 4/3 n^3), all eigenvalues come
// from the tridiagonal form (dsterf, O(n^2)), and only the p wanted
// eigenvectors are computed (dstemr, O(n p)) and back-transformed (dormtr,
// 2 n^2 p); eigen() would back-transform all n (2 n^3).
//
// dsyevr -- what eigen(symmetric = TRUE) calls -- takes the same route for a
// full spectrum; for a subset it uses bisection + inverse iteration
// (dstebz/dstein), whose O(n p^2) reorthogonalisation of clustered eigenvalues
// costs as much as the full back-transformation once p ~ n/2. dstemr's MRRR
// subset mode keeps the vectors orthogonal without it.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif
#include <vector>
#include <memory>
#include <algorithm>
#include <cfloat>
#include <cmath>
using namespace Rcpp;

// dstemr (LAPACK >= 3.1) is not declared in R_ext/Lapack.h but is exported by
// R's bundled LAPACK and by any external LAPACK R can be built against.
extern "C" void F77_NAME(dstemr)(const char* jobz, const char* range,
                                 const int* n, double* d, double* e,
                                 const double* vl, const double* vu,
                                 const int* il, const int* iu, int* m,
                                 double* w, double* z, const int* ldz,
                                 const int* nzc, int* isuppz, int* tryrac,
                                 double* work, const int* lwork, int* iwork,
                                 const int* liwork, int* info FCLEN FCLEN);

// Recursive lower Cholesky factorisation of the n x n block at `a` (leading
// dimension lda), the algorithm of LAPACK's dpotrf2: factor the leading half,
// solve for the off-diagonal block (dtrsm), downdate the trailing half (dsyrk),
// factor it. Every update is a level-3 BLAS call on a block that halves at
// each level, which keeps the work cache-sized under reference BLAS and hands
// it whole to an optimised one; it is ~10% faster than the blocked dpotrf.
// Implemented here rather than called because dpotrf2 needs LAPACK >= 3.6,
// above the 3.2 floor R accepts for an external LAPACK, and a missing symbol
// would stop the package loading; dtrsm and dsyrk are in every BLAS. Only the
// lower triangle is read or written. Returns false at the first pivot that is
// not positive (or NaN), before any work to its right.
static bool CholeskyLower(double* a, int n, int lda) {
  if (n == 1) {
    if (!(a[0] > 0.0)) return false;
    a[0] = std::sqrt(a[0]);
    return true;
  }
  const int n1 = n / 2, n2 = n - n1;
  if (!CholeskyLower(a, n1, lda)) return false;
  const double one = 1.0, minusOne = -1.0;
  const char right = 'R', lower = 'L', trans = 'T', noTrans = 'N';
  double* a21 = a + n1;
  double* a22 = a21 + static_cast<size_t>(n1) * lda;
  F77_CALL(dtrsm)(&right, &lower, &trans, &noTrans, &n2, &n1, &one, a, &lda,
                  a21, &lda FCONE FCONE FCONE FCONE);
  F77_CALL(dsyrk)(&lower, &noTrans, &n2, &n1, &minusOne, a21, &lda, &one, a22,
                  &lda FCONE FCONE);
  return CholeskyLower(a22, n2, lda);
}

// Positive-definiteness certificate: TRUE when the Cholesky factorisation of
// A + delta I succeeds, delta = min(4 n eps ||A||_inf, tol). Every eigenvalue
// of A is then above -delta up to the factorisation's own backward error
// (~n eps ||A||), the resolution of any dense eigen-solver as well, so a
// negative one is round-off, below `tol`, and the clip would remove only
// round-off components (negMass is 0 by its definition); tol = 0 is the plain
// positive-definite test. A is taken as the symmetric matrix defined by its
// lower triangle, which is all that is copied, normed and factorised; the norm
// is accumulated column by column in the same pass. The pass/fail verdict is
// the only output; no factor bits are kept.
// [[Rcpp::export]]
bool CholCertificate_cpp(const NumericMatrix& A, double tol) {
  const int n = A.nrow();
  if (A.ncol() != n) stop("`A` must be square");
  if (n == 0) return true;
  const double* src = A.begin();
  std::unique_ptr<double[]> a(new double[static_cast<size_t>(n) * n]);
  std::vector<double> rowSum(n, 0.0);
  for (int j = 0; j < n; ++j) {
    const size_t off = static_cast<size_t>(j) * n;
    const double* sj = src + off;
    double* aj = a.get() + off;
    double below = 0.0;
    aj[j] = sj[j];
    for (int i = j + 1; i < n; ++i) {
      const double v = sj[i];
      aj[i] = v;
      const double av = std::fabs(v);
      rowSum[i] += av;
      below += av;
    }
    rowSum[j] += std::fabs(sj[j]) + below;
  }
  const double normInf = *std::max_element(rowSum.begin(), rowSum.end());
  const double delta = std::min(4.0 * n * DBL_EPSILON * normInf, tol);
  for (int i = 0; i < n; ++i) a[static_cast<size_t>(i) * n + i] += delta;
  return CholeskyLower(a.get(), n, n);
}

// Number of leading positive eigenvalues (of `vals`, ascending) needed to hold
// a fraction `keep` of the positive eigen-mass -- the truncate rule of
// .MaxEntropyPrepare(): the first r for which cumsum(top r) / sum(pos) >= keep.
// Accumulates in long double as R's sum()/cumsum() do, so r matches the R rule
// bit-for-bit. keep > 1 keeps every positive eigenvalue; keep <= 0 keeps one.
static int TopByMass(const NumericVector& vals, double keep) {
  const int n = vals.size();
  int nPos = 0;
  long double total = 0.0L;
  for (int i = n - 1; i >= 0 && vals[i] > 0.0; --i) { ++nPos; total += vals[i]; }
  if (nPos == 0) return 0;
  const double sumPos = static_cast<double>(total);
  long double cum = 0.0L;
  for (int r = 1; r <= nPos; ++r) {
    cum += vals[n - r];
    if (static_cast<double>(cum) / sumPos >= keep) return r;
  }
  return nPos;
}

// Eigenvalues of a symmetric matrix, plus a chosen subset of its eigenvectors.
//   mode 0: eigenvalues only.
//   mode 1: eigenvectors of the smaller side of zero -- the negative
//           eigenvalues if there are at most as many as non-negative ones
//           (side = -1), else the positive ones (side = +1). side = 0 (no
//           vectors) when no eigenvalue is negative.
//   mode 2: eigenvectors of the top r eigenvalues holding `keep` of the
//           positive eigen-mass (side = +1; r may be 0) -- or, when those are
//           the majority (r > n / 2), of the n - r eigenvalues the truncation
//           DROPS (side = -1): the caller subtracts them instead.
// Returns list(values = all n eigenvalues ascending, side, pvalues = the
// eigenvalues of the returned vectors (ascending), vectors = n x p matrix).
// Only the lower triangle of A is read: pass a symmetric matrix.
// [[Rcpp::export]]
List SymEigenPartial_cpp(const NumericMatrix& A, int mode, double keep) {
  const int n = A.nrow();
  if (A.ncol() != n) stop("`A` must be square");
  if (mode < 0 || mode > 2) stop("`mode` must be 0, 1 or 2");
  NumericVector values(n);
  if (n == 0) {
    return List::create(_["values"] = values, _["side"] = 0,
                        _["pvalues"] = NumericVector(0),
                        _["vectors"] = NumericMatrix(0, 0));
  }
  std::vector<double> a(A.begin(), A.end());        // LAPACK overwrites its input
  std::vector<double> d(n), e(n), tau(n);
  int info = 0, lwork = -1;
  double wq = 0.0;
  const char uplo = 'L';
  F77_CALL(dsytrd)(&uplo, &n, a.data(), &n, d.data(), e.data(), tau.data(),
                   &wq, &lwork, &info FCONE);
  lwork = static_cast<int>(wq);
  std::vector<double> work(std::max(lwork, 1));
  F77_CALL(dsytrd)(&uplo, &n, a.data(), &n, d.data(), e.data(), tau.data(),
                   work.data(), &lwork, &info FCONE);
  if (info != 0) stop("LAPACK dsytrd failed (info = %d)", info);  // # nocov
  std::copy(d.begin(), d.end(), values.begin());
  std::vector<double> e2(e);
  F77_CALL(dsterf)(&n, values.begin(), e2.data(), &info);
  if (info != 0) stop("LAPACK dsterf failed (info = %d)", info);  // # nocov

  int il = 1, iu = 0, side = 0;
  if (mode == 1) {
    int m = 0;
    for (int i = 0; i < n; ++i) if (values[i] < 0.0) ++m;
    if (m > 0) {
      if (m <= n - m) { il = 1; iu = m; side = -1; }
      else            { il = m + 1; iu = n; side = 1; }
    }
  } else if (mode == 2) {
    const int r = TopByMass(values, keep);
    if (r <= n - r) { il = n - r + 1; iu = n; side = 1; }
    else            { il = 1; iu = n - r; side = -1; }
  }
  const int p = iu - il + 1;
  if (p <= 0) {
    return List::create(_["values"] = values, _["side"] = side,
                        _["pvalues"] = NumericVector(0),
                        _["vectors"] = NumericMatrix(n, 0));
  }

  // Eigenpairs il..iu of the tridiagonal form (MRRR), then back-transform.
  const char jobz = 'V', range = 'I';
  const double vl = 0.0, vu = 0.0;
  int mOut = 0, nzc = p, tryrac = 1, liwork = -1, iwq = 0;
  std::vector<double> w(n);
  std::vector<int> isuppz(2 * p);
  NumericMatrix Z(n, p);
  lwork = -1;
  F77_CALL(dstemr)(&jobz, &range, &n, d.data(), e.data(), &vl, &vu, &il, &iu,
                   &mOut, w.data(), Z.begin(), &n, &nzc, isuppz.data(), &tryrac,
                   &wq, &lwork, &iwq, &liwork, &info FCONE FCONE);
  lwork = static_cast<int>(wq);
  liwork = iwq;
  std::vector<double> work2(std::max(lwork, 1));
  std::vector<int> iwork(std::max(liwork, 1));
  F77_CALL(dstemr)(&jobz, &range, &n, d.data(), e.data(), &vl, &vu, &il, &iu,
                   &mOut, w.data(), Z.begin(), &n, &nzc, isuppz.data(), &tryrac,
                   work2.data(), &lwork, iwork.data(), &liwork, &info FCONE FCONE);
  if (info != 0 || mOut != p) {                                     // # nocov start
    stop("LAPACK dstemr failed (info = %d, %d of %d eigenpairs)", info, mOut, p);
  }                                                                 // # nocov end
  const char sideL = 'L', trans = 'N';
  lwork = -1;
  F77_CALL(dormtr)(&sideL, &uplo, &trans, &n, &p, a.data(), &n, tau.data(),
                   Z.begin(), &n, &wq, &lwork, &info FCONE FCONE FCONE);
  lwork = static_cast<int>(wq);
  std::vector<double> work3(std::max(lwork, 1));
  F77_CALL(dormtr)(&sideL, &uplo, &trans, &n, &p, a.data(), &n, tau.data(),
                   Z.begin(), &n, work3.data(), &lwork, &info FCONE FCONE FCONE);
  if (info != 0) stop("LAPACK dormtr failed (info = %d)", info);  // # nocov
  NumericVector pv(p);
  std::copy(w.begin(), w.begin() + p, pv.begin());
  return List::create(_["values"] = values, _["side"] = side,
                      _["pvalues"] = pv, _["vectors"] = Z);
}

// base + V diag(lam) V^T by symmetric rank-p updates (dsyrk) on a copy of
// `base` -- or on zeros when base is NULL -- mirrored to full storage so the
// result is symmetric to the bit. The columns with lam >= 0 go in one
// (added) update and those with lam < 0 in a second (subtracted) one, so a
// signed `lam` -- the eigenvalues a truncation drops, of either sign -- is
// handled without squaring away its signs.
// [[Rcpp::export]]
NumericMatrix RankUpdate_cpp(Nullable<NumericMatrix> base, const NumericMatrix& V,
                             const NumericVector& lam) {
  const int n = V.nrow(), p = V.ncol();
  if (lam.size() != p) stop("`lam` must have one entry per column of `V`");
  NumericMatrix out(n, n);
  if (base.isNotNull()) {
    NumericMatrix b(base);
    if (b.nrow() != n || b.ncol() != n) stop("`base` must be nrow(V) x nrow(V)");
    std::copy(b.begin(), b.end(), out.begin());
  }
  if (p == 0 || n == 0) return out;
  std::vector<double> Vpos, Vneg;
  Vpos.reserve(static_cast<size_t>(n) * p);
  for (int j = 0; j < p; ++j) {
    std::vector<double>& dst = lam[j] < 0.0 ? Vneg : Vpos;
    const double s = std::sqrt(std::fabs(lam[j]));
    const double* col = &V(0, j);
    for (int i = 0; i < n; ++i) dst.push_back(col[i] * s);
  }
  const char uplo = 'L', trans = 'N';
  double beta = base.isNotNull() ? 1.0 : 0.0;
  const int pPos = static_cast<int>(Vpos.size() / n), pNeg = p - pPos;
  if (pPos > 0) {
    const double alpha = 1.0;
    F77_CALL(dsyrk)(&uplo, &trans, &n, &pPos, &alpha, Vpos.data(), &n, &beta,
                    out.begin(), &n FCONE FCONE);
    beta = 1.0;
  }
  if (pNeg > 0) {
    const double alpha = -1.0;
    F77_CALL(dsyrk)(&uplo, &trans, &n, &pNeg, &alpha, Vneg.data(), &n, &beta,
                    out.begin(), &n FCONE FCONE);
  }
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) out(j, i) = out(i, j);
  }
  return out;
}

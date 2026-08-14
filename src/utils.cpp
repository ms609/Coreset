#include <Rcpp.h>
#include <cstdint>
#include <cstring>

// TRUE iff every element of a double or integer vector/matrix is finite --
// no NA, NaN, or +/-Inf. Replaces `anyNA(x) || any(!is.finite(x))` in the
// R-side input guards: that idiom allocates two full-size logical
// intermediates (~2 * 4 bytes per cell), which for a 6000 x 6000 distance
// matrix costs ~200 ms -- 20x the greedy kernel it guards. This scan
// allocates nothing.
//
// A double is non-finite exactly when its biased exponent bits are all ones
// (NA_real_ is a NaN payload, so the same test catches it). The bitwise
// form is branch-free and integer-only, so the loop vectorises without any
// floating-point reassociation concern, and the OR-reduction is
// order-independent -- the parallel result is identical at every thread
// count by construction.
//
// [[Rcpp::export]]
bool AllFinite_cpp(SEXP x, int n_threads = 1) {
  if (TYPEOF(x) == REALSXP) {
    const double* p = REAL(x);
    const R_xlen_t n = XLENGTH(x);
    const uint64_t expMask = 0x7FF0000000000000ULL;
    uint64_t bad = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(|:bad) \
    if(n_threads > 1 && n >= 1048576) num_threads(n_threads) schedule(static)
#endif
    for (R_xlen_t i = 0; i < n; i++) {
      uint64_t bits;
      std::memcpy(&bits, &p[i], sizeof bits);
      bad |= (uint64_t)((bits & expMask) == expMask);
    }
    return bad == 0;
  } else if (TYPEOF(x) == INTSXP) {
    // Integer storage: the only non-finite value is NA_integer_
    // (is.finite() is TRUE for every other int), so match that exactly.
    const int* p = INTEGER(x);
    const R_xlen_t n = XLENGTH(x);
    for (R_xlen_t i = 0; i < n; i++) {
      if (p[i] == NA_INTEGER) return false;
    }
    return true;
  }
  Rcpp::stop("`x` must have double or integer storage");
}

// How far from symmetric is the square matrix `d`? Returns the largest
// |d(i, j) - d(j, i)| over the strict upper triangle, scaled by
// max(1, |d(i, j)|, |d(j, i)|) so the same tolerance serves distances of any
// magnitude; exactly 0 iff the two triangles hold identical bits, and
// infinite if `d` is not square.
//
// Zero and merely-small are separate verdicts to the caller because the
// solvers read whichever of d(i, j) and d(j, i) is the cheaper memory access,
// and which one that is depends on the subset size: a matrix symmetric only
// to rounding would answer differently for different `k` unless the two
// triangles are reconciled first.
//
// Compared tile by tile, so a tile and its transpose are both resident when
// their entries meet; the naive i-outer sweep walks one side of each pair
// with stride n and misses on nearly every read. NaN compares false against
// everything and so leaves the running maximum untouched -- the callers scan
// for finiteness first, and report that instead.
//
// [[Rcpp::export]]
double SymmetryDeviation_cpp(Rcpp::NumericMatrix d) {
  const int n = d.nrow();
  if (n != d.ncol()) {
    return R_PosInf;
  }
  const double* p = REAL(d);
  const int tile = 64;
  double worst = 0;
  for (int j0 = 0; j0 < n; j0 += tile) {
    const int j1 = std::min(j0 + tile, n);
    for (int i0 = j0; i0 < n; i0 += tile) {
      const int i1 = std::min(i0 + tile, n);
      for (int j = j0; j < j1; ++j) {
        const int iFrom = (i0 > j + 1) ? i0 : j + 1;
        for (int i = iFrom; i < i1; ++i) {
          const double a = p[(std::size_t)i + (std::size_t)j * n];
          const double b = p[(std::size_t)j + (std::size_t)i * n];
          if (a != b) {
            const double scale = std::max(1.0, std::max(std::fabs(a),
                                                        std::fabs(b)));
            const double dev = std::fabs(a - b) / scale;
            if (dev > worst) {
              worst = dev;
            }
          }
        }
      }
    }
  }
  return worst;
}

// `d` with each pair of off-diagonal entries replaced by their mean --
// elementwise identical to `(d + t(d)) / 2`, but from one allocation rather
// than three, which at n = 6435 is the difference between 0.3 GB and 1 GB.
// The mean of two equal doubles is exactly that double, so an exactly
// symmetric `d` is returned bit-for-bit.
//
// [[Rcpp::export]]
Rcpp::NumericMatrix Symmetrised_cpp(Rcpp::NumericMatrix d) {
  const int n = d.nrow();
  Rcpp::NumericMatrix out = Rcpp::clone(d);
  const double* p = REAL(d);
  double* q = REAL(out);
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) {
      const std::size_t ij = (std::size_t)i + (std::size_t)j * n;
      const std::size_t ji = (std::size_t)j + (std::size_t)i * n;
      const double mean = (p[ij] + p[ji]) / 2;
      q[ij] = mean;
      q[ji] = mean;
    }
  }
  return out;
}

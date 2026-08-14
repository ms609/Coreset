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

// Is the square matrix `d` exactly symmetric? Compared tile by tile, so a
// tile and its transpose are both resident when their entries meet; the naive
// i-outer sweep walks one side of each pair with stride n and misses on
// nearly every read.
//
// Exact equality, not a tolerance: the solvers read whichever of d(i, j) and
// d(j, i) is the cheaper memory access, and which one that is depends on the
// subset size. A matrix symmetric only to rounding would then answer
// differently for different `k`, so it is rejected rather than silently
// resolved one way.
//
// [[Rcpp::export]]
bool IsExactlySymmetric_cpp(Rcpp::NumericMatrix d) {
  const int n = d.nrow();
  if (n != d.ncol()) {
    return false;
  }
  const double* p = REAL(d);
  const int tile = 64;
  for (int j0 = 0; j0 < n; j0 += tile) {
    const int j1 = std::min(j0 + tile, n);
    for (int i0 = j0; i0 < n; i0 += tile) {
      const int i1 = std::min(i0 + tile, n);
      for (int j = j0; j < j1; ++j) {
        const int iFrom = (i0 > j + 1) ? i0 : j + 1;
        for (int i = iFrom; i < i1; ++i) {
          if (p[(std::size_t)i + (std::size_t)j * n] !=
              p[(std::size_t)j + (std::size_t)i * n]) {
            return false;
          }
        }
      }
    }
  }
  return true;
}

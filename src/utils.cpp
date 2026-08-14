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

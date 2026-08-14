# utils.R

# NULL-coalescing operator: `x` if not NULL, else (lazily) `y`.
# @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

# Worker-thread count for the OpenMP kernels, from the standard `"mc.cores"`
# option. Anything unusable collapses to the default of 1.
# @noRd
.NThreads <- function() {
  n <- suppressWarnings(as.integer(getOption("mc.cores", 1L)))
  if (length(n) != 1L || is.na(n) || n < 1L) 1L else n
}

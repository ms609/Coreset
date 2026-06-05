# utils.R

# NULL-coalescing operator: `x` if not NULL, else (lazily) `y`.
# @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

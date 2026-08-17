# DropAdd tabu search from a distance-column oracle

The pure-R counterpart of `DropAdd_cpp()`, substituting an on-demand
`colFn(i)` call for the matrix-column read `dmat[, i]`.

## Usage

``` r
.DropAddFromColumn(
  colFn,
  N,
  k,
  first,
  plateau = 5000L,
  maxSeconds = Inf,
  maxIter = .Machine$integer.max,
  trace = FALSE
)
```

## Arguments

- colFn:

  Column oracle; see
  [`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md).

- N:

  Integer element count.

- k:

  Integer subset size (`2 <= k <= N`).

- first:

  Integer 1-based index of the seed element.

- plateau:

  Integer: stop after this many consecutive non-improving iterations.

- maxSeconds:

  Numeric wall-clock budget. Checked once per iteration (an iteration is
  dominated by two oracle calls, so
  [`proc.time()`](https://rdrr.io/r/base/proc.time.html) is free by
  comparison and the 1024-iteration stride of the C++ kernel is
  unnecessary); the check is skipped entirely when the budget is
  infinite, keeping the default path deterministic.

- maxIter:

  Integer cap on main-loop iterations.

- trace:

  Logical: also return the dropped/added index sequences, for the
  trajectory-identity tests.

## Value

`.DropAddFromColumn()` returns a list with the same shape as
`DropAdd_cpp()`'s: `indices` (1-based, in FIFO buffer order),
`objective`, `secondary`, `iters`, and – when `trace` is `TRUE` –
`drops` and `adds`.

# Gonzalez maximin from a distance-column oracle

Implements the distance-column oracle path of
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
(dispatched there when `d` is a function); see that function's
*Distance-column oracle* section for the user-facing contract. At each
greedy step the distances from the newly selected element to all `N`
elements are obtained from `colFn`, and a running nearest-distance
vector is maintained, so the `N x N` distance matrix is never
materialised: `O(N * n)` oracle calls and `O(N)` memory.

## Usage

``` r
.GonzalezColumn(
  colFn,
  N,
  n,
  first = NULL,
  progress = getOption("MaxMin.progress", interactive())
)
```

## Arguments

- colFn:

  A function that, when passed an index `i`, must return a vector of
  distances from element `i` to either (i) every element in turn,
  including `i`; or (ii) every other element. See
  [`.DistColumn()`](https://ms609.github.io/MaxMin/reference/dot-DistColumn.md).

- N:

  Integer: the total number of elements.

- n:

  Integer: number of elements to select. If `n > N`, all `N` indices are
  returned in Gonzalez (farthest-first) order.

- first:

  Integer index of the first selected element, or `NULL` (default) to
  use a deterministic peripheral seed computed from two oracle sweeps:
  the element furthest from element 1, then the element furthest from
  that (a diameter-endpoint approximation).

- progress:

  Logical; show a progress bar during greedy selection.

## Value

Integer vector of length `min(n, N)` of selected indices.

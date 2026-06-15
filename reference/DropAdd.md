# DropAdd Tabu Search for the Max-Min Diversity Problem

`DropAdd()` selects a maximally-dispersed subset of `k` points using the
DropAdd tabu search algorithm, which comprises a greedy construction
followed by a first-in, first-out drop-add tabu search, with streamlined
neighbour-evaluation tricks (algorithms 1–4 in Porumbel et al. 2011) .

## Usage

``` r
DropAdd(k, d = NULL, plateau = 5000L, maxSeconds = Inf, points = NULL)
```

## Arguments

- k:

  Integer; subset size, \\2 \le k \le N\\.

- d:

  A `dist` object or square symmetric numeric matrix.

- plateau:

  Integer; stop after this many consecutive drop-add iterations do not
  improve the score.

- maxSeconds:

  Numeric: terminate search after this many seconds have elapsed.

- points:

  A numeric \\N \times \mathrm{dim}\\ coordinate matrix (or an object
  coercible to one via `as.matrix`). Must be complete (no `NA`). Ignored
  if `d` specified. Avoids creating an \\N \times N\\ distance matrix,
  enabling use at \\N \ge 46340\\).

## Value

`DropAdd()` returns an integer vector of length `k` containing the
1-based selected indices **sorted ascending** (unlike
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md),
which returns farthest-first order), with attributes:

- score:

  numeric(1), achieved MaxMin objective \\\min\_{i \ne j \in S}
  d\_{ij}\\.

- secondary:

  numeric(1), achieved sum of pairwise distances over \\S\\
  (upper-triangle sum).

- time_s:

  numeric(1), wall-clock seconds spent.

- iters:

  integer(1), main-loop iterations executed (excluding the construction
  phase).

The vector has class `"MaxMinSelection"` and prints as a one-line
summary (see
[`print.MaxMinSelection()`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md));
it is otherwise an ordinary integer vector.

## Progress bar

In interactive sessions, status messages are shown. To toggle, set
`options("MaxMin.progress" = FALSE)` (or `TRUE`).

## References

Porumbel D, Hao J, Glover F (2011). “A simple and effective algorithm
for the MaxMin diversity problem.” *Annals of Operations Research*,
**186**, 275–293.
[doi:10.1007/s10479-011-0898-z](https://doi.org/10.1007/s10479-011-0898-z)
.

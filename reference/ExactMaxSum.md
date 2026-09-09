# Exact Maximum Diversity Problem (max-sum) solution

`ExactMaxSum()` finds the optimal solution to the Max-Sum Diversity
Problem (the "maximum diversity problem"): it selects the `k` points
that maximizes the total pairwise distance between points. As the
problem is NP-hard it is feasible only for small sets.

## Usage

``` r
ExactMaxSum(k, d, maxSeconds = 30, warmStart = NULL)
```

## Arguments

- k:

  Integer: target subset size, between 2 and `nrow(d)`.

- d:

  `dist` object or a square symmetric numeric distance matrix.

- maxSeconds:

  Numeric: search terminates after this many seconds have elapsed,
  returning largest threshold proven feasible.

- warmStart:

  Optional integer vector giving indices of a candidate subset to add to
  the heuristic warm-start pool.

## Value

`ExactMaxSum()` returns an integer vector of length `k`, sorted
ascending, with class `"MaxSumSelection"`, carrying attributes:

- score:

  Achieved total pairwise distance within the selection. When `proven`
  is `TRUE` this is the optimum; otherwise a lower bound.

- proven:

  Logical: `TRUE` if optimality was certified.

- seconds, N, k:

  Wall-clock seconds elapsed; instance size; target size.

## Details

The solver uses per-node integer-program linearisation (Kuo et al. 1993)
, starting from a multi-start 1-swap local search, whose result is
returned when optimality cannot be proven within `maxSeconds`.

## References

Kuo C, Glover F, Dhir KS (1993). “Analyzing and modeling the maximum
diversity problem by zero-one programming.” *Decision Sciences*,
**24**(6), 1171–1185.
[doi:10.1111/j.1540-5915.1993.tb00509.x](https://doi.org/10.1111/j.1540-5915.1993.tb00509.x)
.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(20), ncol = 2)
ExactMaxSum(3L, dist(pts))
#> 3 elements (1 3 4) by exact MILP, proven optimal, total distance = 9.388
```

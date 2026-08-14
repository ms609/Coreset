# Constructive phase of DropAdd from a distance-column oracle

Algorithm 1 of Porumbel et al. (2011) against a column oracle: the
counterpart of `.DropAddConstruct()`that never materialises the distance
matrix.

## Usage

``` r
.DropAddConstructColumn(colFn, N, k, first)
```

## Arguments

- colFn:

  Column oracle; see
  [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md).

- N:

  Integer element count.

- k:

  Integer subset size (`2 <= k <= N`).

- first:

  Integer 1-based index of the seed element.

## Value

`.DropAddConstructColumn()` returns a `new.env(parent = emptyenv())`
holding `S` (the selection, in add order), `inS`, `minDist`, `sumDist`
and `minDistCount`, matching the record set of `.DropAddConstruct()`.

## Seed deviation

The matrix kernel seeds at the max-row-sum point, \\\mathrm{argmax}\_x
\sum_y d(x, y)\\. Pass `DropAdd(seed=)` to override.

## References

Porumbel D, Hao J, Glover F (2011). “A simple and effective algorithm
for the MaxMin diversity problem.” *Annals of Operations Research*,
**186**, 275–293.
[doi:10.1007/s10479-011-0898-z](https://doi.org/10.1007/s10479-011-0898-z)
.

# Choose the next element to add to a DropAdd selection

The `argmax` over the unselected elements of `(minDist, sumDist)` taken
lexicographically, ties broken to the smallest index.

## Usage

``` r
.DropAddPick(st, exclude = 0L)
```

## Arguments

- st:

  Record environment; see
  [`.DropAddConstructColumn()`](https://ms609.github.io/MaxMin/reference/dot-DropAddConstructColumn.md).

- exclude:

  Integer index barred from selection this iteration (`0L` for none).
  Used for the just-dropped `x#`, which Porumbel et al. (2011, p. 281)
  exclude from `Add X(k)` for one iteration – the tabu rule.

## Value

`.DropAddPick()` returns a single integer index.

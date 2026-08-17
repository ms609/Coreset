# Max-Mean Dispersion Problem solver

`MaxMean()` selects a maximally dispersed subset of elements from a
pairwise distance matrix, maximising the *max-mean* objective: \$\$f(S)
= \frac{\displaystyle\sum\_{i \< j,\\ i,j \in S} d\_{ij}}{\|S\|}\$\$ The
number of elements in the subset \\\|S\| \ge 2\\ is chosen so as to
maximise the mean dispersion.

## Usage

``` r
MaxMean(d, maxSeconds = 0.1, maxIter = 1000, useRL = TRUE)
```

## Arguments

- d:

  A `dist` object or square numeric matrix of pairwise distances; values
  may be negative, and an asymmetric matrix is symmetrized to
  \\(d\_{ij} + d\_{ji})/2\\ before solving.

- maxSeconds:

  Numeric: wall-clock time budget, in seconds.

- maxIter:

  Numeric: cap on the total tabu-search iterations across restarts.

- useRL:

  Logical: if `TRUE` (the default) a \\Q\\-learning layer guides
  initial-solution construction across restarts; if `FALSE`, each
  restart starts from a random solution.

## Value

`MaxMean()` returns an integer vector of selected 1-based indices
(sorted ascending) with class `"MaxMeanSelection"` and attributes:

- score:

  numeric, achieved objective \\\sum\_{i\<j \in S} d\_{ij} / \|S\|\\.

- size:

  integer, number of selected elements \\\|S\|\\.

- time_s:

  numeric, wall-clock seconds spent.

- iters:

  numeric, total tabu-search iterations across restarts. Stored as a
  double, not an integer, because a long run can exceed the 32-bit
  integer range.

The vector has class `"MaxMeanSelection"` and prints as a one-line
summary (see
[`print.MaxMeanSelection()`](https://ms609.github.io/Coreset/reference/print.MaxMeanSelection.md));
it is otherwise an ordinary integer vector that indexes the distance
matrix directly.

## Details

`MaxMean()` implements the reinforcement-learning tabu search algorithm
of (Nijimbere et al. 2020) . An initial solution is constructed randomly
for the first restart and via \\Q\\-learning thereafter; each initial
solution is then refined by a tabu search using one-flip moves (adding
or removing one element per step). Restarts continue until either the
`maxSeconds` or `maxIter` budget is reached.

The reinforcement-learning and tabu hyperparameters are fixed at the
tuned values reported by (Nijimbere et al. 2020) (greedy factor
\\\epsilon = 0.7\\, learning rate \\\alpha = 0.5\\, discount \\\gamma =
0.5\\, maximum tabu tenure \\120\\, search depth 50 000).

## Progress bar

In interactive sessions, status messages are shown. To toggle, set
`options("Coreset.progress" = FALSE)` (or `TRUE`).

## References

Nijimbere D, Zhao S, Gu X, Esangbedo MO, Dominique N (2020). “Tabu
search guided by reinforcement learning for the max-mean dispersion
problem.” *Journal of Industrial & Management Optimization*, **17**,
3223–3254.
[doi:10.3934/jimo.2020115](https://doi.org/10.3934/jimo.2020115) .

## See also

[`MeanDist()`](https://ms609.github.io/Coreset/reference/MeanDist.md) to
score an arbitrary selection under this objective;
[`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md),
[`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md) and
[`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md) for
fixed-cardinality max-min solvers.

## Examples

``` r
# The max-mean problem is defined for signed dissimilarities; with these the
# optimal subset has an interior size, chosen to maximise mean dispersion.
set.seed(1)
x <- matrix(runif(100, -5, 5), 10)
d <- (x + t(x)) / 2          # symmetric, signed
selection <- MaxMean(d)
selection                    # 5 of the 10 elements: {2, 5, 6, 8, 10}
#> 5 elements (2 5 6 8 10) selected by MaxMean RLTS, f = 2.714
MeanDist(d, selection)       # 2.714 — equals attr(selection, "score")
#> [1] 2.714473
```

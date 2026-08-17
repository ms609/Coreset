# Near-optimal discrete k-centre solver

`KCentre()` selects \\k\\ elements (centres) so as to minimize the
largest distance from any point to its nearest centre (the covering
radius), using the Critical Dominating Set heuristic (CDSh) (García-Díaz
et al. 2017; García-Díaz et al. 2019) .

## Usage

``` r
KCentre(k, d, nstart = 1L, effort = 1L)

KCenter(k, d, nstart = 1L, effort = 1L)
```

## Arguments

- k:

  Integer specifying maximum number of centres to identify, from 1 to
  `nrow(d)`.

- d:

  `dist` object or a square symmetric numeric distance matrix.

- nstart:

  Integer specifying how many deterministic peripheral seeds to try.

- effort:

  Integer: if `> 0`, run a parallel
  [`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md)
  search with `effort` random seeds, returning the best of all results.

## Value

`KCentre()` returns an integer vector of length \\\le k\\ specifying the
chosen centres in ascending order. The achieved covering radius is
attached as attribute `radius`. The vector has class
`"KCentreSelection"` and prints as a one-line summary.

## Details

On the benchmark instances of García-Díaz et al. (2019) , the CDS
heuristic reaches roughly 1-3.5% of the optimum at \\O(N^2 \log N)\\,
far tighter than
[`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md)
(typically tens of per cent above optimum).

Despite this good performance in practice, the CDSh is a
3-approximation. To guard against occasional cases where a better
candidate is missed, `KCentre()` runs by default an exhaustive search of
a small candidate grid (for `n` up to ~150); and an additional
[`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md)
pass (controlled via the `effort` argument). These safeguards ensure
that `KCentre()` always returns at least a 2-approximation.

## References

García-Díaz J, Menchaca-Méndez R, Menchaca-Méndez R, Pomares Hernández
S, Pérez-Sansalvador JC, Lakouari N (2019). “Approximation algorithms
for the vertex \\k\\-center problem: survey and experimental
evaluation.” *IEEE Access*, **7**, 109228–109245.
[doi:10.1109/ACCESS.2019.2933875](https://doi.org/10.1109/ACCESS.2019.2933875)
.  
  
García-Díaz J, Sánchez-Hernández J, Menchaca-Méndez R, Menchaca-Méndez R
(2017). “When a worse approximation factor gives better performance: a
3-approximation algorithm for the vertex \\k\\-center problem.” *Journal
of Heuristics*, **23**(5), 349–366.
[doi:10.1007/s10732-017-9345-x](https://doi.org/10.1007/s10732-017-9345-x)
.

## See also

[`ExactKCentre()`](https://ms609.github.io/Coreset/reference/ExactKCentre.md)
for the proven optimum;
[`KCentreRadius()`](https://ms609.github.io/Coreset/reference/KCentreRadius.md)
for a selection's score;
[`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md)
for the González (1985) 2-approximation baseline.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(120), ncol = 2)
d <- dist(pts)
centres <- KCentre(5L, d)
KCentreRadius(d, centres)
#> [1] 1.22535
# Results will beat a Gonzalez 2-approximation:
KCentreRadius(d, FarFirst(5L, d, nSeeds = 1))
#> [1] 1.361541
```

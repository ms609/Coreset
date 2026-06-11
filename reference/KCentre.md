# Near-optimal discrete k-centre by CDSh

Chooses `k` centres minimising the covering radius (the largest distance
from any point to its nearest centre) with the CDSh heuristic of
García-Díaz et al. (2019) (see also (García-Díaz et al. 2017) ). CDSh
binary-searches the achieved distinct distances; at each trial radius it
runs a fixed-`k` farthest-point construction in which every centre is
the highest-degree neighbour (within the trial radius) of the currently
worst-covered vertex, and accepts the radius when the construction
covers all points within it. On the benchmark instances of García-Díaz
et al. (2019) it reaches roughly 1-3.5% of the optimum at \\O(N^2 \log
N)\\, far tighter than the Gonzalez 2-approximation that
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
gives for this objective (typically tens of per cent above optimum).

## Usage

``` r
KCentre(d, k, nstart = 1L, effort = 1L, seeds = NULL)

KCenter(d, k, nstart = 1L, effort = 1L, seeds = NULL)
```

## Arguments

- d:

  A `dist` object or a square symmetric numeric distance matrix.

- k:

  Integer number of centres, `1 <= k <= nrow(d)`.

- nstart:

  Integer; how many deterministic peripheral seeds to try, keeping the
  lowest-radius result. Default `1`. Ignored if `seeds` is supplied.

- effort:

  Integer; controls the Gonzalez floor that keeps the result no worse
  than the 2-approximation. `0` disables it (raw CDSh, fastest, no
  guarantee); `1` (default) runs one deterministic peripheral Gonzalez
  pass; `> 1` runs a distinct-seed Gonzalez restart with
  `nseeds = effort` (a tighter floor that draws on the session RNG –
  call [`set.seed()`](https://rdrr.io/r/base/Random.html) to reproduce).

- seeds:

  Optional integer vector of explicit 1-based seed vertices for the
  construction's first critical vertex (overrides `nstart`).

## Value

Integer vector of length `<= k` (ascending): the chosen centres. The
achieved covering radius is attached as attribute `radius`. The vector
has class `"KCentreSelection"` and prints as a one-line summary; it is
otherwise an ordinary integer vector and indexes a matrix or coordinate
set directly.

## Details

The achieved covering radius is not monotone in the trial radius, so the
binary search can occasionally miss the best candidate. Two safeguards
keep the result robust: for a small candidate grid (`n` up to ~150)
every radius is scanned exhaustively, and the result is floored against
a Gonzalez pass (the `effort` argument, on by default), so `KCentre()`
is **never worse than the 2-approximation**. For a tighter result at
larger `n` raise `nstart` or `effort`; for the proven optimum on a small
instance use
[`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md).

The construction is otherwise fully deterministic: where the reference
seeds its first critical vertex at random, `KCentre()` uses
deterministic peripheral anchors (`nstart` of them, the best kept). Like
[`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
this is a distance-matrix method, \\O(N^2)\\ in memory; for the covering
radius of an existing selection at larger `N`,
[`KCentreRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
has a matrix-free path.

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

[`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
for the proven optimum;
[`KCentreRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
to score a selection;
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) for
the Gonzalez 2-approximation baseline.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(120), ncol = 2)
d <- dist(pts)
centres <- KCentre(d, 5L)
KCentreRadius(d, centres)
#> [1] 1.22535
# CDSh covers at least as tightly as the Gonzalez 2-approximation:
KCentreRadius(d, FarFirst(d, 5L, method = "peripheral"))
#> [1] 1.361541
```

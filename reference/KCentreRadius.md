# Covering radius of a centre set (k-centre objective)

Returns the covering radius of a set of centres: the largest distance
from any of the `N` points to its nearest centre, \\R = \max_p \min\_{c
\in \mathrm{idx}} d(p, c)\\. This is the min-max *k*-centre objective
(González 1985) , the quantity
[`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md) and
[`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
minimise. Lower is better.

## Usage

``` r
KCentreRadius(d = NULL, idx, points = NULL)

KCenterRadius(d = NULL, idx, points = NULL)
```

## Arguments

- d:

  Pairwise distance matrix or `dist` object. Ignored when `points` is
  supplied.

- idx:

  Integer vector of centre indices (`>= 1`).

- points:

  Optional `N x dim` numeric coordinate matrix. When supplied the
  per-point nearest-centre distances are recomputed from coordinates one
  centre column at a time, never materialising the `N x N` matrix (`d`
  is then unused), so the score scales to large `N`. For Euclidean data
  the result is identical to the matrix path.

## Value

Numeric scalar: the covering radius (`0` when the centres include every
point).

## Details

Unlike
[`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md) – the
minimum *pairwise* distance *within* a selection (the MMDP objective) –
the covering radius is taken over *all* `N` points and measures how well
the centres cover the data.

## References

González TF (1985). “Clustering to minimize the maximum intercluster
distance.” *Theoretical Computer Science*, **38**, 293–306.
[doi:10.1016/0304-3975(85)90224-5](https://doi.org/10.1016/0304-3975%2885%2990224-5)
.

## See also

[`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md) and
[`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
(which minimise this);
[`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md) for
the complementary MMDP objective.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
d <- dist(pts)
centres <- KCentre(d, 4L)
KCentreRadius(d, centres)
#> [1] 1.242483
```

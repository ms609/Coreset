# Covering radius of a set of centres

`KCentreRadius()` computes the covering radius of a set of centres: the
largest distance from any of the `N` points to its nearest centre, \\R =
\max_p \min\_{c \in \mathrm{idx}} d(p, c)\\. This is the min-max
*k*-centre objective (González 1985) minimized by
[`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md) and
[`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md).

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

  `N x dim` numeric coordinate matrix. When supplied, the per-point
  nearest-centre distances are computed from coordinates one column at a
  time. This supports sets with large `N`, in which the full `N x N`
  matrix would overflow available memory.

## Value

`KCentreRadius()` returns a numeric denoting the covering radius.

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
centres <- KCentre(4L, d)
KCentreRadius(d, centres)
#> [1] 1.242483
```

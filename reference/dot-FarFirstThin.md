# Run a solver on a farthest-first coreset and map indices back

Implements the composable-coreset path of
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) /
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md)
(dispatched there when their `maxCandidates` cap binds). It builds an
`m`-point coreset with
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md),
restricts the problem to those `m` points, runs the supplied solver on
the restriction, and maps the returned indices back to the original
numbering.

## Usage

``` r
.FarFirstThin(k, m, d = NULL, points = NULL, RunOnSubset, label)
```

## Arguments

- k:

  Integer: target subset size.

- m:

  Integer: coreset size (`k <= m < N`, as returned by
  [`.ResolveCap()`](https://ms609.github.io/MaxMin/reference/dot-ResolveCap.md)).

- d:

  Square distance matrix of the full problem, or `NULL` on the
  coordinate path.

- points:

  `N x dim` coordinate matrix of the full problem, or `NULL` on the
  distance-matrix path. Exactly one of `d` / `points` is non-`NULL`.

- RunOnSubset:

  Function `function(d, points)` that runs the downstream solver on the
  restricted problem and returns its `MaxMinSelection`. It is called
  with the restricted distance matrix (`d = d[core, core]`) on the
  matrix path, or the restricted coordinates (`points = points[core, ]`)
  on the coordinate path; the unused argument is `NULL`.

- label:

  Character naming the calling solver, used in the thinning warning
  (e.g. `"DropAdd"`).

## Value

`.FarFirstThin()` returns the solver's `MaxMinSelection` with its
indices mapped to original-space row indices (sorted ascending).

## Details

Only the integer index values returned by the solver need remapping. The
result is sorted ascending.

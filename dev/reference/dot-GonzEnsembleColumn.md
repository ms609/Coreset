# Seed ensemble from a distance-column oracle

The distance-column counterpart of
[`.GonzEnsemble()`](https://ms609.github.io/Coreset/dev/reference/dot-GonzEnsemble.md),
over the anchors an oracle can reach: start indices, the peripheral
seed, and `nSeeds` distinct random-furthest seeds, drawn exactly as the
matrix path draws them. Each distinct seed runs one
[`.GonzalezColumn()`](https://ms609.github.io/Coreset/dev/reference/dot-GonzalezColumn.md)
pass, serially.

## Usage

``` r
.GonzEnsembleColumn(colFn, N, k, anchors, nSeeds, progress = FALSE)
```

## Arguments

- colFn:

  A function that, when passed an index `i`, must return a vector of
  distances from element `i` to either (i) every element in turn,
  including `i`; or (ii) every other element. See
  [`.DistColumn()`](https://ms609.github.io/Coreset/dev/reference/dot-DistColumn.md).

- N:

  Integer: the total number of elements.

- k:

  Integer: number of elements to select. If `k > N`, all `N` indices are
  returned in Gonzalez (farthest-first) order.

- anchors:

  Any of `"random_furthest"`, `"peripheral"` and start indices as
  strings (see
  [`.NormaliseStrategy()`](https://ms609.github.io/Coreset/dev/reference/dot-NormaliseStrategy.md)).

- nSeeds:

  Integer number of distinct random-furthest seeds.

- progress:

  Logical; show a progress bar during greedy selection.

## Value

`.GonzEnsembleColumn()` returns an integer vector of selected indices
with the attributes of
[`.ResolveEnsemble()`](https://ms609.github.io/Coreset/dev/reference/dot-ResolveEnsemble.md).

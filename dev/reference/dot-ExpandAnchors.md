# Expand ensemble anchor names into labelled seed specs

Maps each anchor name to a `list(label, s1)` spec. The
`"random_furthest"` token expands to one spec per element of `rfSeeds`
(already-resolved seed indices, labelled `random_furthest1`, ...); an
empty vector contributes none.

## Usage

``` r
.ExpandAnchors(anchors, rfSeeds, anchorSeed)
```

## Arguments

- anchors:

  Character vector of (de-duplicated) anchor names.

- rfSeeds:

  Integer vector of already-resolved furthest-point seed indices for the
  `"random_furthest"` token.

- anchorSeed:

  Function mapping a deterministic anchor name to an integer seed index.

## Value

`.ExpandAnchors()` returns a list of `list(label, s1)` specs.

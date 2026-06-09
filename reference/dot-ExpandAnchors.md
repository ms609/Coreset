# Expand ensemble anchor names into labelled seed specs

Maps each anchor name to a `list(label, s1)` spec. The
`"random_furthest"` token expands to one spec per element of `pivots`,
each seeded at the point furthest from that pivot (labelled
`random_furthest1`, ...); an empty `pivots` contributes none.

## Usage

``` r
.ExpandAnchors(anchors, pivots, anchor_seed, rf_seed)
```

## Arguments

- anchors:

  Character vector of (de-duplicated) anchor names.

- pivots:

  Integer vector of pivot indices the `"random_furthest"` token expands
  over (one start per pivot).

- anchor_seed:

  Function mapping a deterministic anchor name to an integer seed index.

- rf_seed:

  Function mapping a pivot index to the furthest-point seed.

## Value

List of `list(label, s1)` specs.

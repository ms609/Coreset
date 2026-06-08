# Expand ensemble anchor names into labelled seed specs

Maps each anchor name to a `list(label, s1, mask)` spec. The
`"random_furthest"` token expands to `n_random` specs, each seeded at
the point furthest from one of `n_random` reproducible random pivots
(labelled `random_furthest1`, ...); with `n_random == 0` it contributes
none.

## Usage

``` r
.ExpandAnchors(anchors, n_random, N, anchor_seed, rf_seed)
```

## Arguments

- anchors:

  Character vector of (de-duplicated) anchor names.

- n_random:

  Integer count the `"random_furthest"` token expands to.

- N:

  Integer element count, for the random pivot draw.

- anchor_seed:

  Function mapping a deterministic anchor name to `list(s1, mask)`.

- rf_seed:

  Function mapping a pivot index to the furthest-point seed.

## Value

List of `list(label, s1, mask)` specs.

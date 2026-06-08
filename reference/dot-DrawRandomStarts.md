# Reproducible random start indices, isolated from ambient RNG

Draws `min(k, N)` distinct indices in `[1, N]` from a fixed internal
seed, the random pivots for the `"random_furthest"` seed. The caller's
RNG kind and `.Random.seed` are saved and restored, so the draw neither
depends on nor perturbs the ambient random stream; the result is
therefore identical across sessions and machines. Restoring the kind
scrambles `.Random.seed`, so the kind is restored before the seed.

## Usage

``` r
.DrawRandomStarts(N, k, seed = .kRandomSeed)
```

## Arguments

- N:

  Integer element count (`>= 1`).

- k:

  Integer number of pivots to draw.

- seed:

  Integer RNG seed; defaults to the package's fixed value.

## Value

Integer vector of `min(k, N)` distinct indices.

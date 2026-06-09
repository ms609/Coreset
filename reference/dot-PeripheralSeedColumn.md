# Deterministic peripheral seed from a column oracle

Two oracle sweeps, no RNG: the element furthest from element 1, then the
element furthest from that. The second is a diameter-endpoint
approximation and a markedly better Gonzalez anchor than an arbitrary
start, at the cost of two of the `O(N * n)` sweeps. The richer
peripheral anchors (diameter, anti-medoid) need `O(N^2)` work and are
unreachable from a column oracle.

## Usage

``` r
.PeripheralSeedColumn(colFn, N)
```

## Arguments

- colFn:

  Column oracle; see
  [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md).

- N:

  Integer element count.

## Value

Integer index of the seed.

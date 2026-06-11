# Stamp the `MaxMinSelection` class onto a solver's index vector

The selection-returning solvers each already attach their score (and any
secondary attributes) via
[`base::structure()`](https://rdrr.io/r/base/structure.html); this adds
the shared S3 class and a `producer` tag the print method reads to name
the algorithm. An empty selection (`length 0`) is left bare: there is
nothing to describe.

## Usage

``` r
.AsMaxMinSelection(x, producer)
```

## Arguments

- x:

  Integer index vector carrying the solver's score attributes.

- producer:

  Character tag naming the solver (`"FarFirst"`, `"DropAdd"`,
  `"Grasp"`).

## Value

`x` with `producer` attribute and `"MaxMinSelection"` class, or `x`
unchanged if it is empty.

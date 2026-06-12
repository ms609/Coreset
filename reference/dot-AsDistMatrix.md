# Coerce distance input to a square matrix, skipping the round-trip when already a matrix.

Coerce distance input to a square matrix, skipping the round-trip when
already a matrix.

## Usage

``` r
.AsDistMatrix(d)
```

## Arguments

- d:

  A `dist` object or a square numeric matrix.

## Value

`.AsDistMatrix()` returns a square numeric matrix.

## Details

Symmetry is not checked; an `O(N^2)` check is intentionally omitted.
Asymmetric matrices are silently accepted, and the algorithm treats
\\d\_{ij}\\ and \\d\_{ji}\\ as independent values.

# Coerce distance input to a square matrix, skipping the round-trip when already a matrix.

Coerce distance input to a square matrix, skipping the round-trip when
already a matrix.

## Usage

``` r
.AsDistMatrix(d, symmetric = TRUE)
```

## Arguments

- d:

  A `dist` object or a square numeric matrix.

- symmetric:

  Logical: reconcile the two triangles of `d`.

## Value

`.AsDistMatrix()` returns a square numeric matrix.

## Details

A `dist` object is symmetric by construction, so it bypasses the check.

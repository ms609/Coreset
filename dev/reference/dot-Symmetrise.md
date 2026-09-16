# Reconcile the two triangles of a distance matrix

Internal helper: returns `d` unchanged when its triangles already agree,
averages them when they differ by no more than `tolerance`, and errors
beyond it.

## Usage

``` r
.Symmetrise(d, dev, tolerance = .SymmetryTolerance())
```

## Arguments

- d:

  A square numeric matrix, known finite.

- dev:

  Numeric: its scaled asymmetry, from `SymmetryScan_cpp()`.

- tolerance:

  Numeric: largest scaled discrepancy to repair.

## Value

`.Symmetrise()` returns an exactly symmetric matrix.

## Details

Averaging costs one `n * n` copy per call: repair `d` yourself if
calling a solver in a loop.

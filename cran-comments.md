## Test environments

* local Windows 10, R-devel
* GitHub Actions: Windows / macOS / Ubuntu (release, devel, oldrel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the CRAN team

* This is a new submission.

* The package name `Coreset` is close in spelling to the existing CRAN package
  `corset` (Arbitrary Bounding of Series and Time Series Objects). The two are
  unrelated: `corset` constrains numerical and time series within given
  boundaries, whereas `Coreset` selects a representative subset of a fixed
  candidate set under an arbitrary distance. The name is the established term
  in computational geometry and machine learning for a small subset that
  preserves a property of the whole, which is what the package computes.

* `Coreset` is likewise distinct from the CRAN package `maximin` (Sun &
  Gramacy), which constructs continuous space-filling designs by generating new
  points in a coordinate region; `Coreset` chooses among points that already
  exist. The package documentation cross-references `maximin` to avoid
  confusion.

* The `highs` package is used only by `ExactKCentre()` and is listed in
  Suggests behind a `requireNamespace()` guard; all examples and tests that
  need it skip when it is absent.

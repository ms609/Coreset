## Test environments

* local Windows 10, R-devel
* GitHub Actions: Windows / macOS / Ubuntu (release, devel, oldrel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the CRAN team

* The package name `MaxMin` is distinct from the existing CRAN package
  `maximin` (Sun & Gramacy): the names differ in spelling, and the two solve
  different problems. `maximin` constructs continuous space-filling designs;
  `MaxMin` selects a subset from a fixed candidate set under an arbitrary
  distance (the discrete p-dispersion / Max-Min Diversity Problem). The README
  and package documentation cross-reference `maximin` to avoid confusion.

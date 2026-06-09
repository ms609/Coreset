# MaxMin: Maximum-Minimum Diversity and Dispersion Subset Selection

Selects a maximally dispersed subset of a fixed candidate set under the
Max-Min Diversity Problem (MMDP, the discrete *p*-dispersion objective):
maximise the minimum pairwise distance within the chosen subset. The
solvers operate on a distance matrix, on Euclidean coordinates (without
materialising the matrix), or on an on-demand distance-column oracle.

## Solvers

- [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md):

  Greedy farthest-first selection from a distance matrix, a coordinate
  matrix, or a distance-column oracle (for spaces with no coordinate
  embedding), with a choice of peripheral seeding strategies and a
  robust ensemble default.

- [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md):

  DropAdd tabu search heuristic.

- [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md):

  Exact node-packing optimum (needs highs).

- [`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md):

  The k-centre objective (minimum pairwise distance).

## Relation to maximin

Not to be confused with the CRAN package maximin (Sun & Gramacy), which
constructs continuous *space-filling designs* — it generates new points
in a coordinate region to maximise the minimum inter-point distance.
`MaxMin` instead *selects a subset* from a *fixed* candidate set under
an arbitrary distance, a combinatorial problem on a different footing.

## See also

Useful links:

- <https://ms609.github.com/MaxMin/>

- Report bugs at <https://github.com/ms609/MaxMin/issues>

## Author

**Maintainer**: Martin R. Smith <martin.smith@durham.ac.uk>
([ORCID](https://orcid.org/0000-0001-5660-1727)) \[copyright holder\]

Authors:

- Martin R. Smith <martin.smith@durham.ac.uk>
  ([ORCID](https://orcid.org/0000-0001-5660-1727)) \[copyright holder\]

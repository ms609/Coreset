# Package index

## Greedy selection

Farthest-first (Gonzalez) selection and its seeding strategies.

- [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) :
  Deterministic Gonzalez furthest-point selection
- [`MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/MaxMinSeed.md)
  : Peripheral seed index for Gonzalez farthest-first selection

## Heuristic and exact solvers

- [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) :
  DropAdd Tabu Search for the Max-Min Diversity Problem
- [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) : GRASP
  with Path Relinking for the Max-Min Diversity Problem
- [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
  : Exact Max-Min Diversity (MMDP) optimum on small instances

## Scoring

- [`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md) :
  Minimum pairwise distance within a selection (T_k = k-centre
  objective)

## Reporting

- [`format(`*`<MaxMinSelection>`*`)`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md)
  [`print(`*`<MaxMinSelection>`*`)`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md)
  [`format(`*`<MaxMinExact>`*`)`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md)
  [`print(`*`<MaxMinExact>`*`)`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md)
  : Format and print MaxMin solver results
- [`summary(`*`<MaxMinSelection>`*`)`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)
  [`summary(`*`<MaxMinExact>`*`)`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)
  : Multi-line summaries of MaxMin solver results

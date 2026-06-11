# Package index

## Greedy selection

Farthest-first (Gonzalez) selection and its seeding strategies.

- [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) :
  Deterministic Gonzalez furthest-point selection
- [`MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/MaxMinSeed.md)
  : Peripheral seed index for Gonzalez farthest-first selection

## MaxMin solvers

- [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) :
  DropAdd Tabu Search for the Max-Min Diversity Problem
- [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) : GRASP
  with Path Relinking for the Max-Min Diversity Problem
- [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
  : Exact Max-Min Diversity (MMDP) optimum on small instances

## k-centre solvers

- [`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md)
  [`KCenter()`](https://ms609.github.io/MaxMin/reference/KCentre.md) :
  Near-optimal discrete k-centre by CDSh
- [`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
  [`ExactKCenter()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
  : Exact discrete k-centre optimum on small instances

## Scoring

- [`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md) :
  Minimum pairwise distance within a selection
- [`KCentreRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
  [`KCenterRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
  : Covering radius of a centre set (k-centre objective)

## Reporting

- [`format(`*`<KCentreSelection>`*`)`](https://ms609.github.io/MaxMin/reference/print.KCentre.md)
  [`print(`*`<KCentreSelection>`*`)`](https://ms609.github.io/MaxMin/reference/print.KCentre.md)
  [`format(`*`<KCentreExact>`*`)`](https://ms609.github.io/MaxMin/reference/print.KCentre.md)
  [`print(`*`<KCentreExact>`*`)`](https://ms609.github.io/MaxMin/reference/print.KCentre.md)
  : Format and print k-centre solver results
- [`format(`*`<MaxMinSelection>`*`)`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md)
  [`print(`*`<MaxMinSelection>`*`)`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md)
  [`format(`*`<MaxMinExact>`*`)`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md)
  [`print(`*`<MaxMinExact>`*`)`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md)
  : Format and print MaxMin solver results
- [`summary(`*`<MaxMinSelection>`*`)`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)
  [`summary(`*`<MaxMinExact>`*`)`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)
  : Multi-line summaries of MaxMin solver results

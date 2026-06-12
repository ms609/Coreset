# Package index

## Greedy selection

Farthest-first (Gonzalez) selection and its seeding strategies.

- [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) :
  Deterministic Gonzalez furthest-point selection
- [`MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/MaxMinSeed.md)
  : Seed to initialize farthest-first selection

## Max-Min diversity problem solvers

- [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) :
  DropAdd Tabu Search for the Max-Min Diversity Problem
- [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) : GRASP
  with Path Relinking for the Max-Min Diversity Problem
- [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
  : Exact Max-Min Diversity Problem optimum on small instances

## k-centre solvers

- [`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md)
  [`KCenter()`](https://ms609.github.io/MaxMin/reference/KCentre.md) :
  Near-optimal discrete k-centre solver
- [`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
  [`ExactKCenter()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
  : Exact discrete k-centre optimum on small instances

## Scoring

- [`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md) :
  Minimum pairwise distance within a selection
- [`KCentreRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
  [`KCenterRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
  : Covering radius of a set of selected centres

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
  : Detailed summaries of MaxMin solver results

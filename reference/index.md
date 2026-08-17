# Package index

## Greedy selection

Farthest-first (Gonzalez) selection and its seeding strategies.

- [`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md)
  : Greedy farthest-first point selection
- [`PickPoint()`](https://ms609.github.io/Coreset/reference/PickPoint.md)
  : Seed to initialize farthest-first selection

## Max-Min diversity problem solvers

- [`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md) :
  DropAdd Tabu Search for the Max-Min Diversity Problem
- [`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md) :
  GRASP with Path Relinking for the Max-Min Diversity Problem
- [`ExactMaxMin()`](https://ms609.github.io/Coreset/reference/ExactMaxMin.md)
  : Exact Max-Min Diversity Problem solution

## Max-Sum diversity problem solver

- [`ExactMaxSum()`](https://ms609.github.io/Coreset/reference/ExactMaxSum.md)
  : Exact Maximum Diversity Problem (max-sum) solution

## Max-Mean dispersion problem solver

- [`MaxMean()`](https://ms609.github.io/Coreset/reference/MaxMean.md) :
  Max-Mean Dispersion Problem solver

## k-centre solvers

- [`KCentre()`](https://ms609.github.io/Coreset/reference/KCentre.md)
  [`KCenter()`](https://ms609.github.io/Coreset/reference/KCentre.md) :
  Near-optimal discrete k-centre solver
- [`ExactKCentre()`](https://ms609.github.io/Coreset/reference/ExactKCentre.md)
  [`ExactKCenter()`](https://ms609.github.io/Coreset/reference/ExactKCentre.md)
  : Exact discrete k-centre optimum on small instances

## Maximum-entropy (maxdet) selection

- [`MaxEntropy()`](https://ms609.github.io/Coreset/reference/MaxEntropy.md)
  : Maximum-entropy (maxdet) subset selection

## Scoring

- [`MinDist()`](https://ms609.github.io/Coreset/reference/MinDist.md) :
  Minimum pairwise distance within a selection
- [`MeanDist()`](https://ms609.github.io/Coreset/reference/MeanDist.md)
  : Mean dispersion of a selection
- [`KCentreRadius()`](https://ms609.github.io/Coreset/reference/KCentreRadius.md)
  [`KCenterRadius()`](https://ms609.github.io/Coreset/reference/KCentreRadius.md)
  : Covering radius of a set of centres

## Reporting

- [`format(`*`<MaxMinSelection>`*`)`](https://ms609.github.io/Coreset/reference/print.Coreset.md)
  [`print(`*`<MaxMinSelection>`*`)`](https://ms609.github.io/Coreset/reference/print.Coreset.md)
  : Format and print Coreset solver results
- [`format(`*`<KCentreSelection>`*`)`](https://ms609.github.io/Coreset/reference/print.KCentre.md)
  [`print(`*`<KCentreSelection>`*`)`](https://ms609.github.io/Coreset/reference/print.KCentre.md)
  [`format(`*`<KCentreExact>`*`)`](https://ms609.github.io/Coreset/reference/print.KCentre.md)
  [`print(`*`<KCentreExact>`*`)`](https://ms609.github.io/Coreset/reference/print.KCentre.md)
  : Format and print k-centre solver results
- [`format(`*`<MaxEntropySelection>`*`)`](https://ms609.github.io/Coreset/reference/print.MaxEntropy.md)
  [`print(`*`<MaxEntropySelection>`*`)`](https://ms609.github.io/Coreset/reference/print.MaxEntropy.md)
  : Format and print maximum-entropy (maxdet) solver results
- [`format(`*`<MaxMeanSelection>`*`)`](https://ms609.github.io/Coreset/reference/print.MaxMeanSelection.md)
  [`print(`*`<MaxMeanSelection>`*`)`](https://ms609.github.io/Coreset/reference/print.MaxMeanSelection.md)
  [`summary(`*`<MaxMeanSelection>`*`)`](https://ms609.github.io/Coreset/reference/print.MaxMeanSelection.md)
  : Format and print Max-Mean solver results
- [`format(`*`<MaxSumSelection>`*`)`](https://ms609.github.io/Coreset/reference/print.MaxSum.md)
  [`print(`*`<MaxSumSelection>`*`)`](https://ms609.github.io/Coreset/reference/print.MaxSum.md)
  : Format and print Max-Sum (maximum diversity) solver results
- [`summary(`*`<MaxMinSelection>`*`)`](https://ms609.github.io/Coreset/reference/summary.Coreset.md)
  : Detailed summaries of Coreset solver results

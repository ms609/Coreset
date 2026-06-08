# Package index

## Greedy selection

Farthest-first (Gonzalez) selection and its seeding strategies.

- [`Gonzalez()`](https://ms609.github.io/MaxMin/reference/Gonzalez.md) :
  Deterministic Gonzalez furthest-point selection
- [`MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/MaxMinSeed.md)
  : Peripheral seed index for Gonzalez farthest-first selection

## Heuristic and exact solvers

- [`DropAddTS()`](https://ms609.github.io/MaxMin/reference/DropAddTS.md)
  : DropAdd Tabu Search for the Max-Min Diversity Problem
- [`DropAddTSPoints()`](https://ms609.github.io/MaxMin/reference/DropAddTSPoints.md)
  : Matrix-free DropAdd Tabu Search for the Max-Min Diversity Problem
- [`GraspPR()`](https://ms609.github.io/MaxMin/reference/GraspPR.md) :
  GRASP with Path Relinking for the Max-Min Diversity Problem
- [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
  : Exact Max-Min Diversity (MMDP) optimum on small instances

## Refinement and scoring

- [`PolishSelection()`](https://ms609.github.io/MaxMin/reference/PolishSelection.md)
  : Local-search polish for a max-min diversity selection
- [`TkScore()`](https://ms609.github.io/MaxMin/reference/TkScore.md) :
  Minimum pairwise distance within a selection (T_k = k-centre
  objective)

# Name the seeding outcome of a [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) selection

A bare single pass is just "Gonzalez farthest-first"; an ensemble pass
(which carries `winning_strategy` / `strategy_results`) additionally
names the winning strategy and how many were tried.

## Usage

``` r
.FarFirstSelectedBy(x)
```

## Arguments

- x:

  A `MaxMinSelection` from
  [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md).

## Value

Length-1 character phrase.

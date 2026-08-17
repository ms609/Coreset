# Promote the kernel's free `t_k` score to the user-facing `score` attribute

The maximin kernels attach the selection's minimum pairwise distance as
a `t_k` attribute (computed during the greedy pass at no extra cost).
The ensemble drivers read it via
[`base::attr()`](https://rdrr.io/r/base/attr.html) into
`strategy_results`; a bare single pass exposes it directly as the
`score` attribute, matching
[`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md) and
[`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md).

## Usage

``` r
.PromoteScore(idx)
```

## Arguments

- idx:

  Integer vector returned by a maximin kernel.

## Value

`.PromoteScore()` returns `idx` with its `t_k` attribute renamed to
`score`.

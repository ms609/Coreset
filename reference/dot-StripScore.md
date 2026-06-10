# Read and detach the kernel's free `t_k` score

The maximin kernels attach the selection's minimum pairwise distance as
a `t_k` attribute (computed during the greedy pass at no extra cost).
The ensemble drivers read it via
[`base::attr()`](https://rdrr.io/r/base/attr.html); bare single passes
strip it with `.StripScore()` so the returned indices carry no
incidental attribute.

## Usage

``` r
.StripScore(idx)
```

## Arguments

- idx:

  Integer vector returned by a maximin kernel.

## Value

`idx` with its `t_k` attribute removed.

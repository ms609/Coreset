# Validate and normalise a `maxCandidates` thinning cap

Shared by
[`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md) and
[`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md). Decides
whether candidate thinning is active and, if so, the intermediate
coreset size `m` to thin to.

## Usage

``` r
.ResolveCap(maxCandidates, n, k)
```

## Arguments

- maxCandidates:

  The user-supplied cap: a positive integer (thin to this many
  candidates when it is below `n`), or `0` / `Inf` to disable thinning.

- n:

  Integer: the number of candidate points in the full problem.

- k:

  Integer: the target subset size.

## Value

`.ResolveCap()` returns `NA_integer_` when thinning is disabled (`0` /
`Inf`) or non-binding (`maxCandidates >= n`); otherwise the integer
coreset size `m` (`k <= m < n`). Errors on a non-integer, negative,
`NA`, or non-scalar cap, or a positive cap below `k`.

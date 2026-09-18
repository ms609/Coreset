# One randomised greedy construction.

One randomised greedy construction.

## Usage

``` r
.GraspConstruct(d, k, alpha, us)
```

## Arguments

- d:

  Square distance matrix.

- k:

  Target subset size.

- alpha:

  RCL threshold parameter (alpha=1 -\> greedy, alpha=0 -\> random).

- us:

  Numeric vector of `k` pre-drawn uniforms. `us[[1]]` picks the first
  vertex and `us[[h]]` step `h`'s RCL member, as `floor(u * size)`
  clamped to `size - 1` — exactly the kernel's mapping, so R and C++
  stay bit-identical. Steps whose RCL is a singleton (or empty; see
  below) leave their uniform unused: consumption is fixed at `k` draws
  either way, which is what lets whole batches be drawn up front and
  constructions run on threads that must not touch R's RNG.

## Value

`.GraspConstruct()` returns an integer vector of length k.

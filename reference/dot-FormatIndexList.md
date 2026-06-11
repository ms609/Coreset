# Format a selected-index list, optionally truncated

Format a selected-index list, optionally truncated

## Usage

``` r
.FormatIndexList(idx, maxShow = 20L)
```

## Arguments

- idx:

  Integer indices in their stored order.

- maxShow:

  Integer: show at most this many before eliding the tail.

## Value

Length-1 character string such as `"6 5 4 3 1 2"`, or
`"1 2 ... (+15 more)"` when `idx` is longer than `maxShow`.

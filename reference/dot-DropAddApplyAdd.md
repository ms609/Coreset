# Fold a newly added element into the streamlined DropAdd records

The ADD pass of Algorithm 3, shared by the construction and the tabu
loop and mirroring the add-pass blocks of `src/dropadd.cpp`. Mutates
`st` in place (the package forbids `<<-`; `st` is a
`new.env(parent = emptyenv())` record bundle). `xNew`'s *own* record is
set by the caller, which knows which peers are current.

## Usage

``` r
.DropAddApplyAdd(st, col, xNew)
```

## Arguments

- st:

  Record environment; see
  [`.DropAddConstructColumn()`](https://ms609.github.io/Coreset/reference/dot-DropAddConstructColumn.md).

- col:

  Self-zeroed distance column of `xNew`; see
  [`.DropAddColumn()`](https://ms609.github.io/Coreset/reference/dot-DropAddColumn.md).

- xNew:

  Integer index of the element just added.

## Value

`.DropAddApplyAdd()` returns `NULL` invisibly, for its effect on `st`.

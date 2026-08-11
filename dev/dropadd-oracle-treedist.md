# DropAdd's distance-column oracle against TreeDist: what the numbers say

Scratch note (not a vignette). Written alongside the `DropAdd(d = <function>)`
oracle path, to settle what claim the user-facing docs are entitled to make:
**faster**, or **lower memory**?

Answer: **lower memory, and a lower distance-evaluation count only against a
naive baseline.** On wall-clock, a bulk all-pairs call plus the compiled matrix
kernel beats every column-oracle variant measured, at every size measured. The
docs say so.

Benchmark: `dev/bench-dropadd-oracle.R`. Metric: `TreeDist::ClusteringInfoDistance`
on random 60-tip trees. MaxMin 0.0.0.9005 installed with
`R CMD INSTALL --no-multiarch --preclean`, run in a fresh `R --vanilla` session
(AGENTS.md: `load_all()` compiles the C++ at `-O0` and would flatter the pure-R
oracle path). Windows 11, R 4.7, TreeDist 2.14.0.

## 1. Does TreeDist expose a columnar primitive?

Yes, and a good one. This was the open question in the briefing, and it resolves
favourably:

- `ClusteringInfoDistance(tree1, tree2)` (and its siblings — `RobinsonFoulds`,
  `TreeDistance`, `MatchingSplitDistance`, …) accept a **single tree against a
  list**, returning a length-`N` vector. That is exactly a `colFn(i)`.
- They accept **precomputed `TreeTools::as.Splits()` objects** in place of
  `phylo` trees, on both sides, and return bit-identical values. So the "good
  closure" of the briefing — extract every tree's splits once in the enclosing
  environment, look them up per call — is directly expressible:

  ```r
  GoodClosure <- function(trees) {
    splits <- as.Splits(trees)                      # once
    function(i) as.numeric(ClusteringInfoDistance(splits[[i]], splits))
  }
  ```

- A column so obtained agrees exactly with the corresponding column of
  `as.matrix(ClusteringInfoDistance(trees))`.

So there is no *interface* obstacle. The obstacle is arithmetic, and it is the
one the briefing did not anticipate.

## 2. The good/bad closure distinction barely matters here

| per column | N = 200 | N = 400 |
|---|---|---|
| good closure (splits precomputed once) | 0.048 s | 0.088 s |
| bad closure (raw `phylo` trees every call) | 0.044 s | 0.088 s |
| bulk all-pairs, amortised per column | 0.009 s | 0.018 s |

The good closure is **not measurably cheaper than the bad one** for this metric.
`ClusteringInfoDistance` already converts its arguments to splits internally, and
that conversion is cheap next to the `N` split-matching problems the call then
solves. The briefing's worry — that a careless caller re-parses trees on every
call and blames `DropAdd()` — is real in principle and worth documenting, but for
TreeDist specifically it is not where the cost is. Precomputing splits is still
the right advice (it costs nothing and other metrics will not be so forgiving);
it is just not the lever.

## 3. The lever is how many columns get requested

This is the finding that determines the claim. The briefing assumed DropAdd
requests `O(n)` or `O(nk)` columns, "never `n^2` of them". That is not right.
The construction is cheap — `k + 1` columns — but the **tabu loop asks for two
columns per iteration and runs until `plateau` consecutive iterations fail to
improve**. The iteration count is not a function of `N`; it is a function of
`plateau`. With the default `plateau = 5000` a run can easily request tens of
thousands of columns from a problem with only `N` distinct ones.

Measured, at `plateau = 200` (already cut 25× below the default just to make the
runs finish):

| N = 200, k = 10 | seconds | oracle calls | pair evaluations | matrix MB |
|---|---|---|---|---|
| bulk all-pairs + matrix kernel | **1.8** | – | 19,900 | 0.3 |
| oracle, good closure | 155.0 | 4,835 | 967,000 | ~0.01 |
| oracle, good closure + memoised | 6.3 | 4,835 (200 distinct) | 40,000 | 0.3 |

| N = 200, k = 25 | seconds | oracle calls | pair evaluations | matrix MB |
|---|---|---|---|---|
| bulk all-pairs + matrix kernel | **1.6** | – | 19,900 | 0.3 |
| oracle, good closure | 62.1 | 1,958 | 391,600 | ~0.01 |
| oracle, good closure + memoised | 6.9 | 1,958 (200 distinct) | 40,000 | 0.3 |

| N = 400, k = 10 | seconds | oracle calls | pair evaluations | matrix MB |
|---|---|---|---|---|
| bulk all-pairs + matrix kernel | **7.2** | – | 79,800 | 1.2 |
| oracle, good closure | 555.8 | 8,314 | 3,325,600 | ~0.02 |
| oracle, good closure + memoised | 26.4 | 8,314 (400 distinct) | 160,000 | 1.2 |

| N = 400, k = 25 | seconds | oracle calls | pair evaluations | matrix MB |
|---|---|---|---|---|
| bulk all-pairs + matrix kernel | **7.3** | – | 79,800 | 1.2 |
| oracle, good closure | 246.2 | 3,567 | 1,426,800 | ~0.02 |
| oracle, good closure + memoised | 27.9 | 3,567 (400 distinct) | 160,000 | 1.2 |

The uncached column oracle performs **18–48× more pairwise distance evaluations
than the full matrix does** across these four cells, because it re-derives the
same columns thousands of times. It is not "less total work at a worse
constant"; it is far more total work. Doubling `N` did not halve the gap: the
oracle's bill scales with the iteration count, which `plateau` sets, not `N`.

Memoising the columns in the closure fixes the redundancy completely — the
oracle collapses to `N` distinct evaluations, which is the *matrix's* work
(×2, since a column pays for both triangles and the bulk routine only computes
one). That leaves a 3.5–4× wall-clock gap to the bulk path — the missing
symmetry plus per-call overhead — and it restores the matrix's memory footprint.

## 4. What the oracle path is actually for

- **Memory.** It never materialises `N × N`. Five length-`N` vectors, whatever
  the closure holds, and nothing else. That is the honest headline, and it is
  the regime where the matrix path simply cannot run: `as.matrix.dist` overflows
  at `N = 46340`, and long before that `N²·8` bytes stops fitting.
- **Metrics with no matrix and no embedding.** The reason this exists: you have
  a distance function and no way to hand `DropAdd()` coordinates.
- **Not speed.** At any `N` where the matrix fits, build the matrix.

The docs (`?DropAdd`, *Distance-column oracle*; `NEWS.md`; the vignette) say
memory and call-count, not speed, and tell callers that `plateau` and
`maxSeconds` set the distance bill directly. The roxygen example shows the
precompute-once **and memoise** closure shape, since §3 shows memoisation is the
one that matters.

## 5. Correctness, incidentally

The benchmark also cross-checks the two paths on a genuinely non-Euclidean
metric: at every `(N, k)` measured, the oracle path returned the **same subset
and the same score** as the matrix path. That is weaker than the trajectory
identity `tests/testthat/test-dropadd-oracle.R` asserts (the benchmark uses each
path's own default seed, which differ), but it is a useful independent check that
nothing about the tree metric — heavy tie structure included — breaks the port.

## 6. Follow-ups not taken in this pass

- **`maxCandidates` thinning on the oracle path.** Currently warns and runs on
  the full problem, because thinning hands the solver an `m × m` coreset matrix.
  It is *not* fundamentally out of reach: `FarFirst()`'s oracle path already
  produces the `m` coreset indices in `O(N·m)` calls, and the `m × m` submatrix
  needs only `m` further columns. Given §3, thinning is arguably the single most
  valuable thing that could be added here — it converts an unbounded per-iteration
  distance bill into a bounded one. Deliberately out of scope; the briefing asked
  for warn-don't-substitute, matching `FarFirst()`.
- **In-package column memoisation.** Left to the caller, because the right cache
  budget depends on `N` (at `N = 1000` the cache is 8 MB; at `N = 50000` it is
  20 GB, i.e. exactly the matrix this path exists to avoid). A caller who knows
  their `N` can make that call; `DropAdd()` cannot.
- **A `plateau` default appropriate to the oracle path.** 5000 is calibrated for
  a compiled kernel reading a matrix. It is the wrong order of magnitude when
  every iteration costs two distance columns. Not changed here, since it would
  make the oracle path's default behaviour differ from the other paths' for
  reasons the API does not express; flagged instead in the docs.

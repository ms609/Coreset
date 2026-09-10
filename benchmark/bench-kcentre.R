source("benchmark/_init.R")

# Exact discrete k-centre: a bisection over the distinct distances whose every
# probe is decided by the exhaustive cover search (src/kcentre_cover.cpp), warm
# started by the CDSh heuristic. Two sizes because the two halves dominate on
# opposite sides of .kCentreExhaustiveMaxCand.

# Below the cutoff KCentre() scans the candidate grid exhaustively, so this
# guards the warm start and the orchestration: the bisection takes 9 probes and
# the search under 200 nodes in total.
d80 <- BenchDist(80L, 5L, seed = 1L)
Benchmark(ExactKCentre(8L, d80))

# Past the cutoff the warm start falls back to a binary search and the solve
# becomes the cover search itself: 13 probes over some 77,000 nodes, of which
# the costliest probe is 29,000. Guards the reduction, the component split and
# the branching.
d200 <- BenchDist(200L, 5L, seed = 2L)
Benchmark(ExactKCentre(8L, d200))

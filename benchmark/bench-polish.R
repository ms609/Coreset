source("benchmark/_init.R")

# Critical-edge 1-swap polish and the T_k objective. Sized (large subset) so the
# timings sit well above CI jitter, where a regression is detectable.
d3000 <- BenchDist(3000L, 8L, seed = 1L)
# Polish an arbitrary 200-subset rather than a near-optimal Gonzalez seed, so
# the 1-swap loop performs many improving moves and the timing carries signal.
Benchmark(PolishSelection(d3000, seq_len(200L)))
Benchmark(TkScore(d3000, seq_len(800L)))

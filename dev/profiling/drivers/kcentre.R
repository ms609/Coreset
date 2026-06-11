# Area 5 — KCentre / CDSh kernel (KCentreCDSh_cpp).
# Matrix-bound O(n^2 log n): binary search over sorted distances, each trial an
# O(n^2) fixed-k dominating-set construction. Clustered Gaussian data gives a
# representative threshold-graph density (the d(i,j) <= r branch distribution
# that drives the inner loops), unlike uniform noise.
# bare: ~? s
suppressMessages(library(MaxMin, lib.loc = "C:/Users/pjjg18/AppData/Local/Temp/maxmin-optlib"))
set.seed(5813)

g <- 12L; per <- 167L                      # n ~ 2004, 12 clusters in 10-D
ctr <- matrix(rnorm(g * 10L, sd = 6), ncol = 10L)
pts <- ctr[rep(seq_len(g), each = per), ] + matrix(rnorm(g * per * 10L), ncol = 10L)
d <- as.matrix(dist(pts))
n <- nrow(d)

reps <- 4L
k <- 20L
t0 <- proc.time()
for (i in seq_len(reps)) {
  res <- KCentre(d, k)
}
el <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("n=%d k=%d reps=%d  elapsed=%.2fs  per-call=%.1fms  radius=%.4f\n",
            n, k, reps, el, 1000 * el / reps, attr(res, "radius")))

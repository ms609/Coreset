o <- readRDS("dev/profiling/.kc-old.rds")
n <- readRDS("dev/profiling/.kc-new.rds")
ok <- TRUE
for (k in names(o$res)) {
  a <- o$res[[k]]; b <- n$res[[k]]
  idx_id <- identical(a$idx, b$idx)
  rad_id <- identical(a$radius, b$radius)
  ok <- ok && idx_id && rad_id
  cat(sprintf("k=%-3s idx_identical=%s radius_identical=%s\n", k, idx_id, rad_id))
}
cat(sprintf("ALL BIT-IDENTICAL: %s   speedup=%.2fx (%.0f -> %.0f ms)\n",
            ok, o$ms / n$ms, o$ms, n$ms))

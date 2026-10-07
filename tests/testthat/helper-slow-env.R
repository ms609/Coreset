# Detects the memcheck/ASan test runs, where every compiled call can be an
# order of magnitude slower than on a plain test box. Tests that race a real
# computation against a short wall-clock budget use this to relax -- never
# skip -- their assertions, so they keep covering the fast path while staying
# green under valgrind/ASan. memcheck/tests.R sets the flag; the shared
# ms609/actions memcheck action runs test_local() itself and never sources
# that script, so a process running under valgrind (its preload library in
# LD_PRELOAD) counts too.
.SlowSearchEnv <- function() {
  identical(Sys.getenv("CORESET_SLOW_TEST_ENV"), "true") ||
    grepl("vgpreload", Sys.getenv("LD_PRELOAD"), fixed = TRUE)
}

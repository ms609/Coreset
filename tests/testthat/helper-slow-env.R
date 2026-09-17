# Detects the memcheck/ASan test runs (flagged by memcheck/tests.R), where
# every compiled call can be an order of magnitude slower than on a plain
# test box. Tests that race a real computation against a short wall-clock
# budget use this to relax -- never skip -- their assertions, so they keep
# covering the fast path while staying green under valgrind/ASan.
.SlowSearchEnv <- function() identical(Sys.getenv("CORESET_SLOW_TEST_ENV"), "true")

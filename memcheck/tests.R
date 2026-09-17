# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < memcheck/thisfile.R
# Flags timing-sensitive tests (see helper-slow-env.R) so they relax their
# wall-clock-budget assertions instead of racing valgrind's slowdown.
Sys.setenv(CORESET_SLOW_TEST_ENV = "true")
testthat::test_local()

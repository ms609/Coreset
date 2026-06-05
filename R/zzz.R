.onUnload <- function(libpath) {
  library.dynam.unload("MaxMin", libpath)
}

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Does R CMD check pass with no NOTEs?",
    "Have the FurthestPoint reproduction tests been re-run against this version?"
  )
}


# Additional tests:
#
# spell_check()
# pkgdown::build_reference_index()
# run_examples()
#
# devtools::check_win_devel(quiet = TRUE); rhub::check_for_cran()
# Check valgrind / ASan results on GitHub Actions
# revdepcheck::revdep_check()
#
# codemeta::write_codemeta()

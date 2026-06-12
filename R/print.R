# print.R
#
# A light S3 layer that makes the solver outputs self-describing at the
# console. The selection-returning solvers ([FarFirst()], [DropAdd()],
# [Grasp()]) attach the `"MaxMinSelection"` class to their integer index
# vector; [ExactMaxMin()] attaches `"MaxMinExact"` to its result list. Both
# remain ordinary objects -- a `MaxMinSelection` still indexes a matrix or
# coordinate set directly, and a `MaxMinExact` is still a plain list with
# `$indices`, `$objective`, ... -- so the class only changes how they print.

#' Stamp the `MaxMinSelection` class onto a solver's index vector
#'
#' The selection-returning solvers each already attach their score (and any
#' secondary attributes) via [base::structure()]; this adds the shared S3 class
#' and a `producer` tag the print method reads to name the algorithm. An empty
#' selection (`length 0`) is left bare: there is nothing to describe.
#' @param x Integer index vector carrying the solver's score attributes.
#' @param producer Character tag naming the solver (`"FarFirst"`, `"DropAdd"`,
#'   `"Grasp"`).
#' @return `.AsMaxMinSelection()` returns `x` with `producer` attribute and
#'   `"MaxMinSelection"` class, or `x` unchanged if it is empty.
#' @keywords internal
.AsMaxMinSelection <- function(x, producer) {
  if (length(x) == 0L) {
    # Return:
    x
  } else {
    attr(x, "producer") <- producer
    class(x) <- "MaxMinSelection"
    # Return:
    x
  }
}

#' Format a selected-index list, optionally truncated
#'
#' @param idx Integer indices in their stored order.
#' @param maxShow Integer: show at most this many before eliding the tail.
#' @return `.FormatIndexList()` returns a length-1 character string such as
#'   `"6 5 4 3 1 2"`, or `"1 2 ... (+15 more)"` when `idx` is longer than `maxShow`.
#' @keywords internal
.FormatIndexList <- function(idx, maxShow = 20L) {
  if (length(idx) > maxShow) {
    head <- paste(idx[seq_len(maxShow)], collapse = " ")
    # Return:
    sprintf("%s ... (+%d more)", head, length(idx) - maxShow)
  } else {
    # Return:
    paste(idx, collapse = " ")
  }
}

#' Compose the one-line selection summary shared by both print methods
#'
#' @param n Integer count of selected elements.
#' @param idx Integer selected indices in stored order.
#' @param by Character phrase naming the algorithm (the "selected by ..." part).
#' @param tk Numeric achieved \eqn{T_k} (the minimum pairwise distance), or `NA`
#'   when there is no pairwise distance to report (a single element, or all
#'   elements selected).
#' @param maxShow Integer index-list truncation threshold; see
#'   [.FormatIndexList()].
#' @return `.MaxMinSummaryLine()` returns a length-1 character string.
#' @keywords internal
.MaxMinSummaryLine <- function(n, idx, by, tk, maxShow = 20L) {
  elements <- sprintf("%d element%s", n, if (n == 1L) "" else "s")
  dist <- if (is.na(tk)) {
    ""
  } else {
    sprintf(", each at distance >= %s", format(signif(tk, 4L)))
  }
  # Return:
  sprintf("%s (%s) selected by %s%s",
          elements, .FormatIndexList(idx, maxShow), by, dist)
}

#' Name the seeding outcome of a [FarFirst()] selection
#'
#' A bare single pass is just "farthest-first"; an ensemble pass
#' (which carries `winning_strategy` / `strategy_results`) additionally names
#' the winning strategy and how many were tried.
#' @param x A `MaxMinSelection` from [FarFirst()].
#' @return `.FarFirstSelectedBy()` returns a length-1 character phrase.
#' @keywords internal
.FarFirstSelectedBy <- function(x) {
  winners <- attr(x, "winning_strategy")
  if (is.null(winners)) {
    # Return:
    "farthest-first"
  } else {
    nStrat <- length(attr(x, "strategy_results"))
    plural <- if (nStrat == 1L) "y" else "ies"
    detail <- if (length(winners) == 1L) {
      sprintf("winner %s", winners)
    } else {
      sprintf("%d tied: %s", length(winners), paste(winners, collapse = ", "))
    }
    # Return:
    sprintf("farthest-first (best of %d strateg%s, %s)",
            nStrat, plural, detail)
  }
}

#' @keywords internal
.MaxMinSelectedBy <- function(x) {
  switch(attr(x, "producer") %||% "",
    FarFirst = .FarFirstSelectedBy(x),
    DropAdd  = "DropAdd tabu search",
    Grasp    = "GRASP with path-relinking",
    "an unrecorded method"
  )
}

#' Format and print MaxMin solver results
#'
#' Terse summaries of the objects returned by the MaxMin solvers.
#'
#' @param x A `MaxMinSelection` or `MaxMinExact` object.
#' @param ... Ignored; present for S3 compatibility.
#' @return
#' `print.MaxMin()` returns `x`, invisibly. It is called for its side-effect
#' of printing `format(x)` to the console.
#' `format.MaxMin()` returns a character string describing a `MaxMinSelection`
#' (from [FarFirst()], [DropAdd()] and [Grasp()]); it reports its size, the
#' selected indices, the algorithm (and if applicable strategy),
#' and the achieved \eqn{T_k}. A `MaxMinExact` (from [ExactMaxMin()]) object
#' additionally states whether optimality was proven.
#'
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' print(FarFirst(5L, dist(pts)))
#' @name print.MaxMin
#' @family reporting functions
#' @export
format.MaxMinSelection <- function(x, ...) {
  idx <- as.integer(x)
  # Return:
  .MaxMinSummaryLine(length(idx), idx, .MaxMinSelectedBy(x), attr(x, "score"))
}

#' @rdname print.MaxMin
#' @export
print.MaxMinSelection <- function(x, ...) {
  cat(format(x, ...), "\n", sep = "")
  # Return:
  invisible(x)
}

#' @rdname print.MaxMin
#' @export
format.MaxMinExact <- function(x, ...) {
  proven <- if (isTRUE(x$proven)) {
    sprintf("exact MILP (%s), proven optimal", x$solver)
  } else {
    sprintf("exact MILP (%s), unproven incumbent", x$solver)
  }
  # Return:
  .MaxMinSummaryLine(length(x$indices), x$indices, proven, x$objective)
}

#' @rdname print.MaxMin
#' @export
print.MaxMinExact <- function(x, ...) {
  cat(format(x, ...), "\n", sep = "")
  # Return:
  invisible(x)
}

#' Format a numeric field for a summary, tolerating `NA`
#' @param v Numeric scalar.
#' @return `.SummaryNum()` returns a length-1 character: `"NA"`, or four significant figures.
#' @keywords internal
.SummaryNum <- function(v) {
  if (length(v) == 0L || is.na(v)) "NA" else format(signif(v, 4L))
}

#' Print a `label: value` detail line under a summary headline
#' @param label,value Character (or coercible) field label and value.
#' @param width Integer column width the labels are padded to.
#' @return `.SummaryField()` returns invisibly `NULL`; called for the side effect.
#' @keywords internal
.SummaryField <- function(label, value, width) {
  cat(sprintf("  %-*s %s\n", width, paste0(label, ":"), value))
}

#' Print the per-strategy \eqn{T_k} table of a [FarFirst()] ensemble
#'
#' One row per strategy tried, ordered best (largest \eqn{T_k}) first, with each
#' tied-best strategy marked `*`. A bare single pass (no `strategy_results`)
#' produces nothing.
#' @param object A `MaxMinSelection` from an ensemble [FarFirst()] call.
#' @return `.SummariseStrategies()` returns invisibly `NULL`; called for the side effect.
#' @keywords internal
.SummariseStrategies <- function(object) {
  sr <- attr(object, "strategy_results")
  if (is.null(sr)) {
    return(invisible(NULL))           # a bare single pass: nothing to tabulate
  }
  winners <- attr(object, "winning_strategy")
  labels  <- names(sr)
  seeds   <- vapply(sr, function(e) as.integer(e$s1), integer(1L))
  tks     <- vapply(sr, function(e) as.numeric(e$t_k), numeric(1L))
  ord     <- order(tks, decreasing = TRUE, na.last = TRUE)
  labels  <- labels[ord]; seeds <- seeds[ord]; tks <- tks[ord]
  mark    <- ifelse(labels %in% winners, "*", " ")
  w       <- max(nchar(labels), nchar("strategy"))

  cat(sprintf("  strategies tried (%d), best marked *:\n", length(sr)))
  cat(sprintf("    %-*s  %5s  %s\n", w, "strategy", "seed", "T_k"))
  for (i in seq_along(labels)) {
    cat(sprintf("  %s %-*s  %5d  %s\n",
                mark[i], w, labels[i], seeds[i], .SummaryNum(tks[i])))
  }
  # Return:
  invisible(NULL)
}

#' Detailed summaries of MaxMin solver results
#'
#' A fuller counterpart to [print.MaxMinSelection()]: the one-line headline,
#' followed by the achieved objective(s), search effort, and -- for a
#' [FarFirst()] ensemble -- the per-strategy \eqn{T_k} table.
#' @param object A `MaxMinSelection` or `MaxMinExact` object.
#' @param ... Ignored; present for S3 compatibility.
#' @return `summary.MaxMin()` returns `object`, invisibly.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' summary(FarFirst(5L, dist(pts)))
#' @name summary.MaxMin
#' @family reporting functions
#' @export
summary.MaxMinSelection <- function(object, ...) {
  cat(format(object), "\n", sep = "")
  switch(attr(object, "producer") %||% "",
    FarFirst = .SummariseStrategies(object),
    DropAdd  = {
      .SummaryField("sum of pairwise distances", .SummaryNum(attr(object, "secondary")), 26L)
      .SummaryField("iterations", attr(object, "iters"), 26L)
      .SummaryField("time", paste(.SummaryNum(attr(object, "time_s")), "s"), 26L)
    },
    Grasp = {
      .SummaryField("refinement iterations", attr(object, "iters"), 26L)
      .SummaryField("path-relinking calls", attr(object, "pr_calls"), 26L)
      .SummaryField("time", paste(.SummaryNum(attr(object, "time_s")), "s"), 26L)
    }
  )
  # Return:
  invisible(object)
}

#' @rdname summary.MaxMin
#' @export
summary.MaxMinExact <- function(object, ...) {
  cat(format(object), "\n", sep = "")
  status <- if (isTRUE(object$proven)) "proven optimal" else "lower bound (unproven)"
  .SummaryField("instance", sprintf("n = %d, k = %d", object$n, object$k), 12L)
  .SummaryField("objective", sprintf("%s (%s)", .SummaryNum(object$objective), status), 12L)
  .SummaryField("solver", object$solver, 12L)
  .SummaryField("time", paste(.SummaryNum(object$time_s), "s"), 12L)
  # Return:
  invisible(object)
}

#' Format and print k-centre solver results
#'
#' Terse summaries of the objects returned by [KCentre()]
#' (`"KCentreSelection"`) and [ExactKCentre()] (`"KCentreExact"`)
#'
#' @param x A `"KCentreSelection"` or `"KCentreExact"` object.
#' @param ... Ignored; present for S3 compatibility.
#' @return
#' `print.KCentre()` returns `x`, invisibly. It is called for its side-effect
#' of printing `format(x)` to the console.
#' `format.KCentre()` returns a character string describing a `KCentreSelection`;
#' it reports the centre count, the chosen indices, the method, and the achieved
#' covering radius (with proof status for the exact solver).
#'
#' @name print.KCentre
#' @family reporting functions
#' @examples
#' set.seed(1)
#' KCentre(4L, dist(matrix(rnorm(60), ncol = 2)))
#' @export
format.KCentreSelection <- function(x, ...) {
  idx <- as.integer(x)
  nc <- length(idx)
  sprintf("%d centre%s (%s) by CDSh, covering radius <= %s",
          nc, if (nc == 1L) "" else "s", .FormatIndexList(idx),
          format(signif(attr(x, "radius"), 4L)))
}

#' @rdname print.KCentre
#' @export
print.KCentreSelection <- function(x, ...) {
  cat(format(x, ...), "\n", sep = "")
  invisible(x)
}

#' @rdname print.KCentre
#' @export
format.KCentreExact <- function(x, ...) {
  nc <- length(x$indices)
  status <- if (isTRUE(x$proven)) {
    sprintf("exact MILP (%s), proven optimal", x$solver)
  } else {
    sprintf("exact MILP (%s), unproven incumbent", x$solver)
  }
  rel <- if (isTRUE(x$proven)) "=" else "<="
  sprintf("%d centre%s (%s) by %s, covering radius %s %s",
          nc, if (nc == 1L) "" else "s", .FormatIndexList(x$indices),
          status, rel, format(signif(x$radius, 4L)))
}

#' @rdname print.KCentre
#' @export
print.KCentreExact <- function(x, ...) {
  cat(format(x, ...), "\n", sep = "")
  invisible(x)
}

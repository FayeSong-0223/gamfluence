#' Methods for gam_influence objects
#'
#' @param x,object A `gam_influence` object.
#' @param n Number of sites to show or return.
#' @param by Column of the site table to rank by (largest first):
#'   `"total"`, `"fixed"` or `"added"`.
#' @param ... Unused.
#' @name gam_influence-methods
NULL

#' @rdname gam_influence-methods
#' @export
print.gam_influence <- function(x, n = 5, ...) {
  s <- x$sites
  cat("gamfluence: fixed vs re-estimated deletion, ", nrow(s), " sites\n", sep = "")
  cat("  model: ", x$info$family, " (", x$info$link, "), ", x$info$n, " rows\n", sep = "")
  bad <- s$status != "ok"
  if (any(bad))
    cat("  status: ", paste(names(table(s$status[bad])), table(s$status[bad]),
                            sep = " x", collapse = "; "), "\n", sep = "")
  cat("\nLargest total change in the other sites' fitted values",
      " (RMS, link scale):\n", sep = "")
  print(format_df(top_units(x, by = "total", n = n)[, c("id", "n_obs", "fixed", "added", "total")]),
        row.names = FALSE)
  if (!is.null(x$hyper)) {
    h <- x$hyper[!is.na(x$hyper$ratio), , drop = FALSE]
    if (nrow(h)) {
      h <- h[order(abs(log(h$ratio)), decreasing = TRUE), , drop = FALSE]
      cat("\nLargest hyperparameter changes when a site is removed:\n")
      print(format_df(utils::head(h[, c("id", "parameter", "full", "deleted", "ratio")], n)),
            row.names = FALSE)
    }
  }
  cat("\nfixed = hyperparameters held; added = what re-estimating them adds.\n",
      "No null reference: compare sites with each other.\n", sep = "")
  invisible(x)
}

#' @rdname gam_influence-methods
#' @export
summary.gam_influence <- function(object, ...) {
  s <- object$sites
  q <- t(vapply(c("fixed", "added", "total"), function(v) {
    z <- s[[v]]
    c(stats::quantile(z, c(0, .5, .9, 1), na.rm = TRUE), n_NA = sum(is.na(z)))
  }, numeric(5)))
  colnames(q) <- c("min", "median", "90%", "max", "NA")
  rng <- NULL
  if (!is.null(object$hyper)) {
    h <- object$hyper
    rng <- do.call(rbind, lapply(split(h, h$parameter), function(d)
      data.frame(parameter = d$parameter[1], type = d$type[1], full = d$full[1],
                 min_ratio = suppressWarnings(min(d$ratio, na.rm = TRUE)),
                 max_ratio = suppressWarnings(max(d$ratio, na.rm = TRUE)))))
    rownames(rng) <- NULL
  }
  out <- list(quantiles = q, hyper = rng, n_sites = nrow(s),
              status = table(s$status), info = object$info)
  class(out) <- "summary.gam_influence"
  out
}

#' @rdname gam_influence-methods
#' @export
print.summary.gam_influence <- function(x, ...) {
  cat("gamfluence summary: ", x$n_sites, " sites\n\n", sep = "")
  print(signif(x$quantiles, 3))
  if (!is.null(x$hyper)) {
    cat("\nHyperparameters: full-data value and range of ratios when a site is removed\n")
    print(format_df(x$hyper), row.names = FALSE)
  }
  cat("\nStatus:\n")
  print(x$status)
  invisible(x)
}

#' @rdname gam_influence-methods
#' @param row.names,optional Accepted for compatibility; ignored.
#' @export
as.data.frame.gam_influence <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$sites
}

#' @rdname gam_influence-methods
#' @export
top_units <- function(x, by = "total", n = 10) {
  if (!inherits(x, "gam_influence")) stop("`x` must be a gam_influence object.", call. = FALSE)
  if (!(by %in% c("fixed", "added", "total"))) stop("`by` must be \"fixed\", \"added\" or \"total\".", call. = FALSE)
  s <- x$sites[order(x$sites[[by]], decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  rownames(s) <- NULL
  utils::head(s, n)
}

format_df <- function(d) {
  for (v in names(d)) if (is.double(d[[v]])) d[[v]] <- signif(d[[v]], 3)
  d
}

#' Compare fixed and re-estimated deletion effects across sites
#'
#' For each site (cluster), the model is refitted without that site twice on
#' the frozen full-data representation: once with the estimated smoothing
#' parameters (and theta) held at their full-fit values, and once with them
#' re-estimated. The changes in the fitted values of the remaining sites are
#' compared.
#'
#' @section What is reported:
#' For each site, three weighted RMS changes in the linear predictor of the
#' *other* sites' rows (the retained rows):
#'
#' * `fixed`: refit without the site, hyperparameters fixed at their full-fit
#'   values. This is what a deletion diagnostic with fixed smoothing
#'   parameters sees.
#' * `total`: refit without the site, hyperparameters re-estimated.
#' * `added`: the difference between the two refits, i.e. what re-estimating
#'   the hyperparameters adds.
#'
#' For every retained observation, fixed + added = total exactly (they are
#' changes in the same fitted value). The three columns are RMS summaries of
#' those changes, and RMS values do not add: `total` is never larger than
#' `fixed + added`, and equals it only when one change is a non-negative
#' multiple of the other (for example when `added` is zero). A smaller
#' `total` therefore only says that the two changes are not perfectly
#' aligned; on its own it does not mean that re-estimation offsets the fixed
#' change. Offsetting shows up as `total`^2 < `fixed`^2 + `added`^2.
#'
#' Magnitudes use the full fit's IRLS (Fisher) weights, so rows with more
#' information count more, and are on the link scale. On a log link a change
#' of d in the linear predictor multiplies the predicted mean by exp(d), a
#' relative change of exp(d) - 1. When the individual log changes are small,
#' exp(d) - 1 is close to d, so an RMS of 0.02 approximates a weighted RMS
#' relative change of about 2% in the other sites' predicted means. The
#' approximation fails for large individual changes (d = 0.2 is +22%,
#' d = -0.2 is -18%), and a small RMS can still hide a few large changes.
#' There is no null reference, so compare sites with each other rather than
#' with a fixed cut-off.
#'
#' The hyperparameters are the smoothing parameters mgcv estimated and, for
#' `nb()` with theta estimated, theta. Smoothing parameters or theta fixed by
#' the user stay fixed in both refits. All results are conditional on the
#' frozen representation: the full fit's model matrix, identifiability
#' constraints and penalties are reused, not rebuilt from the reduced data.
#'
#' @section Status:
#' `status` is `"ok"` or lists what happened: `"rank-deficient"` (the site was
#' the only one informing some coefficients; they are dropped from its refits,
#' which leaves the reported changes valid), `"lambda held"` (a smooth had no
#' data left, so its smoothing parameter was held; see `note`),
#' `"two optima"` (see below), `"not converged"` or `"error"`. Sites whose
#' refits failed have `NA` magnitudes and a warning is issued.
#'
#' The REML criterion can have more than one optimum, typically an interior
#' one and a flat plateau where a smoothing parameter is very large (for a
#' thin plate smooth, a nearly straight line). The
#' re-estimated refit is therefore searched for twice, from mgcv's default
#' starting values and from the full-fit values, and the fit with the lower
#' REML score is used. When the two searches end at fits that differ by more
#' than 10% of the larger total change, the status says `"two optima"` and
#' `note` gives the other fit's `total` and how much worse its REML score was.
#' A small difference in REML score means the data barely prefer one fit over
#' the other, so treat that site's `added` and `total` with care. Two searches
#' still do not guarantee the best optimum overall.
#'
#' @param model A fitted [mgcv::gam()] object within the scope described in
#'   the README; anything the scope guard recognises as outside it is an error.
#' @param cluster The sites: a column name (in the model frame or in `data`),
#'   or a vector with one value per row of `data` or per row used by the
#'   model. Rows dropped by `na.action` or `subset` are matched by row name,
#'   and the match is checked against the data values.
#' @param data The data frame used to fit `model`. Only needed when `cluster`
#'   is not a variable of the model; if `NULL`, the data named in the model
#'   call is used, with a message.
#' @param which Optional subset of site ids to evaluate.
#' @param progress Show a text progress bar.
#'
#' @return An object of class `gam_influence`, a list with
#'   * `sites`: one row per site: `id`, `n_obs`, `fixed`, `added`, `total`,
#'     `status`, `note`;
#'   * `hyper`: one row per site and hyperparameter: `id`, `parameter`,
#'     `type` (`"sd"` for random-effect terms, reported as the standard
#'     deviation; `"lambda"` for other smooths; `"theta"`), `full`, `deleted`
#'     (re-estimated without the site) and `ratio`. A very large lambda
#'     ratio means the smooth became close to a straight line; the exact
#'     value is then not meaningful, because the REML score is flat there;
#'   * `info`: family, number of rows, and the version numbers used.
#'
#' @examples
#' library(mgcv)
#' set.seed(1)
#' d <- data.frame(site = factor(rep(1:12, each = 6)), x = runif(72))
#' d$y <- sin(2 * pi * d$x) + rnorm(12, 0, 0.4)[d$site] + rnorm(72, 0, 0.5)
#' m <- gam(y ~ s(x) + s(site, bs = "re"), data = d, method = "REML")
#' infl <- gam_influence(m, cluster = "site", progress = FALSE)
#' infl
#' top_units(infl, by = "added", n = 3)
#' @export
gam_influence <- function(model, cluster, data = NULL, which = NULL,
                          progress = interactive()) {
  if (missing(cluster))
    stop("`cluster` is required: the column (or vector) that defines the sites.",
         call. = FALSE)
  st    <- frozen_setup(model)
  units <- resolve_units(model, data, cluster, which, env = parent.frame())

  n_units  <- length(units$id)
  all_rows <- seq_along(st$y)
  est  <- st$map$estimated
  re_j <- which(st$map$is_re)
  sd0  <- vapply(seq_along(st$S), function(j)
    if (st$map$is_re[j]) re_sd(st$lambda, st$phi_reml, st, j) else NA_real_, numeric(1))

  site_rows  <- vector("list", n_units)
  hyper_rows <- vector("list", n_units)
  pb <- if (progress && n_units > 1) utils::txtProgressBar(0, n_units, style = 3) else NULL
  for (k in seq_len(n_units)) {
    keep  <- setdiff(all_rows, units$rows[[k]])
    plan  <- NULL
    warns <- character(0)
    res <- tryCatch(
      withCallingHandlers({
          plan <- deletion_plan(st, keep)
          compare_fits(st, keep, plan)
        },
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }),
      error = function(e) e)
    if (is.null(plan)) plan <- no_plan(st)
    held <- st$map$term[plan$hold & est]

    if (inherits(res, "error")) {
      site_rows[[k]] <- data.frame(id = units$id[k], n_obs = length(units$rows[[k]]),
        fixed = NA_real_, added = NA_real_, total = NA_real_, status = "error",
        note = conditionMessage(res), stringsAsFactors = FALSE)
      next
    }
    ok <- res$A$converged && res$T$converged
    w  <- st$w[keep]
    # two REML optima that give noticeably different fits (more than 10% of
    # the larger total change apart): the better one is used, the other noted
    two_opt <- FALSE
    if (ok && !is.null(res$alt)) {
      alt_total <- wrms(res$alt$total, w)
      two_opt <- wrms(res$total - res$alt$total, w) >
                 0.1 * max(wrms(res$total, w), alt_total)
    }
    flags <- c(if (plan$rank_deficient) "rank-deficient",
               if (length(held)) "lambda held",
               if (two_opt) "two optima",
               if (!ok) "not converged")
    notes <- c(if (length(held)) paste("held:", paste(held, collapse = ", ")),
               if (two_opt) sprintf(paste("two REML optima: the better one is used (REML score",
                                          "lower by %.2g); the other gives total %.3g"),
                                    res$alt$reml_gap, alt_total),
               if (length(warns)) paste("warning:", warns[1]))
    site_rows[[k]] <- data.frame(id = units$id[k], n_obs = length(units$rows[[k]]),
      fixed = wrms(res$fixed, w), added = wrms(res$added, w), total = wrms(res$total, w),
      status = if (length(flags)) paste(flags, collapse = "; ") else "ok",
      note = if (length(notes)) paste(notes, collapse = "; ") else NA_character_,
      stringsAsFactors = FALSE)
    if (!ok) site_rows[[k]][c("fixed", "added", "total")] <- NA_real_

    # hyperparameters: random-effect SDs, other estimated lambdas, theta
    hj <- which(est | st$map$is_re)
    if (!length(hj) && is.null(st$psi)) {
      if (!is.null(pb)) utils::setTxtProgressBar(pb, k)
      next
    }
    full <- deleted <- numeric(0)
    for (j in hj) {
      if (st$map$is_re[j]) {
        full    <- c(full, sd0[j])
        deleted <- c(deleted, if (plan$hold[j]) NA_real_ else re_sd(res$T$lambda, res$T$phi, st, j))
      } else {
        full    <- c(full, st$lambda[j])
        deleted <- c(deleted, if (plan$hold[j]) NA_real_ else res$T$lambda[j])
      }
    }
    h <- if (length(hj)) {
      data.frame(id = units$id[k], parameter = st$map$term[hj],
                 type = ifelse(st$map$is_re[hj], "sd", "lambda"),
                 full = full, deleted = deleted, stringsAsFactors = FALSE)
    } else {                                # e.g. nb() with only user-fixed lambda: theta alone
      data.frame(id = character(0), parameter = character(0), type = character(0),
                 full = numeric(0), deleted = numeric(0), stringsAsFactors = FALSE)
    }
    if (!is.null(st$psi))
      h <- rbind(h, data.frame(id = units$id[k], parameter = "theta", type = "theta",
                               full = st$psi, deleted = res$T$psi, stringsAsFactors = FALSE))
    h$ratio <- h$deleted / h$full
    if (!ok) h$deleted <- h$ratio <- NA_real_
    hyper_rows[[k]] <- h
    if (!is.null(pb)) utils::setTxtProgressBar(pb, k)
  }
  if (!is.null(pb)) close(pb)

  sites <- do.call(rbind, site_rows); rownames(sites) <- NULL
  n_bad <- sum(grepl("error|not converged", sites$status))
  if (n_bad > 0)
    warning(n_bad, " of ", n_units, " sites had a refit that failed or did not ",
            "converge; see `status` and `note`.", call. = FALSE)
  hyper <- do.call(rbind, Filter(Negate(is.null), hyper_rows))
  if (!is.null(hyper)) rownames(hyper) <- NULL

  out <- list(sites = sites, hyper = hyper,
              info = list(family = model$family$family, link = model$family$link,
                          n = length(st$y), n_sites = n_units,
                          dropped_at_setup = st$dropped_at_setup,
                          gamfluence = as.character(utils::packageVersion("gamfluence")),
                          mgcv = as.character(utils::packageVersion("mgcv"))))
  class(out) <- "gam_influence"
  out
}

# Two refits per site, compared with the original fit.
#
#   F0  all rows,      full-fit hyperparameters           (the original fit)
#   A   site deleted,  hyperparameters held at F0's values ("fixed")
#   T   site deleted,  hyperparameters re-estimated        ("re-estimated")
#
# Each is a well-defined refit that answers its own question:
#   fixed = A - F0  what a deletion diagnostic with fixed smoothing
#                   parameters (and theta) would see
#   total = T - F0  what refitting without the site gives
#   added = T - A   what letting the hyperparameters re-tune adds
# fixed + added = total exactly. This is a comparison of two refits, not a
# decomposition into causes, so there is no ordering to choose.
#
# The hyperparameters are the smoothing parameters mgcv estimated and, for
# nb() with theta estimated, theta; user-fixed ones stay fixed in both refits.
# Vectors are on the retained rows. `alt` is set when the two REML searches
# for T (see reestimate()) ended at different optima: the other optimum's
# total change and how much worse its REML score was.
compare_fits <- function(st, keep, plan = no_plan(st)) {
  A <- refit(st, keep, lambda = st$lambda, psi = st$psi, plan = plan)
  Tk <- if (!any(st$map$estimated & !plan$hold) && is.null(st$psi)) {
    A                                   # nothing to re-estimate: T is A
  } else {
    reestimate(st, keep, plan)
  }
  e0 <- st$eta0[keep]
  list(fixed = A$eta_all[keep] - e0,
       added = Tk$eta_all[keep] - A$eta_all[keep],
       total = Tk$eta_all[keep] - e0,
       alt   = if (is.null(Tk$alt)) NULL else
                 list(total = Tk$alt$eta_all[keep] - e0, reml_gap = Tk$alt$reml_gap),
       A = A, T = Tk)
}

# The REML criterion can have more than one local optimum: for example an
# interior optimum for a smoothing parameter next to a flat plateau towards
# infinity (the smooth close to linear). mgcv's search from its default
# starting values can stop on the plateau even when the interior optimum is
# better, and a search started at the full-fit values can stay near the old
# optimum when the plateau has become better. T is therefore searched for
# twice, from mgcv's defaults and from the full-fit values, and the converged
# fit with the lower REML score is kept. The other one is returned in `alt`.
# The default search is kept unless the other is better by more than a
# numerical tolerance: on a flat plateau the two can differ by noise only.
# Warnings from the second search are not passed on (started at an optimum,
# mgcv often reports a step failure); it only competes if it converged.
reestimate <- function(st, keep, plan = no_plan(st)) {
  a <- refit(st, keep, plan = plan)
  b <- tryCatch(suppressWarnings(refit(st, keep, plan = plan, start = "full")),
                error = function(e) NULL)
  ok <- Filter(function(f) !is.null(f) && isTRUE(f$converged), list(a, b))
  if (length(ok) < 2L) return(if (length(ok)) ok[[1]] else a)
  s <- vapply(ok, function(f) f$fit$gcv.ubre, numeric(1))
  i <- if (s[2] < s[1] - 1e-6 * (1 + abs(s[1]))) 2L else 1L
  best <- ok[[i]]
  best$alt <- list(eta_all = ok[[3L - i]]$eta_all, reml_gap = abs(s[1] - s[2]))
  best
}

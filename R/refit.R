# Refits on the frozen representation, and what to do when a deletion leaves
# some coefficients or smoothing parameters without information.

# deletion_plan(): inspect the retained rows before refitting.
#
# * Coefficients that are not identifiable from the retained rows plus the
#   penalties (a direction v with X_kept v = 0 and S v = 0, e.g. a factor level
#   whose only data are in the deleted unit) are dropped from the A and T
#   refits. Fitted values on the retained rows do not depend on how such a
#   direction is resolved, so the reported changes stay valid. (Predictions on
#   the deleted rows and coefficient values would depend on it; gamfluence
#   reports neither.)
# * A penalty whose columns are all zero on the retained rows (e.g. the
#   by = factor smooth of a level that only occurs in the deleted unit) has an
#   unidentified smoothing parameter. It is held at its full-fit value in T
#   and reported as held.
# * A random-effect level with no retained rows is NOT a problem: its column is
#   zero but penalised, so it stays in (feasibility check 2.4).
deletion_plan <- function(st, keep, tol = 1e-7) {
  p  <- ncol(st$X)
  Xk <- st$X[keep, , drop = FALSE] * sqrt(st$pw[keep])
  Z  <- rbind(Xk, st$S_half)
  cn <- sqrt(colSums(Z^2)); cn[cn == 0] <- 1
  Zs <- sweep(Z, 2, cn, "/")
  q  <- qr(Zs, tol = tol)
  drop <- integer(0); affected <- integer(0)
  if (q$rank < p) {
    drop <- sort(q$pivot[(q$rank + 1L):p])
    sv <- svd(Zs, nu = 0)
    V0 <- sv$v[, (q$rank + 1L):p, drop = FALSE]
    affected <- which(apply(abs(V0), 1, max) > 1e-6)
  }
  hold <- vapply(seq_along(st$S), function(j) {
    cj <- setdiff(st$pen_cols[[j]], drop)
    length(cj) == 0L || all(Xk[, cj, drop = FALSE] == 0)   # weighted: zero-weight rows carry no data
  }, logical(1))
  # The remaining columns of a held penalty carry no retained data; with the
  # penalty held their coefficients are exactly zero, so they are left out of
  # the refit (mgcv's starting values fail on a penalised block with no data).
  nodata <- unlist(lapply(which(hold), function(j) setdiff(st$pen_cols[[j]], drop)))
  list(drop = sort(unique(c(drop, nodata))), drop_rank = drop, affected = affected,
       hold = hold, rank_deficient = length(drop) > 0L)
}

no_plan <- function(st) {
  list(drop = integer(0), drop_rank = integer(0), affected = integer(0),
       hold = rep(FALSE, length(st$S)), rank_deficient = FALSE)
}

# refit(): fit the frozen model to rows `keep`.
#   lambda = NULL : estimate the estimated smoothing parameters by REML
#                   (user-fixed and plan-held ones stay at the full-fit value)
#   lambda = vec  : hold every smoothing parameter at vec
#   psi    = NULL : re-estimate theta (nb() models whose theta was estimated)
#   psi    = num  : hold theta
#   start  = "default" : mgcv's own starting values for the REML search
#   start  = "full"    : start the search at the full-fit lambda (and theta)
# phi is never passed. Returns coefficients on the full frozen column set
# (dropped columns set to 0), predictions on ALL rows, and hyperparameters in
# the frozen penalty order.
refit <- function(st, keep, lambda = NULL, psi = NULL, plan = no_plan(st),
                  start = c("default", "full")) {
  start <- match.arg(start)
  warm  <- is.null(lambda) && start == "full"
  p    <- ncol(st$X)
  cols <- setdiff(seq_len(p), plan$drop)
  lam_in <- if (is.null(lambda)) {
    ifelse(st$map$estimated & !plan$hold, -1, st$lambda)
  } else lambda

  S    <- lapply(st$S, function(M) M[cols, cols, drop = FALSE])
  live <- vapply(S, function(M) any(M != 0), logical(1))

  family <- st$model$family                         # user-fixed theta, link kept
  if (!is.null(st$psi)) {                           # theta estimated in the full fit
    # nb() and negbin() read `link` with substitute(), so it must arrive as a
    # literal string: do.call() passes the value, not the expression st$link
    # (mgcv >= 1.9-4 rejects the expression).
    family <- if (!is.null(psi)) {
      do.call(mgcv::negbin, list(theta = psi, link = st$link))
    } else {                                        # a negative theta is a starting value
      do.call(mgcv::nb, list(theta = if (warm) -st$psi else NULL, link = st$link))
    }
  }

  dd <- list(y = st$y[keep], X = st$X[keep, cols, drop = FALSE],
             off = st$off[keep], pw = st$pw[keep])
  fit <- if (any(live)) {
    pp <- S[live]
    pp$sp <- lam_in[live]
    # warm start: begin the REML search at the full-fit values (in.out takes the
    # free smoothing parameters only, plus the scale)
    in_out <- if (warm) list(sp = st$lambda[live][lam_in[live] < 0],
                             scale = if (isTRUE(st$model$scale.estimated)) st$phi_reml else 1)
    mgcv::gam(y ~ X - 1 + offset(off), data = dd, weights = pw, family = family,
              method = "REML", paraPen = list(X = pp), in.out = in_out)
  } else {
    mgcv::gam(y ~ X - 1 + offset(off), data = dd, weights = pw, family = family,
              method = "REML")
  }

  beta <- numeric(p)
  beta[cols] <- stats::coef(fit)
  lam_out <- lam_in
  if (any(live)) {
    fs <- if (length(fit$full.sp)) fit$full.sp else fit$sp
    lam_out[live] <- unname(fs)
  }
  unset <- lam_out < 0                               # penalty dropped entirely
  lam_out[unset] <- st$lambda[unset]

  outer <- fit$outer.info$conv
  list(fit = fit, beta = beta,
       eta_all = drop(st$X %*% beta) + st$off,
       lambda  = unname(lam_out),
       # held theta: the value passed in (negbin()'s getTheta() takes no argument)
       psi     = if (is.null(st$psi)) NULL else if (!is.null(psi)) psi else
                   fit$family$getTheta(TRUE),
       phi     = fit$reml.scale,
       inner_converged = isTRUE(fit$converged),
       outer_status    = if (is.null(outer)) NA_character_ else as.character(outer),
       converged = isTRUE(fit$converged) &&
                   (is.null(outer) || identical(outer, "full convergence")))
}

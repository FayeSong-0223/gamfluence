# The frozen representation.
#
# gamfluence defines deletion influence conditional on the full-data
# representation: the constrained model matrix, the identifiability
# constraints and the penalty matrices (with mgcv's scaling) are taken from
# the full fit and never rebuilt. Only rows and hyperparameters change.
#
# Hyperparameters:
#   lambda  smoothing parameters mgcv estimated; user-fixed ones are held at
#           the user's value in every fit (map$estimated == FALSE).
#   psi     theta of nb(), only when mgcv estimated it (n.theta > 0).
#   phi     scale; never passed to a refit (beta-hat given lambda does not
#           depend on it). phi_reml (reml.scale) gives variance components;
#           phi_p (sig2) and edf are kept for the tests' Cook's distance check.

frozen_setup <- function(model, check = TRUE) {
  sc <- check_scope(model)

  # Model matrix on the rows mgcv used. With na.action = na.exclude,
  # predict() would pad the excluded rows back in as NA, so drop na.action.
  m_rows <- model
  m_rows$na.action <- NULL
  X <- stats::predict(m_rows, type = "lpmatrix")
  if (nrow(X) != length(model$y))
    stop("the model matrix and the response have different numbers of rows.",
         call. = FALSE)
  p  <- ncol(X)
  sp_all <- if (length(model$full.sp)) model$full.sp else model$sp

  S <- list(); cols <- list()
  term <- character(); sp_index <- integer(); estimated <- logical(); is_re <- logical()
  for (sm in model$smooth) {
    i <- sm$first.para:sm$last.para
    M <- matrix(0, p, p)
    M[i, i] <- sm$S[[1]]
    S[[length(S) + 1L]] <- M
    cols[[length(cols) + 1L]] <- i
    term      <- c(term, sm$label)
    sp_index  <- c(sp_index, sm$first.sp)
    estimated <- c(estimated, is.null(sm$sp) || sm$sp[1] < 0)
    is_re     <- c(is_re, inherits(sm, "random.effect"))
  }
  map <- data.frame(term = term, sp_index = sp_index, estimated = estimated,
                    is_re = is_re, stringsAsFactors = FALSE)

  st <- list(model = model, X = X, S = S, pen_cols = cols, map = map,
             S_half = make_S_half(S, p), cols = seq_len(p),
             lambda   = unname(sp_all[map$sp_index]),
             psi      = if (sc$theta_estimated) model$family$getTheta(TRUE) else NULL,
             is_nb    = sc$is_nb,
             link     = model$family$link,
             phi_reml = model$reml.scale,
             phi_p    = model$sig2,
             edf      = sum(model$edf),
             off      = if (is.null(model$offset)) rep(0, length(model$y)) else as.numeric(model$offset),
             pw       = model$prior.weights,
             w        = model$weights,     # mgcv's IRLS (Fisher) weights of the full fit
             y        = model$y,
             eta0     = model$linear.predictors,
             dropped_at_setup = character(0))

  # The full model must be identifiable under the frozen representation, so
  # that any rank deficiency found after a deletion is caused by the deletion.
  # If mgcv resolved a collinearity by setting coefficients to exactly zero,
  # those columns are removed once, here.
  all_rows <- seq_along(st$y)
  if (deletion_plan(st, all_rows)$rank_deficient) {
    zero <- which(stats::coef(model) == 0)
    if (length(zero)) {
      st$dropped_at_setup <- colnames(X)[zero]
      st <- drop_columns(st, zero)
    }
    if (!length(zero) || deletion_plan(st, all_rows)$rank_deficient)
      stop("the fitted model is not identifiable (its design is rank deficient); ",
           "gamfluence needs an identifiable model.", call. = FALSE)
  }

  if (check) self_check(st)
  st
}

# A square root of the summed penalty (all lambda = 1). The rank of
# rbind(X_kept, S_half) does not depend on the (positive) lambda values.
make_S_half <- function(S, p) {
  if (!length(S)) return(matrix(0, 0, p))
  e <- eigen(Reduce(`+`, S), symmetric = TRUE)
  pos <- e$values > max(e$values, 0) * 1e-10
  if (!any(pos)) return(matrix(0, 0, p))
  t(e$vectors[, pos, drop = FALSE] %*% diag(sqrt(e$values[pos]), sum(pos), sum(pos)))
}

drop_columns <- function(st, drop) {
  keep <- setdiff(seq_len(ncol(st$X)), drop)
  st$X <- st$X[, keep, drop = FALSE]
  st$S <- lapply(st$S, function(M) M[keep, keep, drop = FALSE])
  st$pen_cols <- lapply(st$pen_cols, function(i) which(keep %in% i))
  st$S_half <- make_S_half(st$S, ncol(st$X))
  st$cols <- st$cols[keep]
  st
}

# The frozen representation must reproduce the fit before any deletion is
# trusted: an error if holding the hyperparameters does not give back the
# fitted linear predictor, a warning if re-estimating them does not.
self_check <- function(st) {
  all_rows <- seq_along(st$y)
  scale <- 1 + max(abs(st$eta0))
  h <- refit(st, all_rows, lambda = st$lambda, psi = st$psi)
  d <- max(abs(h$eta_all - st$eta0))
  if (!is.finite(d) || d > 1e-6 * scale)
    stop("the frozen representation does not reproduce the fitted model ",
         "(max difference in the linear predictor ", signif(d, 3), "); this model ",
         "is probably outside what gamfluence supports.", call. = FALSE)
  if (any(st$map$estimated) || !is.null(st$psi)) {
    r <- reestimate(st, all_rows)
    d <- max(abs(r$eta_all - st$eta0))
    if (!is.finite(d) || d > 1e-5 * scale)
      warning("re-estimating the hyperparameters on all rows does not reproduce ",
              "the fit (max difference ", signif(d, 3), "); the `added` changes ",
              "may include optimiser differences.", call. = FALSE)
  }
  invisible(TRUE)
}

# Variance component of an i.i.d. random-effect term:
# sd = sqrt(phi_REML / (lambda * s)), where the term's penalty is S = s I.
re_sd <- function(lambda, phi, st, j) {
  sm <- st$model$smooth[[j]]
  Sj <- sm$S[[1]]
  s  <- Sj[1, 1]
  if (max(abs(Sj - s * diag(nrow(Sj)))) > 1e-10)
    stop("random-effect penalty is not a multiple of the identity: ", sm$label,
         call. = FALSE)
  unname(sqrt(phi / (lambda[j] * s)))
}

# Weighted RMS of a change in the linear predictor.
wrms <- function(delta, w) sqrt(sum(w * delta^2) / sum(w))

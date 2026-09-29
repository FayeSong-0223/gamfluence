# An independent REML / Laplace-approximate REML criterion for the frozen
# representation, written from the formulas (Wood 2011) without calling mgcv's
# fitting code. Used by T24 to check that the frozen-basis refit's smoothing
# parameters are optimal on the reduced data.
#
# Scope: gaussian (identity), poisson (log), binomial (logit); penalties on
# disjoint column blocks (true for everything in the v0.1 scope).

reml_problem <- function(st, keep) {
  pen_rank <- vapply(st$S, function(M) {
    ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    sum(ev > max(ev) * 1e-10)
  }, numeric(1))
  list(X = st$X[keep, , drop = FALSE], y = st$y[keep], pw = st$pw[keep],
       off = st$off[keep], S = st$S, rank = pen_rank, est = st$map$estimated,
       lambda0 = st$lambda,
       family = if (st$model$family$family == "gaussian") "gaussian" else st$model$family$family)
}

# penalised fit for given lambda: returns beta, eta, W, H = X'WX + S_lambda
pirls <- function(pr, lambda, beta = NULL) {
  X <- pr$X; y <- pr$y; pw <- pr$pw; off <- pr$off
  Sl <- Reduce(`+`, Map(`*`, lambda, pr$S))
  if (pr$family == "gaussian") {
    H <- crossprod(X, pw * X) + Sl
    R <- chol(H)
    beta <- backsolve(R, forwardsolve(t(R), crossprod(X, pw * (y - off))))
    return(list(beta = drop(beta), eta = drop(X %*% beta) + off, W = pw, H = H, Sl = Sl))
  }
  if (is.null(beta)) beta <- numeric(ncol(X))
  for (it in 1:200) {
    eta <- drop(X %*% beta) + off
    mu  <- if (pr$family == "poisson") exp(eta) else plogis(eta)
    v   <- if (pr$family == "poisson") mu else mu * (1 - mu)   # canonical: dmu/deta = V(mu)
    W   <- pw * v
    z   <- eta - off + (y - mu) / v
    H   <- crossprod(X, W * X) + Sl
    R   <- chol(H)
    new <- drop(backsolve(R, forwardsolve(t(R), crossprod(X, W * z))))
    if (max(abs(new - beta)) < 1e-11 * (1 + max(abs(new)))) { beta <- new; break }
    beta <- new
  }
  eta <- drop(X %*% beta) + off
  mu  <- if (pr$family == "poisson") exp(eta) else plogis(eta)
  W   <- pw * (if (pr$family == "poisson") mu else mu * (1 - mu))
  list(beta = beta, eta = eta, W = W, H = crossprod(X, W * X) + Sl, Sl = Sl)
}

# criterion to MINIMISE, as a function of the full lambda vector
reml_value <- function(pr, lambda, fit = pirls(pr, lambda)) {
  ldH <- 2 * sum(log(diag(chol(fit$H))))
  ldS <- sum(pr$rank * log(lambda))              # log|S_lambda|_+ up to a constant
  pen <- sum(fit$beta * drop(fit$Sl %*% fit$beta))
  if (pr$family == "gaussian") {
    Mp <- ncol(pr$X) - sum(pr$rank)
    dp <- sum(pr$pw * (pr$y - fit$eta)^2) + pen
    phi <- dp / (length(pr$y) - Mp)
    return((length(pr$y) - Mp) / 2 * log(phi) + ldH / 2 - ldS / 2)
  }
  ll <- if (pr$family == "poisson") {
    sum(pr$pw * (pr$y * fit$eta - exp(fit$eta)))
  } else {
    mu <- plogis(fit$eta)
    sum(pr$pw * (pr$y * log(mu) + (1 - pr$y) * log(1 - mu)))
  }
  -ll + pen / 2 + ldH / 2 - ldS / 2
}

# optimise over log lambda for the estimated smoothing parameters.
# log lambda is confined to +-12 around the full-fit value: REML is flat as a
# smooth approaches its null space (lambda -> Inf), so the optimum can sit on
# that boundary; fitted values there are indistinguishable.
reml_optimise <- function(pr, start, width = 12) {
  centre <- log(pr$lambda0[pr$est])
  full <- function(rho) { l <- pr$lambda0; l[pr$est] <- exp(rho); l }
  f <- function(rho) {
    rho <- pmin(pmax(rho, centre - width), centre + width)
    v <- tryCatch(reml_value(pr, full(rho)), error = function(e) Inf)
    if (is.finite(v)) v else 1e100
  }
  o <- stats::optim(start, f, method = "BFGS",
                    control = list(reltol = 1e-14, maxit = 1000, ndeps = rep(1e-5, length(start))))
  if (length(start) > 1)
    o <- stats::optim(o$par, f, method = "Nelder-Mead", control = list(reltol = 1e-15, maxit = 5000))
  rho <- pmin(pmax(o$par, centre - width), centre + width)
  list(lambda = full(rho), value = f(rho))
}

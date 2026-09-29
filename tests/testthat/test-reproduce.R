# The frozen representation reproduces what mgcv fitted (full data).
# These are machinery checks: they say nothing about deletion.

all_rows <- function(st) seq_along(st$y)

test_that("T1: hyperparameters held -> the original fit is reproduced", {
  for (nm in names(reef_setups)) {
    st <- reef_setups[[nm]]
    f  <- gamfluence:::refit(st, all_rows(st), lambda = st$lambda, psi = st$psi)
    expect_lt(mx(f$eta_all, st$eta0), 1e-8)
    expect_true(f$converged, info = nm)
  }
})

test_that("T2: REML re-estimation on all rows recovers lambda, theta, phi_REML", {
  for (nm in names(reef_setups)) {
    st <- reef_setups[[nm]]
    f  <- gamfluence:::refit(st, all_rows(st))
    expect_lt(mx(log(f$lambda), log(st$lambda)), 1e-5)
    expect_lt(abs(log(f$phi / st$phi_reml)), 1e-6)
    if (!is.null(st$psi)) expect_lt(abs(log(f$psi / st$psi)), 1e-6)
  }
})

test_that("T3: lambda is mapped through full.sp and first.sp when some are user-fixed", {
  m1 <- reef_models$fixed_sp_in_s; s1 <- reef_setups$fixed_sp_in_s
  m2 <- reef_models$fixed_sp_gam_arg; s2 <- reef_setups$fixed_sp_gam_arg
  expect_length(m1$sp, 3); expect_length(s1$S, 4)
  expect_equal(s1$map$estimated, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(s1$lambda[1], 0.5)
  expect_length(m2$sp, 1); expect_length(s2$S, 2)
  expect_equal(s2$map$estimated, c(FALSE, TRUE))
  expect_equal(s2$lambda[1], 0.5)
  # the user-fixed value stays fixed in a re-estimating refit
  f <- gamfluence:::refit(s1, drop_site(reef, "3"))
  expect_equal(f$lambda[1], 0.5)
})

test_that("T4: full-fit variance component equals gam.vcomp", {
  vc <- function(m) { utils::capture.output(v <- gam.vcomp(m, rescale = FALSE)); if (is.list(v)) v$vc else v }
  for (nm in c("gaussian_w_by", "fixed_sp_in_s", "fixed_sp_gam_arg", "poisson_off", "binomial_cbind")) {
    st <- reef_setups[[nm]]; j <- which(st$map$term == "s(site)")
    expect_equal(gamfluence:::re_sd(st$lambda, st$phi_reml, st, j),
                 unname(vc(st$model)["s(site)", "std.dev"]), tolerance = 1e-6, info = nm)
  }
  # sig2 is not interchangeable with reml.scale when a smoothing parameter is fixed
  st <- reef_setups$fixed_sp_in_s
  expect_gt(abs(st$phi_p - st$phi_reml) / st$phi_reml, 0.01)
})

test_that("T5: beta given lambda does not depend on phi; it does depend on theta", {
  st <- reef_setups$gaussian_w_by
  pp <- st$S; pp$sp <- st$lambda
  dd <- list(y = st$y, X = st$X, pw = st$pw)
  a <- gam(y ~ X - 1, data = dd, weights = pw, method = "REML", paraPen = list(X = pp))
  b <- gam(y ~ X - 1, data = dd, weights = pw, method = "REML", paraPen = list(X = pp),
           scale = 3 * st$phi_reml)
  expect_lt(mx(coef(a), coef(b)), 1e-10)
  sc <- reef_setups$negbin_nb; keep <- drop_site(reef, "7")
  held  <- gamfluence:::refit(sc, keep, lambda = sc$lambda, psi = sc$psi)
  moved <- gamfluence:::refit(sc, keep, lambda = sc$lambda, psi = NULL)
  expect_gt(mx(held$eta_all[keep], moved$eta_all[keep]), 1e-4)
})

test_that("T7: theta is re-estimated only when mgcv estimated it; the link is kept", {
  keep <- drop_site(reef, "7")
  for (nm in c("negbin_theta3", "nb_theta3")) {
    st <- reef_setups[[nm]]
    expect_null(st$psi)
    f <- gamfluence:::refit(st, keep)
    th <- tryCatch(f$fit$family$getTheta(TRUE), error = function(e) f$fit$family$getTheta())
    expect_equal(th, 3, info = nm)
  }
  st <- reef_setups$nb_sqrt_link
  f <- gamfluence:::refit(st, keep)
  expect_identical(f$fit$family$link, "sqrt")
  expect_false(isTRUE(all.equal(f$psi, st$psi)))
})

test_that("T8: prior weights reach the refits", {
  st <- reef_setups$gaussian_w_by
  s0 <- st; s0$pw <- rep(1, length(st$pw))
  expect_lt(mx(gamfluence:::refit(st, all_rows(st), st$lambda)$eta_all, st$eta0), 1e-8)
  expect_gt(mx(gamfluence:::refit(s0, all_rows(st), st$lambda)$eta_all, st$eta0), 1e-3)
})

test_that("the weights used are mgcv's Fisher weights of the full fit", {
  st <- reef_setups$poisson_off; mu <- exp(st$eta0)
  expect_equal(unname(st$w), unname(st$pw * mu), tolerance = 1e-8)
  st <- reef_setups$binomial_cbind; mu <- plogis(st$eta0)
  expect_equal(unname(st$w), unname(st$pw * mu * (1 - mu)), tolerance = 1e-8)
  st <- reef_setups$negbin_nb; mu <- exp(st$eta0)
  expect_equal(unname(st$w), unname(st$pw * mu / (1 + mu / st$psi)), tolerance = 1e-6)
})

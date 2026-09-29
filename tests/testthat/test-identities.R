# Regression tests for identities that hold by construction.
# They are not evidence that the method is right; they catch implementation
# slips (a fit swapped, a sign flipped, a column mislabelled) when the code
# changes.

test_that("fixed + added = total, and the fixed refit keeps the full-fit hyperparameters", {
  for (nm in c("gaussian_w_by", "negbin_nb", "poisson_off", "binomial_cbind", "fixed_sp_in_s")) {
    st <- reef_setups[[nm]]
    keep <- drop_site(reef, "4")
    r <- gamfluence:::compare_fits(st, keep, gamfluence:::deletion_plan(st, keep))
    expect_lt(mx(r$fixed + r$added, r$total), 1e-12)
    expect_equal(r$A$lambda, st$lambda)
    if (!is.null(st$psi)) expect_equal(r$A$psi, st$psi, tolerance = 1e-10)
  }
})

test_that("nothing to re-estimate: the added change is exactly zero", {
  set.seed(3)
  d <- data.frame(site = factor(rep(1:10, each = 4)), z = runif(40))
  d$y <- sin(2 * pi * d$z) + rnorm(40, 0, 0.3)
  m <- gam(y ~ s(z, k = 6), data = d, method = "REML", sp = 0.1)   # user-fixed lambda only
  i <- gam_influence(m, cluster = d$site, progress = FALSE)
  expect_true(all(i$sites$added == 0))
  expect_equal(i$sites$fixed, i$sites$total)
  expect_null(i$hyper)
})

test_that("reported magnitudes match the refit vectors", {
  m  <- reef_models$poisson_off
  st <- reef_setups$poisson_off
  sites <- c("2", "5")
  i <- gam_influence(m, cluster = "site", which = sites, progress = FALSE)
  for (k in seq_along(sites)) {
    keep <- drop_site(reef, sites[k])
    r <- gamfluence:::compare_fits(st, keep)
    w <- st$w[keep]
    expect_equal(i$sites$fixed[k], gamfluence:::wrms(r$fixed, w))
    expect_equal(i$sites$added[k], gamfluence:::wrms(r$added, w))
    expect_equal(i$sites$total[k], gamfluence:::wrms(r$total, w))
  }
})

test_that("hyperparameter table: full values match the fit; deleted values match refit T", {
  m  <- reef_models$negbin_nb
  st <- reef_setups$negbin_nb
  i  <- gam_influence(m, cluster = "site", which = "3", progress = FALSE)
  h  <- i$hyper
  j  <- which(st$map$term == "s(site)")
  keep <- drop_site(reef, "3")
  r <- gamfluence:::compare_fits(st, keep)
  expect_equal(h$full[h$parameter == "s(site)"], gamfluence:::re_sd(st$lambda, st$phi_reml, st, j))
  expect_equal(h$deleted[h$parameter == "s(site)"], gamfluence:::re_sd(r$T$lambda, r$T$phi, st, j))
  expect_equal(h$full[h$parameter == "s(x)"], st$lambda[st$map$term == "s(x)"])
  expect_equal(h$deleted[h$parameter == "theta"], r$T$psi)
  expect_equal(h$ratio, h$deleted / h$full)
})

# T24: an independent reference for the re-estimating deletion fit (T) in
# models WITH covariate smooths. The REML (Gaussian) or Laplace REML
# (Poisson, binomial) criterion is coded from the formulas in helper-reml.R,
# evaluated on the frozen representation of the reduced data, and optimised
# with stats::optim. The refit's smoothing parameters must be optimal under
# that criterion, and the fitted values must agree.
#
# Tolerances: for the Gaussian model one by = factor smooth is nearly linear,
# so REML is almost flat as its lambda -> Inf; the reference optimiser runs to
# its bound there, mgcv stops earlier, and the criterion values differ by
# ~1e-5 with fitted values equal to ~1e-5. The GLM criteria have no flat
# direction here and agree to ~1e-12.

t24_check <- function(nm, sites, tol_v, tol_eta, tol_lambda) {
  st <- reef_setups[[nm]]
  for (s in sites) {
    keep <- if (s == "none") seq_along(st$y) else drop_site(reef, s)
    pr <- reml_problem(st, keep)
    fz <- if (s == "none") list(lambda = st$lambda, eta_all = st$eta0) else gamfluence:::refit(st, keep)
    o  <- reml_optimise(pr, log(st$lambda[pr$est]) + 0.7)
    gap <- reml_value(pr, fz$lambda) - o$value
    expect_lt(gap, tol_v)
    expect_gt(gap, -tol_v)                      # the reference found nothing much better or worse
    expect_lt(mx(pirls(pr, o$lambda)$eta, fz$eta_all[keep]), tol_eta)
    if (!is.null(tol_lambda))
      expect_lt(mx(log(fz$lambda[pr$est]), log(o$lambda[pr$est])), tol_lambda)
  }
}

test_that("T24: the criterion reproduces mgcv's full-data optimum", {
  t24_check("gaussian_w_by", "none", tol_v = 2e-4, tol_eta = 1e-4, tol_lambda = NULL)
  t24_check("poisson_off",   "none", tol_v = 1e-8, tol_eta = 1e-5, tol_lambda = 1e-4)
  t24_check("binomial_cbind", "none", tol_v = 1e-8, tol_eta = 1e-5, tol_lambda = 1e-4)
})

test_that("T24: deletion refits are optimal under the independent criterion", {
  sites <- c("2", "9", "17")
  t24_check("gaussian_w_by",  sites, tol_v = 2e-4, tol_eta = 1e-4, tol_lambda = NULL)
  # the well-determined Gaussian smoothing parameters (s(x) and the site term)
  st <- reef_setups$gaussian_w_by; wd <- match(c("s(x)", "s(site)"), st$map$term)
  for (s in sites) {
    keep <- drop_site(reef, s); pr <- reml_problem(st, keep)
    o <- reml_optimise(pr, log(st$lambda[pr$est]) + 0.7)
    expect_lt(mx(log(gamfluence:::refit(st, keep)$lambda[wd]), log(o$lambda[wd])), 1e-4)
  }
  # a mix of user-fixed and estimated smoothing parameters
  t24_check("fixed_sp_gam_arg", sites, tol_v = 1e-8, tol_eta = 1e-5, tol_lambda = 1e-4)
  t24_check("fixed_sp_in_s",    sites, tol_v = 2e-4, tol_eta = 1e-4, tol_lambda = NULL)
  t24_check("poisson_off",    sites, tol_v = 1e-8, tol_eta = 1e-5, tol_lambda = 1e-4)
  t24_check("binomial_cbind", sites, tol_v = 1e-8, tol_eta = 1e-5, tol_lambda = 1e-4)
})

test_that("T24 has power: the full-fit lambda is detectably suboptimal after deletion", {
  st <- reef_setups$binomial_cbind
  gaps <- vapply(c("2", "9", "17"), function(s) {
    pr <- reml_problem(st, drop_site(reef, s))
    o  <- reml_optimise(pr, log(st$lambda[pr$est]) + 0.7)
    reml_value(pr, st$lambda) - o$value
  }, numeric(1))
  expect_true(all(gaps > 1e-3))
})

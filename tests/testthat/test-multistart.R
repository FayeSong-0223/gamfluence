# The re-estimated refit is searched for from two starting points and the
# better REML fit is kept (found on real data: a year smooth whose REML score
# had an interior optimum and a plateau towards a straight line).

fake_fit <- function(score, conv = TRUE) {
  list(fit = list(gcv.ubre = score), converged = conv, eta_all = rep(score, 3))
}

test_that("reestimate keeps the converged search with the lower REML score", {
  st <- gamfluence:::frozen_setup(reef_models$poisson_off)
  run <- function(a, b) with_mocked_bindings(
    gamfluence:::reestimate(st, 1:3),
    refit = function(st, keep, lambda = NULL, psi = NULL, plan, start = "default") {
      f <- if (start == "full") b else a
      if (is.null(f)) stop("search failed") else f
    })
  r <- run(fake_fit(10), fake_fit(9))                 # the start at the full fit wins
  expect_equal(r$fit$gcv.ubre, 9)
  expect_equal(r$alt$reml_gap, 1)
  expect_equal(r$alt$eta_all, rep(10, 3))
  r <- run(fake_fit(8), fake_fit(9))                  # mgcv's default start wins
  expect_equal(r$fit$gcv.ubre, 8)
  expect_equal(r$alt$eta_all, rep(9, 3))
  r <- run(fake_fit(8, conv = FALSE), fake_fit(9))    # only converged fits compete
  expect_equal(r$fit$gcv.ubre, 9)
  expect_null(r$alt)
  r <- run(fake_fit(5), NULL)                         # a failed second search is ignored
  expect_equal(r$fit$gcv.ubre, 5)
  expect_null(r$alt)
  r <- run(fake_fit(10), fake_fit(10 - 1e-9))         # a difference within noise keeps the default
  expect_equal(r$fit$gcv.ubre, 10)
})

test_that("a search started at the full-fit values reproduces the fit", {
  for (nm in c("gaussian_w_by", "negbin_nb", "poisson_off", "binomial_cbind", "nb_theta3")) {
    m  <- reef_models[[nm]]
    st <- gamfluence:::frozen_setup(m)
    r  <- gamfluence:::refit(st, seq_along(st$y), start = "full")
    expect_true(r$converged, info = nm)
    expect_lt(max(abs(r$eta_all - st$eta0)), 1e-5 * (1 + max(abs(st$eta0))))
    if (!is.null(st$psi)) expect_equal(r$psi, st$psi, tolerance = 1e-4, info = nm)
  }
})

test_that("two different optima are flagged, and near-identical ones are not", {
  m <- reef_models$poisson_off
  with_alt <- function(shift) function(st, keep, plan = gamfluence:::no_plan(st)) {
    r <- gamfluence:::refit(st, keep, plan = plan)
    r$alt <- list(eta_all = r$eta_all + shift, reml_gap = 0.03)
    r
  }
  i <- with_mocked_bindings(
    gam_influence(m, cluster = "site", which = c("1", "2"), progress = FALSE),
    reestimate = with_alt(0.5))
  expect_true(all(grepl("two optima", i$sites$status)))
  expect_true(all(grepl("two REML optima", i$sites$note)))
  expect_true(all(grepl("lower by 0.03", i$sites$note)))
  i <- with_mocked_bindings(
    gam_influence(m, cluster = "site", which = c("1", "2"), progress = FALSE),
    reestimate = with_alt(1e-9))
  expect_true(all(i$sites$status == "ok"))
})

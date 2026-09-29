# Structure of the returned object and the methods.

test_that("gam_influence returns the documented pieces", {
  i <- gam_influence(reef_models$negbin_nb, cluster = "site", which = c("1", "2", "3"),
                     progress = FALSE)
  expect_s3_class(i, "gam_influence")
  expect_named(i, c("sites", "hyper", "info"))
  expect_named(i$sites, c("id", "n_obs", "fixed", "added", "total", "status", "note"))
  expect_equal(nrow(i$sites), 3)
  expect_true(all(i$sites$status == "ok"))
  expect_named(i$hyper, c("id", "parameter", "type", "full", "deleted", "ratio"))
  expect_equal(sort(unique(i$hyper$parameter)), sort(c("s(x)", "s(site)", "theta")))
  expect_equal(unique(i$hyper$type[i$hyper$parameter == "s(site)"]), "sd")
  expect_true(all(i$sites$total >= 0))
})

test_that("theta appears only when it was estimated; user-fixed lambda never appears", {
  i <- gam_influence(reef_models$nb_theta3, cluster = "site", which = "1", progress = FALSE)
  expect_false("theta" %in% i$hyper$parameter)
  i <- gam_influence(reef_models$fixed_sp_in_s, cluster = "site", which = "1", progress = FALSE)
  expect_false("s(x)" %in% i$hyper$parameter)          # s(x, sp = 0.5) is user-fixed
})

test_that("print, summary, as.data.frame and top_units work", {
  i <- gam_influence(reef_models$poisson_off, cluster = "site", which = as.character(1:5),
                     progress = FALSE)
  expect_output(print(i), "fixed vs re-estimated")
  expect_output(print(i), "No null reference")
  expect_output(print(summary(i)), "gamfluence summary")
  expect_identical(as.data.frame(i), i$sites)
  top <- top_units(i, by = "added", n = 2)
  expect_equal(nrow(top), 2)
  expect_true(top$added[1] >= top$added[2])
  expect_error(top_units(i, by = "nope"), "must be")
})

test_that("a failing refit is recorded, not fatal", {
  m <- reef_models$poisson_off
  expect_warning(
    i <- with_mocked_bindings(
      gam_influence(m, cluster = "site", which = c("1", "2"), progress = FALSE),
      compare_fits = function(st, keep, plan) stop("simulated failure")),
    "did not converge")
  expect_identical(unique(i$sites$status), "error")
  expect_identical(unique(i$sites$note), "simulated failure")
  expect_true(all(is.na(i$sites$total)))
})

test_that("cluster is required and must have more than one level", {
  m <- reef_models$poisson_off
  expect_error(gam_influence(m, progress = FALSE), "`cluster` is required")
  expect_error(gam_influence(m, cluster = rep("a", nrow(m$model)), progress = FALSE), "only one level")
})

test_that("theta-only hyperparameter table (nb() with user-fixed lambda or no smooths)", {
  for (f in list(yc ~ mgmt + s(x, sp = 0.5) + offset(log(area)), yc ~ mgmt + x + offset(log(area)))) {
    m <- gam(f, data = reef, family = nb(), method = "REML")
    i <- gam_influence(m, cluster = "site", data = reef, which = c("1", "2"), progress = FALSE)
    expect_true(all(i$sites$status == "ok"))
    expect_identical(i$hyper$parameter, c("theta", "theta"))
    expect_equal(i$hyper$ratio, i$hyper$deleted / i$hyper$full)
    expect_false(anyNA(i$hyper$deleted))
  }
})

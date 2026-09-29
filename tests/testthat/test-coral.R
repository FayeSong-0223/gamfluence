# The shipped coral trout data reproduce the study's Model A, and a real
# negative binomial model with estimated theta runs through gam_influence().

coral_form <- count ~ REGION * NTR + EXPOSURE + depth + s(YEAR, by = REGION, k = 5) +
  s(rugosity, k = 5) + s(LHC, k = 5) + s(kd490, k = 5) + s(maxDHW, k = 5) +
  s(Cyclone, k = 5) + s(SITE, bs = "re")

test_that("coral_trout has the documented shape and reproduces Model A", {
  expect_equal(dim(coral_trout), c(467L, 12L))
  expect_equal(nlevels(coral_trout$SITE), 71L)
  expect_false(anyNA(coral_trout))
  expect_true(all(coral_trout$count >= 0 & coral_trout$count == round(coral_trout$count)))
  m <- gam(coral_form, family = nb(), data = coral_trout, method = "REML")
  expect_equal(m$family$getTheta(TRUE), 5.96007, tolerance = 1e-4)   # the study's log: 5.96
  expect_equal(as.numeric(m$gcv.ubre), 1197.39, tolerance = 1e-5)     # the study's log: 1197.4
})

test_that("gam_influence runs on the coral trout model", {
  skip_on_cran()                                                      # about 10 seconds
  m <- gam(coral_form, family = nb(), data = coral_trout, method = "REML")
  i <- gam_influence(m, cluster = "SITE", which = c("OE4", "B2"), progress = FALSE)
  expect_setequal(i$sites$id, c("OE4", "B2"))
  expect_false(any(grepl("error|not converged", i$sites$status)))
  expect_true("theta" %in% i$hyper$parameter)
  expect_true(all(i$sites$total > 0))
})

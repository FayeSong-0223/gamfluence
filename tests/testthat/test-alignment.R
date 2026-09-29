# Row alignment: units must be defined on the rows mgcv used, not on the rows
# of the data frame, when na.action or subset removed rows.

align_data <- function() {
  d <- make_reef(n_site = 12, n_year = 6, keep_frac = 1, seed = 11)
  d$transect <- factor(paste0("T", (as.integer(d$site) + d$year) %% 5))  # not in the formula
  d$x[c(3, 17, 40)] <- NA                                               # dropped by na.action
  d
}

test_that("NA rows and subset: cluster from `data` matches a fit on the clean rows", {
  d  <- align_data()
  f  <- yg ~ s(x) + s(site, bs = "re")
  m1 <- gam(f, data = d, weights = w, method = "REML", subset = year != 2)
  clean <- d[!is.na(d$x) & d$year != 2, ]
  m2 <- gam(f, data = clean, weights = w, method = "REML")
  expect_equal(unname(m1$linear.predictors), unname(m2$linear.predictors), tolerance = 1e-10)

  i1 <- gam_influence(m1, data = d, cluster = "transect", progress = FALSE)
  i2 <- gam_influence(m2, data = clean, cluster = "transect", progress = FALSE)
  expect_identical(i1$sites$id, i2$sites$id)
  expect_identical(i1$sites$n_obs, i2$sites$n_obs)
  expect_equal(i1$sites$total, i2$sites$total, tolerance = 1e-6)
  expect_equal(i1$sites$added, i2$sites$added, tolerance = 1e-6)
})

test_that("cluster as a vector: data-length and model-length both align", {
  d  <- align_data()
  m  <- gam(yg ~ s(x) + s(site, bs = "re"), data = d, weights = w, method = "REML")
  by_name <- gam_influence(m, data = d, cluster = "transect", progress = FALSE)
  by_data_vec <- gam_influence(m, data = d, cluster = d$transect, progress = FALSE)
  rows <- as.integer(rownames(m$model))
  by_model_vec <- gam_influence(m, cluster = d$transect[rows], progress = FALSE)
  expect_equal(by_name$sites, by_data_vec$sites)
  expect_equal(by_name$sites, by_model_vec$sites)
  expect_error(gam_influence(m, data = d, cluster = d$transect[-1], progress = FALSE),
               "must match either")
})

test_that("cluster variable inside the model frame needs no data", {
  m <- reef_models$poisson_off
  expect_no_error(gam_influence(m, cluster = "site", which = c("1", "2"), progress = FALSE))
})

test_that("modified or re-sorted data is detected", {
  d <- align_data()
  m <- gam(yg ~ s(x) + s(site, bs = "re"), data = d, weights = w, method = "REML")
  shuffled <- d[sample(nrow(d)), ]            # row names kept: still aligns
  expect_equal(gam_influence(m, data = shuffled, cluster = "transect", progress = FALSE)$sites,
               gam_influence(m, data = d, cluster = "transect", progress = FALSE)$sites)
  renumbered <- shuffled; rownames(renumbered) <- NULL
  expect_error(gam_influence(m, data = renumbered, cluster = "transect", progress = FALSE),
               "does not match")
  changed <- d; changed$yg <- changed$yg + 1
  expect_error(gam_influence(m, data = changed, cluster = "transect", progress = FALSE),
               "does not match")
  expect_error(gam_influence(m, data = d[1:10, ], cluster = "transect", progress = FALSE),
               "cannot all be found")
})

test_that("data recovered from the call, with a message", {
  dd <- align_data()
  m <- gam(yg ~ s(x) + s(site, bs = "re"), data = dd, weights = w, method = "REML")
  expect_message(gam_influence(m, cluster = "transect", which = "T1", progress = FALSE),
                 "from the model call")
})

test_that("missing cluster values and bad `which` are errors", {
  d <- align_data()
  d$transect[5] <- NA
  m <- gam(yg ~ s(x) + s(site, bs = "re"), data = d, weights = w, method = "REML")
  expect_error(gam_influence(m, data = d, cluster = "transect", progress = FALSE), "missing")
  expect_error(gam_influence(m, cluster = "site", which = "nope", progress = FALSE), "unknown site")
})

test_that("transformed terms are checked too: re-sorted data with reset row names is caught", {
  d <- align_data()
  d$ypos <- d$yc + 1
  m <- gam(log(ypos) ~ s(sqrt(x)) + log(area), data = d, method = "REML")
  expect_length(intersect(names(m$model), names(d)), 0)   # no plain column to compare
  ok <- gam_influence(m, data = d, cluster = "transect", progress = FALSE)
  resorted <- d[order(d$yg), ]; rownames(resorted) <- NULL
  expect_error(gam_influence(m, data = resorted, cluster = "transect", progress = FALSE),
               "does not match")
  kept <- d[order(d$yg), ]                                 # row names kept: fine
  expect_equal(gam_influence(m, data = kept, cluster = "transect", progress = FALSE)$sites,
               ok$sites)
})

test_that("na.action = na.exclude gives the same result as na.omit", {
  d <- align_data()
  m1 <- gam(yg ~ s(x) + s(site, bs = "re"), data = d, weights = w, method = "REML",
            na.action = na.omit)
  m2 <- gam(yg ~ s(x) + s(site, bs = "re"), data = d, weights = w, method = "REML",
            na.action = na.exclude)
  i1 <- gam_influence(m1, data = d, cluster = "transect", progress = FALSE)
  i2 <- gam_influence(m2, data = d, cluster = "transect", progress = FALSE)
  expect_equal(i1$sites, i2$sites)
})

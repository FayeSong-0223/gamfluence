# T6: the scope guard accepts the in-scope models and rejects every tested
# out-of-scope case before any refit.

test_that("all in-scope test models pass the guard", {
  for (nm in names(reef_models))
    expect_no_error(gamfluence:::check_scope(reef_models[[nm]]), message = nm)
})

test_that("out-of-scope models are rejected with a reason", {
  d  <- reef
  d6 <- transform(d, yl = log1p(yc))
  dz <- list(yg = d$yg, Z = cbind(d$x, d$x^2), year = d$year)
  cases <- list(
    "te(x, year)"          = list(quote(gam(yg ~ te(x, year), data = d, method = "REML")), "smooth class"),
    "te(x)"                = list(quote(gam(yg ~ te(x), data = d, method = "REML")), "smooth class"),
    "ti()"                 = list(quote(gam(yg ~ ti(x) + ti(year, k = 5), data = d, method = "REML")), "smooth class"),
    "t2()"                 = list(quote(gam(yg ~ t2(x, year, k = 4), data = d, method = "REML")), "smooth class"),
    "fs"                   = list(quote(gam(yg ~ s(x, site, bs = "fs", k = 4), data = d, method = "REML")), "smooth class"),
    "select = TRUE"        = list(quote(gam(yg ~ s(x), data = d, method = "REML", select = TRUE)), "more than one penalty"),
    "id ="                 = list(quote(gam(yg ~ s(x, id = 1, k = 5) + s(year, id = 1, k = 5), data = d, method = "REML")), "id ="),
    "bs = 'cr'"            = list(quote(gam(yg ~ s(x, bs = "cr"), data = d, method = "REML")), "basis"),
    "numeric by"           = list(quote(gam(yg ~ s(year, by = x, k = 5), data = d, method = "REML")), "numeric by"),
    "fx = TRUE"            = list(quote(gam(yg ~ s(x, k = 5, fx = TRUE), data = d, method = "REML")), "unpenalised"),
    "Tweedie"              = list(quote(gam(yc ~ s(x), data = d, family = tw(), method = "REML")), "family"),
    "Gamma"                = list(quote(gam(I(yc + 1) ~ s(x), data = d, family = Gamma(link = "log"), method = "REML")), "family"),
    "gaussian log link"    = list(quote(gam(I(yc + 1) ~ s(x), data = d, family = gaussian(link = "log"), method = "REML")), "link"),
    "ML"                   = list(quote(gam(yg ~ s(x), data = d, method = "ML")), "REML"),
    "GCV.Cp"               = list(quote(gam(yg ~ s(x), data = d)), "REML"),
    "gamma = 1.4"          = list(quote(gam(yg ~ s(x), data = d, method = "REML", gamma = 1.4)), "gamma"),
    "Gaussian scale = 1"   = list(quote(gam(yg ~ s(x), data = d, method = "REML", scale = 1)), "scale"),
    "Poisson scale = -1"   = list(quote(gam(yp ~ s(x), data = d, family = poisson, method = "REML", scale = -1)), "scale"),
    "min.sp"               = list(quote(gam(yg ~ s(x), data = d, method = "REML", min.sp = 1e-3)), "min.sp"),
    "H"                    = list(quote(gam(yg ~ s(x), data = d, method = "REML", H = diag(0.1, 10))), "H"),
    "paraPen"              = list(quote(gam(yg ~ Z + s(year, k = 5), data = dz, method = "REML",
                                            paraPen = list(Z = list(diag(2))))), "paraPen"),
    "bam()"                = list(quote(bam(yg ~ s(x), data = d, method = "REML")), "bam"),
    "gamm()$gam"           = list(quote(gamm(yl ~ s(x), random = list(site = ~1), data = d6,
                                             method = "REML")$gam), "gamm"),
    "lm()"                 = list(quote(lm(yg ~ x, data = d)), "mgcv::gam"),
    "glm()"                = list(quote(glm(yp ~ x, data = d, family = poisson)), "mgcv::gam"),
    "optimizer = 'efs'"    = list(quote(gam(yg ~ s(x), data = d, method = "REML", optimizer = "efs")), "optimizer"),
    "2-D thin plate"       = list(quote(gam(yg ~ s(x, year, k = 10), data = d, method = "REML")), "one-dimensional"),
    "thin plate m = 1"     = list(quote(gam(yg ~ s(x, m = 1), data = d, method = "REML")), "one-dimensional"),
    "re random slope"      = list(quote(gam(yg ~ s(x) + s(site, bs = "re") + s(site, year, bs = "re"),
                                            data = d, method = "REML")), "single factor"),
    "re two factors"       = list(quote(gam(yg ~ s(x) + s(site, region, bs = "re"), data = d,
                                            method = "REML")), "single factor"),
    "re numeric"           = list(quote(gam(yg ~ s(x) + s(year, bs = "re"), data = d,
                                            method = "REML")), "single factor"))
  expect_length(cases, 31)
  for (nm in names(cases)) {
    fit <- suppressWarnings(eval(cases[[nm]][[1]]))
    expect_error(gamfluence:::frozen_setup(fit), cases[[nm]][[2]], info = nm)
    expect_error(gam_influence(fit, cluster = "site", data = d, progress = FALSE),
                 cases[[nm]][[2]], info = nm)
  }
})

test_that("binomial response formats: 0/1, factor, proportion + weights, cbind", {
  d <- reef
  d$yf <- factor(ifelse(d$yb == 1, "yes", "no"), levels = c("no", "yes"))
  d$prop <- d$succ / d$ntrial
  m01 <- gam(yb ~ s(x) + s(site, bs = "re"), data = d, family = binomial, method = "REML")
  mf  <- gam(yf ~ s(x) + s(site, bs = "re"), data = d, family = binomial, method = "REML")
  mc  <- gam(cbind(succ, fail) ~ s(x) + s(site, bs = "re"), data = d, family = binomial, method = "REML")
  mp  <- suppressWarnings(gam(prop ~ s(x) + s(site, bs = "re"), data = d, weights = ntrial,
                              family = binomial, method = "REML"))
  for (m in list(m01, mf, mc, mp)) expect_no_error(gamfluence:::check_scope(m))
  # factor response is the same model as 0/1
  expect_equal(unname(m01$linear.predictors), unname(mf$linear.predictors), tolerance = 1e-10)
  # cbind and proportion + weights are the same model
  expect_equal(unname(mc$linear.predictors), unname(mp$linear.predictors), tolerance = 1e-6)
  expect_equal(range(mc$prior.weights), range(d$ntrial))
})

test_that("response guard rejects impossible responses", {
  m <- reef_models$gaussian_w_by
  bad <- m; bad$model <- m$model[-1, ]
  expect_error(gamfluence:::check_response(bad, "gaussian", FALSE), "different numbers of rows")
  pm <- reef_models$poisson_off
  bad <- pm; bad$y[1] <- -1
  expect_error(gamfluence:::check_response(bad, "poisson", FALSE), "negative")
  # a two-column response is only accepted for binomial (guard logic, on a modified object)
  mat <- m; mat$model[[1]] <- cbind(m$model[[1]], m$model[[1]])
  expect_error(gamfluence:::check_response(mat, "gaussian", FALSE), "matrix responses")
})

test_that("binomial response forms give identical influence results", {
  d <- reef
  d$yf <- factor(ifelse(d$yb == 1, "yes", "no"), levels = c("no", "yes"))
  d$yl <- d$yb == 1
  d$prop <- d$succ / d$ntrial
  run <- function(m) gam_influence(m, cluster = "site", which = c("2", "7"), progress = FALSE)$sites
  f01 <- run(gam(yb ~ s(x) + s(site, bs = "re"), data = d, family = binomial, method = "REML"))
  ffa <- run(gam(yf ~ s(x) + s(site, bs = "re"), data = d, family = binomial, method = "REML"))
  flg <- run(gam(yl ~ s(x) + s(site, bs = "re"), data = d, family = binomial, method = "REML"))
  fcb <- run(gam(cbind(succ, fail) ~ s(x) + s(site, bs = "re"), data = d, family = binomial, method = "REML"))
  fpr <- run(suppressWarnings(gam(prop ~ s(x) + s(site, bs = "re"), data = d, weights = ntrial,
                                  family = binomial, method = "REML")))
  for (v in c("fixed", "added", "total")) {
    expect_equal(ffa[[v]], f01[[v]], tolerance = 1e-8)
    expect_equal(flg[[v]], f01[[v]], tolerance = 1e-8)
    expect_equal(fpr[[v]], fcb[[v]], tolerance = 1e-4)
  }
})

# Rank deficiency after deletion.

rank_data <- function() {
  d <- make_reef(n_site = 16, n_year = 6, keep_frac = 1, seed = 21)
  # a factor level carried by one site only
  d$zone <- factor(ifelse(d$site == "3", "rare", ifelse(as.integer(d$site) %% 2 == 0, "a", "b")))
  # a by = factor level carried by one site only
  d$grp <- factor(ifelse(d$site == "5", "solo", "main"))
  d
}

test_that("a factor level held by one cluster: detected, retained rows still exact", {
  d <- rank_data()
  m <- gam(yg ~ zone + s(site, bs = "re"), data = d, weights = w, method = "REML")
  i <- gam_influence(m, cluster = "site", which = c("3", "4"), progress = FALSE)$sites
  expect_identical(i$status, c("rank-deficient", "ok"))
  expect_false(anyNA(i[c("fixed", "added", "total")]))

  # independent reference: ordinary gam() on the reduced data, level dropped
  st <- gamfluence:::frozen_setup(m)
  keep <- which(d$site != "3")
  plan <- gamfluence:::deletion_plan(st, keep)
  expect_identical(colnames(st$X)[plan$drop], "zonerare")
  nv <- gam(yg ~ zone + s(site, bs = "re"), data = droplevels(d[keep, ]), weights = w, method = "REML")
  tk <- gamfluence:::refit(st, keep, plan = plan)
  expect_lt(mx(nv$linear.predictors, tk$eta_all[keep]), 1e-5)
  ak <- gamfluence:::refit(st, keep, lambda = st$lambda, plan = plan)
  expect_true(ak$converged && tk$converged)
})

test_that("a by = factor smooth whose level is one cluster: lambda held, fit consistent", {
  d <- rank_data()
  m <- gam(yg ~ grp + s(x, by = grp) + s(site, bs = "re"), data = d, weights = w, method = "REML")
  st <- gamfluence:::frozen_setup(m)
  keep <- which(d$site != "5")
  plan <- gamfluence:::deletion_plan(st, keep)
  solo <- which(st$map$term == "s(x):grpsolo")
  expect_true(plan$hold[solo])
  expect_false(any(plan$hold[-solo]))
  expect_true(plan$rank_deficient)

  i <- gam_influence(m, cluster = "site", which = c("5", "6"), progress = FALSE)
  expect_identical(i$sites$status[1], "rank-deficient; lambda held")
  # site 6: nothing to drop or hold (its REML score has a second optimum, a
  # plateau for the solo level's smooth, which may be flagged as "two optima")
  expect_false(grepl("rank-deficient|lambda held|error|not converged", i$sites$status[2]))
  expect_identical(i$sites$note[1], "held: s(x):grpsolo")
  h5 <- i$hyper[i$hyper$id == "5" & i$hyper$parameter == "s(x):grpsolo", ]
  expect_true(is.na(h5$deleted) && is.na(h5$ratio))

  # consistency: same as refitting without the solo level's columns at all
  tk <- gamfluence:::refit(st, keep, plan = plan)
  gone <- c(which(colnames(st$X) == "grpsolo"), st$pen_cols[[solo]])
  Sr <- lapply(st$S[-solo], function(M) M[-gone, -gone]); Sr$sp <- rep(-1, length(Sr))
  ref <- gam(y ~ X - 1, data = list(y = st$y[keep], X = st$X[keep, -gone], pw = st$pw[keep]),
             weights = pw, method = "REML", paraPen = list(X = Sr))
  expect_lt(mx(ref$linear.predictors, tk$eta_all[keep]), 1e-5)
  expect_lt(mx(log(ref$sp), log(tk$lambda[-solo])), 1e-3)
})

test_that("ordinary deletions are not flagged", {
  st <- reef_setups$gaussian_w_by
  for (s in c("1", "9", "20"))
    expect_false(gamfluence:::deletion_plan(st, drop_site(reef, s))$rank_deficient)
})

test_that("a collinear full design: mgcv's zeroed column is removed once; ill-conditioning is warned about", {
  set.seed(4)
  d <- data.frame(site = factor(rep(1:15, each = 8)), year = rep(2001:2020, length.out = 120))
  d$y <- 0.02 * (d$year - 2010) + rnorm(15, 0, 0.3)[d$site] + rnorm(120, 0, 0.4)
  m <- gam(y ~ year + I(year^2) + I(year^3) + s(site, bs = "re"), data = d, method = "REML")
  expect_identical(names(coef(m))[coef(m) == 0], "(Intercept)")
  expect_warning(st <- gamfluence:::frozen_setup(m), "does not reproduce")
  expect_identical(st$dropped_at_setup, "(Intercept)")
  i <- suppressWarnings(gam_influence(m, cluster = "site", progress = FALSE))$sites
  expect_true(all(i$status == "ok"))           # deletions add no new deficiency
})

test_that("rows with zero prior weight carry no data for the held-lambda check", {
  d <- rank_data()
  extra <- d[d$site == "6", ][1, ]; extra$grp <- "solo"; extra$w <- 0
  d <- rbind(d, extra); rownames(d) <- NULL
  m <- gam(yg ~ grp + s(x, by = grp) + s(site, bs = "re"), data = d, weights = w, method = "REML")
  i <- gam_influence(m, cluster = "site", which = c("5", "6"), progress = FALSE)$sites
  expect_false(any(grepl("error|not converged", i$status)))
  expect_identical(i$note[1], "held: s(x):grpsolo")
  expect_false(grepl("held", i$note[2]))              # NA, or a "two REML optima" note
})

test_that("no variance-component ratio is reported for an unidentified random-effect term", {
  d <- rank_data()
  m <- gam(yg ~ grp + s(site, by = grp, bs = "re"), data = d, weights = w, method = "REML")
  i <- gam_influence(m, cluster = "site", which = c("5", "6"), progress = FALSE)
  expect_identical(i$sites$note, c("held: s(site):grpsolo", NA))
  held <- i$hyper$id == "5" & i$hyper$parameter == "s(site):grpsolo"
  expect_true(is.na(i$hyper$ratio[held]))
  expect_false(anyNA(i$hyper$ratio[!held]))
})

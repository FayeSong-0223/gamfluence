# Deletion checked against independent references, where one exists.

test_that("T12: Gaussian closed form equals the lambda-held deletion refit", {
  st <- reef_setups$gaussian_w_by
  for (s in c("3", "7", "12")) {
    k <- which(reef$site == s); keep <- setdiff(seq_along(st$y), k)
    Sl <- Reduce(`+`, Map(`*`, st$lambda, st$S)); W <- st$pw
    A  <- solve(crossprod(st$X, W * st$X) + Sl)
    beta <- drop(A %*% crossprod(st$X, W * st$y))
    e  <- st$y - drop(st$X %*% beta); Xk <- st$X[k, , drop = FALSE]
    Hkk <- Xk %*% A %*% t(W[k] * Xk)
    bdel <- beta - drop(A %*% t(W[k] * Xk) %*% solve(diag(length(k)) - Hkk, e[k]))
    expect_lt(mx(beta, coef(reef_models$gaussian_w_by)), 1e-8)
    expect_lt(mx(bdel, gamfluence:::refit(st, keep, st$lambda)$beta), 1e-8)
  }
})

test_that("T13: a deleted random-effect level can keep its zero column (consistency)", {
  st <- reef_setups$gaussian_w_by
  keep <- drop_site(reef, "7")
  plan <- gamfluence:::deletion_plan(st, keep)
  expect_false(plan$rank_deficient)            # the zero re column is penalised: kept
  j  <- which(st$map$term == "s(site)")
  cd <- st$model$smooth[[j]]$first.para - 1 + which(levels(reef$site) == "7")
  Sr <- lapply(st$S, function(M) M[-cd, -cd]); Sr$sp <- rep(-1, length(Sr))
  rm_ <- gam(y ~ X - 1, data = list(y = st$y[keep], X = st$X[keep, -cd], pw = st$pw[keep]),
             weights = pw, method = "REML", paraPen = list(X = Sr))
  tk <- gamfluence:::refit(st, keep)
  expect_lt(mx(rm_$linear.predictors, tk$eta_all[keep]), 1e-5)   # REML optimiser tolerance
})

test_that("T14: invariant to row order and to renaming and reordering levels", {
  d <- droplevels(reef[reef$site %in% as.character(c(1:4, 11:14)), ])
  f <- yg ~ mgmt + s(x) + s(year, by = region, k = 4) + s(site, bs = "re")
  m1 <- gam(f, data = d, weights = w, method = "REML")
  d2 <- d[sample(nrow(d)), ]
  old <- levels(d2$site); new <- paste0("S", sample(100:999, length(old)))
  d2$site <- factor(new[match(as.character(d2$site), old)], levels = sample(new))
  m2 <- gam(f, data = d2, weights = w, method = "REML")
  i1 <- gam_influence(m1, cluster = "site", progress = FALSE)$sites
  i2 <- gam_influence(m2, cluster = "site", progress = FALSE)$sites
  i2 <- i2[match(new, i2$id), ]                # new[i] is old[i] renamed
  i1 <- i1[match(old, i1$id), ]
  for (v in c("fixed", "added", "total"))
    expect_lt(mx(i1[[v]], i2[[v]]), 1e-6)
})

test_that("T16: re-estimating deletion = ordinary gam() refit when the basis is data-independent", {
  # parametric + random-effect models: mgcv builds the same basis from the
  # reduced data, so the naive refit is an independent reference for fit T
  specs <- list(
    gaussian_w = list(f = yg ~ mgmt + s(site, bs = "re"), fam = function() gaussian(), wt = TRUE),
    poisson    = list(f = yp ~ mgmt + s(site, bs = "re") + offset(log(area)), fam = function() poisson(), wt = FALSE),
    negbin     = list(f = yc ~ mgmt + s(site, bs = "re") + offset(log(area)), fam = function() nb(), wt = FALSE),
    binom_cbind = list(f = cbind(succ, fail) ~ mgmt + s(site, bs = "re"), fam = function() binomial(), wt = FALSE))
  fit <- function(z, dat) if (z$wt) gam(z$f, data = dat, weights = w, family = z$fam(), method = "REML") else
                                    gam(z$f, data = dat, family = z$fam(), method = "REML")
  for (nm in names(specs)) {
    z <- specs[[nm]]; st <- gamfluence:::frozen_setup(fit(z, reef))
    for (s in c("1", "6", "15")) {
      keep <- drop_site(reef, s)
      nv <- fit(z, reef[keep, ]); fz <- gamfluence:::refit(st, keep)
      expect_lt(mx(nv$linear.predictors, fz$eta_all[keep]), 1e-4)
      sd_nv <- sqrt(nv$reml.scale / (nv$sp[1] * nv$smooth[[1]]$S[[1]][1, 1]))
      expect_lt(abs(sd_nv - gamfluence:::re_sd(fz$lambda, fz$phi, st, 1)), 1e-4)
      if (!is.null(st$psi)) expect_lt(abs(log(nv$family$getTheta(TRUE) / fz$psi)), 1e-4)
    }
  }
})

test_that("T15: every test model: all refits converge for a set of deletions", {
  for (nm in names(reef_models)) {
    i <- gam_influence(reef_models[[nm]], cluster = "site", which = as.character(1:6),
                       progress = FALSE)
    expect_true(all(i$sites$status == "ok"), info = nm)
  }
})

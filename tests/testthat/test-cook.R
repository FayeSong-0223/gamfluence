# Classical references for the fixed-hyperparameter refit (engine tests; the
# package does not report Cook's distance). Cook's distance for one row is the
# all-rows change of the fixed refit, scaled by EDF x scale.

cook_data <- function() {
  set.seed(7)
  d <- data.frame(x1 = rnorm(60), x2 = rnorm(60), z = runif(60),
                  f = factor(sample(letters[1:3], 60, TRUE)))
  d$y  <- 1 + d$x1 - 0.5 * d$x2 + (d$f == "b") + rnorm(60)
  d$y2 <- 1 + d$x1 + sin(2 * pi * d$z) + rnorm(60, 0, 0.5)
  d$cnt <- rpois(60, exp(0.5 + 0.4 * d$x1))
  d
}

engine_cook <- function(m) {
  st <- gamfluence:::frozen_setup(m)
  n <- length(st$y)
  vapply(seq_len(n), function(i) {
    a <- gamfluence:::refit(st, setdiff(seq_len(n), i), lambda = st$lambda)
    sum(st$w * (a$eta_all - st$eta0)^2) / (st$edf * st$phi_p)
  }, numeric(1))
}

test_that("T9: unpenalised Gaussian: the fixed refit gives stats::cooks.distance(lm)", {
  d <- cook_data()
  m <- gam(y ~ x1 + x2 + f, data = d, method = "REML")
  expect_lt(mx(engine_cook(m), unname(cooks.distance(lm(y ~ x1 + x2 + f, data = d)))), 1e-10)
})

test_that("T10: penalised fixed refit -> classical Cook's as the held lambda -> 0", {
  d <- cook_data()
  err <- sapply(c(1e-1, 1e-3, 1e-6), function(lam) {
    m <- gam(y2 ~ x1 + s(z, k = 6), data = d, method = "REML", sp = lam)
    X <- predict(m, type = "lpmatrix")
    mx(engine_cook(m), unname(cooks.distance(lm(d$y2 ~ X - 1))))
  })
  expect_true(all(diff(err) < 0))
  expect_lt(err[3], 1e-4)
})

test_that("T11: unpenalised Poisson: coefficient changes equal exact glm refits", {
  d  <- cook_data()
  m  <- gam(cnt ~ x1 + x2 + f, data = d, family = poisson, method = "REML")
  g  <- glm(cnt ~ x1 + x2 + f, data = d, family = poisson)
  st <- gamfluence:::frozen_setup(m)
  dc <- sapply(1:60, function(i) {
    b <- gamfluence:::refit(st, setdiff(1:60, i))$beta
    mx(coef(m) - b, coef(g) - coef(update(g, data = d[-i, ])))
  })
  expect_lt(max(dc), 1e-8)
})

# Shared simulated data: reef-survey-like, unbalanced, with prior weights,
# an area offset, a management factor and a region-specific trend.
suppressMessages(library(mgcv))
# Make mgcv's functions visible however the tests are run (under
# pkgload::load_all() the helper environment does not see attached packages).
for (f in c("gam", "bam", "gamm", "nb", "negbin", "tw", "gam.vcomp"))
  assign(f, getExportedValue("mgcv", f))
rm(f)

make_reef <- function(n_site = 20, n_year = 8, keep_frac = 0.85, seed = 42) {
  set.seed(seed)
  d <- expand.grid(site = factor(seq_len(n_site)), year = seq_len(n_year))
  d$region <- factor(ifelse(as.integer(d$site) <= n_site / 2, "A", "B"))
  d$mgmt   <- factor(ifelse(as.integer(d$site) %% 3 == 0, "NTR", "fished"))
  d$x      <- runif(nrow(d))
  d$area   <- runif(nrow(d), 0.5, 1.5)
  d$w      <- runif(nrow(d), 0.5, 2)
  d$ntrial <- sample(5:15, nrow(d), TRUE)
  eta <- 1 + 0.5 * (d$mgmt == "NTR") + sin(2 * pi * d$x) +
         0.1 * d$year * (d$region == "A") + rnorm(n_site, 0, 0.4)[d$site]
  d$yg <- eta + rnorm(nrow(d), 0, 0.5) / sqrt(d$w)
  d$yc <- rnbinom(nrow(d), mu = d$area * exp(eta), size = 3)
  d$yp <- rpois(nrow(d), d$area * exp(eta - 0.5))
  d$succ <- rbinom(nrow(d), d$ntrial, plogis(eta - 1.5))
  d$fail <- d$ntrial - d$succ
  d$yb <- rbinom(nrow(d), 1, plogis(eta - 1.5))
  d <- d[sort(sample(nrow(d), round(keep_frac * nrow(d)))), ]
  rownames(d) <- NULL
  d
}

reef <- make_reef()

fit_models <- function(d) {
  list(
    gaussian_w_by = gam(yg ~ mgmt + s(x) + s(year, by = region, k = 5) + s(site, bs = "re"),
                        data = d, weights = w, method = "REML"),
    negbin_nb     = gam(yc ~ mgmt + s(x) + s(site, bs = "re") + offset(log(area)),
                        data = d, family = nb(), method = "REML"),
    poisson_off   = gam(yp ~ mgmt + s(x) + s(site, bs = "re"), offset = log(area),
                        data = d, family = poisson, method = "REML"),
    binomial_01   = gam(yb ~ mgmt + s(x) + s(site, bs = "re"), data = d,
                        family = binomial, method = "REML"),
    binomial_cbind = gam(cbind(succ, fail) ~ mgmt + s(x) + s(site, bs = "re"), data = d,
                         family = binomial, method = "REML"),
    fixed_sp_in_s = gam(yg ~ mgmt + s(x, sp = 0.5) + s(year, by = region, k = 5) +
                          s(site, bs = "re"), data = d, weights = w, method = "REML"),
    fixed_sp_gam_arg = gam(yg ~ mgmt + s(x) + s(site, bs = "re"), data = d, weights = w,
                           method = "REML", sp = c(0.5, -1)),
    negbin_theta3 = gam(yc ~ mgmt + s(x) + s(site, bs = "re") + offset(log(area)), data = d,
                        family = negbin(3), method = "REML"),
    nb_theta3     = gam(yc ~ mgmt + s(x) + s(site, bs = "re") + offset(log(area)), data = d,
                        family = nb(theta = 3), method = "REML"),
    nb_sqrt_link  = gam(yp ~ mgmt + s(x) + s(site, bs = "re"), data = d,
                        family = nb(link = "sqrt"), method = "REML")
  )
}

reef_models <- fit_models(reef)
reef_setups <- lapply(reef_models, gamfluence:::frozen_setup)

mx <- function(a, b) max(abs(a - b))
drop_site <- function(d, s) which(d$site != s)

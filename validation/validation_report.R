# validation_report.R: the numbers behind the test suite.
# The tests assert thresholds; this script prints the actual values quoted in
# the spec. Run from the package root after installing the package:
#   Rscript validation/validation_report.R
suppressMessages(library(gamfluence))
source("tests/testthat/helper-data.R")
source("tests/testthat/helper-reml.R")
say <- function(...) cat(sprintf(...), "\n", sep = "")
fz  <- gamfluence:::frozen_setup; rf <- gamfluence:::refit
dp  <- gamfluence:::deletion_plan; cf <- gamfluence:::compare_fits
say("gamfluence %s, mgcv %s, R %s", packageVersion("gamfluence"), packageVersion("mgcv"), getRversion())
say("data: %d rows, %d sites, %d years (85%% of site-years kept)", nrow(reef), nlevels(reef$site), 8)

cat("\n[V1] Reproduction of the full fit (10 test models)\n")
for (nm in names(reef_setups)) { st <- reef_setups[[nm]]; a <- seq_along(st$y)
  h <- rf(st, a, st$lambda, st$psi); t <- rf(st, a)
  say("  %-17s held: max|d eta| %.1e   re-estimated: max|d log lambda| %.1e", nm,
      mx(h$eta_all, st$eta0), mx(log(t$lambda), log(st$lambda))) }

tab <- function(x) { t <- table(x); paste(names(t), t, collapse = ", ") }
cat("\n[V2] Status and time, all 20 site deletions per model\n")
cat("     (held refit, plus the re-estimated refit searched from two starts)\n")
for (nm in names(reef_models)) {
  t0 <- proc.time()[["elapsed"]]
  i <- gam_influence(reef_models[[nm]], cluster = "site", progress = FALSE)
  say("  %-17s %.2f s per site (n = %d)   status: %s", nm,
      (proc.time()[["elapsed"]] - t0) / 20, nrow(reef), tab(i$sites$status))
  two <- i$sites[grepl("two optima", i$sites$status), ]
  for (k in seq_len(nrow(two))) say("      site %s: %s", two$id[k], two$note[k]) }

cat("\n[V3] Identity (regression test): fixed + added - total, largest violation over all 20 sites\n")
for (nm in c("gaussian_w_by", "negbin_nb", "binomial_cbind")) { st <- reef_setups[[nm]]
  v <- sapply(levels(reef$site), function(s) { k <- drop_site(reef, s); r <- cf(st, k, dp(st, k))
    max(abs(r$fixed + r$added - r$total)) })
  say("  %-15s %.1e", nm, max(v)) }

cat("\n[V4] Independent reference for fit T, data-independent basis (all 20 sites)\n")
specs <- list(
  gaussian_w  = list(f = yg ~ mgmt + s(site, bs = "re"), fam = function() gaussian(), wt = TRUE),
  poisson     = list(f = yp ~ mgmt + s(site, bs = "re") + offset(log(area)), fam = function() poisson(), wt = FALSE),
  negbin      = list(f = yc ~ mgmt + s(site, bs = "re") + offset(log(area)), fam = function() nb(), wt = FALSE),
  binom_cbind = list(f = cbind(succ, fail) ~ mgmt + s(site, bs = "re"), fam = function() binomial(), wt = FALSE))
fitz <- function(z, dat) if (z$wt) gam(z$f, data = dat, weights = w, family = z$fam(), method = "REML") else
                                   gam(z$f, data = dat, family = z$fam(), method = "REML")
for (nm in names(specs)) { z <- specs[[nm]]; st <- fz(fitz(z, reef))
  r <- sapply(levels(reef$site), function(s) { k <- drop_site(reef, s); nv <- fitz(z, reef[k, ]); f <- rf(st, k)
    c(mx(nv$linear.predictors, f$eta_all[k]),
      abs(sqrt(nv$reml.scale / (nv$sp[1] * nv$smooth[[1]]$S[[1]][1, 1])) - gamfluence:::re_sd(f$lambda, f$phi, st, 1))) })
  say("  %-11s max|d eta| %.1e   max|d site SD| %.1e", nm, max(r[1, ]), max(r[2, ])) }

cat("\n[V5] T24: independent REML / Laplace-REML optimiser on the frozen X, models with covariate smooths\n")
say("  gap = criterion at the refit's lambda minus the reference optimum (>= 0 means the reference did no worse)")
for (nm in c("gaussian_w_by", "fixed_sp_in_s", "fixed_sp_gam_arg", "poisson_off", "binomial_cbind")) {
  st <- reef_setups[[nm]]
  for (s in c("none", "2", "6", "9", "13", "17")) {
    k  <- if (s == "none") seq_along(st$y) else drop_site(reef, s)
    pr <- reml_problem(st, k)
    f  <- if (s == "none") list(lambda = st$lambda, eta_all = st$eta0) else rf(st, k)
    o  <- reml_optimise(pr, log(st$lambda[pr$est]) + 0.7)
    dl <- abs(log(f$lambda / o$lambda))[pr$est]
    say("  %-16s deleted %-4s gap %9.1e   full-fit lambda gap %8.1e   max|d eta| %.1e   |d log lambda| (estimated only) %s",
        nm, s, reml_value(pr, f$lambda) - o$value,
        if (s == "none") NA else reml_value(pr, st$lambda) - o$value,
        mx(pirls(pr, o$lambda)$eta, f$eta_all[k]), paste(sprintf("%.0e", dl), collapse = " ")) } }

cat("  gaussian_w_by / fixed_sp_in_s: the 2nd and 3rd estimated lambda are the by = region smooths\n")
cat("\n[V6] Rank deficiency after deletion\n")
rd <- make_reef(n_site = 16, n_year = 6, keep_frac = 1, seed = 21)
rd$zone <- factor(ifelse(rd$site == "3", "rare", ifelse(as.integer(rd$site) %% 2 == 0, "a", "b")))
rd$grp  <- factor(ifelse(rd$site == "5", "solo", "main"))
m1 <- gam(yg ~ zone + s(site, bs = "re"), data = rd, weights = w, method = "REML"); s1 <- fz(m1)
k  <- which(rd$site != "3"); p1 <- dp(s1, k)
nv <- gam(yg ~ zone + s(site, bs = "re"), data = droplevels(rd[k, ]), weights = w, method = "REML")
say("  factor level only in site 3: dropped %s; retained-row eta vs ordinary gam() on reduced data %.1e",
    paste(colnames(s1$X)[p1$drop], collapse = ","), mx(nv$linear.predictors, rf(s1, k, plan = p1)$eta_all[k]))
m2 <- gam(yg ~ grp + s(x, by = grp) + s(site, bs = "re"), data = rd, weights = w, method = "REML"); s2 <- fz(m2)
k  <- which(rd$site != "5"); p2 <- dp(s2, k); solo <- which(s2$map$term == "s(x):grpsolo")
t2 <- rf(s2, k, plan = p2)
gone <- c(which(colnames(s2$X) == "grpsolo"), s2$pen_cols[[solo]])
Sr <- lapply(s2$S[-solo], function(M) M[-gone, -gone]); Sr$sp <- rep(-1, length(Sr))
ref <- gam(y ~ X - 1, data = list(y = s2$y[k], X = s2$X[k, -gone], pw = s2$pw[k]), weights = pw,
           method = "REML", paraPen = list(X = Sr))
say("  by = factor level only in site 5: held %s; rank-deficient %s; eta vs refit without that level %.1e; other log lambda %.1e",
    s2$map$term[p2$hold], p2$rank_deficient, mx(ref$linear.predictors, t2$eta_all[k]), mx(log(ref$sp), log(t2$lambda[-solo])))
i2 <- gam_influence(m2, cluster = "site", progress = FALSE)
say("  gam_influence on that model: status of site 5: %s; other 15 sites: %s",
    i2$sites$status[i2$sites$id == "5"], tab(i2$sites$status[i2$sites$id != "5"]))

cat("\n")

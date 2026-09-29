# planted_site.R: does re-estimation reveal influence that a fixed-hyperparameter
# deletion misses? A small, known-truth check run through the package itself.
#   Rscript validation/planted_site.R      (about 5 minutes)
#
# Design (the same as the feasibility check, with the same random draws):
#   y = sin(2 pi x) + site effect + noise; site SD 0.4, residual SD 0.5;
#   40 sites x 8 observations; model y ~ s(x) + s(site, bs = "re"), REML.
#   Null: 20 datasets, every site deleted in turn (800 deletions).
#   Planted: 20 datasets in which site 1's effect is set to +2.0 (5 SD);
#   site 1 is compared with the null 99% lines.
suppressMessages({ library(mgcv); library(gamfluence) })
say <- function(...) cat(sprintf(...), "\n", sep = "")
cp  <- function(x, n) { b <- binom.test(x, n)$conf.int
  sprintf("%d/%d (95%% CI %.0f-%.0f%%)", x, n, 100 * b[1], 100 * b[2]) }
q3  <- function(v) paste(signif(quantile(v, c(.5, .95, .99)), 3), collapse = " / ")

sim_data <- function(plant = c("none", "level")) {
  plant <- match.arg(plant)
  dd <- data.frame(site = factor(rep(1:40, each = 8)), x = runif(320))
  b  <- rnorm(40, 0, 0.4)
  if (plant == "level") b[1] <- 2.0
  dd$y <- sin(2 * pi * dd$x) + b[dd$site] + rnorm(320, 0, 0.5)
  dd
}
run <- function(dd, which = NULL) {
  m <- gam(y ~ s(x) + s(site, bs = "re"), data = dd, method = "REML")
  i <- gam_influence(m, cluster = "site", which = which, progress = FALSE)
  h <- i$hyper
  cbind(i$sites[c("fixed", "added", "total")],
        ok = i$sites$status == "ok",
        sd_ratio = h$ratio[h$parameter == "s(site)"],
        lambda_x_ratio = h$ratio[h$parameter == "s(x)"])
}
meas <- c("fixed", "added", "total")
say("gamfluence %s, mgcv %s, R %s", packageVersion("gamfluence"), packageVersion("mgcv"), getRversion())

set.seed(2026); R_null <- 20
null <- do.call(rbind, lapply(seq_len(R_null), function(r) cbind(rep = r, run(sim_data("none")))))
say("\n[Null] %d datasets x 40 sites = %d deletions; status ok: %d", R_null, nrow(null), sum(null$ok))
for (v in meas) say("  %-6s median / 95%% / 99%%: %s", v, q3(null[[v]]))
say("  site-SD ratio median %.3f; s(x) lambda ratio median %.3f",
    median(null$sd_ratio), median(null$lambda_x_ratio))

line <- sapply(meas, function(v) unname(quantile(null[[v]], .99)))
set.seed(99)
boot <- replicate(2000, { rr <- sample(R_null, R_null, TRUE)
  idx <- unlist(lapply(rr, function(r) which(null$rep == r)))
  sapply(meas, function(v) unname(quantile(null[[v]][idx], .99))) })
say("\n[Null 99%% lines] each rests on about 8 of 800 values")
for (v in meas) {
  fa <- sum(tapply(null[[v]] > line[v], null$rep, any))
  say("  %-6s %.4f   dataset-bootstrap 90%% interval %.4f-%.4f   null datasets with >= 1 site above: %s",
      v, line[v], quantile(boot[v, ], .05), quantile(boot[v, ], .95), cp(fa, R_null))
}

R_p <- 20
pl <- do.call(rbind, lapply(seq_len(R_p), function(r) run(sim_data("level"), which = "1")))
say("\n[Planted] site 1 at +5 SD, %d datasets; status ok: %d", R_p, sum(pl$ok))
for (v in meas) say("  %-6s median %.4f   above its null 99%% line: %s",
                    v, median(pl[[v]]), cp(sum(pl[[v]] > line[v]), R_p))
say("  added above its line while fixed is not: %s",
    cp(sum(pl$added > line["added"] & pl$fixed <= line["fixed"]), R_p))
say("  site-SD ratio median %.3f; s(x) lambda ratio median %.3f",
    median(pl$sd_ratio), median(pl$lambda_x_ratio))
cat("\n")

# data-raw/coral_trout_influence.R: the saved gam_influence() result that the
# vignette loads, because deleting all 71 sites takes about 3 minutes.
# Run from the package root, with the package installed:
#   Rscript data-raw/coral_trout_influence.R
library(mgcv)
library(gamfluence)
m <- gam(count ~ REGION * NTR + EXPOSURE + depth + s(YEAR, by = REGION, k = 5) +
           s(rugosity, k = 5) + s(LHC, k = 5) + s(kd490, k = 5) + s(maxDHW, k = 5) +
           s(Cyclone, k = 5) + s(SITE, bs = "re"),
         family = nb(), data = coral_trout, method = "REML")
infl <- gam_influence(m, cluster = "SITE", progress = FALSE)
saveRDS(infl, "inst/extdata/coral_trout_influence.rds", compress = "xz")

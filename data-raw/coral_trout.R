# data-raw/coral_trout.R: builds data/coral_trout.rda from the analysis dataset
# of the coral trout project (Song 2026, https://doi.org/10.5281/zenodo.22933039).
# That file is written by the project's coral_trout_stage1.R from the AIMS
# extract; the project's DATA.md says where the AIMS files come from.
# Run from the package root:
#   Rscript data-raw/coral_trout.R path/to/analysis_dataset.csv
args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args)) args[1] else
  "~/Desktop/Coral trout Study/Stage1/outputs/analysis_dataset.csv"
d <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)

coral_trout <- data.frame(
  SITE     = factor(d$SITE),
  REGION   = factor(d$REGION, levels = c("Palm", "Whitsunday")),
  YEAR     = as.integer(d$YEAR),
  NTR      = factor(d$NTR, levels = c("Fished", "NTR 1987", "NTR 2004")),
  EXPOSURE = factor(d$EXPOSURE, levels = c("Sheltered", "Semi-Exposed", "Exposed")),
  depth    = d[["Corrected depth"]],
  rugosity = d$rugosity,
  LHC      = d[["LHC_%"]],
  kd490    = d$kd490,
  maxDHW   = d$maxDHW,
  Cyclone  = d$Cyclone,
  count    = as.integer(round(d$count))
)
stopifnot(nrow(coral_trout) == 467, nlevels(coral_trout$SITE) == 71,
          !anyNA(coral_trout), all(coral_trout$count == d$count))
save(coral_trout, file = "data/coral_trout.rda", compress = "xz")

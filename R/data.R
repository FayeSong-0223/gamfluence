#' Coral trout counts on inshore reefs of the Great Barrier Reef
#'
#' Counts of coral trout (*Plectropomus leopardus*) at 71 inshore reef sites in
#' the Palm and Whitsunday island groups, surveyed six or seven times between
#' 2007 and 2018 (467 site-surveys), with each site's management zone, wave
#' exposure, and habitat and environmental covariates. It is the analysis
#' dataset of the coral trout study cited below and is used in the vignette.
#'
#' @format A data frame with 467 rows and 12 columns:
#' \describe{
#'   \item{SITE}{Site, a factor with 71 levels.}
#'   \item{REGION}{Island group: `Palm` or `Whitsunday`.}
#'   \item{YEAR}{Survey year.}
#'   \item{NTR}{Management zone: `Fished`, or a no-take reserve zoned in 1987
#'     (`NTR 1987`) or 2004 (`NTR 2004`). Constant within a site.}
#'   \item{EXPOSURE}{Wave exposure: `Sheltered`, `Semi-Exposed` or `Exposed`.
#'     Constant within a site.}
#'   \item{depth}{Survey depth, the `Corrected depth` column of the source.
#'     Constant within a site.}
#'   \item{rugosity}{Reef structural complexity.}
#'   \item{LHC}{Live hard coral cover (percent).}
#'   \item{kd490}{Water clarity: diffuse light attenuation at 490 nm (higher
#'     is more turbid).}
#'   \item{maxDHW}{Thermal stress: maximum degree heating weeks.}
#'   \item{Cyclone}{Cyclone exposure, as given in the source.}
#'   \item{count}{Number of *P. leopardus* counted in the survey. The source
#'     gives densities; they are an exact multiple of the counts with a
#'     year-specific constant, so the counts were recovered from them.}
#' }
#'
#' @source Derived from Australian Institute of Marine Science (AIMS). (2022).
#'   Spatio-temporal dynamics of coral reef fish assemblages on inshore reefs
#'   of the Great Barrier Reef.
#'   <https://apps.aims.gov.au/metadata/view/814a0be3-ed85-4a43-87b7-59ece6eb6a05>,
#'   accessed 02-Sep-2026. Licensed under Creative Commons Attribution 3.0
#'   Australia (CC BY 3.0 AU). Based on Australian Institute of Marine Science
#'   data. The data were restricted to *P. leopardus* in the Palm and
#'   Whitsunday island groups from 2007 to 2018, and densities were converted
#'   back to counts. See `LICENSE.note` for the licence of this dataset.
#'
#' @references Song, Z. (2026). *What explains variation in Plectropomus
#'   leopardus density* (Version v1.0). Zenodo.
#'   \doi{10.5281/zenodo.22933039}
#'
#' @examples
#' str(coral_trout)
#' table(coral_trout$REGION, coral_trout$NTR)
"coral_trout"

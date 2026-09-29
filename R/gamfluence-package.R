#' gamfluence: deletion influence diagnostics for penalised GAMs
#'
#' Site (cluster) deletion diagnostics for models fitted with [mgcv::gam()]:
#' each site is removed and the model refitted with the estimated
#' hyperparameters fixed and re-estimated, on the frozen full-data
#' representation. The main function is [gam_influence()].
#'
#' @keywords internal
"_PACKAGE"

# `pw` and `off` are columns of the data list passed to mgcv::gam() in refit().
utils::globalVariables(c("pw", "off"))

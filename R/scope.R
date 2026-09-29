# Scope and response-format guards.
#
# v0.1 supports a deliberately narrow class of models: the ones the feasibility
# checks and the test suite cover. Everything the guard checks outside that
# class stops with an error before any refit. The guard cannot catch every
# possible departure; the cases it is tested against are in test-scope.R.

supported_links <- function(family_name, is_nb) {
  if (is_nb) return(c("log", "sqrt"))
  switch(family_name,
         gaussian = "identity",
         poisson  = "log",
         binomial = "logit",
         character(0))
}

# Returns list(is_nb, theta_estimated) or stops.
check_scope <- function(model) {
  if (!inherits(model, "gam"))
    stop("`model` must be a fitted mgcv::gam() object.", call. = FALSE)
  if (inherits(model, "bam"))
    stop("bam() models are not supported.", call. = FALSE)
  if (isTRUE(grepl("^lme", model$method)))
    stop("gamm() models are not supported.", call. = FALSE)
  if (is.null(model$coefficients) || is.null(model$linear.predictors))
    stop("`model` has not been fitted.", call. = FALSE)
  if (!identical(model$method, "REML"))
    stop("only method = \"REML\" is supported in v0.1 (model used \"",
         model$method, "\").", call. = FALSE)

  fam   <- model$family$family
  is_nb <- grepl("^Negative Binomial", fam)
  if (!(fam %in% c("gaussian", "poisson", "binomial") || is_nb))
    stop("family not supported in v0.1: ", fam, call. = FALSE)
  if (fam == "gaussian" && !isTRUE(model$scale.estimated))
    stop("a user-supplied scale is not supported.", call. = FALSE)
  if (fam != "gaussian" && isTRUE(model$scale.estimated))
    stop("an estimated scale is not supported for this family.", call. = FALSE)
  if (!(model$family$link %in% supported_links(fam, is_nb)))
    stop("link not supported in v0.1 for this family: ", model$family$link,
         call. = FALSE)

  # The refits use mgcv's default outer Newton optimiser; a different optimiser
  # can stop at a slightly different lambda and pollute the hyperparameter channel.
  if (!identical(model$optimizer, c("outer", "newton")))
    stop("only the default optimizer (outer Newton) is supported; model used \"",
         paste(model$optimizer, collapse = " "), "\".", call. = FALSE)

  # Arguments that must be left at their defaults. H (a fixed extra penalty)
  # and paraPen (extra penalties) are not in the frozen penalty list.
  for (a in c("gamma", "scale", "min.sp", "H", "paraPen")) {
    v <- model$call[[a]]
    if (is.null(v)) next
    val <- tryCatch(eval(v, environment(stats::formula(model))),
                    error = function(e) NA)
    at_default <- (a == "gamma" && isTRUE(all.equal(val, 1))) ||
                  (a == "scale" && isTRUE(all.equal(val, 0)))
    if (!at_default)
      stop("argument not supported in v0.1: ", a, call. = FALSE)
  }

  for (sm in model$smooth) {
    if (inherits(sm, c("tensor.smooth", "t2.smooth", "fs.interaction")))
      stop("smooth class not supported in v0.1: ", sm$label, call. = FALSE)
    if (!(class(sm)[1] %in% c("tprs.smooth", "random.effect")))
      stop("basis not supported in v0.1 (default thin plate and bs = \"re\" only): ",
           sm$label, call. = FALSE)
    if (inherits(sm, "tprs.smooth") && (sm$dim != 1 || !isTRUE(sm$p.order[1] == 0)))
      stop("only one-dimensional thin plate smooths with the default penalty order ",
           "are supported in v0.1: ", sm$label, call. = FALSE)
    if (inherits(sm, "random.effect") &&
        (length(sm$term) != 1 || !is.factor(model$model[[sm$term[1]]])))
      stop("random effects must be s(g, bs = \"re\") with a single factor g ",
           "(random slopes are not supported in v0.1): ", sm$label, call. = FALSE)
    if (!is.null(sm$id))
      stop("smooths linked by id = are not supported.", call. = FALSE)
    if (!identical(sm$by, "NA") && is.null(sm$by.level))
      stop("numeric by = variables are not supported: ", sm$label, call. = FALSE)
    if (length(sm$S) == 0)
      stop("unpenalised smooths (fx = TRUE) are not supported: ", sm$label,
           call. = FALSE)
    if (length(sm$S) != 1)
      stop("smooths with more than one penalty are not supported in v0.1: ",
           sm$label, call. = FALSE)
  }

  check_response(model, fam, is_nb)

  theta_est <- is_nb && inherits(model$family, "extended.family") &&
               isTRUE(model$family$n.theta > 0)
  list(is_nb = is_nb, theta_estimated = theta_est)
}

# The refits use model$y and model$prior.weights. That is correct for a vector
# response and for binomial responses given as 0/1, a factor, proportions with
# weights, or cbind(successes, failures) (mgcv stores the proportion and puts
# the trials in the prior weights). Anything else is rejected.
check_response <- function(model, fam, is_nb) {
  resp <- stats::model.response(model$model)
  nc <- NCOL(resp)
  if (fam == "binomial") {
    if (nc > 2)
      stop("binomial response with more than two columns is not supported.",
           call. = FALSE)
    if (any(model$y < 0 | model$y > 1))
      stop("binomial response outside [0, 1] after mgcv's processing.",
           call. = FALSE)
  } else if (nc != 1) {
    stop("matrix responses are only supported for binomial models ",
         "(cbind(successes, failures)).", call. = FALSE)
  }
  if ((fam == "poisson" || is_nb) && any(model$y < 0))
    stop("negative counts in the response.", call. = FALSE)
  if (length(model$y) != nrow(model$model) ||
      length(model$prior.weights) != length(model$y))
    stop("response, weights and model frame have different numbers of rows.",
         call. = FALSE)
  invisible(TRUE)
}

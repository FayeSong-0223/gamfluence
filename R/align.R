# Row alignment: deletion units are defined on the rows mgcv actually used.
#
# Every quantity in the refits (model matrix, response, weights, offset) is
# indexed by the rows of the model frame, model$model. Those can differ from
# the rows of the user's data frame: na.action drops incomplete rows and
# `subset` drops others. Units are therefore always mapped onto model-frame
# rows, using the row names mgcv keeps from the original data, and the mapping
# is checked against the data values rather than trusted.

get_model_data <- function(model, data, env = parent.frame()) {
  if (!is.null(data)) {
    if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
    return(data)
  }
  expr <- model$call$data
  if (is.null(expr))
    stop("`data` was not supplied and the model call has no `data` argument.",
         call. = FALSE)
  # look where the caller is first, then where the model formula was made
  out <- NULL
  for (e in list(env, environment(stats::formula(model)), globalenv())) {
    if (is.null(e)) next
    out <- tryCatch(eval(expr, e), error = function(err) NULL)
    if (is.data.frame(out)) break
  }
  if (!is.data.frame(out))
    stop("`data` was not supplied and could not be recovered from the model call; ",
         "pass the data frame used to fit the model.", call. = FALSE)
  message("gamfluence: using `", deparse(expr),
          "` from the model call as `data`.")
  out
}

# Find the model-frame rows inside `data` and check that the values agree.
# Every column of the model frame is re-evaluated on the matched rows of
# `data` (so transformed terms such as log(y) or offset(log(area)) are checked
# too, as are the prior weights) and must equal what mgcv stored. If nothing
# can be re-evaluated, the match cannot be verified and that is an error.
align_rows <- function(model, data) {
  rn <- rownames(model$model)
  idx <- match(rn, rownames(data))
  if (anyNA(idx))
    stop("the rows used to fit the model cannot all be found in `data` by row ",
         "name; pass the data frame that was used to fit the model.", call. = FALSE)
  sub <- data[idx, , drop = FALSE]
  env <- environment(stats::formula(model))
  checked <- 0L
  for (v in names(model$model)) {
    expr <- if (v == "(weights)") model$call$weights else
            if (v == "(offset)") model$call$offset else
            tryCatch(str2lang(v), error = function(e) NULL)
    if (is.null(expr)) next
    b <- tryCatch(eval(expr, sub, env), error = function(e) NULL)
    if (is.null(b)) next
    if (!same_values(model$model[[v]], b))
      stop("`", v, "` computed from `data` does not match the rows the model was ",
           "fitted to. Was `data` modified or re-sorted without its row names ",
           "after fitting?", call. = FALSE)
    checked <- checked + 1L
  }
  if (checked == 0L)
    stop("could not check that `data` matches the rows the model was fitted to; ",
         "pass `cluster` as a vector with one value per row used by the model.",
         call. = FALSE)
  idx
}

same_values <- function(a, b) {
  if (is.factor(a) || is.character(a) || is.logical(a) ||
      is.factor(b) || is.character(b) || is.logical(b))
    return(identical(as.character(a), as.character(b)))
  a <- as.numeric(a); b <- as.numeric(b)
  length(a) == length(b) && isTRUE(all.equal(a, b, check.attributes = FALSE))
}

# Returns list(id, rows): rows[[k]] are the model-frame rows deleted for site k.
resolve_units <- function(model, data, cluster, which, env = parent.frame()) {
  n <- nrow(model$model)
  if (is.character(cluster) && length(cluster) == 1L) {
    if (cluster %in% names(model$model)) {
      g <- model$model[[cluster]]
    } else {
      data <- get_model_data(model, data, env)
      if (!(cluster %in% names(data)))
        stop("`cluster` column \"", cluster, "\" is in neither the model frame ",
             "nor `data`.", call. = FALSE)
      g <- data[[cluster]][align_rows(model, data)]
    }
  } else {
    if (length(cluster) == n) {
      g <- cluster
    } else {
      data <- get_model_data(model, data, env)
      if (length(cluster) != nrow(data))
        stop("`cluster` has length ", length(cluster), "; it must match either ",
             "the ", n, " rows used by the model or the ", nrow(data),
             " rows of `data`.", call. = FALSE)
      g <- cluster[align_rows(model, data)]
    }
  }
  if (anyNA(g))
    stop("`cluster` is missing for some rows used by the model.", call. = FALSE)
  g <- droplevels(factor(g))
  rows <- split(seq_len(n), g)
  ids <- names(rows)

  if (!is.null(which)) {
    which <- as.character(which)
    bad <- setdiff(which, ids)
    if (length(bad))
      stop("unknown site(s) in `which`: ", paste(utils::head(bad, 5), collapse = ", "),
           call. = FALSE)
    keep <- ids %in% which
    ids <- ids[keep]; rows <- rows[keep]
  }
  if (length(ids) == 1L && length(levels(g)) == 1L)
    stop("`cluster` has only one level; there is nothing to compare.", call. = FALSE)
  list(id = ids, rows = unname(rows))
}

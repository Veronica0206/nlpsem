#' Inspect the Specification and Provenance of a Fitted Model
#'
#' @description Returns the resolved specification, fitting settings and provenance
#' recorded by the estimation functions. No caller environment or duplicate data
#' frame is retained. The original S4 slots remain unchanged, so saved models from
#' earlier versions remain readable.
#'
#' @param object An object of class \code{myMxOutput} returned by an estimation function.
#' @return A list containing \code{schema_version}, \code{call}, \code{function_name},
#' \code{specification}, \code{fitting}, \code{data}, \code{versions}, and
#' \code{diagnostics}, or \code{NULL} for older objects without saved metadata.
#' @details
#' The specification contains resolved model arguments, with canonical curve names.
#' The recorded call uses \code{.data} and, when explicitly supplied, \code{.starts}
#' placeholders; it is a provenance record, not an automatically replayable call.
#' The fitting section contains the resolved initial values, optimizer, thread count, retry settings,
#' RNG kind and the pre-initialization RNG state if one existed. Recording metadata
#' neither sets a seed nor creates an RNG state. Unused argument expressions are
#' not evaluated; their fields are NULL and their names appear in
#' \code{fitting$unevaluated_arguments}. Literal unused defaults may be recorded.
#' Reproducibility also depends on the
#' data, software versions, optimizer and numerical platform.
#'
#' The data section records original row count/names, retained and dropped row indices,
#' and, for manifest groups, the original label-to-code mapping and original row indices
#' in each group. Raw observations remain in \code{object@mxOutput}.
#' Diagnostics report the optimizer status without reclassifying a nonconverged fit
#' as successful. Modifying the underlying OpenMx object afterwards does not update
#' this stored provenance snapshot.
#'
#' @examples
#' # Historical objects remain valid and return NULL when metadata was not recorded.
#' load(system.file("extdata", "getLGCM_examples.RData", package = "nlpsem"))
#' getSpecification(Math_LGCM_BLS_f)
#'
#' @seealso \code{\link{getLGCM}}, \code{\link{getMGroup}}, \code{\link{getMIX}}
#' @export
getSpecification <- function(object) {
  if (!methods::is(object, "myMxOutput")) {
    stop("object must be a myMxOutput fitted model.", call. = FALSE)
  }
  metadata <- attr(object, "nlpsem", exact = TRUE)
  if (is.null(metadata)) return(NULL)
  if (!is.list(metadata) || !identical(metadata$schema_version, 1L)) {
    stop("Unsupported or malformed nlpsem specification metadata.", call. = FALSE)
  }
  metadata
}

.capture_fit_metadata <- function(call, function_name, frame, dat) {
  # Snapshot before initialization, without forcing any fitting argument promises.
  call[[1L]] <- as.name(function_name)
  call$dat <- as.name(".data")
  if ("starts" %in% names(call)) call$starts <- as.name(".starts")
  list(schema_version = 1L, call = call, function_name = function_name,
       specification = list(),
       fitting = list(starts_source = "user",
                      optimizer = OpenMx::mxOption(NULL, "Default optimizer"),
                      number_of_threads = OpenMx::mxOption(NULL, "Number of Threads"),
                      RNG_kind = RNGkind(),
                      RNG_state = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL),
       data = list(original_n = nrow(dat), retained_rows = seq_len(nrow(dat)),
                    dropped_rows = integer(), row_names = rownames(dat),
                    group_mapping = NULL, group_rows = NULL),
       versions = list(R = R.version.string, OpenMx = as.character(utils::packageVersion("OpenMx")),
                       nlpsem = as.character(utils::packageVersion("nlpsem"))), diagnostics = list())
}

.resolve_fit_metadata <- function(metadata, frame) {
  unevaluated <- character()
  value <- function(name) {
    if (!exists(name, envir = frame, inherits = FALSE)) return(NULL)
    if (!rlang::env_binding_are_lazy(frame, name)) {
      return(get(name, envir = frame, inherits = FALSE))
    }
    # Literal defaults and do.call() values are safe to copy without forcing.
    # Expressions for unused options must remain unevaluated: even reading an
    # unused retry control can otherwise consume randomness or raise an error.
    expression <- substitute(ARG, list(ARG = as.name(name)))
    expression <- do.call(substitute, list(expression, frame))
    if (is.null(expression) || is.atomic(expression)) return(expression)
    unevaluated <<- c(unevaluated, name)
    NULL
  }
  function_name <- metadata$function_name
  specification_keys <- c("t_var", "y_var", "curveFun", "intrinsic", "records", "growth_TIC",
                          "TVC", "decompose", "x_type", "x_var", "m_var", "grp_var", "cluster_TIC")
  spec <- stats::setNames(lapply(specification_keys, value), specification_keys)
  family <- switch(function_name, getLGCM = "LGCM", getLCSM = "LCSM", getTVCmodel = "TVC",
                    getMGM = "MGM", getMediation = "MED", value("sub_Model"))
  spec$sub_Model <- family
  spec$y_model <- if (family %in% c("LGCM", "LCSM")) family else if (family == "MED") "LGCM" else value("y_model")
  spec$container <- if (function_name == "getMGroup") "multigroup" else if (function_name == "getMIX") "mixture" else "single"
  spec["nClass"] <- list(NULL)
  curve_names <- c("linear" = "LIN", "LIN" = "LIN", "quadratic" = "QUAD", "QUAD" = "QUAD",
                    "negative exponential" = "EXP", "EXP" = "EXP", "Jenss-Bayley" = "JB", "JB" = "JB",
                    "bilinear spline" = "BLS", "BLS" = "BLS", "nonparametric" = "NonP", "NonP" = "NonP")
  spec$curveFun <- unname(curve_names[[as.character(spec$curveFun)]])
  fitting_keys <- c("tries", "OKStatus", "jitterD", "loc", "scale", "res_scale", "res_cor",
                    "prop_starts", "paramOut", "names")
  fitting <- stats::setNames(lapply(fitting_keys, value), fitting_keys)
  metadata$specification <- spec
  metadata$fitting[names(fitting)] <- fitting
  metadata$fitting$unevaluated_arguments <- unique(unevaluated)
  metadata
}

.attach_fit_metadata <- function(object, metadata, starts, nClass = NULL, group_recode = NULL,
                                 frame) {
  metadata <- .resolve_fit_metadata(metadata, frame)
  metadata$specification["nClass"] <- list(nClass)
  metadata$fitting["starts"] <- list(starts)
  if (!is.null(group_recode)) {
    codes <- group_recode$codes
    metadata$data$retained_rows <- which(!is.na(codes))
    metadata$data$dropped_rows <- which(is.na(codes))
    metadata$data$group_mapping <- data.frame(code = seq_len(group_recode$nClass),
                                              label = as.character(group_recode$levels),
                                              stringsAsFactors = FALSE)
    metadata$data$group_rows <- lapply(seq_len(group_recode$nClass), function(k) which(codes == k))
    names(metadata$data$group_rows) <- as.character(seq_len(group_recode$nClass))
  }
  mx <- object@mxOutput
  status <- mx@output$status$code
  likelihood <- mx@output$Minus2LogLikelihood
  metadata$diagnostics <- list(status = if (is.null(status)) NA_integer_ else as.integer(status),
                                minus2LL = if (is.null(likelihood)) NA_real_ else as.numeric(likelihood),
                                n_parameters = length(OpenMx::omxGetParameters(mx, free = TRUE)),
                                n_observations = length(metadata$data$retained_rows))
  attr(object, "nlpsem") <- metadata
  object
}

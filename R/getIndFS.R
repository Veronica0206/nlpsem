#' @title Derive Individual Factor Scores for Each Latent Variable Included in Model
#'
#' @description This function computes individual factor scores for each latent variable in a given model. It supports three
#' types of factor scores: maximum likelihood, weighted maximum likelihood, and regression.
#'
#' @param model A fitted mxModel object. Specifically, this should be the \code{mxOutput} slot from the result returned by
#' a single-process or multivariate model with one RAM or LISREL expectation
#' and its own raw data. Multiple-group and mixture parent models are not supported.
#' @param FS_type A string specifying the type of factor scores to compute. Supported options include \code{"ML"} (for
#' Maximum Likelihood), \code{"WeightedML"} (for Weighted Maximum Likelihood), and \code{"Regression"}. Default is \code{"Regression"}.
#'
#' @return An object of class \code{FSOutput} with two slots:
#' \itemize{
#'   \item \code{scores_est}: Contains the factor score estimates.
#'   \item \code{scores_se}: Contains the standard errors of the factor score estimates.
#' }
#' The content of these slots can be printed using the \code{printTable()} method for S4 objects.
#'
#' @details Before scoring, the observed data stored in the model are reduced to the variables the model actually uses
#' (manifest variables and definition variables), and the model is given a name without spaces or punctuation. Both
#' steps are required by \code{OpenMx::mxFactorScores()}: columns that are not part of the model make the
#' \code{"ML"} and \code{"WeightedML"} scoring runs fail, and the descriptive model names used by this package
#' are not valid identifiers for the scoring model.
#'
#' @references
#' \itemize{
#'   \item {Estabrook, R. & Neale, M. C. (2013). A Comparison of Factor Score Estimation Methods in the Presence of Missing
#'   Data: Reliability and an Application to Nicotine Dependence. Multivariate Behavioral Research, 48, 1-27. \doi{10.1080/00273171.2012.730072}}
#'   \item {Priestley, M. & Subba Rao, T. (1975). The Estimation of Factor Scores and Kalman Filtering For Discrete Parameter
#'   Stationary Processes. International Journal of Control, 21, 971-975. \doi{10.1080/00207177508922050}}
#' }
#'
#' @seealso \code{\link{getEstimateStats}}
#'
#' @export
#'
#' @examples
#' # Fitted bilinear spline model (fixed knot) shipped with the package
#' load(system.file("extdata", "getLGCM_examples.RData", package = "nlpsem"))
#' \donttest{
#' # Regression scores for the 500 persons (individual measurement occasions are definition
#' # variables, so the scores are computed person by person; about 15 seconds)
#' FS <- getIndFS(model = Math_LGCM_BLS_r@mxOutput, FS_type = "Regression")
#' head(FS@scores_est)
#' # Maximum-likelihood scores re-optimize one model per person
#' FS_ML <- getIndFS(model = Math_LGCM_BLS_r@mxOutput, FS_type = "ML")
#' head(FS_ML@scores_est)
#' }
#'
#' @importFrom OpenMx mxFactorScores mxRename mxModel mxData
#' @importFrom methods new
#'
getIndFS <- function(model, FS_type = "Regression"){
  if (!is.character(FS_type) || length(FS_type) != 1L || is.na(FS_type) ||
      !FS_type %in% c("ML", "WeightedML", "Regression")) {
    stop("FS_type must be one of 'ML', 'WeightedML', or 'Regression'.", call. = FALSE)
  }
  if (!methods::is(model, "MxModel")) {
    stop("model must be a fitted OpenMx model (fit@mxOutput).", call. = FALSE)
  }
  if (length(model@submodels) ||
      !(methods::is(model@expectation, "MxExpectationRAM") ||
        methods::is(model@expectation, "MxExpectationLISREL")) ||
      is.null(model$data) || !identical(model$data$type, "raw")) {
    stop("getIndFS supports a single RAM or LISREL model with its own raw data; multiple-group and mixture parent models are not supported.", call. = FALSE)
  }
  if (is.null(model@output$estimate)) {
    stop("model must be fitted before factor scoring.", call. = FALSE)
  }
  out <- mxFactorScores(model = .prepare_FS_model(model), type = FS_type)
  # Keep matrix dimensions for a single person or latent variable.
  estimates <- matrix(out[, , 1], nrow = dim(out)[1], ncol = dim(out)[2],
                      dimnames = dimnames(out)[1:2])
  errors <- matrix(out[, , 2], nrow = dim(out)[1], ncol = dim(out)[2],
                   dimnames = dimnames(out)[1:2])
  FSOut <- new("FSOutput", scores_est = estimates, scores_se = errors)
  return(FSOut)
}

## Variables referenced by a fitted model: manifest variables plus definition variables (`data.*` labels or
## algebra references).
.model_variables <- function(model){
  labels <- unlist(lapply(model@matrices, function(mat) mat$labels))
  formulas <- vapply(model@algebras, function(alg) paste(deparse(alg$formula), collapse = " "), character(1))
  text <- c(labels, formulas)
  defvars <- unlist(regmatches(text, gregexpr("data\\.[A-Za-z0-9_.]+", text)))
  unique(c(model@manifestVars, sub("^data\\.", "", defvars)))
}

## Reduce the stored raw data to the variables the model uses and give the model a scoring-safe name.
.prepare_FS_model <- function(model){
  observed <- model$data$observed
  if (!is.null(observed)){
    used <- intersect(.model_variables(model), colnames(observed))
    if (length(used) > 0 && length(used) < ncol(observed)){
      model <- mxModel(model, mxData(observed = observed[, used, drop = FALSE], type = "raw"))
    }
  }
  safe_name <- gsub("[^A-Za-z0-9_]", "", model@name)
  if (!nzchar(safe_name)) safe_name <- "nlpsemModel"
  if (!identical(safe_name, model@name)) model <- mxRename(model, safe_name)
  model
}

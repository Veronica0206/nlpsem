#' @title Summarize Model Fit Statistics for Fitted Models
#'
#' @description This function summarizes the model fit statistics for a list of fitted models. The summary includes the
#' number of parameters, estimated likelihood (-2ll), AIC, BIC, and other relevant statistics.
#'
#' @param model_list A list of fitted mxModel objects. Specifically, each element of the list should be the \code{mxOutput}
#' slot from the result returned by one of the estimation functions provided by this package.
#' @param HetModels A logical flag indicating whether a mixture model or a multiple group model is included in the list.
#' If set to \code{TRUE}, the function can also be used for the enumeration process, allowing the determination of the
#' optimal number of latent classes based on model fit statistics such as BIC. The default value is \code{FALSE}.
#'
#' @return A data frame summarizing model fit statistics (number of parameters, estimated likelihood, AIC, BIC, etc.)
#' for each model. A plain \code{data.frame} is returned for both \code{HetModels = TRUE} and \code{HetModels = FALSE}.
#'
#' @details Residual columns are combined by name, with NA where a model does not
#' contain that parameter. For heterogeneous models, each class has one proportion
#' column regardless of the number of longitudinal outcomes. Mixture proportions
#' count maximum-posterior assignments; multiple-group proportions count observations.
#' A summary table does not establish that its models are statistically comparable:
#' compare information criteria only for appropriate models fitted to the same data
#' and observed variables.
#'
#' @seealso \code{\link{getLGCM}}, \code{\link{getMIX}}, \code{\link{getMGroup}}
#'
#' @export
#'
#' @examples
#' # Fitted models shipped with the package: one-, two- and three-class bilinear spline models
#' load(system.file("extdata", "getMIX_examples.RData", package = "nlpsem"))
#' # Single model
#' getSummary(model_list = list(Math_BLS_LGCM1@mxOutput), HetModels = FALSE)
#' # Enumeration of latent classes
#' getSummary(
#'   model_list = list(Math_BLS_LGCM1@mxOutput, Math_BLS_LGCM2@mxOutput, Math_BLS_LGCM3@mxOutput),
#'   HetModels = TRUE
#' )
#'
#' @importFrom dplyr left_join
#' @importFrom tidyr separate spread pivot_wider
#' @importFrom stringr str_replace_all
#' @importFrom stats AIC BIC
#' @importFrom methods new
#'
getSummary <- function(model_list, HetModels = FALSE){
  if (!is.list(model_list) || !length(model_list) ||
      any(!vapply(model_list, methods::is, logical(1), class2 = "MxModel"))) {
    stop("'model_list' must be a nonempty list of fitted mxModel objects.", call. = FALSE)
  }
  if (!is.logical(HetModels) || length(HetModels) != 1L || is.na(HetModels)) {
    stop("'HetModels' must be TRUE or FALSE.", call. = FALSE)
  }
  rows <- lapply(seq_along(model_list), function(i) {
    model <- model_list[[i]]
    params <- omxGetParameters(model)
    likelihood <- model@output$Minus2LogLikelihood
    if (!is.numeric(likelihood) || length(likelihood) != 1L || !is.finite(likelihood)) {
      stop("Model ", i, " has no finite fitted -2 log-likelihood.", call. = FALSE)
    }
    out <- data.frame(Model = paste0("Model", i), No_Params = length(params),
                      `-2ll` = round(likelihood, 2), AIC = round(AIC(model), 2),
                      BIC = round(BIC(model), 2), check.names = FALSE)
    if (!HetModels) {
      residual_names <- grep("_res|_RES", names(params), value = TRUE)
      for (nm in residual_names) out[[nm]] <- round(unname(params[nm]), 4)
      return(out)
    }
    nClass <- length(model@submodels)
    if (!nClass) {
      proportion <- "100%"
      residual_names <- grep("_residuals$", names(params), value = TRUE)
      residual_columns <- sub("_residuals$", "_res_c1", residual_names)
    } else {
      if (is.null(model$data$observed)) {
        counts <- vapply(model@submodels, function(submodel) {
          observed <- submodel$data$observed
          if (is.null(observed)) return(0L)
          nrow(observed)
        }, integer(1))
        if (!sum(counts)) stop("Model ", i, " has no group observations.", call. = FALSE)
      } else {
        membership <- getPosterior(model, nClass = nClass)@membership
        counts <- tabulate(membership, nbins = nClass)
      }
      proportion <- paste0(counts / sum(counts) * 100, "%")
      residual_names <- grep("_residuals$", names(params), value = TRUE)
      classes <- sub("^(c[0-9]+).*", "\\1", residual_names)
      if (any(!classes %in% paste0("c", seq_len(nClass)))) {
        stop("Model ", i, " has residual labels without valid class prefixes.", call. = FALSE)
      }
      residual_columns <- paste0(sub("_residuals$", "_res_", sub("^c[0-9]+", "", residual_names)), classes)
    }
    for (j in seq_along(residual_names)) out[[residual_columns[j]]] <- round(unname(params[residual_names[j]]), 4)
    for (k in seq_along(proportion)) out[[paste0("%Class", k)]] <- proportion[k]
    out
  })
  # Union schemas rather than relying on positional rbind; every model retains
  # exactly one row, and absent residuals or classes are explicitly missing.
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    for (nm in setdiff(columns, names(row))) row[[nm]] <- NA
    row[columns]
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

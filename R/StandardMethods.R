# Define S4 Class
## S4 Class for the output structure for estimate functions, which include a slot for MxModel object
## and a slot for data frame of estimates when applicable.
#' @title S4 Class for optimized MxModel and point estimates with standard errors (when applicable)
#' @description S4 Class for the output structure for estimate functions.
#' @slot mxOutput An object of class "MxModel".
#' @slot Estimates A data frame of estimates.
#' @export
setClass("myMxOutput",
         slots = list(mxOutput = "MxModel", Estimates = "data.frame"))
#'
## S4 Class for the output structure output for the \code{getEstimateStats()} function, which contains
## slots for Wald, Likelihood, and Bootstrap estimates
#' @title S4 Class for p values and confidence intervals (when specified).
#' @description S4 Class for the output structure for the \code{getEstimateStats()} function.
#' @slot wald A data frame for p values and Wald confidence intervals (when specified).
#' @slot likelihood A data frame for Likelihood confidence intervals (when specified).
#' @slot bootstrap A data frame for Bootstrap confidence intervals (when specified).
#' @export
setClass("StatsOutput",
         slots = list(wald = "data.frame", likelihood = "data.frame", bootstrap = "data.frame"))
#'
## S4 Class for the output structure output for the \code{getIndFS()} function, which contains slots
## for estimated factor scores and their standard errors.
#' @title S4 Class for estimated factor scores and their standard errors.
#' @description S4 Class for the output structure for the \code{getIndFS()} function.
#' @slot scores_est A matrix of estimated factor scores.
#' @slot scores_se A matrix of standard errors of estimated factor scores.
#' @export
setClass("FSOutput",
         slots = list(scores_est = "matrix", scores_se = "matrix"))
#'
## S4 Class for the output structure output for the \code{getLatentKappa()}, which contains slots for
## kappa value, confidence interval, and judgment.
#' @title S4 Class for kappa statistic with confidence interval and judgment.
#' @description S4 Class for the output structure for the \code{getLatentKappa()} function.
#' @slot kappa_value A character vector for the kappa statistic with 95\% CI.
#' @slot judgment A character vector for the judgement for agreement.
#' @export
setClass("KappaOutput",
         slots = list(kappa_value = "character", judgment = "character"))
#'
## S4 Class for the output structure output for the \code{getPosterior()}, which contains slots for
## posterior probabilities, membership, entropy, and accuracy.
#' @title S4 Class for posterior probabilities, membership, entropy, and accuracy (when applicable)
#' @description S4 Class for the output structure for the \code{getPosterior()} function.
#' @slot prob A matrix of posterior probabilities.
#' @slot membership A numeric vector for membership.
#' @slot entropy A numeric value for entropy.
#' @slot accuracy A numeric value for accuracy.
#' @export
setClass("postOutput",
         slots = list(prob = "matrix", membership = "numeric", entropy = "numeric", accuracy = "numeric"))
#'
## S4 Class for the output structure output for the \code{getFigure()}, which contains slots for figures.
#' @title S4 Class for displaying figures
#' @description S4 Class to hold the figures output from the getFigure function.
#' @slot figures A list of lists containing figures for each specified sub_model and y_model (when
#' applicable).
#' @export
setClass("figOutput", slots = list(figures = "list"))
#'
########################################################################################################################
#'
# Define S4 Generic
## Define Generic for printing model summary of MxModel object. This Generic is only useful for an
## MxModel object.
#' @title S4 Generic for summarizing an optimized MxModel.
#' @description Generic function for printing model summary of MxModel object.
#' @param object An object of the appropriate class.
#' @return Invisibly returns the summary of the MxModel object. Called primarily for its
#' side effect of printing the model summary.
#' @examples
#' \donttest{
#' mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' baseT <- RMS_dat0$T1
#' RMS_dat0$T1 <- RMS_dat0$T1 - baseT
#' RMS_dat0$T2 <- RMS_dat0$T2 - baseT
#' RMS_dat0$T3 <- RMS_dat0$T3 - baseT
#' RMS_dat0$T4 <- RMS_dat0$T4 - baseT
#' RMS_dat0$T5 <- RMS_dat0$T5 - baseT
#' RMS_dat0$T6 <- RMS_dat0$T6 - baseT
#' RMS_dat0$T7 <- RMS_dat0$T7 - baseT
#' RMS_dat0$T8 <- RMS_dat0$T8 - baseT
#' RMS_dat0$T9 <- RMS_dat0$T9 - baseT
#' BLS_LGCM <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#'   records = 1:9
#' )
#' ModelSummary(BLS_LGCM)
#' }
#' @export
setGeneric("ModelSummary", function(object) {
  standardGeneric("ModelSummary")
})
#'
## Define Generic for printing output that are tables
#' @title S4 Generic for displaying output in a table format.
#' @description Generic function for printing output as tables.
#' @param object An object of the appropriate class.
#' @return Called primarily for its side effect of printing formatted output. Returns
#' the object invisibly.
#' @examples
#' \donttest{
#' mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' baseT <- RMS_dat0$T1
#' RMS_dat0$T1 <- RMS_dat0$T1 - baseT
#' RMS_dat0$T2 <- RMS_dat0$T2 - baseT
#' RMS_dat0$T3 <- RMS_dat0$T3 - baseT
#' RMS_dat0$T4 <- RMS_dat0$T4 - baseT
#' RMS_dat0$T5 <- RMS_dat0$T5 - baseT
#' RMS_dat0$T6 <- RMS_dat0$T6 - baseT
#' RMS_dat0$T7 <- RMS_dat0$T7 - baseT
#' RMS_dat0$T8 <- RMS_dat0$T8 - baseT
#' RMS_dat0$T9 <- RMS_dat0$T9 - baseT
#' BLS_LGCM <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#'   records = 1:9, paramOut = TRUE
#' )
#' printTable(BLS_LGCM)
#' }
#' @export
setGeneric("printTable", function(object){
  standardGeneric("printTable")
})
#'
########################################################################################################################
#'
# Define S4 Methods
## Define Method for printing model summary of MxModel object
#' @title S4 Method for summarizing an optimized MxModel.
#' @description Method for printing model summary of MxModel object.
#' @param object An object of class "myMxOutput".
#' @return Invisibly returns the summary of the MxModel object. Called primarily for its
#' side effect of printing the model summary.
#' @export
setMethod("ModelSummary", "myMxOutput", function(object){
  ### Print the summary of MxModel
  print(summary(object@mxOutput))
})
#'
## Define Method for printing point estimates and standard errors
#' @title S4 Method for printing point estimates with standard errors
#' @description Method for printing point estimates and standard errors.
#' @param object An object of class "myMxOutput".
#' @return Called for its side effect of printing point estimates and standard errors.
#' Returns the object invisibly.
#' @export
setMethod("printTable", "myMxOutput", function(object) {
  # Print the Estimates if they exist
  if (nrow(object@Estimates) > 0) {
    cat("\nPoint Estimates and Standard Errors:\n")
    print(object@Estimates)
  }
  else {
    message("To obtain Estimates, please set 'paramOut = TRUE' and specify parameters via 'names' argument!")
  }
})
## Define Method for printing p values and confidence intervals (when applicable)
#' @title S4 Method for printing p values and confidence intervals (when applicable)
#' @description Method for printing p values and confidence intervals.
#' @param object An object of class "StatsOutput".
#' @return Called for its side effect of printing p values and confidence intervals.
#' Returns the object invisibly.
#' @export
setMethod("printTable", "StatsOutput", function(object) {
  cat("Wald Estimates:\n")
  print(object@wald)
  if (nrow(object@likelihood) > 0) {
    cat("\nLikelihood Estimates:\n")
    print(object@likelihood)
  }
  if (nrow(object@bootstrap) > 0) {
    cat("\nBootstrap Estimates:\n")
    print(object@bootstrap)
  }
})
## Define Method for printing estimated factor scores and their standard errors
#' @title S4 Method for printing estimated factor scores and their standard errors
#' @description Method for printing estimated factor scores and their standard errors.
#' @param object An object of class "FSOutput".
#' @return Called for its side effect of printing factor scores and their standard errors.
#' Returns the object invisibly.
#' @export
setMethod("printTable", "FSOutput", function(object) {
  cat("\nEstimated factor scores:\n")
  print(object@scores_est)
  cat("\nStandard errors of estimated factor scores:\n")
  print(object@scores_se)
})
## Define Method for printing kappa statistic with 95% CI and judgement for agreement.
#' @title S4 Method for printing kappa statistic with 95\% CI and judgement for agreement.
#' @description Method for printing kappa statistic with 95\% CI and judgement for agreement.
#' @param object An object of class "KappaOutput".
#' @return Called for its side effect of printing kappa statistic with confidence interval
#' and agreement judgment. Returns the object invisibly.
#' @export
setMethod("printTable", "KappaOutput", function(object) {
  cat("\nKappa statistic (95% CI):\n")
  print(object@kappa_value)
  cat("\nJudgement for agreement:\n")
  print(object@judgment)
})
## Define Method for printing posterior probabilities, membership, entropy, and accuracy.
#' @title S4 Method for printing posterior probabilities, membership, entropy, and accuracy.
#' @description Method for printing posterior probabilities, membership, entropy, and accuracy.
#' @param object An object of class "postOutput".
#' @return Called for its side effect of printing posterior probabilities, membership,
#' entropy, and accuracy. Returns the object invisibly.
#' @export
setMethod("printTable", "postOutput", function(object) {
  cat("Posterior probabilities:\n")
  print(object@prob)
  cat("\nMembership:\n")
  print(object@membership)
  cat("\nEntropy:\n")
  print(object@entropy)
  if (length(object@accuracy) > 0) {
    cat("\nAccuracy:\n")
    print(object@accuracy)
  }
})
## Define Method for displaying figures.
#' @title S4 Method for displaying figures.
#' @description Method to display a summary of the figOutput object when printed.
#' @param object An object of class "figOutput".
#' @return Called for its side effect of displaying figures. Returns the object invisibly.
#' @export
#' @importFrom methods show setClass setGeneric setMethod
setMethod("show", "figOutput", function(object) {
  cat("figOutput Object\n")
  cat("--------------------\n")
  cat("Trajectories:", length(object@figures), "\n")
  for (j in seq_along(object@figures)) {
    if (length(object@figures) > 1) {
      cat("\nTrajectory", j, ":\n")
      prefix <- "  "
    } else {
      prefix <- ""
    }
    for (k in seq_along(object@figures[[j]])) {
      cat(prefix, "Figure ", k, ":\n", sep = "")
      print(object@figures[[j]][[k]])
    }
  }
})



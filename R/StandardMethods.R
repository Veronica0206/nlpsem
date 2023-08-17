#' @title Standard Methods (S4) for the package
#'
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
#' @slot kappa_value A character vector for the kappa statistic with $95\%$ CI.
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
#' @export
setGeneric("ModelSummary", function(object) {
  standardGeneric("ModelSummary")
})
#'
## Define Generic for printing output that are tables
#' @title S4 Generic for displaying output in a table format.
#' @description Generic function for printing output that are tables.
#' @param object An object of the appropriate class.
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
#' @export
setMethod("printTable", "myMxOutput", function(object) {
  # Print the Estimates if they exist
  if (!is.null(object@Estimates)) {
    cat("\nPoint Estimates and Standard Errors:\n")
    print(object@Estimates)
  }
  else if (is.null(object@Estimates)) {
    cat("\nTo obtain Estimates, please set 'paramOut = T' and specify parameters via 'names' argument!\n")
  }
})
## Define Method for printing p values and confidence intervals (when applicable)
#' @title S4 Method for printing p values and confidence intervals (when applicable)
#' @description Method for printing p values and confidence intervals.
#' @param object An object of class "StatsOutput".
#' @export
setMethod("printTable", "StatsOutput", function(object) {
  cat("Wald Estimates:\n")
  print(object@wald)
  if (!is.null(object@likelihood)) {
    cat("\nLikelihood Estimates:\n")
    print(object@likelihood)
  }
  if (!is.null(object@bootstrap)) {
    cat("\nBootstrap Estimates:\n")
    print(object@bootstrap)
  }
})
## Define Method for printing estimated factor scores and their standard errors
#' @title S4 Method for printing estimated factor scores and their standard errors
#' @description Method for printing estimated factor scores and their standard errors.
#' @param object An object of class "FSOutput".
#' @export
setMethod("printTable", "FSOutput", function(object) {
  cat("\nEstimated factor scores:\n")
  print(object@scores_est)
  cat("\nStandard errors of estimated factor scores:\n")
  print(object@scores_se)
})
## Define Method for printing kappa statistic with 95% CI and judgement for agreement.
#' @title S4 Method for printing kappa statistic with $95\%$ CI and judgement for agreement.
#' @description Method for printing kappa statistic with $95\%$ CI and judgement for agreement.
#' @param object An object of class "KappaOutput".
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
#' @export
setMethod("printTable", "postOutput", function(object) {
  cat("Posterior probabilities:\n")
  print(object@prob)
  cat("\nMembership:\n")
  print(object@membership)
  cat("\nEntropy:\n")
  print(object@entropy)
  if (!is.null(object@accuracy)) {
    cat("\nAccuracys:\n")
    print(object@accuracy)
  }
})
## Define Method for displaying figures.
#' @title S4 Method for displaying figures.
#' @description Method to display a summary of the figOutput object when printed.
#' @param object An object of class "figOutput".
#' @export
#' @importFrom methods show
setMethod("show", "figOutput", function(object) {
  cat("figOutput Object\n")
  cat("--------------------\n")
  cat("Trajectories:", length(object@figures), "\n")
  # Check the length of object@figures
  if (length(object@figures) == 1) {
    # Check the length of the first element of object@figures
    if (length(object@figures[[1]]) == 1) {
      cat("Figure 1:\n")
      print(object@figures[[1]][[1]])
    } else { # Assuming length is 2
      cat("Figure 1:\n")
      print(object@figures[[1]][[1]])
      cat("Figure 2:\n")
      print(object@figures[[1]][[2]])
    }
  } else { # If length of object@figures is not 1
    for (j in seq_along(object@figures)) {
      cat("\nTrajectory", j, ":\n")
      # Check the length of the j-th element of object@figures
      if (length(object@figures[[j]]) == 1) {
        cat("  Figure 1:\n")
        print(object@figures[[j]][[1]])
      } else { # Assuming length is 2
        cat("  Figure 1:\n")
        print(object@figures[[j]][[1]])
        cat("  Figure 2:\n")
        print(object@figures[[j]][[2]])
      }
    }
  }
})



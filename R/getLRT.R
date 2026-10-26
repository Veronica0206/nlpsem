#' @title Perform Bootstrap Likelihood Ratio Test for Comparing Full and Reduced Models
#'
#' @description Wraps OpenMx model comparison for a fitted full model and a
#' corresponding reduced model, with optional parametric bootstrapping. A returned
#' p-value is meaningful only when the comparison satisfies the required assumptions.
#'
#' @param full A fitted mxModel object for the full model. Specifically, this should be the \code{mxOutput} slot from the result
#' returned by one of the estimation functions provided by this package.
#' @param reduced A fitted mxModel object for the reduced model. Specifically, this should be the \code{mxOutput} slot from the
#' result returned by one of the estimation functions provided by this package.
#' @param boot A logical flag indicating whether to perform bootstrapping for the comparison. Default is \code{FALSE}.
#' @param rep An integer specifying the number of bootstrap replications if \code{boot} is \code{TRUE}. Default is \code{NA}.
#'
#' @return A data frame containing the number of free parameters, estimated likelihood (-2ll), degrees of freedom, differences
#' in log-likelihood and degrees of freedom, p-values, AIC, and BIC for both the full and reduced models.
#'
#' @details Both models must use the same observations and observed variables,
#' and the reduced model must be nested in the full model. Neither this wrapper
#' nor OpenMx establishes mathematical nesting automatically. Inspect
#' \code{getDiagnostics()} and check solution stability before interpreting results.
#' The ordinary chi-square reference requires the usual likelihood regularity
#' conditions; boundary parameters, unidentified parameters, and latent class
#' enumeration can violate these conditions.
#'
#' In particular, random-knot and fixed-knot bilinear spline fits use different
#' smooth/absolute-value loading formulas, and latent covariances are not generally
#' constrained to be positive semidefinite. Their automatic comparison is not a
#' validated significance test for a random-knot effect. Establish the null model,
#' admissible parameter space and calibration for the intended analysis.
#'
#' \code{boot = TRUE} delegates to \code{OpenMx::mxCompare()}. It does not
#' repair invalid models or non-nesting. Use enough successful replications to
#' control Monte Carlo error, and inspect bootstrap fitting behavior. A reported
#' zero bootstrap p-value is bounded by the simulation resolution, not evidence
#' of a probability exactly equal to zero.
#'
#' @seealso \code{\link{getSummary}}
#'
#' @export
#'
#' @examples
#' # Historical comparison output: inspect the assumptions in Details before
#' # interpreting the p-value for these nonlinear models.
#' load(system.file("extdata", "getLGCM_examples.RData", package = "nlpsem"))
#' comparison <- getLRT(full = Math_LGCM_BLS_f@mxOutput,
#'                      reduced = Math_LGCM_BLS_r@mxOutput)
#' comparison[, c("# of Free Param", "-2loglik", "AIC", "BIC")]
#'
#' @importFrom OpenMx mxCompare
#'
getLRT <- function(full, reduced, boot = FALSE, rep = NA){
  if (!is.logical(boot) || length(boot) != 1L || is.na(boot)) {
    stop("boot must be a single TRUE or FALSE.", call. = FALSE)
  }
  if (boot && (!is.numeric(rep) || length(rep) != 1 || !is.finite(rep) ||
               rep < 1 || rep != floor(rep) || rep > .Machine$integer.max)){
    stop("Please specify bootstrap replications as a positive integer!")
  }
  LRT_out <- mxCompare(full, reduced, boot = boot, replications = rep)[, c(3:5, 7:9, 6)]
  rownames(LRT_out) <- c("Full Model", "Reduced Model")
  colnames(LRT_out) <- c("# of Free Param", "-2loglik", "Degree of Freedom", "Diff in loglik",
                         "Diff in DoF", "p.values", "AIC")
  # p.values stays numeric (0.4 replaced small values by the string "<0.0001", which turned the column into text)
  LRT_out$BIC <- c(BIC(full), BIC(reduced))
  return(LRT_out)
}

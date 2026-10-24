#' @title Perform Bootstrap Likelihood Ratio Test for Comparing Full and Reduced Models
#'
#' @description This function performs the likelihood ratio test (LRT) to compare a full model (an intrinsically nonlinear
#' longitudinal model) with a corresponding parsimonious alternative (a non-intrinsically nonlinear longitudinal model).
#' It also provides an option to perform bootstrapping for the comparison.
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
#' @seealso \code{\link{getSummary}}
#'
#' @export
#'
#' @examples
#' # Fitted bilinear spline models shipped with the package: random knot (full) versus
#' # fixed knot (reduced)
#' load(system.file("extdata", "getLGCM_examples.RData", package = "nlpsem"))
#' comparison <- getLRT(full = Math_LGCM_BLS_f@mxOutput,
#'                      reduced = Math_LGCM_BLS_r@mxOutput)
#' stopifnot(comparison["Reduced Model", "Diff in DoF"] == 4,
#'           is.finite(comparison["Reduced Model", "p.values"]))
#' comparison
#'
#' @importFrom OpenMx mxCompare
#'
getLRT <- function(full, reduced, boot = FALSE, rep = NA){
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

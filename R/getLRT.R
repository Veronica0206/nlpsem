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
#' @export
#'
#' @examples
#'
#' mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' # Load ECLS-K (2011) data
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' # Re-baseline the data so that the estimated initial status is for the starting point of the study
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
#'
#' \donttest{
#' # Fit bilinear spline growth model with random knot (intrinsically nonlinear model)
#' BLS_LGCM_f <- getLGCM(dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#'                       intrinsic = TRUE, records = 1:9, res_scale = 0.1)
#' # Fit bilinear spline growth model with fix knot (non-intrinsically nonlinear model)
#' BLS_LGCM_r <- getLGCM(dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#'                       intrinsic = FALSE, records = 1:9, res_scale = 0.1)
#' # Likelihood ratio test
#' getLRT(full = BLS_LGCM_f@mxOutput, reduced = BLS_LGCM_r@mxOutput, boot = FALSE, rep = NA)
#' }
#'
#' @importFrom OpenMx mxCompare
#'
getLRT <- function(full, reduced, boot = FALSE, rep = NA){
  if (I(boot & is.na(rep))){
    stop("Please specify bootstrap replications!")
  }
  LRT_out <- mxCompare(full, reduced, boot = boot, replications = rep)[, c(3:5, 7:9, 6)]
  rownames(LRT_out) <- c("Full Model", "Reduced Model")
  colnames(LRT_out) <- c("# of Free Param", "-2loglik", "Degree of Freedom", "Diff in loglik",
                         "Diff in DoF", "p.values", "AIC")
  LRT_out$p.values[2] <- ifelse(LRT_out$p.values[2] < 0.0001, "<0.0001", round(LRT_out$p.values[2], 4))
  LRT_out$BIC <- c(BIC(full), BIC(reduced))
  return(LRT_out)
}

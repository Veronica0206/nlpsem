#' @title Derive Individual Factor Scores for Each Latent Variable Included in Model
#'
#' @description This function computes individual factor scores for each latent variable in a given model. It supports three
#' types of factor scores: maximum likelihood, weighted maximum likelihood, and regression.
#'
#' @param model A fitted mxModel object. This is the output from one of the estimation functions in this
#' package.
#' @param FS_type A string specifying the type of factor scores to compute. Supported options include \code{"ML"} (for
#' Maximum Likelihood), \code{"WeightedML"} (for Weighted Maximum Likelihood), and \code{"Regression"}. Default is \code{"Regression"}.
#'
#' @return A list containing two elements: scores_est, the factor score estimates, and scores_se, the standard errors of the
#' factor score estimates.
#'
#' @references
#' \itemize{
#'   \item {Estabrook, R. & Neale, M. C. (2013). A Comparison of Factor Score Estimation Methods in the Presence of Missing
#'   Data: Reliability and an Application to Nicotine Dependence. Multivariate Behavioral Research, 48, 1-27. \doi{10.1080/00273171.2012.730072}}
#'   \item {Priestley, M. & Subba Rao, T. (1975). The Estimation of Factor Scores and Kalman Filtering For Discrete Parameter
#'   Stationary Processes. International Journal of Control, 21, 971-975. \doi{10.1080/00207177508922050}}
#' }
#'
#' @export
#'
#' @examples
#' OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
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
#' # Standardized time-invariant covariates
#' RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#' RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
#' \donttest{
#' # Fit bilinear spline latent growth curve model (fixed knots)
#' LIN_LGCM <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "linear",
#'   intrinsic = FALSE, records = 1:9, growth_TIC = NULL, res_scale = 0.1
#' )
#' getIndFS(model = LIN_LGCM, FS_type = "Regression")
#' # Fit bilinear spline latent growth curve model (random knots) with time-invariant covariates for
#' # mathematics development
#' ## Fit the model
#' BLS_LGCM.TIC_f <- getLGCM(dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS",
#'                           intrinsic = TRUE, records = 1:9, growth_TIC = c("ex1", "ex2"),
#'                           res_scale = 0.1)
#' ## Output point estimate and standard errors
#' getIndFS(model = BLS_LGCM.TIC_f, FS_type = "Regression")
#' }
#'
#' @importFrom OpenMx mxFactorScores
#'
getIndFS <- function(model, FS_type = "Regression"){
  out <- mxFactorScores(model = model, type = FS_type)
  return(list(scores_est = out[, , 1], scores_se <- out[, , 2]))
}

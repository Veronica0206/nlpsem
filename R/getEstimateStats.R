#' @title Calculate p-Values and Confidence Intervals of Parameters for a Fitted Model
#'
#' @description This function calculates p-values and confidence intervals (CIs) of parameters for a given model. It supports
#' different types of CIs, including Wald CIs, likelihood-based CIs, bootstrap CIs, or all three.
#'
#' @param model A fitted mxModel object. Specifically, this should be the \code{mxOutput} slot from the result returned by
#' one of the estimation functions provided by this package. The default value is \code{NULL}. Providing this parameter is
#' essential when generating likelihood-based and bootstrap confidence intervals (CIs).
#' @param est_in The \code{Estimates} slot from the result returned by one of the estimation functions provided by this
#' package, which contains a dataframe with point estimates and standard errors.
#' @param p_values A logical flag indicating whether to calculate p-values. Default is \code{TRUE}.
#' @param CI A logical flag indicating whether to compute confidence intervals. Default is \code{TRUE}.
#' @param CI_type A string specifying the type of confidence interval to compute. Supported options include
#' \code{"Wald"}, \code{"likelihood"}, \code{"bootstrap"}, or \code{"all"}. Default is \code{"Wald"}.
#' @param rep An integer specifying the number of replications for bootstrap. This is applicable if \code{CI_type} is
#' \code{"bootstrap"} or \code{"all"}. Default is \code{NA}.
#' @param conf.level A numeric value representing the confidence level for confidence interval calculation. Default is
#' \code{0.95}.
#'
#' @return An object of class \code{StatsOutput} with potential slots:
#' \itemize{
#'   \item \code{wald}: Contains a data frame with, point estimates, standard errors, p-values, and Wald confidence intervals
#'   (when specified).
#'   \item \code{likelihood}: Contains a data frame with likelihood-based confidence intervals (when specified).
#'   \item \code{bootstrap}: Contains a data frame with bootstrap confidence intervals (when specified).
#' }
#' The content of these slots can be printed using the \code{printTable()} method for S4 objects.
#'
#' @references
#' \itemize{
#'   \item {Casella, G. & Berger, R.L. (2002). Statistical Inference (2nd ed.). Duxbury Press.}
#'   \item {Madansky, A. (1965). Approximate Confidence Limits for the Reliability of Series and Parallel Systems.
#'   Technometrics, 7(4), 495-503. Taylor & Francis, Ltd. \url{https://www.jstor.org/stable/1266390}}
#'   \item {Matthews, D. E. (1988). Likelihood-Based Confidence Intervals for Functions of Many Parameters.
#'   Biometrika, 75(1), 139-144. Oxford University Press. \url{https://www.jstor.org/stable/2336444}}
#'   \item {Efron, B. & Tibshirani, R. J. (1994). An Introduction to the Bootstrap. CRC press.}
#' }
#'
#' @seealso \code{\link{getLGCM}}, \code{\link{getIndFS}}, \code{\link{getPosterior}}
#'
#' @export
#'
#' @examples
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
#' # Standardized time-invariant covariates
#' RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#' RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
#' \donttest{
#' # Fit bilinear spline latent growth curve model (fixed knots)
#' BLS_LGCM_r <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#'   records = 1:9, paramOut = TRUE)
#' ## Generate P value and Wald confidence intervals
#' getEstimateStats(
#'   est_in = BLS_LGCM_r@Estimates, CI_type = "Wald"
#'   )
#' }
#' \dontrun{
#' # Fit bilinear spline latent growth curve model (random knots) with time-invariant covariates for
#' # mathematics development
#' ## Fit the model
#' set.seed(20191029)
#' BLS_LGCM.TIC_f <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = TRUE, records = 1:9,
#'   growth_TIC = c("ex1", "ex2"), tries = 20, paramOut = TRUE
#'   )
#' ## Change optimizer to "SLSQP" for getting likelihood-based confidence interval
#' mxOption(model = NULL, key = "Default optimizer", "SLSQP", reset = FALSE)
#' ## Generate P value and all three types of confidence intervals
#' getEstimateStats(
#'   model = BLS_LGCM.TIC_f@mxOutput, est_in = BLS_LGCM.TIC_f@Estimates, CI_type = "all", rep = 20
#'   )
#' }
#'
#' @importFrom OpenMx omxGetParameters mxTryHard mxModel mxCI mxBootstrap
#' @importFrom dplyr %>% select left_join
#' @importFrom stats quantile
#' @importFrom methods new
#'
getEstimateStats <- function(model = NULL, est_in, p_values = TRUE, CI = TRUE, CI_type = "Wald", rep = NA,
                             conf.level = 0.95){
  est_out_wald <- est_in[, -1]
  rownames(est_out_wald) <- est_in$Name
  if (p_values){
    se <- est_in[, 3]
    z_ratio <- ifelse(se == 0, NA_real_, est_in[, 2] / se)
    p_value <- round((1 - pnorm(abs(z_ratio), lower.tail = TRUE)) * 2, 4)
    est_out_wald$p.value <- ifelse(is.na(p_value), NA_character_,
                                   ifelse(p_value < 0.0001, "<0.0001",
                                          ifelse(p_value > 0.9999, ">0.9999", p_value)))
  }
  # Helper: compute Wald CIs
  .add_wald_CI <- function(est_out, est_in, conf.level) {
    z <- qnorm(1 - (1 - conf.level)/2)
    est_out$wald_lbound <- round(est_in[, 2] - z * est_in[, 3], 4)
    est_out$wald_ubound <- round(est_in[, 2] + z * est_in[, 3], 4)
    est_out
  }
  # Helper: compute likelihood-based CIs
  .get_likelihood_CI <- function(model) {
    free_para <- omxGetParameters(model = model)
    model_CI <- mxTryHard(mxModel(model = model, mxCI(names(free_para))), intervals = TRUE)
    est_out_lik <- as.data.frame(round(model_CI$output$confidenceIntervals, 4))[, c(2, 1, 3)]
    names(est_out_lik) <- c("Estimate", "lik_lbound", "lik_ubound")
    est_out_lik
  }
  # Helper: compute bootstrap CIs
  .get_bootstrap_CI <- function(model, rep) {
    boot <- mxBootstrap(model, replications = rep)
    boot_raw <- boot$compute$output$raw %>%
      select(-seed, -fit, -statusCode)
    boot_CI <- as.data.frame(t(sapply(boot_raw, function(x) quantile(x, c(0.025, 0.975)))))
    boot_CI$name <- row.names(boot_CI)
    boot_est <- summary(boot)$parameters
    est_out_boot <- data.frame(round(boot_est %>% select(name, Estimate) %>% left_join(boot_CI, by = c("name" = "name")) %>%
                                       select(-name), 4))
    names(est_out_boot) <- c("Estimate", "boot_lbound", "boot_ubound")
    rownames(est_out_boot) <- row.names(boot_CI)
    est_out_boot
  }
  if (CI){
    if (CI_type == "Wald"){
      est_out_wald <- .add_wald_CI(est_out_wald, est_in, conf.level)
      statOut <- new("StatsOutput", wald = est_out_wald)
    }
    else if (CI_type == "likelihood"){
      if (is.null(model)) stop("A fitted mxModel object must be provided via 'model' for likelihood-based CIs.")
      est_out_lik <- .get_likelihood_CI(model)
      statOut <- new("StatsOutput", likelihood = est_out_lik)
    }
    else if (CI_type == "bootstrap"){
      if (is.null(model)) stop("A fitted mxModel object must be provided via 'model' for bootstrap CIs.")
      est_out_boot <- .get_bootstrap_CI(model, rep)
      statOut <- new("StatsOutput", bootstrap = est_out_boot)
    }
    else if (CI_type == "all"){
      if (is.null(model)) stop("A fitted mxModel object must be provided via 'model' for likelihood-based and bootstrap CIs.")
      est_out_wald <- .add_wald_CI(est_out_wald, est_in, conf.level)
      est_out_lik <- .get_likelihood_CI(model)
      est_out_boot <- .get_bootstrap_CI(model, rep)
      statOut <- new("StatsOutput", wald = est_out_wald, likelihood = est_out_lik, bootstrap = est_out_boot)
    }
  }
  else {
    statOut <- new("StatsOutput", wald = est_out_wald)
  }
  return(statOut)
}

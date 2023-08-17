#' @title Calculate p-Values and Confidence Intervals of Parameters for a Fitted Model
#'
#' @description This function calculates p-values and confidence intervals (CIs) of parameters for a given model.
#' It supports different types of CIs, including Wald CIs, likelihood-based CIs, bootstrap CIs, or all three.
#'
#' @param model A fitted mxModel object. This is the output from one of the estimation functions in this package.
#' Default is \code{NULL}. This is only required when generating likelihood-based and bootstrap CIs.
#' @param est_in A data frame containing input estimates.
#' @param p_values A logical flag indicating whether to calculate p-values. Default is \code{TRUE}.
#' @param CI A logical flag indicating whether to compute confidence intervals. Default is \code{TRUE}.
#' @param CI_type A string specifying the type of confidence interval to compute. Supported options include
#' \code{"Wald"}, \code{"likelihood"}, \code{"bootstrap"}, or \code{"all"}. Default is \code{"Wald"}.
#' @param rep An integer specifying the number of replications for bootstrap. This is applicable if \code{CI_type} is
#' \code{"bootstrap"} or \code{"all"}. Default is \code{NA}.
#' @param conf.level A numeric value representing the confidence level for confidence interval calculation. Default is
#' \code{0.95}.
#'
#' @return A data frame with calculated statistics (p-value, confidence intervals) added to the input estimates.
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
#' # Fit bilinear spline latent growth curve model (fixed knots)
#' \donttest{
#' paraBLS_LGCM.r <- c(
#'   "mueta0", "mueta1", "mueta2", "knot",
#'   paste0("psi", c("00", "01", "02", "11", "12", "22")),
#'   "residuals"
#'   )
#' BLS_LGCM_r <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#'   records = 1:9, res_scale = 0.1, paramOut = TRUE, names = paraBLS_LGCM.r)
#' ## Output point estimate and standard errors
#' getEstimateStats(
#'   est_in = BLS_LGCM_r@Estimates, CI_type = "Wald"
#'   )
#' # Fit bilinear spline latent growth curve model (random knots) with time-invariant covariates for
#' # mathematics development
#' ## Define parameter names
#' paraBLS.TIC_LGCM.f <- c(
#'   "alpha0", "alpha1", "alpha2", "alphag",
#'   paste0("psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")), "residuals",
#'   paste0("beta1", c(0:2, "g")), paste0("beta2", c(0:2, "g")), paste0("mux", 1:2),
#'   paste0("phi", c("11", "12", "22")), "mueta0", "mueta1", "mueta2", "mu_knot"
#'   )
#' ## Fit the model
#' BLS_LGCM.TIC_f <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = TRUE, records = 1:9,
#'   growth_TIC = c("ex1", "ex2"), res_scale = 0.1, paramOut = TRUE, names = paraBLS.TIC_LGCM.f
#'   )
#' ## Output point estimate and standard errors
#' getEstimateStats(
#'   model = BLS_LGCM.TIC_f@mxOutput, est_in = BLS_LGCM.TIC_f@Estimates, CI_type = "all", rep = 1000
#'   )
#' }
#' @importFrom OpenMx omxGetParameters mxTryHard mxModel mxCI mxBootstrap
#' @importFrom dplyr %>% select summarize_all
#' @importFrom stats quantile
#' @importFrom methods new
#'
getEstimateStats <- function(model = NULL, est_in, p_values = TRUE, CI = TRUE, CI_type = "Wald", rep = NA,
                             conf.level = 0.95){
  est_out_wald <- est_in[, -1]
  rownames(est_out_wald) <- est_in$Name
  if (p_values){
    p_value <- round((1 - pnorm(abs(est_in[, 2]/est_in[, 3]), lower.tail = T)) * 2, 4)
    est_out_wald$p.value <- ifelse(p_value < 0.0001, "<0.0001", ifelse(p_value > 0.9999, ">0.9999", p_value))
  }
  if (CI){
    if (CI_type == "Wald"){
      est_out_wald$wald_lbound <- round(est_in[, 2] - qnorm(1 - (1 - conf.level)/2) * est_in[, 3], 4)
      est_out_wald$wald_ubound <- round(est_in[, 2] + qnorm(1 - (1 - conf.level)/2) * est_in[, 3], 4)
      statOut <- new("StatsOutput", wald = est_out_wald)
    }
    else if (CI_type == "likelihood"){
      free_para <- omxGetParameters(model = model)
      model_CI <- mxTryHard(mxModel(model = model, mxCI(names(free_para))), intervals = TRUE)
      est_out_lik <- as.data.frame(round(model_CI$output$confidenceIntervals, 4))[, c(2, 1, 3)]
      names(est_out_lik) <- c("Estimate", "lik_lbound", "lik_ubound")
      statOut <- new("StatsOutput", likelihood = est_out_lik)
    }
    else if (CI_type == "bootstrap"){
      boot <- mxBootstrap(model, replications = rep)
      boot_CI  <- as.data.frame(t(boot$compute$output$raw %>%
                                    select(-seed, -fit, -statusCode) %>%
                                    summarize_all(function(x) quantile(x, c(0.025, 0.975)))))
      boot_CI$name <- row.names(boot_CI)
      boot_est <- summary(boot)$parameters
      est_out_boot <- data.frame(round(boot_est %>% select(name, Estimate) %>% left_join(boot_CI, by = c("name" = "name")) %>%
                                         select(-name), 4))
      names(est_out_boot) <- c("Estimate", "boot_lbound", "boot_ubound")
      rownames(est_out_boot) <- row.names(boot_CI)
      statOut <- new("StatsOutput", bootstrap = est_out_boot)
    }
    else if (CI_type == "all"){
      est_out_wald$wald_lbound <- round(est_in[, 2] - qnorm(1 - (1 - conf.level)/2) * est_in[, 3], 4)
      est_out_wald$wald_ubound <- round(est_in[, 2] + qnorm(1 - (1 - conf.level)/2) * est_in[, 3], 4)
      free_para <- omxGetParameters(model = model)
      model_CI <- mxTryHard(mxModel(model = model, mxCI(names(free_para))), intervals = TRUE)
      est_out_lik <- as.data.frame(round(model_CI$output$confidenceIntervals, 4))[, c(2, 1, 3)]
      names(est_out_lik) <- c("Estimate", "lik_lbound", "lik_ubound")
      boot <- mxBootstrap(model, replications = rep)
      boot_CI  <- as.data.frame(t(boot$compute$output$raw %>%
                                    select(-seed, -fit, -statusCode) %>%
                                    summarize_all(function(x) quantile(x, c(0.025, 0.975)))))
      boot_CI$name <- row.names(boot_CI)
      boot_est <- summary(boot)$parameters
      est_out_boot <- data.frame(round(boot_est %>% select(name, Estimate) %>% left_join(boot_CI, by = c("name" = "name")) %>%
                                         select(-name), 4))
      rownames(est_out_boot) <- row.names(boot_CI)
      names(est_out_boot) <- c("Estimate", "boot_lbound", "boot_ubound")
      statOut <- new("StatsOutput", wald = est_out_wald, likelihood = est_out_lik, bootstrap = est_out_boot)
    }
  }
  else if (!CI){
    statOut <- new("StatsOutput", wald = est_out_wald)
  }
  return(statOut)
}


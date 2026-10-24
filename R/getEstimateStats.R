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
#' \code{"Wald"}, \code{"likelihood"}, \code{"bootstrap"}, or \code{"all"}. Default is \code{"Wald"}. An
#' unrecognized value triggers an informative error.
#' @param rep An integer specifying the number of replications for bootstrap. A positive integer must be
#' supplied when \code{CI_type} is \code{"bootstrap"} or \code{"all"}. Default is \code{NA}.
#' @param conf.level A numeric value representing the confidence level for confidence interval calculation. It applies
#' to Wald, likelihood-based, and bootstrap confidence intervals. Default is \code{0.95}.
#'
#' @return An object of class \code{StatsOutput} with potential slots:
#' \itemize{
#'   \item \code{wald}: Contains a data frame with, point estimates, standard errors, p-values, and Wald confidence intervals
#'   (when specified).
#'   \item \code{likelihood}: Contains a data frame with likelihood-based confidence intervals (when specified).
#'   \item \code{bootstrap}: Contains a data frame with bootstrap confidence intervals (when specified).
#' }
#' The content of these slots can be printed using the \code{printTable()} method for S4 objects.
#' Bootstrap intervals use only replications with an acceptable OpenMx status and finite
#' parameter estimates. Discarded replications trigger a warning; fewer than two usable
#' replications trigger an error. The bootstrap data frame has a \code{bootstrap_diagnostics}
#' attribute recording requested, attempted, successful and discarded counts and status frequencies.
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
#' # Fitted bilinear spline latent growth curve models shipped with the package
#' load(system.file("extdata", "getLGCM_examples.RData", package = "nlpsem"))
#' # p values and Wald confidence intervals
#' stats_wald <- getEstimateStats(est_in = Math_LGCM_TIC_BLS_f@Estimates, CI_type = "Wald")
#' printTable(stats_wald)
#' \dontrun{
#' # Likelihood-based and bootstrap intervals re-optimize the model; this takes several minutes
#' mxOption(model = NULL, key = "Default optimizer", "SLSQP", reset = FALSE)
#' getEstimateStats(
#'   model = Math_LGCM_TIC_BLS_f@mxOutput, est_in = Math_LGCM_TIC_BLS_f@Estimates,
#'   CI_type = "all", rep = 20
#' )
#' }
#'
#' @importFrom OpenMx omxGetParameters mxTryHard mxModel mxCI mxBootstrap
#' @importFrom dplyr %>% select left_join
#' @importFrom stats quantile
#' @importFrom methods new
#'
getEstimateStats <- function(model = NULL, est_in, p_values = TRUE, CI = TRUE, CI_type = "Wald", rep = NA,
                             conf.level = 0.95){
  CI_type <- match.arg(CI_type, c("Wald", "likelihood", "bootstrap", "all"))
  if (CI) .ci_probs(conf.level)
  if (CI && CI_type %in% c("bootstrap", "all")){
    if (!is.numeric(rep) || length(rep) != 1 || !is.finite(rep) ||
        rep < 1 || rep != floor(rep) || rep > .Machine$integer.max){
      stop("A positive integer must be provided via 'rep' when bootstrap CIs are requested (e.g., rep = 1000).")
    }
  }
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
  .get_likelihood_CI <- function(model, conf.level) {
    free_para <- omxGetParameters(model = model)
    model_CI <- mxTryHard(mxModel(model = model, mxCI(names(free_para), interval = conf.level)), intervals = TRUE)
    est_out_lik <- as.data.frame(round(model_CI$output$confidenceIntervals, 4))[, c(2, 1, 3)]
    names(est_out_lik) <- c("Estimate", "lik_lbound", "lik_ubound")
    est_out_lik
  }
  # Helper: compute bootstrap CIs
  .get_bootstrap_CI <- function(model, rep, conf.level) {
    boot <- mxBootstrap(model, replications = rep)
    .bootstrap_percentile_CI(boot$compute$output$raw, boot$compute$OK,
                             omxGetParameters(model), rep, conf.level)
  }
  if (CI){
    if (CI_type == "Wald"){
      est_out_wald <- .add_wald_CI(est_out_wald, est_in, conf.level)
      statOut <- new("StatsOutput", wald = est_out_wald)
    }
    else if (CI_type == "likelihood"){
      if (is.null(model)) stop("A fitted mxModel object must be provided via 'model' for likelihood-based CIs.")
      est_out_lik <- .get_likelihood_CI(model, conf.level)
      statOut <- new("StatsOutput", likelihood = est_out_lik)
    }
    else if (CI_type == "bootstrap"){
      if (is.null(model)) stop("A fitted mxModel object must be provided via 'model' for bootstrap CIs.")
      est_out_boot <- .get_bootstrap_CI(model, rep, conf.level)
      statOut <- new("StatsOutput", bootstrap = est_out_boot)
    }
    else if (CI_type == "all"){
      if (is.null(model)) stop("A fitted mxModel object must be provided via 'model' for likelihood-based and bootstrap CIs.")
      est_out_wald <- .add_wald_CI(est_out_wald, est_in, conf.level)
      est_out_lik <- .get_likelihood_CI(model, conf.level)
      est_out_boot <- .get_bootstrap_CI(model, rep, conf.level)
      statOut <- new("StatsOutput", wald = est_out_wald, likelihood = est_out_lik, bootstrap = est_out_boot)
    }
  }
  else {
    statOut <- new("StatsOutput", wald = est_out_wald)
  }
  return(statOut)
}

## Lower and upper tail probabilities for a two-sided interval at the requested confidence level.
.ci_probs <- function(conf.level){
  if (!is.numeric(conf.level) || length(conf.level) != 1 || !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1){
    stop("conf.level must be a single number strictly between 0 and 1.")
  }
  c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
}

# Match by parameter name, and use the same accepted-status set as the OpenMx
# bootstrap compute plan. Never let a failed fit or non-finite estimate enter a quantile.
.bootstrap_percentile_CI <- function(raw, accepted, estimates, requested, conf.level) {
  probs <- .ci_probs(conf.level)
  pars <- names(estimates)
  if (!is.data.frame(raw) || !length(pars) ||
      !all(c("statusCode", pars) %in% names(raw)) || !length(accepted)) {
    stop("Bootstrap results are missing parameter estimates or convergence statuses.")
  }
  values <- as.matrix(raw[, pars, drop = FALSE])
  usable <- as.character(raw$statusCode) %in% as.character(accepted) &
    rowSums(!is.finite(values)) == 0
  diagnostics <- list(requested = requested, attempted = nrow(raw),
                      successful = sum(usable), discarded = sum(!usable),
                      status_counts = table(raw$statusCode, useNA = "ifany"))
  if (sum(usable) < 2L) {
    stop(sprintf("Only %d of %d bootstrap replications were usable; at least two are needed for intervals.",
                 sum(usable), nrow(raw)))
  }
  if (any(!usable)) {
    warning(sprintf("Discarded %d of %d bootstrap replications with unacceptable status or non-finite estimates; intervals use %d successful replications.",
                    sum(!usable), nrow(raw), sum(usable)), call. = FALSE)
  }
  bounds <- vapply(seq_along(pars), function(j) quantile(values[usable, j], probs), numeric(2))
  out <- data.frame(Estimate = round(unname(estimates), 4),
                    boot_lbound = round(bounds[1L, ], 4),
                    boot_ubound = round(bounds[2L, ], 4), row.names = pars)
  attr(out, "bootstrap_diagnostics") <- diagnostics
  out
}

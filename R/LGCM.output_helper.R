#' @title Extract Point Estimates and Standard Errors of Latent Growth Curve Model with Time-invariant Covariates (If Any)
#'
#' @description This function computes and returns a data frame containing point estimates and standard errors for the parameters
#' of a latent growth curve model with time-invariant covariates (if any)
#'
#' @param model An object representing a fitted latent growth curve model.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing to the
#' variability of growth factors if any. It takes the value passed from \code{getLGCM()}.
#' @param names A character vector specifying parameter names. It takes the value passed from \code{getLGCM()}.
#'
#' @return A data frame containing the point estimates and standard errors for parameters of a latent growth curve
#' model with time-invariant covariates (if any).
#'
#' @keywords internal
#'
#' @importFrom OpenMx mxEval mxSE
#'
getLGCM.output <- function(model, growth_TIC, names){
  if (!is.null(growth_TIC)){
    model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                         model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                         mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                         model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model)), 4)
    model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                        model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                        mxSE(beta, model), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                        model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model)), 4)
  }
  else if (is.null(growth_TIC)){
    model.est <- round(c(mxEval(Y_mean0, model), mxEval(Y_psi0, model)[row(mxEval(Y_psi0, model)) >= col(mxEval(Y_psi0, model))],
                         model@output$estimate[grep("Y_residuals", names(model@output$estimate))]), 4)
    model.se <- round(c(mxSE(Y_mean0, model), mxSE(Y_psi0, model)[row(mxSE(Y_psi0, model)) >= col(mxSE(Y_psi0, model))],
                        model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))]), 4)
  }
  estimate_out <- data.frame(Name = names, Estimate = model.est, SE = model.se)
  return(estimate_out)
}

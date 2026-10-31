#' @title Extract Point Estimates And Standard Errors of Latent Change Score Model with Time-invariant Covariates (If Any)
#'
#' @description This function computes and returns a data frame containing point estimates and standard errors for the parameters
#' of a latent change score model with time-invariant covariates (if any).
#'
#' @param model An object representing a fitted latent change score model.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for latent change score
#' models include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}), \code{"Jenss-Bayley"}
#' (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the value passed from \code{getLCSM()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing to the
#' variability of growth factors if any. It takes the value passed from \code{getLCSM()}.
#' @param names A character vector specifying parameter names for the core model parameters (growth factor means/intercepts,
#' variance-covariance elements, residuals, and TIC-related parameters if any). If the vector does not include names for the
#' derived change score quantities (interval-specific slopes, interval changes, change from baseline), these are auto-generated
#' with suffixes \code{slp_val_est}, \code{slp_var_est}, \code{chg_inv_val_est}, \code{chg_inv_var_est}, \code{chg_bl_val_est},
#' and \code{chg_bl_var_est}. It takes the value passed from \code{getLCSM()}.
#'
#' @return A data frame containing the point estimates and standard errors for parameters of a latent change score
#' model with time-invariant covariates (if any).
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom OpenMx mxEval mxSE diag2vec
#'
getLCSM.output <- function(model, curveFun, growth_TIC, names){
  if (curveFun %in% c("nonparametric", "NonP")){
    if (!is.null(growth_TIC)){
      model.est <- c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                           model@output$estimate[grep("Y_rel_rate", names(model@output$estimate))],
                           model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                           mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                           model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model),
                           mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                           mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                           mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model)))
      model.se <- c(mxSE(Y_alpha0, model, silent = TRUE), mxSE(Y_psi_r, model, silent = TRUE)[row(mxSE(Y_psi_r, model, silent = TRUE)) >= col(mxSE(Y_psi_r, model, silent = TRUE))],
                          model@output$standardErrors[, 1][grep("Y_rel_rate", names(model@output$standardErrors[, 1]))],
                          model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                          mxSE(beta, model, silent = TRUE), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                          model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model, silent = TRUE),
                          mxSE(Yslp_m, model, silent = TRUE), diag(mxSE(Yslp_v, model, silent = TRUE)),
                          mxSE(Ychg_inv_m, model, silent = TRUE), mxSE(Ychg_inv_v, model, silent = TRUE),
                          mxSE(Ychg_bl_m, model, silent = TRUE), diag(mxSE(Ychg_bl_v, model, silent = TRUE)))
    }
    else if (is.null(growth_TIC)){
      model.est <- c(mxEval(Y_mean0, model), mxEval(Y_psi0, model)[row(mxEval(Y_psi0, model)) >= col(mxEval(Y_psi0, model))],
                           model@output$estimate[grep("Y_rel_rate", names(model@output$estimate))],
                           model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                           mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                           mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                           mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model)))
      model.se <- c(mxSE(Y_mean0, model, silent = TRUE), mxSE(Y_psi0, model, silent = TRUE)[row(mxSE(Y_psi0, model, silent = TRUE)) >= col(mxSE(Y_psi0, model, silent = TRUE))],
                          model@output$standardErrors[, 1][grep("Y_rel_rate", names(model@output$standardErrors[, 1]))],
                          model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                          mxSE(Yslp_m, model, silent = TRUE), diag(mxSE(Yslp_v, model, silent = TRUE)),
                          mxSE(Ychg_inv_m, model, silent = TRUE), mxSE(Ychg_inv_v, model, silent = TRUE),
                          mxSE(Ychg_bl_m, model, silent = TRUE), diag(mxSE(Ychg_bl_v, model, silent = TRUE)))
    }
  }
  else{
    if (!is.null(growth_TIC)){
      model.est <- c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                           model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                           mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                           model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model),
                           mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                           mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                           mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model)))
      model.se <- c(mxSE(Y_alpha0, model, silent = TRUE), mxSE(Y_psi_r, model, silent = TRUE)[row(mxSE(Y_psi_r, model, silent = TRUE)) >= col(mxSE(Y_psi_r, model, silent = TRUE))],
                          model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                          mxSE(beta, model, silent = TRUE), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                          model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model, silent = TRUE),
                          mxSE(Yslp_m, model, silent = TRUE), diag(mxSE(Yslp_v, model, silent = TRUE)),
                          mxSE(Ychg_inv_m, model, silent = TRUE), mxSE(Ychg_inv_v, model, silent = TRUE),
                          mxSE(Ychg_bl_m, model, silent = TRUE), diag(mxSE(Ychg_bl_v, model, silent = TRUE)))
    }
    else if (is.null(growth_TIC)){
      model.est <- c(mxEval(Y_mean0, model), mxEval(Y_psi0, model)[row(mxEval(Y_psi0, model)) >= col(mxEval(Y_psi0, model))],
                           model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                           mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                           mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                           mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model)))
      model.se <- c(mxSE(Y_mean0, model, silent = TRUE), mxSE(Y_psi0, model, silent = TRUE)[row(mxSE(Y_psi0, model, silent = TRUE)) >= col(mxSE(Y_psi0, model, silent = TRUE))],
                          model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                          mxSE(Yslp_m, model, silent = TRUE), diag(mxSE(Yslp_v, model, silent = TRUE)),
                          mxSE(Ychg_inv_m, model, silent = TRUE), mxSE(Ychg_inv_v, model, silent = TRUE),
                          mxSE(Ychg_bl_m, model, silent = TRUE), diag(mxSE(Ychg_bl_v, model, silent = TRUE)))
    }
  }
  # Accept either the complete labels or the exact core-model label count.
  nIntervals <- length(mxEval(Yslp_m, model))
  derived_names <- unlist(lapply(c("slp_val_est", "slp_var_est", "chg_inv_val_est",
                                   "chg_inv_var_est", "chg_bl_val_est", "chg_bl_var_est"),
                                 function(prefix) paste0(prefix, seq_len(nIntervals))))
  names <- .validate_output_names(names, length(model.est), derived_names)
  estimate_out <- data.frame(Name = names, Estimate = model.est, SE = model.se)
  return(estimate_out)
}

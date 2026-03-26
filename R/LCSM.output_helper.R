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
      model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                           model@output$estimate[grep("Y_rel_rate", names(model@output$estimate))],
                           model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                           mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                           model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model),
                           mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                           mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                           mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
      model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                          model@output$standardErrors[, 1][grep("Y_rel_rate", names(model@output$standardErrors[, 1]))],
                          model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                          mxSE(beta, model), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                          model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model),
                          mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                          mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                          mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
    }
    else if (is.null(growth_TIC)){
      model.est <- round(c(mxEval(Y_mean0, model), mxEval(Y_psi0, model)[row(mxEval(Y_psi0, model)) >= col(mxEval(Y_psi0, model))],
                           model@output$estimate[grep("Y_rel_rate", names(model@output$estimate))],
                           model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                           mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                           mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                           mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
      model.se <- round(c(mxSE(Y_mean0, model), mxSE(Y_psi0, model)[row(mxSE(Y_psi0, model)) >= col(mxSE(Y_psi0, model))],
                          model@output$standardErrors[, 1][grep("Y_rel_rate", names(model@output$standardErrors[, 1]))],
                          model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                          mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                          mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                          mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
    }
  }
  else{
    if (!is.null(growth_TIC)){
      model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                           model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                           mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                           model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model),
                           mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                           mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                           mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
      model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                          model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                          mxSE(beta, model), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                          model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model),
                          mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                          mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                          mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
    }
    else if (is.null(growth_TIC)){
      model.est <- round(c(mxEval(Y_mean0, model), mxEval(Y_psi0, model)[row(mxEval(Y_psi0, model)) >= col(mxEval(Y_psi0, model))],
                           model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                           mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                           mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                           mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
      model.se <- round(c(mxSE(Y_mean0, model), mxSE(Y_psi0, model)[row(mxSE(Y_psi0, model)) >= col(mxSE(Y_psi0, model))],
                          model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                          mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                          mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                          mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
    }
  }
  # Auto-generate derived quantity names if the user-supplied names vector is shorter than model.est
  nDerived <- length(model.est) - length(names)
  if (nDerived > 0){
    nIntervals <- nDerived %/% 6
    derived_names <- c(paste0("slp_val_est", seq_len(nIntervals)),
                       paste0("slp_var_est", seq_len(nIntervals)),
                       paste0("chg_inv_val_est", seq_len(nIntervals)),
                       paste0("chg_inv_var_est", seq_len(nIntervals)),
                       paste0("chg_bl_val_est", seq_len(nIntervals)),
                       paste0("chg_bl_var_est", seq_len(nIntervals)))
    names <- c(names, derived_names)
  }
  # Ensure names and estimates have matching length
  if (length(names) != length(model.est)){
    warning("Number of parameter names (", length(names), ") does not match number of estimates (",
            length(model.est), "). Please check the 'names' argument. ",
            "Auto-padding with generic labels.")
    if (length(names) < length(model.est)){
      names <- c(names, paste0("param_", seq(length(names) + 1, length(model.est))))
    } else {
      names <- names[seq_len(length(model.est))]
    }
  }
  estimate_out <- data.frame(Name = names, Estimate = model.est, SE = model.se)
  return(estimate_out)
}

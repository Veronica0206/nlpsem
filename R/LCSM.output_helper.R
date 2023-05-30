#' Extract Point Estimates And Standard Errors of Latent Change Score Model with Time-invariant Covariates (If Any)
#'
#' This function computes and returns a data frame containing point estimates and standard errors for the parameters
#' of a latent change score model with time-invariant covariates (if any)
#'
#' @param model An object representing a fitted latent change score model.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for latent change score
#' models include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}), \code{"Jenss-Bayley"}
#' (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the value passed from \code{getLCSM()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing to the
#' variability of growth factors if any. It takes the value passed from \code{getLCSM()}.
#' @param names A character vector specifying parameter names. It takes the value passed from \code{getLCSM()}.
#'
#' @return A data frame containing the point estimates and standard errors for parameters of a latent change  score
#' model with time-invariant covariates (if any).
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
  estimate_out <- data.frame(Name = names, Estimate = model.est, SE = model.se)
  return(estimate_out)
}

#' @title Extract Point Estimates And Standard Errors of Latent Growth Curve Model Or Latent Change Score Model
#' with a Time-varying Covariate and Time-invariant Covariates (If Any)
#'
#' @description This function computes and returns a data frame containing point estimates and standard errors for the
#' parameters of a latent growth curve model or a latent change score model with a time-varying covariate
#' and time-invariant covariates (if any)
#'
#' @param model An object representing a fitted latent growth curve model or a latent change score model with a TVC.
#' It takes the value passed from \code{getTVCmodel()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the
#' value passed from \code{getTVCmodel()}.
#' @param records A numeric vector specifying the indices of the observed study waves. It takes the value passed from
#' \code{getTVCmodel()}.
#' @param y_model A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' It takes the value passed from \code{getTVCmodel()}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). It takes the
#' value passed from \code{getTVCmodel()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) that account for the
#' variability of growth factors, if any. It takes the value passed from \code{getTVCmodel()}.
#' @param names A character vector specifying parameter names. It takes the value passed from \code{getTVCmodel()}.
#'
#' @return A data frame containing the point estimates and standard errors for parameters of a latent growth
#' curve model or a latent change score model with a time-varying covariate and time-invariant covariates (if
#' any).
#'
#' @importFrom OpenMx mxEval mxSE diag2vec
#'
getTVC.output <- function(model, curveFun, records, y_model, decompose, growth_TIC, names){
  if (decompose == 0){
    if (y_model == "LGCM"){
      if (!is.null(growth_TIC)){
        model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                             model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                             model@output$estimate[grep("TVC_m", names(model@output$estimate))],
                             model@output$estimate[grep("TVC_v", names(model@output$estimate))],
                             mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                             model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model),
                             model@output$estimate[grep("kappa", names(model@output$estimate))]), 4)
        model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                            model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("TVC_m", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("TVC_v", names(model@output$standardErrors[, 1]))],
                            mxSE(beta, model), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model),
                            model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))]), 4)
      }
      else if (is.null(growth_TIC)){
        model.est <- round(c(mxEval(Y_mean0, model), mxEval(Y_psi0, model)[row(mxEval(Y_psi0, model)) >= col(mxEval(Y_psi0, model))],
                             model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                             model@output$estimate[grep("TVC_m", names(model@output$estimate))],
                             model@output$estimate[grep("TVC_v", names(model@output$estimate))],
                             model@output$estimate[grep("kappa", names(model@output$estimate))]), 4)
        model.se <- round(c(mxSE(Y_mean0, model), mxSE(Y_psi0, model)[row(mxSE(Y_psi0, model)) >= col(mxSE(Y_psi0, model))],
                            model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("TVC_m", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("TVC_v", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))]), 4)
      }
    }
    else if (y_model == "LCSM"){
      if (curveFun %in% c("nonparametric", "NonP")){
        if (!is.null(growth_TIC)){
          model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                               model@output$estimate[grep("Y_rel_rate", names(model@output$estimate))],
                               model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                               model@output$estimate[grep("TVC_m", names(model@output$estimate))],
                               model@output$estimate[grep("TVC_v", names(model@output$estimate))],
                               mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                               model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model),
                               model@output$estimate[grep("kappa", names(model@output$estimate))],
                               mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                               mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                               mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
          model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                              model@output$standardErrors[, 1][grep("Y_rel_rate", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("TVC_m", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("TVC_v", names(model@output$standardErrors[, 1]))],
                              mxSE(beta, model), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model),
                              model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))],
                              mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                              mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                              mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
        }
        else if (is.null(growth_TIC)){
          model.est <- round(c(mxEval(Y_mean0, model), mxEval(Y_psi0, model)[row(mxEval(Y_psi0, model)) >= col(mxEval(Y_psi0, model))],
                               model@output$estimate[grep("Y_rel_rate", names(model@output$estimate))],
                               model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                               model@output$estimate[grep("TVC_m", names(model@output$estimate))],
                               model@output$estimate[grep("TVC_v", names(model@output$estimate))],
                               model@output$estimate[grep("kappa", names(model@output$estimate))],
                               mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                               mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                               mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
          model.se <- round(c(mxSE(Y_mean0, model), mxSE(Y_psi0, model)[row(mxSE(Y_psi0, model)) >= col(mxSE(Y_psi0, model))],
                              model@output$standardErrors[, 1][grep("Y_rel_rate", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("TVC_m", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("TVC_v", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))],
                              mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                              mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                              mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
        }
      }
      else{
        if (!is.null(growth_TIC)){
          model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                               model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                               model@output$estimate[grep("TVC_m", names(model@output$estimate))],
                               model@output$estimate[grep("TVC_v", names(model@output$estimate))],
                               mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                               model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model),
                               model@output$estimate[grep("kappa", names(model@output$estimate))],
                               mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                               mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                               mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
          model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                              model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("TVC_m", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("TVC_v", names(model@output$standardErrors[, 1]))],
                              mxSE(beta, model), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model),
                              model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))],
                              mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                              mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                              mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
        }
        else if (is.null(growth_TIC)){
          model.est <- round(c(mxEval(Y_mean0, model), mxEval(Y_psi0, model)[row(mxEval(Y_psi0, model)) >= col(mxEval(Y_psi0, model))],
                               model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                               model@output$estimate[grep("TVC_m", names(model@output$estimate))],
                               model@output$estimate[grep("TVC_v", names(model@output$estimate))],
                               model@output$estimate[grep("kappa", names(model@output$estimate))],
                               mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                               mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                               mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
          model.se <- round(c(mxSE(Y_mean0, model), mxSE(Y_psi0, model)[row(mxSE(Y_psi0, model)) >= col(mxSE(Y_psi0, model))],
                              model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("TVC_m", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("TVC_v", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))],
                              mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                              mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                              mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
        }
      }
    }
  }
  else if (decompose != 0){
    X_slp_m.est <- X_slp_m.se <- rep(0, length(records) - 1)
    for (j in records[-1]){
      X_slp_m.est[j - 1] <- mxEvalByName(paste0("X_abs_rate", j), model)
      X_slp_m.se[j - 1] <- mxSE(paste0("X_abs_rate", j), model, forceName = T)
    }
    if (y_model == "LGCM"){
      if (!is.null(growth_TIC)){
        model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                             model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                             mxEval(X_mean0, model), mxEval(X_psi0, model)[row(mxEval(X_psi0, model)) >= col(mxEval(X_psi0, model))],
                             model@output$estimate[grep("X_rel_rate", names(model@output$estimate))], X_slp_m.est,
                             model@output$estimate[grep("X_residuals", names(model@output$estimate))],
                             mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                             model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model),
                             model@output$estimate[grep("covBL", names(model@output$estimate))],
                             model@output$estimate[grep("kappa", names(model@output$estimate))],
                             model@output$estimate[grep("Cov_XYres", names(model@output$estimate))]), 4)
        model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                            model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                            mxSE(X_mean0, model), mxSE(X_psi0, model)[row(mxSE(X_psi0, model)) >= col(mxSE(X_psi0, model))],
                            model@output$standardErrors[, 1][grep("X_rel_rate", names(model@output$standardErrors[, 1]))], X_slp_m.se,
                            model@output$standardErrors[, 1][grep("X_residuals", names(model@output$standardErrors[, 1]))],
                            mxSE(beta, model), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model),
                            model@output$standardErrors[, 1][grep("covBL", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("Cov_XYres", names(model@output$standardErrors[, 1]))]), 4)
      }
      else if (is.null(growth_TIC)){
        model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                             model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                             mxEval(X_mean0, model), mxEval(X_psi0, model)[row(mxEval(X_psi0, model)) >= col(mxEval(X_psi0, model))],
                             model@output$estimate[grep("X_rel_rate", names(model@output$estimate))], X_slp_m.est,
                             model@output$estimate[grep("X_residuals", names(model@output$estimate))],
                             mxEval(beta, model), mxEval(Y_mean0, model),
                             model@output$estimate[grep("covBL", names(model@output$estimate))],
                             model@output$estimate[grep("kappa", names(model@output$estimate))],
                             model@output$estimate[grep("Cov_XYres", names(model@output$estimate))]), 4)
        model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                            model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                            mxSE(X_mean0, model), mxSE(X_psi0, model)[row(mxSE(X_psi0, model)) >= col(mxSE(X_psi0, model))],
                            model@output$standardErrors[, 1][grep("X_rel_rate", names(model@output$standardErrors[, 1]))], X_slp_m.se,
                            model@output$standardErrors[, 1][grep("X_residuals", names(model@output$standardErrors[, 1]))],
                            mxSE(beta, model), mxSE(Y_mean0, model),
                            model@output$standardErrors[, 1][grep("covBL", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))],
                            model@output$standardErrors[, 1][grep("Cov_XYres", names(model@output$standardErrors[, 1]))]), 4)
      }
    }
    else if (y_model == "LCSM"){
      if (curveFun %in% c("nonparametric", "NonP")){
        if (!is.null(growth_TIC)){
          model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                               model@output$estimate[grep("Y_rel_rate", names(model@output$estimate))],
                               model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                               mxEval(X_mean0, model), mxEval(X_psi0, model)[row(mxEval(X_psi0, model)) >= col(mxEval(X_psi0, model))],
                               model@output$estimate[grep("X_rel_rate", names(model@output$estimate))], X_slp_m.est,
                               model@output$estimate[grep("X_residuals", names(model@output$estimate))],
                               mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                               model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model),
                               model@output$estimate[grep("covBL", names(model@output$estimate))],
                               model@output$estimate[grep("kappa", names(model@output$estimate))],
                               model@output$estimate[grep("Cov_XYres", names(model@output$estimate))],
                               mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                               mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                               mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
          model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                              model@output$standardErrors[, 1][grep("Y_rel_rate", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                              mxSE(X_mean0, model), mxSE(X_psi0, model)[row(mxSE(X_psi0, model)) >= col(mxSE(X_psi0, model))],
                              model@output$standardErrors[, 1][grep("X_rel_rate", names(model@output$standardErrors[, 1]))], X_slp_m.se,
                              model@output$standardErrors[, 1][grep("X_residuals", names(model@output$standardErrors[, 1]))],
                              mxSE(beta, model), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model),
                              model@output$standardErrors[, 1][grep("covBL", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("Cov_XYres", names(model@output$standardErrors[, 1]))],
                              mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                              mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                              mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
        }
        else if (is.null(growth_TIC)){
          model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                               model@output$estimate[grep("Y_rel_rate", names(model@output$estimate))],
                               model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                               mxEval(X_mean0, model), mxEval(X_psi0, model)[row(mxEval(X_psi0, model)) >= col(mxEval(X_psi0, model))],
                               model@output$estimate[grep("X_rel_rate", names(model@output$estimate))], X_slp_m.est,
                               model@output$estimate[grep("X_residuals", names(model@output$estimate))],
                               mxEval(beta, model), mxEval(Y_mean0, model),
                               model@output$estimate[grep("covBL", names(model@output$estimate))],
                               model@output$estimate[grep("kappa", names(model@output$estimate))],
                               model@output$estimate[grep("Cov_XYres", names(model@output$estimate))],
                               mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                               mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                               mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
          model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                              model@output$standardErrors[, 1][grep("Y_rel_rate", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                              mxSE(X_mean0, model), mxSE(X_psi0, model)[row(mxSE(X_psi0, model)) >= col(mxSE(X_psi0, model))],
                              model@output$standardErrors[, 1][grep("X_rel_rate", names(model@output$standardErrors[, 1]))], X_slp_m.se,
                              model@output$standardErrors[, 1][grep("X_residuals", names(model@output$standardErrors[, 1]))],
                              mxSE(beta, model), mxSE(Y_mean0, model),
                              model@output$standardErrors[, 1][grep("covBL", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("Cov_XYres", names(model@output$standardErrors[, 1]))],
                              mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                              mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                              mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
        }
      }
      else{
        if (!is.null(growth_TIC)){
          model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                               model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                               mxEval(X_mean0, model), mxEval(X_psi0, model)[row(mxEval(X_psi0, model)) >= col(mxEval(X_psi0, model))],
                               model@output$estimate[grep("X_rel_rate", names(model@output$estimate))], X_slp_m.est,
                               model@output$estimate[grep("X_residuals", names(model@output$estimate))],
                               mxEval(beta, model), model@output$estimate[grep("mux", names(model@output$estimate))],
                               model@output$estimate[grep("phi", names(model@output$estimate))], mxEval(Y_mean0, model),
                               model@output$estimate[grep("covBL", names(model@output$estimate))],
                               model@output$estimate[grep("kappa", names(model@output$estimate))],
                               model@output$estimate[grep("Cov_XYres", names(model@output$estimate))],
                               mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                               mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                               mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
          model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                              model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                              mxSE(X_mean0, model), mxSE(X_psi0, model)[row(mxSE(X_psi0, model)) >= col(mxSE(X_psi0, model))],
                              model@output$standardErrors[, 1][grep("X_rel_rate", names(model@output$standardErrors[, 1]))], X_slp_m.se,
                              model@output$standardErrors[, 1][grep("X_residuals", names(model@output$standardErrors[, 1]))],
                              mxSE(beta, model), model@output$standardErrors[, 1][grep("mux", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("phi", names(model@output$standardErrors[, 1]))], mxSE(Y_mean0, model),
                              model@output$standardErrors[, 1][grep("covBL", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("Cov_XYres", names(model@output$standardErrors[, 1]))],
                              mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                              mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                              mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
        }
        else if (is.null(growth_TIC)){
          model.est <- round(c(mxEval(Y_alpha0, model), mxEval(Y_psi_r, model)[row(mxEval(Y_psi_r, model)) >= col(mxEval(Y_psi_r, model))],
                               model@output$estimate[grep("Y_residuals", names(model@output$estimate))],
                               mxEval(X_mean0, model), mxEval(X_psi0, model)[row(mxEval(X_psi0, model)) >= col(mxEval(X_psi0, model))],
                               model@output$estimate[grep("X_rel_rate", names(model@output$estimate))], X_slp_m.est,
                               model@output$estimate[grep("X_residuals", names(model@output$estimate))],
                               mxEval(beta, model), mxEval(Y_mean0, model),
                               model@output$estimate[grep("covBL", names(model@output$estimate))],
                               model@output$estimate[grep("kappa", names(model@output$estimate))],
                               model@output$estimate[grep("Cov_XYres", names(model@output$estimate))],
                               mxEval(Yslp_m, model), diag(mxEval(Yslp_v, model)),
                               mxEval(Ychg_inv_m, model), mxEval(Ychg_inv_v, model),
                               mxEval(Ychg_bl_m, model), diag(mxEval(Ychg_bl_v, model))), 4)
          model.se <- round(c(mxSE(Y_alpha0, model), mxSE(Y_psi_r, model)[row(mxSE(Y_psi_r, model)) >= col(mxSE(Y_psi_r, model))],
                              model@output$standardErrors[, 1][grep("Y_residuals", names(model@output$standardErrors[, 1]))],
                              mxSE(X_mean0, model), mxSE(X_psi0, model)[row(mxSE(X_psi0, model)) >= col(mxSE(X_psi0, model))],
                              model@output$standardErrors[, 1][grep("X_rel_rate", names(model@output$standardErrors[, 1]))], X_slp_m.se,
                              model@output$standardErrors[, 1][grep("X_residuals", names(model@output$standardErrors[, 1]))],
                              mxSE(beta, model), mxSE(Y_mean0, model),
                              model@output$standardErrors[, 1][grep("covBL", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("kappa", names(model@output$standardErrors[, 1]))],
                              model@output$standardErrors[, 1][grep("Cov_XYres", names(model@output$standardErrors[, 1]))],
                              mxSE(Yslp_m, model), diag(mxSE(Yslp_v, model)),
                              mxSE(Ychg_inv_m, model), mxSE(Ychg_inv_v, model),
                              mxSE(Ychg_bl_m, model), diag(mxSE(Ychg_bl_v, model))), 4)
        }
      }
    }
  }
  estimate_out <- data.frame(Name = names, Estimate = model.est, SE = model.se)
  return(estimate_out)
}



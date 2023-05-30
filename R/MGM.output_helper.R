#' @title Extract Point Estimates And Standard Errors of Multivariate Latent Growth Curve Models Or Multivariate Latent Change
#' Score Models
#'
#' @description This function computes and returns a data frame containing point estimates and standard errors for the parameters of
#' a multivariate latent growth curve model or a multivariate latent change score model.
#'
#' @param model An object representing a fitted multivariate latent growth curve model or latent change score model.
#' @param y_var A vector of strings, with each element representing the prefix for column names corresponding to a
#' particular outcome variable at each study wave. It takes the value passed from \code{getMGM()}.
#' @param records A list of numeric vectors, with each vector specifying the indices of the observed study waves for
#' the corresponding outcome variable. It takes the value passed from \code{getMGM()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the
#' value passed from \code{getMGM()}.
#' @param y_model A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' It takes the value passed from \code{getMGM()}.
#' @param names A character vector specifying parameter names. It takes the value passed from \code{getMGM()}.
#'
#' @return A data frame containing the point estimates and standard errors for parameters of a multivariate
#' latent growth curve model or a multivariate latent change score model.
#'
#' @importFrom OpenMx mxEval mxEvalByName mxSE diag2vec
#'
getMGM.output <- function(model, y_var, records, curveFun, y_model, names){
  mean_est <- mean_se <- psi_est <- psi_se <- res_est <- res_se <-
    psi_btw_est <- psi_btw_se <- res_btw_est <- res_btw_se <- outcome_est <- outcome_se <- btw_est <- btw_se <- list()
  for (traj in 1:length(y_var)){
    mean_est[[length(mean_est) + 1]] <- mxEvalByName(paste0(y_var[traj], "_mean0"), model = model)
    mean_se[[length(mean_se) + 1]] <- mxSE(paste0(y_var[traj], "_mean0"), model, forceName = T)
    psi_est[[length(psi_est) + 1]] <- mxEvalByName(paste0(y_var[traj], "_psi0"), model = model)
    psi_se[[length(psi_se) + 1]] <- mxSE(paste0(y_var[traj], "_psi0"), model, forceName = T)
    res_est[[length(res_est) + 1]] <- model@output$estimate[grep(paste0(y_var[traj], "_residuals"), names(model@output$estimate))]
    res_se[[length(res_se) + 1]] <- model@output$standardErrors[, 1][grep(paste0(y_var[traj], "_residuals"),
                                                                          names(model@output$standardErrors[, 1]))]
    outcome_est[[length(outcome_est) + 1]] <- c(unlist(mean_est[[traj]]), psi_est[[traj]][row(psi_est[[traj]]) >= col(psi_est[[traj]])],
                                                res_est[[traj]])
    outcome_se[[length(outcome_se) + 1]] <- c(unlist(mean_se[[traj]]), psi_se[[traj]][row(psi_se[[traj]]) >= col(psi_se[[traj]])], res_se[[traj]])

  }
  for (traj_i in 1:(length(y_var) - 1)){
    for (traj_j in traj_i:(length(y_var) - 1)){
      psi_btw_est[[traj_i + traj_j - 1]] <- mxEvalByName(paste0(y_var[traj_i], y_var[traj_j + 1],
                                                                "_psi"), model = model)
      psi_btw_se[[traj_i + traj_j - 1]] <- mxSE(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"),
                                                model, forceName = T)
      res_btw_est[[length(res_btw_est) + 1]] <- model@output$estimate[grep(paste0(y_var[traj_i], y_var[traj_j + 1], "_RES"), names(model@output$estimate))]
      res_btw_se[[length(res_btw_se) + 1]] <- model@output$standardErrors[, 1][grep(paste0(y_var[traj_i], y_var[traj_j + 1], "_RES"),
                                                                                    names(model@output$standardErrors[, 1]))]
      btw_est[[length(btw_est) + 1]] <- c(unlist(c(psi_btw_est)), unlist(res_btw_est))
      btw_se[[length(btw_se) + 1]] <- c(unlist(c(psi_btw_se)), unlist(res_btw_se))
    }
  }
  if (y_model == "LGCM"){
    model.est <- round(c(unlist(outcome_est), unlist(btw_est)), 4)
    model.se <- round(c(unlist(outcome_se), unlist(btw_se)), 4)
  }
  if (y_model == "LCSM"){
    slp_m.est_L <- slp_v.est_L <- slp_cov.est_L <- chg_inv_m.est_L <- chg_inv_v.est_L <- chg_bl_m.est_L <- chg_bl_v.est_L <-
      slp_m.se_L <- slp_v.se_L <- slp_cov.se_L <- chg_inv_m.se_L <- chg_inv_v.se_L <- chg_bl_m.se_L <- chg_bl_v.se_L <-
      rel_m_est_L <- rel_m_se_L <- list()
    for (traj in 1:length(y_var)){
      slp_m.est_L[[traj]] <- mxEvalByName(paste0(y_var[traj], "slp_m"), model = model)
      slp_v.est_L[[traj]] <- diag(mxEvalByName(paste0(y_var[traj], "slp_v"), model = model))
      chg_inv_m.est_L[[traj]] <- mxEvalByName(paste0(y_var[traj], "chg_inv_m"), model = model)
      chg_inv_v.est_L[[traj]] <- mxEvalByName(paste0(y_var[traj], "chg_inv_v"), model = model)
      chg_bl_m.est_L[[traj]] <- mxEvalByName(paste0(y_var[traj], "chg_bl_m"), model = model)
      chg_bl_v.est_L[[traj]] <- diag(mxEvalByName(paste0(y_var[traj], "chg_bl_v"), model = model))
      slp_m.se_L[[traj]] <- mxSE(paste0(y_var[traj], "slp_m"), model, forceName = T)
      slp_v.se_L[[traj]] <- diag(mxSE(paste0(y_var[traj], "slp_v"), model, forceName = T))
      chg_inv_m.se_L[[traj]] <- mxSE(paste0(y_var[traj], "chg_inv_m"), model, forceName = T)
      chg_inv_v.se_L[[traj]] <- mxSE(paste0(y_var[traj], "chg_inv_v"), model, forceName = T)
      chg_bl_m.se_L[[traj]] <- mxSE(paste0(y_var[traj], "chg_bl_m"), model, forceName = T)
      chg_bl_v.se_L[[traj]] <- diag(mxSE(paste0(y_var[traj], "chg_bl_v"), model, forceName = T))
      if (curveFun %in% c("nonparametric", "NonP")){
        rel_m_est_L[[traj]] <- model@output$estimate[grep(paste0(y_var[traj], "_rel_rate"),
                                                          names(model@output$estimate))]
        rel_m_se_L[[traj]] <- model@output$standardErrors[, 1][grep(paste0(y_var[traj], "_rel_rate"),
                                                                    names(model@output$standardErrors[, 1]))]
      }
    }
    if (curveFun %in% c("nonparametric", "NonP")){
      model.est <- round(c(unlist(outcome_est), unlist(rel_m_est_L), unlist(btw_est), unlist(slp_m.est_L), unlist(slp_v.est_L),
                           unlist(chg_inv_m.est_L), unlist(chg_inv_v.est_L), unlist(chg_bl_m.est_L), unlist(chg_bl_v.est_L)), 4)
      model.se <- round(c(unlist(outcome_se), unlist(rel_m_se_L), unlist(btw_se), unlist(slp_m.se_L), unlist(slp_v.se_L),
                          unlist(chg_inv_m.se_L), unlist(chg_inv_v.se_L), unlist(chg_bl_m.se_L), unlist(chg_bl_v.se_L)), 4)
    }
    else{
      model.est <- round(c(unlist(outcome_est), unlist(btw_est), unlist(slp_m.est_L), unlist(slp_v.est_L), unlist(chg_inv_m.est_L),
                           unlist(chg_inv_v.est_L), unlist(chg_bl_m.est_L), unlist(chg_bl_v.est_L)), 4)
      model.se <- round(c(unlist(outcome_se), unlist(btw_se), unlist(slp_m.se_L), unlist(slp_v.se_L), unlist(chg_inv_m.se_L),
                          unlist(chg_inv_v.se_L), unlist(chg_bl_m.se_L), unlist(chg_bl_v.se_L)), 4)
    }
  }
  estimate_out <- data.frame(Name = names, Estimate = model.est, SE = model.se)
  return(estimate_out)
}






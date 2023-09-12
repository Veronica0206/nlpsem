#' @title Extract Point Estimates And Standard Errors of Longitudinal Mediation Model
#'
#' @description This function computes and returns a data frame containing point estimates and standard errors for the parameters of
#' a longitudinal mediation model.
#'
#' @param model An object representing a fitted longitudinal mediation model.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getMediation()}.
#' @param m_var A string specifying the prefix of the column names corresponding to the mediator variable at each study wave.
#' It takes the value passed from \code{getMediation()}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. It takes the value passed from \code{getMediation()}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. It takes the value passed from
#' \code{getMediation()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options include: \code{"linear"}
#' (or \code{"LIN"}), and \code{"bilinear spline"} (or \code{"BLS"}). It takes the value passed from \code{getMediation()}.
#' @param names A character vector specifying parameter names. It takes the value passed from \code{getMediation()}.
#'
#' @return A data frame containing the point estimates and standard errors for parameters of a longitudinal mediation model.
#'
#' @keywords internal
#'
#' @importFrom OpenMx mxEval mxSE
#'
getMED.output <- function(model, y_var, m_var, x_type, x_var, curveFun, names){
  if (x_type == "baseline"){
    if (curveFun %in% c("linear", "LIN")){
      model.est <- round(c(model@output$estimate[grep("muX", names(model@output$estimate))],
                           model@output$estimate[grep("phi11", names(model@output$estimate))],
                           model$M_alpha$result, c(model$M_psi_r$result)[c(1, 2, 4)],
                           model$Y_alpha$result, c(model$Y_psi_r$result)[c(1, 2, 4)],
                           model$beta_xm$result, model$beta_xy$result,
                           c(model$beta_my$result)[c(1, 2, 4)],
                           model$M_mean$result, model$Y_mean$result,
                           model$mediator$result, model$total$result,
                           model@output$estimate[names(model@output$estimate) == paste0(y_var, "_residuals")],
                           model@output$estimate[names(model@output$estimate) == paste0(m_var, "_residuals")],
                           model@output$estimate[names(model@output$estimate) == paste0(y_var, m_var, "_RES")]), 4)
      model.se <- round(c(model@output$standardErrors[, 1][grep("muX", names(model@output$standardErrors[, 1]))],
                          model@output$standardErrors[, 1][grep("phi11", names(model@output$standardErrors[, 1]))],
                          mxSE(M_alpha, model), c(mxSE(M_psi_r, model))[c(1, 2, 4)],
                          mxSE(Y_alpha, model), c(mxSE(Y_psi_r, model))[c(1, 2, 4)],
                          mxSE(beta_xm, model), mxSE(beta_xy, model),
                          c(mxSE(beta_my, model))[c(1, 2, 4)],
                          mxSE(M_mean, model), mxSE(Y_mean, model),
                          mxSE(mediator, model), mxSE(total, model),
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(y_var, "_residuals")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(m_var, "_residuals")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(y_var, m_var, "_RES")]), 4)
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      model.est <- round(c(model@output$estimate[grep("muX", names(model@output$estimate))],
                           model@output$estimate[grep("phi11", names(model@output$estimate))],
                           model$M_alpha$result, c(model$M_psi_r$result)[c(1:3, 5:6, 9)],
                           model$Y_alpha$result, c(model$Y_psi_r$result)[c(1:3, 5:6, 9)],
                           model$beta_xm$result, model$beta_xy$result,
                           c(model$beta_my$result)[c(1:3, 5:6, 9)],
                           model$M_mean$result, model$Y_mean$result,
                           model$mediator$result, model$total$result,
                           model@output$estimate[names(model@output$estimate) == paste0(y_var, "_residuals")],
                           model@output$estimate[names(model@output$estimate) == paste0(m_var, "_residuals")],
                           model@output$estimate[names(model@output$estimate) == paste0(y_var, m_var, "_RES")]), 4)
      model.se <- round(c(model@output$standardErrors[, 1][grep("muX", names(model@output$standardErrors[, 1]))],
                          model@output$standardErrors[, 1][grep("phi11", names(model@output$standardErrors[, 1]))],
                          mxSE(M_alpha, model), c(mxSE(M_psi_r, model))[c(1:3, 5:6, 9)],
                          mxSE(Y_alpha, model), c(mxSE(Y_psi_r, model))[c(1:3, 5:6, 9)],
                          mxSE(beta_xm, model), mxSE(beta_xy, model),
                          c(mxSE(beta_my, model))[c(1:3, 5:6, 9)],
                          mxSE(M_mean, model), mxSE(Y_mean, model),
                          mxSE(mediator, model), mxSE(total, model),
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(m_var, "_residuals")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(y_var, "_residuals")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(y_var, m_var, "_RES")]), 4)
    }
  }
  if (x_type == "longitudinal"){
    if (curveFun %in% c("linear", "LIN")){
      model.est <- round(c(model$X_mean$result, c(model$X_psi0$result)[c(1, 2, 4)],
                           model$M_alpha$result, c(model$M_psi_r$result)[c(1, 2, 4)],
                           model$Y_alpha$result, c(model$Y_psi_r$result)[c(1, 2, 4)],
                           c(model$beta_xm$result)[c(1, 2, 4)], c(model$beta_xy$result)[c(1, 2, 4)],
                           c(model$beta_my$result)[c(1, 2, 4)], model$M_mean$result, model$Y_mean$result,
                           model$mediator$result, model$total$result[c(1, 2, 4)],
                           model@output$estimate[names(model@output$estimate) == paste0(x_var, "_residuals")],
                           model@output$estimate[names(model@output$estimate) == paste0(y_var, "_residuals")],
                           model@output$estimate[names(model@output$estimate) == paste0(m_var, "_residuals")],
                           model@output$estimate[names(model@output$estimate) == paste0(m_var, x_var, "_RES")],
                           model@output$estimate[names(model@output$estimate) == paste0(y_var, x_var, "_RES")],
                           model@output$estimate[names(model@output$estimate) == paste0(y_var, m_var, "_RES")]), 4)
      model.se <- round(c(mxSE(X_mean, model), c(mxSE(X_psi0, model))[c(1, 2, 4)],
                          mxSE(M_alpha, model), c(mxSE(M_psi_r, model))[c(1, 2, 4)],
                          mxSE(Y_alpha, model), c(mxSE(Y_psi_r, model))[c(1, 2, 4)],
                          c(mxSE(beta_xm, model))[c(1, 2, 4)], c(mxSE(beta_xy, model))[c(1, 2, 4)],
                          c(mxSE(beta_my, model))[c(1, 2, 4)], mxSE(M_mean, model), mxSE(Y_mean, model),
                          mxSE(mediator, model), mxSE(total, model)[c(1, 2, 4)],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(x_var, "_residuals")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(m_var, "_residuals")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(y_var, "_residuals")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(m_var, x_var, "_RES")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(y_var, x_var, "_RES")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(y_var, m_var, "_RES")]), 4)
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      model.est <- round(c(model$X_mean$result, c(model$X_psi0$result)[c(1:3, 5:6, 9)],
                           model$M_alpha$result, c(model$M_psi_r$result)[c(1:3, 5:6, 9)],
                           model$Y_alpha$result, c(model$Y_psi_r$result)[c(1:3, 5:6, 9)],
                           c(model$beta_xm$result)[c(1:3, 5:6, 9)], c(model$beta_xy$result)[c(1:3, 5:6, 9)],
                           c(model$beta_my$result)[c(1:3, 5:6, 9)], model$M_mean$result, model$Y_mean$result,
                           model$mediator$result, model$total$result[c(1:3, 5:6, 9)],
                           model@output$estimate[names(model@output$estimate) == paste0(x_var, "_residuals")],
                           model@output$estimate[names(model@output$estimate) == paste0(y_var, "_residuals")],
                           model@output$estimate[names(model@output$estimate) == paste0(m_var, "_residuals")],
                           model@output$estimate[names(model@output$estimate) == paste0(m_var, x_var, "_RES")],
                           model@output$estimate[names(model@output$estimate) == paste0(y_var, x_var, "_RES")],
                           model@output$estimate[names(model@output$estimate) == paste0(y_var, m_var, "_RES")]), 4)
      model.se <- round(c(mxSE(X_mean, model), c(mxSE(X_psi0, model))[c(1:3, 5:6, 9)],
                          mxSE(M_alpha, model), c(mxSE(M_psi_r, model))[c(1:3, 5:6, 9)],
                          mxSE(Y_alpha, model), c(mxSE(Y_psi_r, model))[c(1:3, 5:6, 9)],
                          c(mxSE(beta_xm, model))[c(1:3, 5:6, 9)], c(mxSE(beta_xy, model))[c(1:3, 5:6, 9)],
                          c(mxSE(beta_my, model))[c(1:3, 5:6, 9)], mxSE(M_mean, model), mxSE(Y_mean, model),
                          mxSE(mediator, model), mxSE(total, model)[c(1:3, 5:6, 9)],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(x_var, "_residuals")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(m_var, "_residuals")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(y_var, "_residuals")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(m_var, x_var, "_RES")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(y_var, x_var, "_RES")],
                          model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0(y_var, m_var, "_RES")]), 4)
    }
  }
  estimate_out <- data.frame(Name = names, Estimate = model.est, SE = model.se)
  return(estimate_out)
}

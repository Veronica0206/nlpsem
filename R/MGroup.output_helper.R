#' @title Extract Point Estimates And Standard Errors of Longitudinal Multiple Group Models
#'
#' @description This function computes and returns a data frame containing point estimates and standard errors for the parameters of
#' a longitudinal multiple group model.
#'
#' @param model An object representing a fitted mixture model.
#' @param nClass An integer specifying the number of latent classes for the mixture model. It takes the value passed from \code{getMGroup()}.
#' @param sub_Model A string that specifies the sub-model for latent classes. Supported sub-models include \code{"LGCM"} (for latent
#' growth curve models), \code{"LCSM"} (for latent change score models), \code{"TVC"} (for latent growth curve models or latent change
#' score models with a time-varying covariate), \code{"MGM"} (for multivariate latent growth curve models or latent change score models),
#' and \code{"MED"} (for longitudinal mediation models). It takes the value passed from \code{getMGroup()}.
#' @param y_var A string defining the prefix of the column names corresponding to the outcome variable for each study wave. This is applicable
#' when \code{sub_Model} is not \code{"MGM"}. For \code{sub_Model} being \code{"MGM"}, \code{y_var} should be a string vector where each element
#' corresponds to the prefix of the column names for each outcome variable across the study waves. It takes the value passed from \code{getMGroup()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model = "LGCM"} include:
#' \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported options for
#' \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the value passed from \code{getMGroup()}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. It takes the value passed from \code{getMGroup()}.
#' @param records A numeric vector denoting the indices of the observed study waves. This applies when \code{sub_Model} is \code{"LGCM"},
#' \code{"LCSM"} or \code{"TVC"}. For \code{sub_Model} being \code{"MGM"} or \code{"MED"}, \code{records} should be a list of numeric vectors,
#' where each vector provides the indices of the observed study waves for each longitudinal process. It takes the value passed from \code{getMGroup()}.
#' @param growth_TIC A string or character vector of column names of time-invariant covariate(s) accounting for the variability
#' of growth factors if any. It takes the value passed from \code{getMGroup()}.
#' @param y_model A string that specifies how to fit longitudinal outcomes. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' It takes the value passed from \code{getMGroup()}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). It takes the value passed
#' from \code{getMGroup()}.
#' @param names A character vector specifying parameter names. It takes the value passed from \code{getMGroup()}.
#'
#' @return A dataframe containing point estimates and standard errors for the parameters of interest for a mixture model.
#'
#' @importFrom OpenMx mxEvalByName mxSE diag2vec
#'
getMGroup.output <- function(model, nClass, sub_Model, y_var, curveFun, x_type, records, growth_TIC,
                             y_model, decompose, names){
  model.est <- model.se <- est <- list()
  for (k in 1:nClass){
    if (sub_Model == "LGCM"){
      if (!is.null(growth_TIC)){
        model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                  mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                    row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                      col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                  model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                  mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                  model@output$estimate[grep(paste0("c", k, "mux"), names(model@output$estimate))],
                                  model@output$estimate[grep(paste0("c", k, "phi"), names(model@output$estimate))],
                                  mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]])), 4)
        model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                 mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                   row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                     col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                 model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                 mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                 model@output$standardErrors[, 1][grep(paste0("c", k, "mux"), names(model@output$standardErrors[, 1]))],
                                 model@output$standardErrors[, 1][grep(paste0("c", k, "phi"), names(model@output$standardErrors[, 1]))],
                                 mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T)), 4)
      }
      else if (is.null(growth_TIC)){
        model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                  mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])[
                                    row(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])) >=
                                      col(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]]))],
                                  model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))]), 4)
        model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                 mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)[
                                   row(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)) >=
                                     col(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T))],
                                 model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))]), 4)
      }
    }
    else if (sub_Model == "LCSM"){
      if (curveFun %in% c("nonparametric", "NonP")){
        if (!is.null(growth_TIC)){
          model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                      row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                        col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                    model@output$estimate[grep(paste0("c", k, "Y_rel_rate"), names(model@output$estimate))],
                                    model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                    mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                    model@output$estimate[grep(paste0("c", k, "mux"), names(model@output$estimate))],
                                    model@output$estimate[grep(paste0("c", k, "phi"), names(model@output$estimate))],
                                    mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                    diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                    mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                    diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
          model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                     row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                       col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "Y_rel_rate"), names(model@output$standardErrors[, 1]))],
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                   mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "mux"), names(model@output$standardErrors[, 1]))],
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "phi"), names(model@output$standardErrors[, 1]))],
                                   mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                   diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                   diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
        }
        else if (is.null(growth_TIC)){
          model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])[
                                      row(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])) >=
                                        col(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]]))],
                                    model@output$estimate[grep(paste0("c", k, "Y_rel_rate"), names(model@output$estimate))],
                                    model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                    mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                    diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                    mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                    diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
          model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)[
                                     row(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)) >=
                                       col(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T))],
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "Y_rel_rate"), names(model@output$standardErrors[, 1]))],
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                   mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                   diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                   diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
        }
      }
      else {
        if (!is.null(growth_TIC)){
          model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                      row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                        col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                    model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                    mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                    model@output$estimate[grep(paste0("c", k, "mux"), names(model@output$estimate))],
                                    model@output$estimate[grep(paste0("c", k, "phi"), names(model@output$estimate))],
                                    mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                    diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                    mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                    diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
          model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                     row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                       col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                   mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "mux"), names(model@output$standardErrors[, 1]))],
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "phi"), names(model@output$standardErrors[, 1]))],
                                   mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                   diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                   diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
        }
        else if (is.null(growth_TIC)){
          model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])[
                                      row(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])) >=
                                        col(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]]))],
                                    model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                    mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                    diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                    mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                    diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
          model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)[
                                     row(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)) >=
                                       col(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T))],
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                   mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                   diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                   diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
        }
      }
    }
    else if (sub_Model == "TVC"){
      if (decompose == 0){
        if (y_model == "LGCM"){
          if (!is.null(growth_TIC)){
            model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                      mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                        row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                          col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                      model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                      model@output$estimate[grep(paste0("c", k, "TVC_m"), names(model@output$estimate))],
                                      model@output$estimate[grep(paste0("c", k, "TVC_v"), names(model@output$estimate))],
                                      mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                      model@output$estimate[grep(paste0("c", k, "mux"), names(model@output$estimate))],
                                      model@output$estimate[grep(paste0("c", k, "phi"), names(model@output$estimate))],
                                      mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                      model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))]), 4)
            model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                     mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                       row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                         col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_m"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_v"), names(model@output$standardErrors[, 1]))],
                                     mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "mux"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "phi"), names(model@output$standardErrors[, 1]))],
                                     mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "covBL"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "Cov_XYres"), names(model@output$standardErrors[, 1]))]), 4)
          }
          else if (is.null(growth_TIC)){
            model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                      mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])[
                                        row(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])) >=
                                          col(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]]))],
                                      model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                      model@output$estimate[grep(paste0("c", k, "TVC_m"), names(model@output$estimate))],
                                      model@output$estimate[grep(paste0("c", k, "TVC_v"), names(model@output$estimate))],
                                      model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))]), 4)
            model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                     mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)[
                                       row(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)) >=
                                         col(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_m"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_v"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))]), 4)
          }
          est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
        }
        else if (y_model == "LCSM"){
          if (curveFun %in% c("nonparametric", "NonP")){
            if (!is.null(growth_TIC)){
              model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "Y_rel_rate"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "TVC_m"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "TVC_v"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                        model@output$estimate[grep(paste0("c", k, "mux"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "phi"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                        model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
              model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_rel_rate"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_m"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_v"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "mux"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "phi"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
            }
            else if (is.null(growth_TIC)){
              model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "Y_rel_rate"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "TVC_m"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "TVC_v"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
              model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_rel_rate"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_m"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_v"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
            }
          }
          else{
            if (!is.null(growth_TIC)){
              model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "TVC_m"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "TVC_v"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                        model@output$estimate[grep(paste0("c", k, "mux"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "phi"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                        model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
              model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_m"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_v"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "mux"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "phi"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
            }
            else if (is.null(growth_TIC)){
              model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "Y_psi0"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "TVC_m"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "TVC_v"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
              model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "Y_psi0"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_m"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "TVC_v"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
            }
          }
          est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
        }
      }
      else if (decompose != 0){
        if (y_model == "LGCM"){
          X_slp_m.est <- X_slp_m.se <- rep(0, length(records) - 1)
          for (j in records[-1]){
            X_slp_m.est[j - 1] <- mxEvalByName(paste0("c", k, "abs_rate", j), model = model@submodels[[k]])
            X_slp_m.se[j - 1] <- mxSE(paste0("Class", k, ".c", k, "abs_rate", j), model = model, forceName = T)
          }
          if (!is.null(growth_TIC)){
            model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                      mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                        row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                          col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                      model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                      mxEvalByName(paste0("c", k, "X_mean0"), model = model@submodels[[k]]),
                                      mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])[
                                        row(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])) >=
                                          col(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]]))],
                                      model@output$estimate[grep(paste0("c", k, "X_rel_rate"), names(model@output$estimate))], X_slp_m.est,
                                      model@output$estimate[grep(paste0("c", k, "X_residuals"), names(model@output$estimate))],
                                      mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                      model@output$estimate[grep(paste0("c", k, "mux"), names(model@output$estimate))],
                                      model@output$estimate[grep(paste0("c", k, "phi"), names(model@output$estimate))],
                                      mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                      model@output$estimate[grep(paste0("c", k, "covBL"), names(model@output$estimate))],
                                      model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))],
                                      model@output$estimate[grep(paste0("c", k, "Cov_XYres"), names(model@output$estimate))]), 4)
            model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                     mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                       row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                         col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                     mxSE(paste0("Class", k, ".c", k, "X_mean0"), model, forceName = T),
                                     mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)[
                                       row(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)) >=
                                         col(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "X_rel_rate"), names(model@output$standardErrors[, 1]))], X_slp_m.se,
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "X_residuals"), names(model@output$standardErrors[, 1]))],
                                     mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "mux"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "phi"), names(model@output$standardErrors[, 1]))],
                                     mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "covBL"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "Cov_XYres"), names(model@output$standardErrors[, 1]))]), 4)
          }
          else if (is.null(growth_TIC)){
            model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                      mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                        row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                          col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                      model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                      mxEvalByName(paste0("c", k, "X_mean0"), model = model@submodels[[k]]),
                                      mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])[
                                        row(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])) >=
                                          col(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]]))],
                                      model@output$estimate[grep(paste0("c", k, "X_rel_rate"), names(model@output$estimate))], X_slp_m.est,
                                      model@output$estimate[grep(paste0("c", k, "X_residuals"), names(model@output$estimate))],
                                      mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                      mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                      model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))],
                                      model@output$estimate[grep(paste0("c", k, "Cov_XYres"), names(model@output$estimate))]), 4)
            model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                     mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                       row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                         col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                     mxSE(paste0("Class", k, ".c", k, "X_mean0"), model, forceName = T),
                                     mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)[
                                       row(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)) >=
                                         col(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "X_rel_rate"), names(model@output$standardErrors[, 1]))], X_slp_m.se,
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "X_residuals"), names(model@output$standardErrors[, 1]))],
                                     mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                     mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                     model@output$standardErrors[, 1][grep(paste0("c", k, "Cov_XYres"), names(model@output$standardErrors[, 1]))]), 4)
          }
          est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
        }
        else if (y_model == "LCSM"){
          X_slp_m.est <- X_slp_m.se <- rep(0, length(records) - 1)
          for (j in records[-1]){
            X_slp_m.est[j - 1] <- mxEvalByName(paste0("c", k, "abs_rate", j), model = model@submodels[[k]])
            X_slp_m.se[j - 1] <- mxSE(paste0("Class", k, ".c", k, "abs_rate", j), model = model, forceName = T)
          }
          if (curveFun %in% c("nonparametric", "NonP")){
            if (!is.null(growth_TIC)){
              model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "Y_rel_rate"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "X_mean0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "X_rel_rate"), names(model@output$estimate))], X_slp_m.est,
                                        model@output$estimate[grep(paste0("c", k, "X_residuals"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                        model@output$estimate[grep(paste0("c", k, "mux"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "phi"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                        model@output$estimate[grep(paste0("c", k, "covBL"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "Cov_XYres"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
              model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_rel_rate"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "X_mean0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "X_rel_rate"), names(model@output$standardErrors[, 1]))], X_slp_m.se,
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "X_residuals"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "mux"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "phi"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "covBL"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Cov_XYres"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
            }
            else if (is.null(growth_TIC)){
              model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "Y_rel_rate"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "X_mean0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "X_rel_rate"), names(model@output$estimate))], X_slp_m.est,
                                        model@output$estimate[grep(paste0("c", k, "X_residuals"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                        model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "Cov_XYres"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
              model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_rel_rate"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "X_mean0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "X_rel_rate"), names(model@output$standardErrors[, 1]))], X_slp_m.se,
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "X_residuals"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Cov_XYres"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
            }
          }
          else{
            if (!is.null(growth_TIC)){
              model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "X_mean0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "X_rel_rate"), names(model@output$estimate))], X_slp_m.est,
                                        model@output$estimate[grep(paste0("c", k, "X_residuals"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                        model@output$estimate[grep(paste0("c", k, "mux"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "phi"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                        model@output$estimate[grep(paste0("c", k, "covBL"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "Cov_XYres"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
              model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "X_mean0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "X_rel_rate"), names(model@output$standardErrors[, 1]))], X_slp_m.se,
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "X_residuals"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "mux"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "phi"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "covBL"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Cov_XYres"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
            }
            else if (is.null(growth_TIC)){
              model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "Y_alpha0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "Y_residuals"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "X_mean0"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])[
                                          row(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]])) >=
                                            col(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]]))],
                                        model@output$estimate[grep(paste0("c", k, "X_rel_rate"), names(model@output$estimate))], X_slp_m.est,
                                        model@output$estimate[grep(paste0("c", k, "X_residuals"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Y_mean0"), model = model@submodels[[k]]),
                                        model@output$estimate[grep(paste0("c", k, "kappa"), names(model@output$estimate))],
                                        model@output$estimate[grep(paste0("c", k, "Cov_XYres"), names(model@output$estimate))],
                                        mxEvalByName(paste0("c", k, "Yslp_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Yslp_v"), model = model@submodels[[k]])),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_m"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_inv_v"), model = model@submodels[[k]]),
                                        mxEvalByName(paste0("c", k, "Ychg_bl_m"), model = model@submodels[[k]]),
                                        diag(mxEvalByName(paste0("c", k, "Ychg_bl_v"), model = model@submodels[[k]]))), 4)
              model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "Y_alpha0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Y_residuals"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "X_mean0"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)[
                                         row(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T)) >=
                                           col(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "X_rel_rate"), names(model@output$standardErrors[, 1]))], X_slp_m.se,
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "X_residuals"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Y_mean0"), model, forceName = T),
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "kappa"), names(model@output$standardErrors[, 1]))],
                                       model@output$standardErrors[, 1][grep(paste0("c", k, "Cov_XYres"), names(model@output$standardErrors[, 1]))],
                                       mxSE(paste0("Class", k, ".c", k, "Yslp_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Yslp_v"), model, forceName = T)),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_m"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_inv_v"), model, forceName = T),
                                       mxSE(paste0("Class", k, ".c", k, "Ychg_bl_m"), model, forceName = T),
                                       diag(mxSE(paste0("Class", k, ".c", k, "Ychg_bl_v"), model, forceName = T))), 4)
            }
          }
          est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
        }
      }
    }
    else if (sub_Model == "MGM"){
      mean_est <- mean_se <- psi_est <- psi_se <- res_est <- res_se <-
        psi_btw_est <- psi_btw_se <- res_btw_est <- res_btw_se <- outcome_est <- outcome_se <- btw_est <- btw_se <- list()
      for (traj in 1:length(y_var)){
        mean_est[[traj]] <- mxEvalByName(paste0("c", k, y_var[traj], "_mean0"), model = model@submodels[[k]])
        mean_se[[traj]] <- mxSE(paste0("Class", k, ".c", k, y_var[traj], "_mean0"), model, forceName = T)
        psi_est[[traj]] <- mxEvalByName(paste0("c", k, y_var[traj], "_psi0"), model = model@submodels[[k]])
        psi_se[[traj]] <- mxSE(paste0("Class", k, ".c", k, y_var[traj], "_psi0"), model, forceName = T)
        res_est[[traj]] <- model@output$estimate[grep(paste0("c", k, y_var[traj], "_residuals"),
                                                      names(model@output$estimate))]
        res_se[[traj]] <- model@output$standardErrors[, 1][grep(paste0("c", k, y_var[traj], "_residuals"),
                                                                names(model@output$standardErrors[, 1]))]
        outcome_est[[traj]] <- c(unlist(mean_est[[traj]]), psi_est[[traj]][row(psi_est[[traj]]) >= col(psi_est[[traj]])],
                                 res_est[[traj]])
        outcome_se[[traj]] <- c(unlist(mean_se[[traj]]), psi_se[[traj]][row(psi_se[[traj]]) >= col(psi_se[[traj]])], res_se[[traj]])

      }
      for (traj_i in 1:(length(y_var) - 1)){
        for (traj_j in traj_i:(length(y_var) - 1)){
          psi_btw_est[[traj_i + traj_j - 1]] <- mxEvalByName(paste0("c", k, y_var[traj_i], y_var[traj_j + 1],
                                                                    "_psi"), model = model@submodels[[k]])
          psi_btw_se[[traj_i + traj_j - 1]] <- mxSE(paste0("Class", k, ".c", k, y_var[traj_i], y_var[traj_j + 1], "_psi"),
                                                    model, forceName = T)
          res_btw_est[[length(res_btw_est) + 1]] <- model@output$estimate[grep(paste0("c", k, y_var[traj_i], y_var[traj_j + 1], "_RES"), names(model@output$estimate))]
          res_btw_se[[length(res_btw_se) + 1]] <- model@output$standardErrors[, 1][grep(paste0("c", k, y_var[traj_i], y_var[traj_j + 1], "_RES"),
                                                                                        names(model@output$standardErrors[, 1]))]
          btw_est[[length(btw_est) + 1]] <- c(unlist(c(psi_btw_est)), unlist(res_btw_est))
          btw_se[[length(btw_se) + 1]] <- c(unlist(c(psi_btw_se)), unlist(res_btw_se))
        }
      }
      if (y_model == "LGCM"){
        model.est[[k]] <- round(c(unlist(outcome_est), unlist(btw_est)), 4)
        model.se[[k]] <- round(c(unlist(outcome_se), unlist(btw_se)), 4)
        est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
      }
      else if (y_model == "LCSM"){
        slp_m.est_L <- slp_v.est_L <- slp_cov.est_L <- chg_inv_m.est_L <- chg_inv_v.est_L <- chg_bl_m.est_L <- chg_bl_v.est_L <-
          slp_m.se_L <- slp_v.se_L <- slp_cov.se_L <- chg_inv_m.se_L <- chg_inv_v.se_L <- chg_bl_m.se_L <- chg_bl_v.se_L <-
          rel_val_est_L <- rel_val_se_L <- list()
        for (traj in 1:length(y_var)){
          slp_m.est_L[[traj]] <- mxEvalByName(paste0("c", k, y_var[traj], "slp_m"), model = model@submodels[[k]])
          slp_v.est_L[[traj]] <- diag(mxEvalByName(paste0("c", k, y_var[traj], "slp_v"), model = model@submodels[[k]]))
          chg_inv_m.est_L[[traj]] <- mxEvalByName(paste0("c", k, y_var[traj], "chg_inv_m"), model = model@submodels[[k]])
          chg_inv_v.est_L[[traj]] <- mxEvalByName(paste0("c", k, y_var[traj], "chg_inv_v"), model = model@submodels[[k]])
          chg_bl_m.est_L[[traj]] <- mxEvalByName(paste0("c", k, y_var[traj], "chg_bl_m"), model = model@submodels[[k]])
          chg_bl_v.est_L[[traj]] <- diag(mxEvalByName(paste0("c", k, y_var[traj], "chg_bl_v"), model = model@submodels[[k]]))
          slp_m.se_L[[traj]] <- mxSE(paste0("Class", k, ".c", k, y_var[traj], "slp_m"), model, forceName = T)
          slp_v.se_L[[traj]] <- diag(mxSE(paste0("Class", k, ".c", k, y_var[traj], "slp_v"), model, forceName = T))
          chg_inv_m.se_L[[traj]] <- mxSE(paste0("Class", k, ".c", k, y_var[traj], "chg_inv_m"), model, forceName = T)
          chg_inv_v.se_L[[traj]] <- mxSE(paste0("Class", k, ".c", k, y_var[traj], "chg_inv_v"), model, forceName = T)
          chg_bl_m.se_L[[traj]] <- mxSE(paste0("Class", k, ".c", k, y_var[traj], "chg_bl_m"), model, forceName = T)
          chg_bl_v.se_L[[traj]] <- diag(mxSE(paste0("Class", k, ".c", k, y_var[traj], "chg_bl_v"), model, forceName = T))
          if (curveFun %in% c("nonparametric", "NonP")){
            rel_val_est_L[[traj]] <- model@output$estimate[grep(paste0("c", k, y_var[traj], "_rel_rate"),
                                                                names(model@output$estimate))]
            rel_val_se_L[[traj]] <- model@output$standardErrors[, 1][grep("c", k, paste0(y_var[traj], "_rel_rate"),
                                                                          names(model@output$standardErrors[, 1]))]
          }
        }
        if (curveFun %in% c("nonparametric", "NonP")){
          model.est[[k]] <- round(c(unlist(outcome_est), unlist(rel_val_est_L), unlist(btw_est), unlist(slp_m.est_L), unlist(slp_v.est_L),
                                    unlist(chg_inv_m.est_L), unlist(chg_inv_v.est_L), unlist(chg_bl_m.est_L), unlist(chg_bl_v.est_L)), 4)
          model.se[[k]] <- round(c(unlist(outcome_se), unlist(rel_val_se_L), unlist(btw_se), unlist(slp_m.se_L), unlist(slp_v.se_L),
                                   unlist(chg_inv_m.se_L), unlist(chg_inv_v.se_L), unlist(chg_bl_m.se_L), unlist(chg_bl_v.se_L)), 4)
        }
        else{
          model.est[[k]] <- round(c(unlist(outcome_est), unlist(btw_est), unlist(slp_m.est_L), unlist(slp_v.est_L), unlist(chg_inv_m.est_L),
                                    unlist(chg_inv_v.est_L), unlist(chg_bl_m.est_L), unlist(chg_bl_v.est_L)), 4)
          model.se[[k]] <- round(c(unlist(outcome_se), unlist(btw_se), unlist(slp_m.se_L), unlist(slp_v.se_L), unlist(chg_inv_m.se_L),
                                   unlist(chg_inv_v.se_L), unlist(chg_bl_m.se_L), unlist(chg_bl_v.se_L)), 4)
          est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
        }
      }
    }
    else if (sub_Model == "MED"){
      if (x_type == "baseline"){
        if (curveFun %in% c("linear", "LIN")){
          model.est[[k]] <- round(c(model@output$estimate[grep(paste0("c", k, "muX"), names(model@output$estimate))],
                                    model@output$estimate[grep(paste0("c", k, "phi11"), names(model@output$estimate))],
                                    mxEvalByName(paste0("c", k, "M_alpha"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "M_psi_r"), model = model@submodels[[k]]))[c(1, 2, 4)],
                                    mxEvalByName(paste0("c", k, "Y_alpha"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))[c(1, 2, 4)],
                                    mxEvalByName(paste0("c", k, "beta_xm"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "beta_xy"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "beta_my"), model = model@submodels[[k]]))[c(1, 2, 4)],
                                    mxEvalByName(paste0("c", k, "M_mean"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Y_mean"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "mediator"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "total"), model = model@submodels[[k]]),
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "Y_residuals")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "M_residuals")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "YM_RES")]), 4)
          model.se[[k]] <- round(c(model@output$standardErrors[, 1][grep(paste0("c", k, "muX"), names(model@output$standardErrors[, 1]))],
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "phi11"), names(model@output$standardErrors[, 1]))],
                                   mxSE(paste0("Class", k, ".c", k, "M_alpha"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "M_psi_r"), model, forceName = T))[c(1, 2, 4)],
                                   mxSE(paste0("Class", k, ".c", k, "Y_alpha"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))[c(1, 2, 4)],
                                   mxSE(paste0("Class", k, ".c", k, "beta_xm"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "beta_xy"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "beta_my"), model, forceName = T))[c(1, 2, 4)],
                                   mxSE(paste0("Class", k, ".c", k, "M_mean"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Y_mean"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "mediator"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "total"), model, forceName = T),
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "Y_residuals")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "M_residuals")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "YM_RES")]), 4)
          est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
        }
        else if (curveFun %in% c("bilinear spline", "BLS")){
          model.est[[k]] <- round(c(model@output$estimate[grep(paste0("c", k, "muX"), names(model@output$estimate))],
                                    model@output$estimate[grep(paste0("c", k, "phi11"), names(model@output$estimate))],
                                    mxEvalByName(paste0("c", k, "M_alpha"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "M_psi_r"), model = model@submodels[[k]]))[c(1:3, 5:6, 9)],
                                    mxEvalByName(paste0("c", k, "Y_alpha"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))[c(1:3, 5:6, 9)],
                                    mxEvalByName(paste0("c", k, "beta_xm"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "beta_xy"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "beta_my"), model = model@submodels[[k]]))[c(1:3, 5:6, 9)],
                                    mxEvalByName(paste0("c", k, "M_mean"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Y_mean"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "mediator"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "total"), model = model@submodels[[k]]),
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "Y_residuals")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "M_residuals")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "YM_RES")]), 4)
          model.se[[k]] <- round(c(model@output$standardErrors[, 1][grep(paste0("c", k, "muX"), names(model@output$standardErrors[, 1]))],
                                   model@output$standardErrors[, 1][grep(paste0("c", k, "phi11"), names(model@output$standardErrors[, 1]))],
                                   mxSE(paste0("Class", k, ".c", k, "M_alpha"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "M_psi_r"), model, forceName = T))[c(1:3, 5:6, 9)],
                                   mxSE(paste0("Class", k, ".c", k, "Y_alpha"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))[c(1:3, 5:6, 9)],
                                   mxSE(paste0("Class", k, ".c", k, "beta_xm"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "beta_xy"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "beta_my"), model, forceName = T))[c(1:3, 5:6, 9)],
                                   mxSE(paste0("Class", k, ".c", k, "M_mean"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Y_mean"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "mediator"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "total"), model, forceName = T),
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "Y_residuals")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "M_residuals")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "YM_RES")]), 4)
          est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
        }
      }
      if (x_type == "longitudinal"){
        if (curveFun %in% c("linear", "LIN")){
          model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "X_mean"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]]))[c(1, 2, 4)],
                                    mxEvalByName(paste0("c", k, "M_alpha"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "M_psi_r"), model = model@submodels[[k]]))[c(1, 2, 4)],
                                    mxEvalByName(paste0("c", k, "Y_alpha"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))[c(1, 2, 4)],
                                    c(mxEvalByName(paste0("c", k, "beta_xm"), model = model@submodels[[k]]))[c(1, 2, 4)],
                                    c(mxEvalByName(paste0("c", k, "beta_xy"), model = model@submodels[[k]]))[c(1, 2, 4)],
                                    c(mxEvalByName(paste0("c", k, "beta_my"), model = model@submodels[[k]]))[c(1, 2, 4)],
                                    mxEvalByName(paste0("c", k, "M_mean"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Y_mean"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "mediator"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "total"), model = model@submodels[[k]]))[c(1, 2, 4)],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "Y_residuals")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "M_residuals")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "X_residuals")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "YM_RES")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "YX_RES")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "MX_RES")]), 4)
          model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "X_mean"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T))[c(1, 2, 4)],
                                   mxSE(paste0("Class", k, ".c", k, "M_alpha"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "M_psi_r"), model, forceName = T))[c(1, 2, 4)],
                                   mxSE(paste0("Class", k, ".c", k, "Y_alpha"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))[c(1, 2, 4)],
                                   c(mxSE(paste0("Class", k, ".c", k, "beta_xm"), model, forceName = T))[c(1, 2, 4)],
                                   c(mxSE(paste0("Class", k, ".c", k, "beta_xy"), model, forceName = T))[c(1, 2, 4)],
                                   c(mxSE(paste0("Class", k, ".c", k, "beta_my"), model, forceName = T))[c(1, 2, 4)],
                                   mxSE(paste0("Class", k, ".c", k, "M_mean"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Y_mean"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "mediator"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "total"), model, forceName = T))[c(1, 2, 4)],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "Y_residuals")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "M_residuals")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "X_residuals")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "YM_RES")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "YX_RES")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "MX_RES")]), 4)
          est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
        }
        else if (curveFun %in% c("bilinear spline", "BLS")){
          model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "X_mean"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "X_psi0"), model = model@submodels[[k]]))[c(1:3, 5:6, 9)],
                                    mxEvalByName(paste0("c", k, "M_alpha"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "M_psi_r"), model = model@submodels[[k]]))[c(1:3, 5:6, 9)],
                                    mxEvalByName(paste0("c", k, "Y_alpha"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "Y_psi_r"), model = model@submodels[[k]]))[c(1:3, 5:6, 9)],
                                    c(mxEvalByName(paste0("c", k, "beta_xm"), model = model@submodels[[k]]))[c(1:3, 5:6, 9)],
                                    c(mxEvalByName(paste0("c", k, "beta_xy"), model = model@submodels[[k]]))[c(1:3, 5:6, 9)],
                                    c(mxEvalByName(paste0("c", k, "beta_my"), model = model@submodels[[k]]))[c(1:3, 5:6, 9)],
                                    mxEvalByName(paste0("c", k, "M_mean"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "Y_mean"), model = model@submodels[[k]]),
                                    mxEvalByName(paste0("c", k, "mediator"), model = model@submodels[[k]]),
                                    c(mxEvalByName(paste0("c", k, "total"), model = model@submodels[[k]]))[c(1:3, 5:6, 9)],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "Y_residuals")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "M_residuals")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "X_residuals")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "YM_RES")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "YX_RES")],
                                    model@output$estimate[names(model@output$estimate) == paste0("c", k, "MX_RES")]), 4)
          model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "X_mean"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "X_psi0"), model, forceName = T))[c(1:3, 5:6, 9)],
                                   mxSE(paste0("Class", k, ".c", k, "M_alpha"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "M_psi_r"), model, forceName = T))[c(1:3, 5:6, 9)],
                                   mxSE(paste0("Class", k, ".c", k, "Y_alpha"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "Y_psi_r"), model, forceName = T))[c(1:3, 5:6, 9)],
                                   c(mxSE(paste0("Class", k, ".c", k, "beta_xm"), model, forceName = T))[c(1:3, 5:6, 9)],
                                   c(mxSE(paste0("Class", k, ".c", k, "beta_xy"), model, forceName = T))[c(1:3, 5:6, 9)],
                                   c(mxSE(paste0("Class", k, ".c", k, "beta_my"), model, forceName = T))[c(1:3, 5:6, 9)],
                                   mxSE(paste0("Class", k, ".c", k, "M_mean"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "Y_mean"), model, forceName = T),
                                   mxSE(paste0("Class", k, ".c", k, "mediator"), model, forceName = T),
                                   c(mxSE(paste0("Class", k, ".c", k, "total"), model, forceName = T))[c(1:3, 5:6, 9)],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "Y_residuals")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "M_residuals")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "X_residuals")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "YM_RES")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "YX_RES")],
                                   model@output$standardErrors[, 1][names(model@output$standardErrors[, 1]) == paste0("c", k, "MX_RES")]), 4)
          est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
        }
      }
    }
    est[[k]] <- data.frame(Name = paste0("c", k, names), Estimate = model.est[[k]], SE = model.se[[k]])
  }
  estimate_out <- rbind(do.call(rbind.data.frame, est))
  return(estimate_out)
}



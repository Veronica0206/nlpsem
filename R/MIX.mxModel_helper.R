#' @title Construct An Object of mxModel for Longitudinal Mixture Models To Be Evaluated
#'
#' @description This function builds up an object of mxModel for a mixture model with user-specified number of latent classes, class-
#' specific models and functional form (including whether intrinsically nonlinear).
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions for each longitudinal process, and time-invariant covariates (TICs) if any. It takes the value
#' passed from \code{getMIX()}.
#' @param nClass An integer specifying the number of latent classes for the mixture model. It takes the value passed from \code{getMIX()}.
#' @param sub_Model A string that specifies the sub-model for latent classes. Supported sub-models include \code{"LGCM"} (for latent
#' growth curve models), \code{"LCSM"} (for latent change score models), \code{"TVC"} (for latent growth curve models or latent change
#' score models with a time-varying covariate), \code{"MGM"} (for multivariate latent growth curve models or latent change score models),
#' and \code{"MED"} (for longitudinal mediation models). It takes the value passed from \code{getMIX()}.
#' @param cluster_TIC A string or character vector representing the column name(s) for time-invariant covariate(s) indicating cluster
#' formations. It takes the value passed from \code{getMIX()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable for each study wave. This applies when
#' \code{sub_Model} is \code{"LGCM"}, \code{"LCSM"} or \code{"TVC"}. For \code{sub_Model} being \code{"MGM"} or \code{"MED"}, \code{t_var} should
#' be a string vector where each element corresponds to the time variable prefix for each respective longitudinal process. It takes the value
#' passed from \code{getMIX()}.
#' @param records A numeric vector denoting the indices of the observed study waves. This applies when \code{sub_Model} is \code{"LGCM"},
#' \code{"LCSM"} or \code{"TVC"}. For \code{sub_Model} being \code{"MGM"} or \code{"MED"}, \code{records} should be a list of numeric vectors,
#' where each vector provides the indices of the observed study waves for each longitudinal process. It takes the value passed from \code{getMIX()}.
#' @param y_var A string defining the prefix of the column names corresponding to the outcome variable for each study wave. This is applicable
#' when \code{sub_Model} is not \code{"MGM"}. For \code{sub_Model} being \code{"MGM"}, \code{y_var} should be a string vector where each element
#' corresponds to the prefix of the column names for each outcome variable across the study waves. It takes the value passed from \code{getMIX()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model = "LGCM"} include:
#' \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported options for
#' \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the value passed from \code{getMIX()}.
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. It takes the value passed
#' from \code{getMIX()}.
#' @param y_model A string that specifies how to fit longitudinal outcomes. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' It takes the value passed from \code{getMIX()}.
#' @param m_var A string that specifies the prefix of the column names corresponding to the mediator variable at each study wave.
#' It takes the value passed from \code{getMIX()}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. It takes the value passed from \code{getMIX()}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. It takes the value passed from \code{getMIX()}.
#' @param TVC A string that specifies the prefix of the column names corresponding to the time-varying covariate at each time
#' point. It takes the value passed from \code{getMIX()}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). It takes the value passed
#' from \code{getMIX()}.
#' @param growth_TIC A string or character vector of column names of time-invariant covariate(s) accounting for the variability
#' of growth factors if any. It takes the value passed from \code{getMIX()}.
#' @param starts A list of initial values for the parameters, either takes the value passed from \code{getMIX()} or derived by the
#' helper function \code{getMIX.initial()}.
#' @param res_cor A list where each element is a (vector of) numeric initial value(s) for residual correlation in each class. It
#' needs to be specified if the sub_Model is \code{"TVC"} (when \code{decompose != 0}), \code{"MGM"}, or \code{"MED"}. It takes the
#' value passed from \code{getMIX()}.
#'
#' @return A pre-optimized mxModel for a longitudinal mixture model.
#'
#' @keywords internal
#'
#' @importFrom OpenMx mxFitFunctionML mxExpectationMixture
#'
getMIX.mxModel <- function(dat, nClass, sub_Model, cluster_TIC, t_var, y_var, curveFun, intrinsic, records,
                           y_model, m_var, x_var, x_type, TVC, decompose, growth_TIC, starts, res_cor){
  if (sub_Model == "LGCM"){
    MIX_submodels <- getsub.LGCM_l(dat = dat, nClass = nClass, t_var = t_var, y_var = y_var, curveFun = curveFun,
                                   intrinsic = intrinsic, records = records, growth_TIC = growth_TIC, starts = starts)
  }
  else if (sub_Model == "LCSM"){
    MIX_submodels <- getsub.LCSM_l(dat = dat, nClass = nClass, t_var = t_var, y_var = y_var, curveFun = curveFun,
                                   intrinsic = intrinsic, records = records, growth_TIC = growth_TIC, starts = starts)
  }
  else if (sub_Model == "TVC"){
    MIX_submodels <- getsub.TVC_l(dat = dat, nClass = nClass, t_var = t_var, y_var = y_var, curveFun = curveFun,
                                  intrinsic = intrinsic, records = records, y_model = y_model, TVC = TVC,
                                  decompose = decompose, growth_TIC = growth_TIC, starts = starts)
  }
  else if (sub_Model == "MGM"){
    MIX_submodels <- getsub.MGM_l(dat = dat, nClass = nClass, t_var = t_var, y_var = y_var, curveFun = curveFun,
                                  intrinsic = intrinsic, records = records, y_model = y_model, starts = starts)
  }
  else if (sub_Model == "MED"){
    MIX_submodels <- getsub.MED_l(dat = dat, nClass = nClass, t_var = t_var, y_var = y_var, curveFun = curveFun,
                                  records = records, m_var = m_var, x_var = x_var, x_type = x_type, starts = starts,
                                  res_cor = res_cor)
  }
  if (!is.null(cluster_TIC)){
    classBeta <- mxMatrix(type = "Full", nrow = nClass, ncol = length(cluster_TIC) + 1,
                          free = rep(c(F, rep(T, nClass - 1)), length(cluster_TIC) + 1), values = starts[[length(starts)]],
                          labels = paste0("beta", rep(1:nClass), rep(0:length(cluster_TIC), each = nClass)),
                          name = "classbeta")
    classPV <- mxMatrix(nrow = length(cluster_TIC) + 1, ncol = 1, labels = c("ONE", paste0("data.", cluster_TIC)),
                        values = 1, name = "weightsV")
    classP <- mxAlgebra(classbeta %*% weightsV, name = "weights")
    algebraObjective <- mxExpectationMixture(paste0("Class", 1:nClass),
                                             weights = "weights", scale = "softmax")
    objective <- mxFitFunctionML()
    model_mx <- mxModel("Mixture Models with Cluster TICs", mxData(observed = dat, type = "raw"), MIX_submodels, classBeta,
                        classPV, classP, algebraObjective, objective)
  }
  else if (is.null(cluster_TIC)){
    classP <- mxMatrix("Full", nClass, 1, free = c(F, rep(T, nClass - 1)), values = starts[[length(starts)]],
                       labels = paste0("w", 1:nClass), name = "weights")
    algebraObjective <- mxExpectationMixture(paste0("Class", 1:nClass), weights = "weights", scale = "softmax")
    objective <- mxFitFunctionML()
    model_mx <- mxModel("Mixture Models without Cluster TICs",
                        mxData(observed = dat, type = "raw"), MIX_submodels, classP, algebraObjective, objective)
  }
  return(model_mx)
}

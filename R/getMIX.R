#' @title Fit a Longitudinal Mixture Model
#'
#' @description This function fits a longitudinal mixture model based on the specified sub-model. Supported submodels include:
#' \itemize{
#'   \item Latent growth curve models,
#'   \item Latent change score models,
#'   \item Latent growth curve models or latent change score models with a time-varying covariate,
#'   \item Multivariate latent growth curve models or multivariate latent change score models,
#'   \item Longitudinal mediation models.
#' } Time-invariant covariates are allowed for the first three submodels.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions for each longitudinal process, and time-invariant covariates (TICs) if any.
#' @param prop_starts A numeric vector of user-specified initial component proportions of latent classes.
#' @param sub_Model A string that specifies the sub-model for latent classes. Supported sub-models include \code{"LGCM"} (for latent
#' growth curve models), \code{"LCSM"} (for latent change score models), \code{"TVC"} (for latent growth curve models or latent change
#' score models with a time-varying covariate), \code{"MGM"} (for multivariate latent growth curve models or latent change score models),
#' and \code{"MED"} (for longitudinal mediation models).
#' @param cluster_TIC A string or character vector representing the column name(s) for time-invariant covariate(s) indicating cluster
#' formations. Default is \code{NULL}, indicating no such time-invariant covariates are present in the model.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable for each study wave. This applies when
#' \code{sub_Model} is \code{"LGCM"}, \code{"LCSM"} or \code{"TVC"}. For \code{sub_Model} being \code{"MGM"} or \code{"MED"}, \code{t_var} should
#' be a string vector where each element corresponds to the time variable prefix for each respective longitudinal process.
#' @param records A numeric vector denoting the indices of the observed study waves. This applies when \code{sub_Model} is \code{"LGCM"},
#' \code{"LCSM"} or \code{"TVC"}. For \code{sub_Model} being \code{"MGM"} or \code{"MED"}, \code{records} should be a list of numeric vectors,
#' where each vector provides the indices of the observed study waves for each longitudinal process.
#' @param y_var A string defining the prefix of the column names corresponding to the outcome variable for each study wave. This is applicable
#' when \code{sub_Model} is not \code{"MGM"}. For \code{sub_Model} being \code{"MGM"}, \code{y_var} should be a string vector where each element
#' corresponds to the prefix of the column names for each outcome variable across the study waves.
#' @param curveFun A string specifying the functional forms of the growth curve(s). Supported options for \code{y_model = "LGCM"} include:
#' \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported options for
#' \code{y_model = "LCSM"} include: \code{"nonparametric"} (or \code{"NonP"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' and \code{"Jenss-Bayley"} (or \code{"JB"}).
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. By default, this is
#' \code{NULL} as it is unnecessary when \code{sub_Model} is \code{"MED"}.
#' @param y_model A string that specifies how to fit longitudinal outcomes. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' By default, this is \code{NULL} as this argument only requires when \code{sub_Model} is \code{"TVC"} or \code{"MGM"}.
#' @param m_var A string that specifies the prefix of the column names corresponding to the mediator variable at each study wave.
#' By default, this is \code{NULL} as this argument only requires when \code{sub_Model} is \code{"MED"}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. By default, this is \code{NULL} as this argument only requires when \code{sub_Model} is \code{"MED"}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. By default, this is \code{NULL} as
#' this argument only requires when \code{sub_Model} is \code{"MED"}.
#' @param TVC A string that specifies the prefix of the column names corresponding to the time-varying covariate at each time
#' point. By default, this is \code{NULL} as this argument only requires when \code{sub_Model} is \code{"TVC"}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). By default, this
#' is \code{NULL} as this argument only requires when \code{sub_Model} is \code{"TVC"}.
#' @param growth_TIC A string or character vector of column names of time-invariant covariate(s) accounting for the variability
#' of growth factors if any. Default is \code{NULL}, indicating no growth TICs present in the model.
#' @param starts A list containing initial values for the parameters. Default is \code{NULL}, indicating no user-specified
#' initial values.
#' @param res_scale A list where each element is a (vector of) numeric scaling factor(s) for residual variance to calculate the
#' corresponding initial value for a latent class, between \code{0} and \code{1} exclusive. By default, this is \code{NULL}, as it
#' is unnecessary when the user specifies the initial values using the \code{starts} argument.
#' @param res_cor A list where each element is a (vector of) numeric initial value(s) for residual correlation in each class. It
#' needs to be specified if the sub_Model is \code{"TVC"} (when \code{decompose != 0}), \code{"MGM"}, or \code{"MED"}. By default,
#' this is \code{NULL}, as it is unnecessary when the user specifies the initial values using the \code{starts} argument.
#' @param tries An integer specifying the number of additional optimization attempts. Default is \code{NULL}.
#' @param OKStatus An integer (vector) specifying acceptable status codes for convergence. Default is \code{0}.
#' @param jitterD A string specifying the distribution for jitter. Supported values are: \code{"runif"} (uniform
#' distribution), \code{"rnorm"} (normal distribution), and \code{"rcauchy"} (Cauchy distribution). Default is \code{"runif"}.
#' @param loc A numeric value representing the location parameter of the jitter distribution. Default is \code{1}.
#' @param scale A numeric value representing the scale parameter of the jitter distribution. Default is \code{0.25}.
#' @param paramOut A logical flag indicating whether to output the parameter estimates and standard errors. Default is \code{FALSE}.
#' @param names A character vector specifying parameter names. Default is \code{NULL}.
#'
#' @return A list containing the fitted latent change score model and, if \code{paramOut = TRUE}, a data frame with parameter
#' estimates and standard errors.
#'
#' @references
#' \itemize{
#'   \item {Liu, J., & Perera, R. A. (2022). Extending Mixture of Experts Model to Investigate Heterogeneity of Trajectories:
#'   When, Where and How to Add Which Covariates. Psychological Methods (Advance online publication).
#'   \doi{10.1037/met0000436}}
#'   \item {Liu, J., & Perera, R. A. (2022). Extending Growth Mixture Model to Assess Heterogeneity in Joint Development with
#'   Piecewise Linear Trajectories in the Framework of Individual Measurement Occasions. Psychological Methods (Advance online
#'   publication). \doi{10.1037/met0000500}}
#'   \item {Liu, J., & Perera, R. A. (2023). "Estimating Rate of Change for Nonlinear Trajectories in the Framework of Individual Measurement
#'   Occasions: A New Perspective on Growth Curves." Behavior Research Methods. In press.}
#'   \item {Liu, J. (2023). Further Exploration of the Effects of Time-varying Covariate in Growth Mixture Models with Nonlinear
#'   Trajectories. arXiv. \url{https://arxiv.org/abs/2301.06014}}
#' }
#'
#' @export
#'
#' @examples
#' OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' baseT <- RMS_dat0$T1
#' RMS_dat0$T1 <- RMS_dat0$T1 - baseT
#' RMS_dat0$T2 <- RMS_dat0$T2 - baseT
#' RMS_dat0$T3 <- RMS_dat0$T3 - baseT
#' RMS_dat0$T4 <- RMS_dat0$T4 - baseT
#' RMS_dat0$T5 <- RMS_dat0$T5 - baseT
#' RMS_dat0$T6 <- RMS_dat0$T6 - baseT
#' RMS_dat0$T7 <- RMS_dat0$T7 - baseT
#' RMS_dat0$T8 <- RMS_dat0$T8 - baseT
#' RMS_dat0$T9 <- RMS_dat0$T9 - baseT
#' RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#' RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
#' RMS_dat0$gx1 <- scale(RMS_dat0$INCOME)
#' RMS_dat0$gx2 <- scale(RMS_dat0$EDU)
#' \donttest{
#' MIX_BLS_LGCM.TIC_r <-  getMIX(
#'   dat = RMS_dat0, prop_starts = c(0.45, 0.55), sub_Model = "LGCM",
#'   cluster_TIC = NULL, y_var = "M", t_var = "T", records = 1:9,
#'   curveFun = "BLS", intrinsic = FALSE, res_scale = list(0.3, 0.3)
#'   )
#' paraBLS.TIC_LGCM.r <- c(
#'   "alpha0", "alpha1", "alpha2", "knot",
#'   paste0("psi", c("00", "01", "02", "11", "12", "22")), "residuals",
#'   paste0("beta1", 0:2), paste0("beta2", 0:2),
#'   paste0("mux", 1:2), paste0("phi", c("11", "12", "22")),
#'   "mueta0", "mueta1", "mueta2"
#' )
#' set.seed(20191029)
#' MIX_BLS_LGCM.TIC_r <-  getMIX(
#'   dat = RMS_dat0, prop_starts = c(0.33, 0.34, 0.33), sub_Model = "LGCM",
#'   cluster_TIC = c("gx1", "gx2"), y_var = "M", t_var = "T", records = 1:9,
#'   curveFun = "BLS", intrinsic = FALSE, res_scale = list(0.3, 0.3, 0.3),
#'   growth_TIC = c("ex1", "ex2"), tries = 10, paramOut = TRUE,
#'   names = paraBLS.TIC_LGCM.r
#' )
#' MIX_BLS_LGCM.TIC_r[[2]]
#' }
#'
#' @importFrom OpenMx mxTryHard mxRun
#'
getMIX <- function(dat, prop_starts, sub_Model, cluster_TIC = NULL, t_var, records, y_var, curveFun, intrinsic = NULL,
                   y_model = NULL, m_var = NULL, x_type = NULL, x_var = NULL, TVC = NULL, decompose = NULL,
                   growth_TIC = NULL, starts = NULL, res_scale = NULL, res_cor = NULL, tries = NULL, OKStatus = 0,
                   jitterD = "runif", loc = 1, scale = 0.25, paramOut = FALSE, names = NULL){
  if (I(paramOut & is.null(names))){
    stop("Please enter the original parameters if want to obtain them!")
  }
  if (any(unlist(res_scale) < 0 | any(unlist(res_scale) > 1))){
    stop("Please enter a value between 0 and 1 (exclusive) for res_scale!")
  }
  if (sub_Model %in% c("LGCM", "LCSM", "TVC", "MGM")){
    if (is.null(intrinsic)){
      stop("Please specify if an intrinsic nonlinear model will be built if the sub_Model is LGCM, LCSM, MGM, or TVC!")
    }
    else if (I(intrinsic & curveFun %in% c("linear", "LIN", "quadratic", "QUAD"))){
      stop("An intrinsic nonlinear function should be one of the following: negative exponential, Jenss-Bayley, or bilinear spline!")
    }
  }
  if (I(sub_Model %in% c("MGM", "MED") | (sub_Model == "TVC" && decompose != 0))){
    if (is.null(res_cor)){
      stop("Please specify res_cor if the sub_Model is TVC (decomposed), MGM, or MED!")
    }
  }
  if (sub_Model %in% c("TVC", "MGM")){
    if (is.null(y_model)){
      stop("Please specify y_model (either LGCM or LCSM) if the sub_Model is TVC or MGM!")
    }
  }
  if (sub_Model %in% c("MED")){
    if (is.null(x_type)){
      stop("Please specify x_type for a longitudinal mediation model!")
    }
    else if (I(is.null(x_var)|is.null(m_var)|is.null(y_var))){
      stop("Please specify predictor, mediator and outcome for a longitudinal mediation model!")
    }
    else if (!I(curveFun %in% c("linear", "LIN", "bilinear spline", "BLS"))){
      stop("Longitudinal mediation model only allows for linear or bilinear spline functional form!")
    }
  }
  nClass <- length(prop_starts)
  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    starts <- getMIX.initial(dat = dat, nClass = nClass, prop_starts = prop_starts, sub_Model = sub_Model,
                             cluster_TIC = cluster_TIC, t_var = t_var, y_var = y_var, curveFun = curveFun,
                             records = records, m_var = m_var, x_var = x_var, x_type = x_type, TVC = TVC,
                             decompose = decompose, growth_TIC = growth_TIC, res_scale = res_scale,
                             res_cor = res_cor)
  }
  ## Build up a latent change score model, with or without TICs
  model_mx <- getMIX.mxModel(dat = dat, nClass = nClass, sub_Model = sub_Model, cluster_TIC = cluster_TIC,
                             t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                             records = records, y_model = y_model, m_var = m_var, x_var = x_var, x_type = x_type,
                             TVC = TVC, decompose = decompose, growth_TIC = growth_TIC, starts = starts,
                             res_cor = res_cor)
  ## Optimize the constructed latent change score model
  if (!is.null(tries)){
    model0 <- mxTryHard(model_mx, extraTries = tries, OKstatuscodes = OKStatus, jitterDistrib = jitterD,
                        loc = loc, scale = scale)
    model <- mxRun(model0)
  }
  else{
    model <- mxRun(model_mx)
  }
  ## Print out the point estimates and standard errors for the parameters of interest
  if(paramOut){
    MIX_output <- getMIX.output(model = model, nClass = nClass, sub_Model = sub_Model, cluster_TIC = cluster_TIC,
                                y_var = y_var, curveFun = curveFun, x_type = x_type, records = records,
                                growth_TIC = growth_TIC, y_model = y_model, decompose = decompose, names = names)
    return(list(model, MIX_output))
  }
  return(model)
}

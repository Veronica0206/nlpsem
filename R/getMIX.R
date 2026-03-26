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
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable for each study wave.
#' This applies when \code{sub_Model} is \code{"LGCM"}, \code{"LCSM"} or \code{"TVC"}. For \code{sub_Model} being \code{"MGM"}
#' or \code{"MED"}, \code{t_var} should be a string vector where each element corresponds to the time variable prefix for each
#' respective longitudinal process.
#' @param records A numeric vector denoting the indices of the observed study waves. This applies when \code{sub_Model} is
#' \code{"LGCM"}, \code{"LCSM"} or \code{"TVC"}. For \code{sub_Model} being \code{"MGM"} or \code{"MED"}, \code{records} should
#' be a list of numeric vectors, where each vector provides the indices of the observed study waves for each longitudinal process.
#' @param y_var A string defining the prefix of the column names corresponding to the outcome variable for each study wave. This
#' is applicable when \code{sub_Model} is not \code{"MGM"}. For \code{sub_Model} being \code{"MGM"}, \code{y_var} should be a
#' string vector where each element corresponds to the prefix of the column names for each outcome variable across the study waves.
#' @param curveFun A string specifying the functional forms of the growth curve(s). Supported options for \code{y_model = "LGCM"}
#' include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported options for \code{y_model =
#' "LCSM"} include: \code{"nonparametric"} (or \code{"NonP"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), and \code{"Jenss-Bayley"} (or \code{"JB"}).
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. By default, this is
#' \code{NULL} as it is unnecessary when \code{sub_Model} is \code{"MED"}.
#' @param y_model A string that specifies how to fit longitudinal outcomes. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' By default, this is \code{NULL} as this argument is only required when \code{sub_Model} is \code{"TVC"} or \code{"MGM"}.
#' @param m_var A string that specifies the prefix of the column names corresponding to the mediator variable at each study wave.
#' By default, this is \code{NULL} as this argument is only required when \code{sub_Model} is \code{"MED"}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. By default, this is \code{NULL} as this argument is only required when \code{sub_Model} is \code{"MED"}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. By default, this is \code{NULL} as
#' this argument is only required when \code{sub_Model} is \code{"MED"}.
#' @param TVC A string that specifies the prefix of the column names corresponding to the time-varying covariate at each time
#' point. By default, this is \code{NULL} as this argument is only required when \code{sub_Model} is \code{"TVC"}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). By default, this
#' is \code{NULL} as this argument is only required when \code{sub_Model} is \code{"TVC"}.
#' @param growth_TIC A string or character vector of column names of time-invariant covariate(s) accounting for the variability
#' of growth factors if any. Default is \code{NULL}, indicating no growth TICs present in the model.
#' @param starts A list containing initial values for the parameters. Default is \code{NULL}, indicating no user-specified
#' initial values.
#' @param res_scale An optional list where each element is a (vector of) numeric scaling factor(s) for residual variance to
#' calculate the corresponding initial value for a latent class, between \code{0} and \code{1} exclusive. Default is \code{NULL},
#' in which case data-driven residual variance estimation is used. If data-driven estimation fails, a heuristic of \code{0.1} is
#' applied as fallback.
#' @param res_cor An optional list where each element is a (vector of) numeric initial value(s) for residual correlation in each
#' class. Applicable when the sub_Model is \code{"TVC"} (when \code{decompose != 0}), \code{"MGM"}, or \code{"MED"}. Default is
#' \code{NULL}, in which case data-driven residual correlation estimation is used. If data-driven estimation fails, a heuristic
#' of \code{0.3} is applied as fallback.
#' @param tries An integer specifying the number of additional optimization attempts. Default is \code{NULL}.
#' @param OKStatus An integer (vector) specifying acceptable status codes for convergence. Default is \code{0}.
#' @param jitterD A string specifying the distribution for jitter. Supported values are: \code{"runif"} (uniform
#' distribution), \code{"rnorm"} (normal distribution), and \code{"rcauchy"} (Cauchy distribution). Default is \code{"runif"}.
#' @param loc A numeric value representing the location parameter of the jitter distribution. Default is \code{1}.
#' @param scale A numeric value representing the scale parameter of the jitter distribution. Default is \code{0.25}.
#' @param paramOut A logical flag indicating whether to output the parameter estimates and standard errors. Default is \code{FALSE}.
#' @param names A character vector specifying parameter names. Default is \code{NULL}, in which case
#' meaningful names are automatically generated based on the model configuration.
#'
#' @return An object of class \code{myMxOutput}. Depending on the \code{paramOut} argument, the object may contain the following slots:
#' \itemize{
#'   \item \code{mxOutput}: This slot contains the fitted longitudinal mixture model. A summary of this model can be obtained using
#'   the \code{ModelSummary()} function.
#'   \item \code{Estimates} (optional): If \code{paramOut = TRUE}, a data frame with parameter estimates and standard errors. The content
#'   of this slot can be printed using the \code{printTable()} method for S4 objects.
#' }
#'
#' @references
#' \itemize{
#'   \item {Liu, J., & Perera, R. A. (2022). Extending Mixture of Experts Model to Investigate Heterogeneity of Trajectories:
#'   When, Where and How to Add Which Covariates. Psychological Methods, 28(1), 152-178. \doi{10.1037/met0000436}}
#'   \item {Liu, J., & Perera, R. A. (2022). Extending Growth Mixture Model to Assess Heterogeneity in Joint Development with
#'   Piecewise Linear Trajectories in the Framework of Individual Measurement Occasions. Psychological Methods, 28(5), 1029-1051.
#'   \doi{10.1037/met0000500}}
#'   \item {Liu, J., & Perera, R. A. (2023). Estimating Rate of Change for Nonlinear Trajectories in the Framework of Individual Measurement
#'   Occasions: A New Perspective on Growth Curves. Behavior Research Methods, 56(3), 1349-1375. \doi{10.3758/s13428-023-02097-2}}
#'   \item {Liu, J. (2023). Further Exploration of the Effects of Time-varying Covariate in Growth Mixture Models with Nonlinear
#'   Trajectories. Behavior Research Methods, 56(4), 2804-2827. \doi{10.3758/s13428-023-02183-5}}
#' }
#'
#' @seealso \code{\link{getLGCM}}, \code{\link{getLCSM}}, \code{\link{getMGroup}}, \code{\link{getFigure}}
#'
#' @export
#'
#' @examples
#' mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' # Re-baseline the data so that the estimated initial status is for the starting point of the study
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
#' # Fit longitudinal mixture group model of bilinear spline functional form with fixed knot
#' # (2 classes)
#' MIX_BLS_LGCM.TIC_r <-  getMIX(
#'   dat = RMS_dat0, prop_starts = c(0.45, 0.55), sub_Model = "LGCM",
#'   cluster_TIC = NULL, y_var = "M", t_var = "T", records = 1:9,
#'   curveFun = "BLS", intrinsic = FALSE
#'   )
#' # Fit longitudinal mixture group model of bilinear spline functional form with fixed knot
#' # (3 classes)
#' set.seed(20191029)
#' MIX_BLS_LGCM.TIC_r <-  getMIX(
#'   dat = RMS_dat0, prop_starts = c(0.33, 0.34, 0.33), sub_Model = "LGCM",
#'   cluster_TIC = c("gx1", "gx2"), y_var = "M", t_var = "T", records = 1:9,
#'   curveFun = "BLS", intrinsic = FALSE,
#'   growth_TIC = c("ex1", "ex2"), tries = 10, paramOut = TRUE
#' )
#' printTable(MIX_BLS_LGCM.TIC_r)
#' }
#'
#' @importFrom OpenMx mxTryHard mxRun
#' @importFrom methods new
#'
getMIX <- function(dat, prop_starts, sub_Model, cluster_TIC = NULL, t_var, records, y_var, curveFun, intrinsic = NULL,
                   y_model = NULL, m_var = NULL, x_type = NULL, x_var = NULL, TVC = NULL, decompose = NULL,
                   growth_TIC = NULL, starts = NULL, res_scale = NULL, res_cor = NULL, tries = NULL, OKStatus = 0,
                   jitterD = "runif", loc = 1, scale = 0.25, paramOut = FALSE, names = NULL){
  dat <- as.data.frame(dat)
  validate_paramOut(paramOut, names)
  validate_res_scale(res_scale, is_list = TRUE)
  validate_res_cor(res_cor)
  validate_curveFun(curveFun)
  validate_decompose(decompose)
  validate_x_type(x_type)
  validate_sub_Model(sub_Model, intrinsic = intrinsic, curveFun = curveFun,
                     res_cor = res_cor, y_model = y_model, x_type = x_type,
                     x_var = x_var, m_var = m_var, y_var = y_var,
                     decompose = decompose)
  validate_columns(dat, t_var = t_var, y_var = y_var, records = records,
                   growth_TIC = growth_TIC, TVC = TVC, m_var = m_var,
                   x_var = x_var, x_type = x_type, cluster_TIC = cluster_TIC)
  nClass <- length(prop_starts)
  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    starts <- getMIX.initial(dat = dat, nClass = nClass, prop_starts = prop_starts, sub_Model = sub_Model,
                             cluster_TIC = cluster_TIC, t_var = t_var, y_var = y_var, curveFun = curveFun,
                             records = records, m_var = m_var, x_var = x_var, x_type = x_type, TVC = TVC,
                             decompose = decompose, growth_TIC = growth_TIC, res_scale = res_scale,
                             res_cor = res_cor, intrinsic = intrinsic)
  }
  ## Build up a mixture model
  model_mx <- getMIX.mxModel(dat = dat, nClass = nClass, sub_Model = sub_Model, cluster_TIC = cluster_TIC,
                             t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                             records = records, y_model = y_model, m_var = m_var, x_var = x_var, x_type = x_type,
                             TVC = TVC, decompose = decompose, growth_TIC = growth_TIC, starts = starts,
                             res_cor = res_cor)
  ## Optimize the constructed mixture model
  if (!is.null(tries)){
    model <- mxTryHard(model_mx, extraTries = tries, OKstatuscodes = OKStatus, jitterDistrib = jitterD,
                        loc = loc, scale = scale)
  }
  else{
    model <- mxRun(model_mx)
  }
  ## Print out the point estimates and standard errors for the parameters of interest
  if(paramOut){
    if (is.null(names)) {
      names <- .auto_names_subModel(sub_Model, curveFun, intrinsic, growth_TIC,
                                     records, y_model, decompose, y_var, m_var, x_type, x_var)
    }
    MIX_output <- getMIX.output(model = model, nClass = nClass, sub_Model = sub_Model, cluster_TIC = cluster_TIC,
                                y_var = y_var, curveFun = curveFun, x_type = x_type, records = records,
                                growth_TIC = growth_TIC, y_model = y_model, decompose = decompose, names = names)
    model <- new("myMxOutput", mxOutput = model, Estimates = MIX_output)
  }
  else{
    model <- new("myMxOutput", mxOutput = model)
  }
  return(model)
}

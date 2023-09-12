#' @title Fit a Latent Growth Curve Model or Latent Change Score Model with Time-varying and Time-invariant Covariates
#'
#' @description This function fits a latent growth curve model or latent change score model with a time-varying covariate
#' and potential time-invariant covariates to the provided data. It manages model setup, optimization, and if requested,
#' outputs parameter estimates and standard errors.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements (for the longitudinal outcome and time-varying covariates), occasions, and time-invariant covariates
#' (TICs) if any.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"nonparametric"} (or \code{"NonP"}), \code{"quadratic"} (or \code{"QUAD"}),
#' \code{"negative exponential"} (or \code{"EXP"}), and \code{"Jenss-Bayley"} (or \code{"JB"}.
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. Default is
#' \code{TRUE}.
#' @param records A numeric vector specifying the indices of the observed study waves.
#' @param y_model A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' @param TVC A string specifying the prefix of the column names corresponding to the time-varying covariate at each study wave.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states).
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) that account for the
#' variability of growth factors, if any. Default is \code{NULL}, indicating no growth TICs present in the model.
#' @param starts A list containing initial values for the parameters. Default is \code{NULL}, indicating no user-specified initial
#' values.
#' @param res_scale A numeric value or numeric vector. For a model with \code{decompose = 0}, it is a numeric value representing
#' the scaling factor used to calculate the initial value for the residual variance of the longitudinal outcome. In cases where
#' \code{decompose != 0}, it is a numeric vector of user-specified scaling factors used to calculate the initial values for the
#' residual variance of both the longitudinal outcome and the time-varying covariate. By default, this is \code{NULL}, as it is
#' unnecessary when the user specifies the initial values using the \code{starts} argument.
#' @param res_cor A numeric value. When \code{decompose != 0}, this represents the user-specified residual correlation between the
#' longitudinal outcome and the time-varying covariate, which is used to calculate the corresponding initial value. If \code{decompose = 0},
#' this should be \code{NULL}. By default, this is \code{NULL}, as it is unnecessary when the user specifies the initial values
#' using the \code{starts} argument.
#' @param tries An integer specifying the number of additional optimization attempts. Default is \code{NULL}.
#' @param OKStatus An integer (vector) specifying acceptable status codes for convergence. Default is \code{0}.
#' @param jitterD A string specifying the distribution for jitter. Supported values are: \code{"runif"} (uniform
#' distribution), \code{"rnorm"} (normal distribution), and \code{"rcauchy"} (Cauchy distribution). Default is \code{"runif"}.
#' @param loc A numeric value representing the location parameter of the jitter distribution. Default is \code{1}.
#' @param scale A numeric value representing the scale parameter of the jitter distribution. Default is \code{0.25}.
#' @param paramOut A logical flag indicating whether to output the parameter estimates and standard errors. Default is \code{FALSE}.
#' @param names A character vector specifying parameter names. Default is \code{NULL}.
#'
#' @return An object of class \code{myMxOutput}. Depending on the \code{paramOut} argument, the object may contain the following slots:
#' \itemize{
#'   \item \code{mxOutput}: This slot contains the fitted latent growth curve model or latent change score model with a time-varying
#'   covariate. A summary of this model can be obtained using the \code{ModelSummary()} function.
#'   \item \code{Estimates} (optional): If \code{paramOut = TRUE}, a data frame with parameter estimates and standard errors. The content
#'   of this slot can be printed using the \code{printTable()} method for S4 objects.
#' }
#'
#' @references
#' \itemize{
#'   \item {Liu, J., & Perera, R. A. (2023). Estimating Rate of Change for Nonlinear Trajectories in the Framework of Individual Measurement
#'   Occasions: A New Perspective on Growth Curves. Behavior Research Methods. \doi{10.3758/s13428-023-02097-2}}
#'   \item {Liu, J. (2022). "Decomposing Impact on Longitudinal Outcome of Time-varying Covariate into Baseline Effect and Temporal Effect."
#'   arXiv. https://arxiv.org/abs/2210.16916}
#' }
#'
#' @export
#'
#' @examples
#' mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' # Re-baseline the data so that the estimated initial status is for the starting point of the study
#' baseT <- RMS_dat0$T1
#' RMS_dat0$T1 <- (RMS_dat0$T1 - baseT)/12
#' RMS_dat0$T2 <- (RMS_dat0$T2 - baseT)/12
#' RMS_dat0$T3 <- (RMS_dat0$T3 - baseT)/12
#' RMS_dat0$T4 <- (RMS_dat0$T4 - baseT)/12
#' RMS_dat0$T5 <- (RMS_dat0$T5 - baseT)/12
#' RMS_dat0$T6 <- (RMS_dat0$T6 - baseT)/12
#' RMS_dat0$T7 <- (RMS_dat0$T7 - baseT)/12
#' RMS_dat0$T8 <- (RMS_dat0$T8 - baseT)/12
#' RMS_dat0$T9 <- (RMS_dat0$T9 - baseT)/12
#' RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#' RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
#'
#' # Standardize reading ability over time with its baseline value
#' BL_mean <- mean(RMS_dat0[, "R1"])
#' BL_var <- var(RMS_dat0[, "R1"])
#' RMS_dat0$Rs1 <- (RMS_dat0$R1 - BL_mean)/sqrt(BL_var)
#' RMS_dat0$Rs2 <- (RMS_dat0$R2 - BL_mean)/sqrt(BL_var)
#' RMS_dat0$Rs3 <- (RMS_dat0$R3 - BL_mean)/sqrt(BL_var)
#' RMS_dat0$Rs4 <- (RMS_dat0$R4 - BL_mean)/sqrt(BL_var)
#' RMS_dat0$Rs5 <- (RMS_dat0$R5 - BL_mean)/sqrt(BL_var)
#' RMS_dat0$Rs6 <- (RMS_dat0$R6 - BL_mean)/sqrt(BL_var)
#' RMS_dat0$Rs7 <- (RMS_dat0$R7 - BL_mean)/sqrt(BL_var)
#' RMS_dat0$Rs8 <- (RMS_dat0$R8 - BL_mean)/sqrt(BL_var)
#' RMS_dat0$Rs9 <- (RMS_dat0$R9 - BL_mean)/sqrt(BL_var)
#'
#' \donttest{
#' # Fit bilinear spline latent growth curve model (fixed knot) with a time-varying
#' # reading ability for mathematics development
#' BLS_TVC_LGCM1 <- getTVCmodel(
#'  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
#'  records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 0,  growth_TIC = NULL,
#'  res_scale = 0.1
#'  )
#'
#' # Fit negative exponential latent growth curve model (random ratio) with a
#' # decomposed time-varying reading ability and time-invariant covariates for
#' # mathematics development
#' paraEXP_LGCM3.f <- c(
#'   "Y_alpha0", "Y_alpha1", "Y_alphag",
#'   paste0("Y_psi", c("00", "01", "0g", "11", "1g", "gg")), "Y_residuals",
#'   "X_mueta0", "X_mueta1", paste0("X_psi", c("00", "01", "11")),
#'   paste0("X_rel_rate", 2:8), paste0("X_abs_rate", 1:8), "X_residuals",
#'   paste0("betaTIC", c(0:1, "g")), paste0("betaTIC", c(0:1, "g")),
#'   paste0("betaTVC", c(0:1, "g")),
#'   "muTIC1", "muTIC2", "phiTIC11", "phiTIC12", "phiTIC22",
#'   "Y_mueta0", "Y_mueta1", "Y_mu_slp_ratio",
#'   "covBL1", "covBL2", "kappa", "Cov_XYres")
#' set.seed(20191029)
#' EXP_TVCslp_LGCM3.f <- getTVCmodel(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "EXP", intrinsic = TRUE,
#'   records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 1,
#'   growth_TIC = c("ex1", "ex2"), res_scale = c(0.1, 0.1),
#'   res_cor = 0.3, tries = 10, paramOut = TRUE, names = paraEXP_LGCM3.f
#' )
#' printTable(EXP_TVCslp_LGCM3.f)
#' }
#'
#' @importFrom OpenMx mxTryHard mxRun
#' @importFrom methods new
#'
getTVCmodel <- function(dat, t_var, y_var, curveFun, intrinsic = TRUE, records, y_model, TVC, decompose,
                        growth_TIC = NULL, starts = NULL, res_scale = NULL, res_cor = NULL, tries = NULL,
                        OKStatus = 0, jitterD = "runif", loc = 1, scale = 0.25, paramOut = FALSE, names = NULL){
  if (I(paramOut & is.null(names))){
    stop("Please enter the original parameters if want to obtain them!")
  }
  if (I(any(res_scale <= 0) | any(res_scale >= 1))){
    stop("Please enter a value between 0 and 1 (exclusive) for res_scale!")
  }
  if (I(intrinsic & curveFun %in% c("linear", "LIN", "quadratic", "QUAD", "nonparametric", "NonP"))){
    stop("An intrinsic nonlinear function should be negative exponential, Jenss-Bayley, or bilinear spline!")
  }

  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    starts <- getTVC.initial(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, records = records,
                             growth_TIC = growth_TIC, TVC = TVC, decompose = decompose, res_scale = res_scale,
                             res_cor = res_cor)
  }
  ## Build up a latent growth curve model or latent change score model with a decomposed TVC (and TICs if any)
  model_mx <- getTVC.mxModel(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                             records = records, growth_TIC = growth_TIC, TVC = TVC, decompose = decompose,
                             y_model = y_model, starts = starts)
  ## Optimize the constructed latent growth curve model or latent change score model with a decomposed TVC
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
    TVC_output <- getTVC.output(model = model, curveFun = curveFun, records = records, y_model = y_model,
                                decompose = decompose, growth_TIC = growth_TIC, names = names)
    model <- new("myMxOutput", mxOutput = model, Estimates = TVC_output)
  }
  else{
    model <- new("myMxOutput", mxOutput = model)
  }
  return(model)
}


#' @title Fit a Latent Growth Curve Model with Time-invariant Covariate (If Any)
#'
#' @description This function fits a latent growth curve model with or without time-invariant covariates to the provided data.
#' It manages model setup, optimization, and if requested, outputs parameter estimates and standard errors.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions, and time-invariant covariates (TICs) if any.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for latent growth curve
#' models are: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}).
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. Default is
#' \code{TRUE}.
#' @param records A numeric vector specifying indices of the study waves.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing to the
#' variability of growth factors if any. Default is \code{NULL}, indicating no growth TICs are included in the model.
#' @param starts A list containing initial values for the parameters. Default is \code{NULL}, indicating no user-specified
#' initial values.
#' @param res_scale A numeric value representing the scaling factor for the initial calculation of the residual variance. This
#' value should be between \code{0} and \code{1}, exclusive. By default, this is \code{NULL}, as it is unnecessary when the
#' user specifies the initial values using the \code{starts} argument.
#' @param tries An integer specifying the number of additional optimization attempts. Default is \code{NULL}.
#' @param OKStatus An integer (vector) specifying acceptable status codes for convergence. Default is \code{0}.
#' @param jitterD A string specifying the distribution for jitter. Supported values are: \code{"runif"} (uniform
#' distribution), \code{"rnorm"} (normal distribution), and \code{"rcauchy"} (Cauchy distribution). Default is \code{"runif"}.
#' @param loc A numeric value representing the location parameter of the jitter distribution. Default is \code{1}.
#' @param scale A numeric value representing the scale parameter of the jitter distribution. Default is \code{0.25}.
#' @param paramOut A logical flag indicating whether to output the parameter estimates and standard errors. Default is \code{FALSE}.
#' @param names A character vector specifying parameter names. Default is \code{NULL}.
#'
#' @return A list containing the fitted latent growth curve model and, if \code{paramOut = TRUE}, a data frame with parameter
#' estimates and standard errors.
#'
#' @references
#' \itemize{
#'   \item {Liu, J., Perera, R. A., Kang, L., Kirkpatrick, R. M., & Sabo, R. T. (2021). "Obtaining Interpretable Parameters from
#'   Reparameterizing Longitudinal Models: Transformation Matrices between Growth Factors in Two Parameter Spaces". Journal of
#'   Educational and Behavioral Statistics. \doi{10.3102/10769986211052009}}
#'   \item {Sterba, S. K. (2014). "Fitting Nonlinear Latent Growth Curve Models With Individually Varying Time Points". Structural
#'   Equation Modeling: A Multidisciplinary Journal, 21(4), 630-647. \doi{10.1080/10705511.2014.919828}}
#' }
#'
#' @export
#'
#' @examples
#' OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' # Load ECLS-K (2011) data
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' # Re-baseline the data so that the estimated initial status is for the
#' # starting point of the study
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
#' # Standardized time-invariant covariates
#' RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#' RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
#'
#' \donttest{
#' # Fit bilinear spline latent growth curve model (fixed knots)
#' BLS_LGCM_r <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#'   intrinsic = FALSE, records = 1:9, growth_TIC = NULL, res_scale = 0.1
#' )
#' # Fit bilinear spline latent growth curve model (random knots) with
#' # time-invariant covariates for mathematics development
#' ## Define parameter names
#' paraBLS.TIC_LGCM.f <- c(
#'   "alpha0", "alpha1", "alpha2", "alphag",
#'   paste0("psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")),
#'   "residuals", paste0("beta1", c(0:2, "g")), paste0("beta2", c(0:2, "g")),
#'   paste0("mux", 1:2), paste0("phi", c("11", "12", "22")),
#'   "mueta0", "mueta1", "mueta2", "mu_knot"
#' )
#' ## Fit the model
#' BLS_LGCM.TIC_f <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#'   intrinsic = TRUE, records = 1:9, growth_TIC = c("ex1", "ex2"), res_scale = 0.1,
#'   paramOut = TRUE, names = paraBLS.TIC_LGCM.f
#' )
#' ## Output point estimate and standard errors
#' BLS_LGCM.TIC_f[[2]]
#' }
#'
#' @importFrom OpenMx mxTryHard mxRun
#'
getLGCM <- function(dat, t_var, y_var, curveFun, intrinsic = TRUE, records, growth_TIC = NULL, starts = NULL,
                    res_scale = NULL, tries = NULL, OKStatus = 0, jitterD = "runif", loc = 1, scale = 0.25,
                    paramOut = FALSE, names = NULL){
  if (paramOut & is.null(names)){
    stop("Please provide the original parameters if you want to obtain them!")
  }
  if (res_scale <= 0 | res_scale >= 1){
    stop("Please enter a value between 0 and 1 (exclusive) for res_scale!")
  }
  if (intrinsic & curveFun %in% c("linear", "LIN", "quadratic", "QUAD")){
    stop("An intrinsic nonlinear function should be negative exponential, Jenss-Bayley, or bilinear spline for a LGCM!")
  }
  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    starts <- getUNI.initial(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, records = records,
                             growth_TIC = growth_TIC, res_scale = res_scale)
  }
  ## Build up a latent growth curve model, with or without TICs
  model_mx <- getLGCM.mxModel(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                              records = records, growth_TIC = growth_TIC, starts = starts)
  ## Optimize the constructed latent growth curve model
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
    LGCM_output <- getLGCM.output(model = model, growth_TIC = growth_TIC, names = names)
    return(list(model, LGCM_output))
  }
  return(model)
}


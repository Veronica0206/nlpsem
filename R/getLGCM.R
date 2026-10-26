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
#' @param res_scale An optional numeric value representing the scaling factor for the initial calculation of the residual
#' variance. This value should be between \code{0} and \code{1}, exclusive. Default is \code{NULL}, in which case data-driven
#' residual variance estimation is used. If data-driven estimation fails, a heuristic of \code{0.1} is applied as fallback.
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
#'   \item \code{mxOutput}: This slot contains the fitted latent growth curve model. A summary of this model can be obtained using the
#'   \code{ModelSummary()} function.
#'   \item \code{Estimates} (optional): If \code{paramOut = TRUE}, a data frame with parameter estimates and standard errors. The content
#'   of this slot can be printed using the \code{printTable()} method for S4 objects.
#' }
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
#' @seealso \code{\link{getLCSM}}, \code{\link{getTVCmodel}}, \code{\link{getFigure}}
#'
#' @details New fits retain their resolved model specification and fitting provenance.
#' Use \code{\link{getSpecification}} to inspect it; standard \code{summary()}, \code{coef()},
#' \code{logLik()} and \code{plot()} methods are available.
#'
#' @section Standardized fitting:
#' Models with a quadratic, negative exponential or Jenss-Bayley functional form are fitted to standardized data and converted back to the original units. Each outcome is
#' centred at the mean and divided by the standard deviation of its first recorded occasion, a
#' time-varying covariate is divided by its first-occasion standard deviation, and time is divided by the
#' mean of its last recorded occasion; time-invariant covariates are not changed. The model starts from the
#' same starting values, mapped exactly to the standardized scale, and is optimized with the same optimizer
#' and retry settings. The estimates (means, variances, covariances, residual variances, path coefficients,
#' growth ratios and, with a time-varying covariate, its effects), the Hessian, the parameter covariance
#' matrix, the standard errors and the -2 log-likelihood are converted back exactly, so the result is the
#' same model in the original units; only the numerical conditioning of the optimization differs.
#' \code{getSpecification()} records the centres, scales and condition numbers under
#' \code{fitting$standardization}. Set \code{options(nlpsem.standardize = FALSE)} to fit on the original scale.
#'
#' @export
#'
#' @examples
#' mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' # Load the ECLS-K (2011) sample and re-baseline the measurement occasions (in months)
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' baseT <- RMS_dat0$T1
#' for (j in 1:9) RMS_dat0[[paste0("T", j)]] <- RMS_dat0[[paste0("T", j)]] - baseT
#' # Standardized time-invariant covariates
#' RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#' RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
#'
#' \donttest{
#' # Bilinear spline latent growth curve model with a fixed knot
#' BLS_LGCM_r <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#'   intrinsic = FALSE, records = 1:9, paramOut = TRUE
#' )
#' printTable(BLS_LGCM_r)
#' }
#' \dontrun{
#' # Intrinsically nonlinear model (random knot) with time-invariant covariates and 20 extra
#' # optimization attempts; this takes several minutes
#' set.seed(20191029)
#' BLS_LGCM.TIC_f <- getLGCM(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
#'   intrinsic = TRUE, records = 1:9, growth_TIC = c("ex1", "ex2"),
#'   tries = 20, paramOut = TRUE
#' )
#' printTable(BLS_LGCM.TIC_f)
#' }
#'
#' @importFrom OpenMx mxTryHard mxRun
#' @importFrom methods new
#'
getLGCM <- function(dat, t_var, y_var, curveFun, intrinsic = TRUE, records, growth_TIC = NULL, starts = NULL,
                    res_scale = NULL, tries = NULL, OKStatus = 0, jitterD = "runif", loc = 1, scale = 0.25,
                    paramOut = FALSE, names = NULL){
  dat <- as.data.frame(dat)
  validate_paramOut(paramOut, names)
  validate_res_scale(res_scale)
  validate_curveFun(curveFun)
  validate_intrinsic(intrinsic, curveFun, model_type = "LGCM")
  validate_columns(dat, t_var = t_var, y_var = y_var, records = records,
                   growth_TIC = growth_TIC)
  if (paramOut && !is.null(names)) {
    expected_names <- .auto_names_LGCM(curveFun, intrinsic, growth_TIC)
    names <- .validate_output_names(names, length(expected_names))
  }
  metadata <- .capture_fit_metadata(match.call(), "getLGCM", environment(), dat)
  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    metadata$fitting$starts_source <- "automatic"
    starts <- getUNI.initial(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, records = records,
                             growth_TIC = growth_TIC, res_scale = res_scale, intrinsic = intrinsic)
  }
  ## Build up a latent growth curve model, with or without TICs
  model_mx <- getLGCM.mxModel(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                              records = records, growth_TIC = growth_TIC, starts = starts)
  ## Optimize the constructed latent growth curve model
  ## Quadratic, negative exponential and Jenss-Bayley models are fitted on the standardized scale and
  ## converted back exactly (see ?getLGCM, Details); other forms are fitted as built.
  fit_std <- .fit_model(model_mx, run = function(model_mx){
    if (!is.null(tries)){
      mxTryHard(model_mx, extraTries = tries, OKstatuscodes = OKStatus, jitterDistrib = jitterD,
                loc = loc, scale = scale)
    }
    else{
      mxRun(model_mx)
    }
  }, plan = .std_plan(dat, "LGCM", curveFun, y_var, records, t_var),
  build = function(plan) getLGCM.mxModel(dat = .std_data(dat, plan), t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                                records = records, growth_TIC = growth_TIC, starts = starts))
  model <- fit_std$model
  metadata$fitting$standardization <- fit_std$standardization
  ## Print out the point estimates and standard errors for the parameters of interest
  if(paramOut){
    if (is.null(names)) names <- .auto_names_LGCM(curveFun, intrinsic, growth_TIC)
    LGCM_output <- getLGCM.output(model = model, growth_TIC = growth_TIC, names = names)
    model <- new("myMxOutput", mxOutput = model, Estimates = LGCM_output)
  }
  else{
    model <- new("myMxOutput", mxOutput = model)
  }
  model <- .attach_fit_metadata(model, metadata, starts, frame = environment())
  return(model)
}

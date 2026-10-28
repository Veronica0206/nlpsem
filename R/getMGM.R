#' @title Fit a Multivariate Latent Growth Curve Model or Multivariate Latent Change Score Model
#'
#' @description This function fits a multivariate latent growth curve model or a multivariate latent change score model
#' with the provided data. It manages model setup, optimization, and if requested, outputs parameter estimates and
#' standard errors.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables
#' with repeated measurements and occasions for multiple longitudinal outcomes.
#' @param t_var A vector of strings, with each element representing the prefix for column names related to the time
#' variable for the corresponding outcome variable at each study wave.
#' @param y_var A vector of strings, with each element representing the prefix for column names corresponding to a
#' particular outcome variable at each study wave. The parameters shared by two outcomes are labelled by joining their
#' prefixes (\code{"Y"} and \code{"Z"} give \code{"YZ"}), so prefixes for which a joined pair equals a prefix or another
#' pair (such as \code{"A"}, \code{"B"} and \code{"AB"}) are rejected; with \code{y_model = "LCSM"}, the prefixes must
#' also differ in more than letter case.
#' @param curveFun A string specifying the functional forms of the growth curve(s). Supported options for \code{y_model = "LGCM"} include:
#' \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported options for
#' \code{y_model = "LCSM"} include: \code{"nonparametric"} (or \code{"NonP"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' and \code{"Jenss-Bayley"} (or \code{"JB"}).
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. Default is
#' \code{TRUE}.
#' @param records A list of numeric vectors, with each vector specifying the indices of the observed study waves for
#' the corresponding outcome variable. Every pair of outcomes must share at least one
#' record ID to estimate its residual covariance. Residual pairs are matched by record
#' ID, independently of the order of the selected columns.
#' @param y_model A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' @param starts A list containing initial values for the parameters. Default is \code{NULL}, indicating no user-specified initial
#' values.
#' @param res_scale An optional numeric vector with each element representing the scaling factor for the initial calculation of
#' the residual variance. These values should be between \code{0} and \code{1}, exclusive. Default is \code{NULL}, in which case
#' data-driven residual variance estimation is used. If data-driven estimation fails, a heuristic of \code{0.1} is applied as
#' fallback.
#' @param res_cor An optional numeric value or vector for user-specified residual correlation between any two longitudinal outcomes
#' to calculate the corresponding initial value. Default is \code{NULL}, in which case data-driven residual correlation estimation
#' is used. If data-driven estimation fails, a heuristic of \code{0.3} is applied as fallback.
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
#'   \item \code{mxOutput}: This slot contains the fitted multivariate latent growth curve model or a multivariate latent change score
#'   model. A summary of this model can be obtained using the \code{ModelSummary()} function.
#'   \item \code{Estimates} (optional): If \code{paramOut = TRUE}, a data frame with parameter estimates and standard errors. The content
#'   of this slot can be printed using the \code{printTable()} method for S4 objects.
#' }
#'
#' @references
#' \itemize{
#'   \item {Liu, J., & Perera, R. A. (2021). "Estimating Knots and Their Association in Parallel Bilinear Spline Growth Curve Models in
#'   the Framework of Individual Measurement Occasions," Psychological Methods, 27(5), 703-729. \doi{10.1037/met0000309}}
#'   \item {Blozis, S. A. (2004). "Structured Latent Curve Models for the Study of Change in Multivariate Repeated Measures," Psychological
#'   Methods, 9(3), 334-353. \doi{10.1037/1082-989X.9.3.334}}
#' }
#'
#' @seealso \code{\link{getLGCM}}, \code{\link{getMGroup}}, \code{\link{getMIX}}, \code{\link{getFigure}}
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
#'
#' \donttest{
#' # Linear multivariate latent growth curve model for reading and mathematics
#' LIN_PLGCM_f <- getMGM(
#'   dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "LIN",
#'   intrinsic = FALSE, records = list(1:9, 1:9), y_model = "LGCM", paramOut = TRUE
#' )
#' printTable(LIN_PLGCM_f)
#' }
#' \dontrun{
#' # Bilinear spline multivariate model with random knots and 20 extra optimization
#' # attempts; this takes several minutes
#' set.seed(20191029)
#' BLS_PLGCM_f <- getMGM(
#'   dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS",
#'   intrinsic = TRUE, records = list(1:9, 1:9), y_model = "LGCM",
#'   tries = 20, paramOut = TRUE
#' )
#' printTable(BLS_PLGCM_f)
#' }
#'
#' @importFrom OpenMx mxTryHard mxRun
#' @importFrom methods new
#'
getMGM <- function(dat, t_var, y_var, curveFun, intrinsic = TRUE, records, y_model, starts = NULL, res_scale = NULL,
                   res_cor = NULL, tries = NULL, OKStatus = 0, jitterD = "runif", loc = 1, scale = 0.25,
                   paramOut = FALSE, names = NULL){
  dat <- as.data.frame(dat)
  validate_paramOut(paramOut, names)
  validate_res_scale(res_scale)
  validate_res_cor(res_cor)
  validate_curveFun(curveFun)
  if (length(t_var) != length(y_var) || length(t_var) != length(records)) {
    stop("Lengths of t_var, y_var, and records must be equal!")
  }
  validate_intrinsic(intrinsic, curveFun, model_type = "general")
  validate_curveFun_ymodel(curveFun, y_model)
  validate_columns(dat, t_var = t_var, y_var = y_var, records = records)
  validate_mgm_names(y_var, y_model)
  if (paramOut && !is.null(names)) {
    expected_names <- .auto_names_MGM(curveFun, intrinsic, y_var, records, y_model)
    names <- .validate_output_names(names, length(expected_names))
  }
  metadata <- .capture_fit_metadata(match.call(), "getMGM", environment(), dat)
  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    metadata$fitting$starts_source <- "automatic"
    starts <- getMULTI.initial(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, records = records,
                               res_scale = res_scale, res_cor = res_cor, intrinsic = intrinsic)
  }
  ## Build up a multivariate latent growth curve model or latent change score model
  ### Obtain manifest and latent variables, paths of the longitudinal outcomes
  model_mx <- getMGM.mxModel(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                             records = records, y_model = y_model, starts = starts)
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
  }, plan = .std_plan(dat, "MGM", curveFun, y_var, records, t_var),
  build = function(plan) getMGM.mxModel(dat = .std_data(dat, plan), t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                               records = records, y_model = y_model, starts = starts))
  model <- fit_std$model
  metadata$fitting$standardization <- fit_std$standardization
  if(paramOut){
    if (is.null(names)) names <- .auto_names_MGM(curveFun, intrinsic, y_var, records, y_model)
    MGM_output <- getMGM.output(model = model, y_var = y_var, records = records, curveFun = curveFun, y_model = y_model,
                                names = names)
    model <- new("myMxOutput", mxOutput = model, Estimates = MGM_output)
  }
  else{
    model <- new("myMxOutput", mxOutput = model)
  }
  model <- .attach_fit_metadata(model, metadata, starts, frame = environment())
  return(model)
}

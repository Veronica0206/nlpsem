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
#' \code{"negative exponential"} (or \code{"EXP"}), and \code{"Jenss-Bayley"} (or \code{"JB"}).
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
#' @param res_scale An optional numeric value or numeric vector. For a model with \code{decompose = 0}, it is a numeric value
#' representing the scaling factor used to calculate the initial value for the residual variance of the longitudinal outcome. In
#' cases where \code{decompose != 0}, it is a numeric vector of scaling factors used to calculate the initial values for the
#' residual variance of both the longitudinal outcome and the time-varying covariate. Default is \code{NULL}, in which case
#' data-driven residual variance estimation is used. If data-driven estimation fails, a heuristic of \code{0.1} is applied as
#' fallback.
#' @param res_cor An optional numeric value. When \code{decompose != 0}, this represents the user-specified residual correlation
#' between the longitudinal outcome and the time-varying covariate, which is used to calculate the corresponding initial value. If
#' \code{decompose = 0}, this should be \code{NULL}. Default is \code{NULL}, in which case data-driven residual correlation
#' estimation is used when \code{decompose != 0}. If data-driven estimation fails, a heuristic of \code{0.3} is applied as fallback.
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
#'   \item \code{mxOutput}: This slot contains the fitted latent growth curve model or latent change score model with a time-varying
#'   covariate. A summary of this model can be obtained using the \code{ModelSummary()} function.
#'   \item \code{Estimates} (optional): If \code{paramOut = TRUE}, a data frame with parameter estimates and standard errors. The content
#'   of this slot can be printed using the \code{printTable()} method for S4 objects.
#' }
#'
#' @references
#' \itemize{
#'   \item {Liu, J., & Perera, R. A. (2023). Estimating Rate of Change for Nonlinear Trajectories in the Framework of Individual Measurement
#'   Occasions: A New Perspective on Growth Curves. Behavior Research Methods, 56(3), 1349-1375. \doi{10.3758/s13428-023-02097-2}}
#'   \item {Liu, J. (2025). Decomposing Impact on Longitudinal Outcome of Time-varying Covariate into Baseline Effect and Temporal Effect.
#'   Journal of Educational and Behavioral Statistics, 50(5), 765-805. \doi{10.3102/10769986241263975}}
#' }
#'
#' @seealso \code{\link{getLGCM}}, \code{\link{getLCSM}}, \code{\link{getFigure}}
#'
#' @details New fits retain their resolved model specification and fitting provenance.
#' A new fit whose \code{\link{getDiagnostics}} reports issues (for example a nonzero optimizer status or an
#' indefinite covariance matrix) is returned with a warning that lists them.
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
#' @section Approximations:
#' With \code{intrinsic = TRUE}, the model is fitted through a first-order Taylor expansion of the growth curve
#' around the mean growth ratio or knot. When the expected change of the TVC over time has a similar shape to the
#' growth curve of the outcome, the TVC effect absorbs part of the error of this expansion: in data generated by
#' \code{getSimData()} with its default values and a negative exponential curve, the TVC effect was 5\% to 10\%
#' too high and the mean of \code{eta1} or of the ratio up to 3\% too low with \code{decompose = 1} or \code{3}
#' (one generated data set of 5,000 persons for each model).
#' A fixed ratio (\code{intrinsic = FALSE}) avoids this. With \code{y_model = "LCSM"}, the midpoint
#' approximation described in \code{\link{getLCSM}} also applies.
#'
#' @export
#'
#' @examples
#' mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' # Load the ECLS-K (2011) sample and re-baseline the measurement occasions (in years)
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' baseT <- RMS_dat0$T1
#' for (j in 1:9) RMS_dat0[[paste0("T", j)]] <- (RMS_dat0[[paste0("T", j)]] - baseT) / 12
#' # Standardized time-invariant covariates
#' RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#' RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
#' # Standardize the time-varying covariate (reading) by its baseline distribution
#' BL_mean <- mean(RMS_dat0$R1)
#' BL_sd <- sd(RMS_dat0$R1)
#' for (j in 1:9) RMS_dat0[[paste0("Rs", j)]] <- (RMS_dat0[[paste0("R", j)]] - BL_mean) / BL_sd
#'
#' \donttest{
#' # Linear latent growth curve model with the time-varying covariate in its original form
#' LIN_TVC_LGCM <- getTVCmodel(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "LIN", intrinsic = FALSE,
#'   records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 0, paramOut = TRUE
#' )
#' printTable(LIN_TVC_LGCM)
#' }
#' \dontrun{
#' # Negative exponential model with a random ratio, a decomposed time-varying covariate
#' # (interval-specific slopes) and time-invariant covariates; this takes several minutes
#' set.seed(20191029)
#' EXP_TVCslp_LGCM3.f <- getTVCmodel(
#'   dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "EXP", intrinsic = TRUE,
#'   records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 1,
#'   growth_TIC = c("ex1", "ex2"), tries = 20, paramOut = TRUE
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
  dat <- as.data.frame(dat)
  validate_paramOut(paramOut, names)
  validate_res_scale(res_scale)
  validate_res_cor(res_cor)
  validate_curveFun(curveFun)
  validate_decompose(decompose)
  validate_intrinsic(intrinsic, curveFun, model_type = "general")
  validate_curveFun_ymodel(curveFun, y_model)
  validate_columns(dat, t_var = t_var, y_var = y_var, records = records,
                   growth_TIC = growth_TIC, TVC = TVC)

  if (paramOut && !is.null(names)) {
    expected_names <- .auto_names_TVC(curveFun, intrinsic, records, y_model, decompose, growth_TIC)
    names <- .validate_output_names(names, length(expected_names))
  }
  metadata <- .capture_fit_metadata(match.call(), "getTVCmodel", environment(), dat)
  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    metadata$fitting$starts_source <- "automatic"
    starts <- getTVC.initial(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, records = records,
                             growth_TIC = growth_TIC, TVC = TVC, decompose = decompose, res_scale = res_scale,
                             res_cor = res_cor, intrinsic = intrinsic)
  }
  ## Build up a latent growth curve model or latent change score model with a decomposed TVC (and TICs if any)
  model_mx <- getTVC.mxModel(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                             records = records, growth_TIC = growth_TIC, TVC = TVC, decompose = decompose,
                             y_model = y_model, starts = starts)
  ## Optimize the constructed latent growth curve model or latent change score model with a decomposed TVC
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
  }, plan = .std_plan(dat, "TVC", curveFun, y_var, records, t_var, TVC, decompose, covariates = growth_TIC),
  build = function(plan) getTVC.mxModel(dat = .std_data(dat, plan), t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                               records = records, growth_TIC = .std_cov(plan, growth_TIC), TVC = TVC, decompose = decompose,
                               y_model = y_model, starts = starts))
  model <- fit_std$model
  metadata$fitting$standardization <- fit_std$standardization
  ## Print out the point estimates and standard errors for the parameters of interest
  if(paramOut){
    if (is.null(names)) names <- .auto_names_TVC(curveFun, intrinsic, records, y_model, decompose, growth_TIC)
    TVC_output <- getTVC.output(model = model, curveFun = curveFun, records = records, y_model = y_model,
                                decompose = decompose, growth_TIC = growth_TIC, names = names)
    model <- new("myMxOutput", mxOutput = model, Estimates = TVC_output)
  }
  else{
    model <- new("myMxOutput", mxOutput = model)
  }
  model <- .attach_fit_metadata(model, metadata, starts, frame = environment())
  model <- .warn_fit_issues(model)
  return(model)
}

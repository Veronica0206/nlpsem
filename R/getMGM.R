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
#' particular outcome variable at each study wave.
#' @param curveFun A string specifying the functional forms of the growth curve(s). Supported options for \code{y_model = "LGCM"} include:
#' \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported options for
#' \code{y_model = "LCSM"} include: \code{"nonparametric"} (or \code{"NonP"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' and \code{"Jenss-Bayley"} (or \code{"JB"}).
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. Default is
#' \code{TRUE}.
#' @param records A list of numeric vectors, with each vector specifying the indices of the observed study waves for
#' the corresponding outcome variable.
#' @param y_model A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' @param starts A list containing initial values for the parameters. Default is \code{NULL}, indicating no user-specified initial
#' values.
#' @param res_scale A numeric vector with each element representing the scaling factor for the initial calculation of the residual
#' variance. These values should be between \code{0} and \code{1}, exclusive. By default, this is \code{NULL}, as it is unnecessary
#' when the user specifies the initial values using the \code{starts} argument.
#' @param res_cor A numeric value or vector for user-specified residual correlation between any two longitudinal outcomes to calculate
#' the corresponding initial value. By default, this is \code{NULL}, as it is unnecessary when the user specifies the initial values
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
#'   \item \code{mxOutput}: This slot contains the fitted multivariate latent growth curve model or a multivariate latent change score
#'   model. A summary of this model can be obtained using the \code{ModelSummary()} function.
#'   \item \code{Estimates} (optional): If \code{paramOut = TRUE}, a data frame with parameter estimates and standard errors. The content
#'   of this slot can be printed using the \code{printTable()} method for S4 objects.
#' }
#'
#' @references
#' \itemize{
#'   \item {Liu, J., & Perera, R. A. (2021). "Estimating Knots and Their Association in Parallel Bilinear Spline Growth Curve Models in
#'   the Framework of Individual Measurement Occasions," Psychological Methods (Advance online publication).
#'   \doi{10.1037/met0000309}}
#'   \item {Blozis, S. A. (2004). "Structured Latent Curve Models for the Study of Change in Multivariate Repeated Measures," Psychological
#'   Methods, 9(3), 334-353. \doi{10.1037/1082-989X.9.3.334}}
#' }
#'
#' @export
#'
#' @examples
#' mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' # Load ECLS-K (2011) data
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
#'
#' \donttest{
#' # Fit linear multivariate latent growth curve model
#' LIN_PLGCM_f <- getMGM(
#'   dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "LIN",
#'   intrinsic = FALSE, records = list(1:9, 1:9), y_model = "LGCM", res_scale = c(0.1, 0.1),
#'   res_cor = 0.3
#'   )
#' # Fit bilinear spline multivariate latent growth curve model (random knots)
#' ## Define parameter names
#' paraBLS_PLGCM.f <- c(
#'   "Y_mueta0", "Y_mueta1", "Y_mueta2", "Y_knot",
#'   paste0("Y_psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")), "Y_res",
#'   "Z_mueta0", "Z_mueta1", "Z_mueta2", "Z_knot",
#'   paste0("Z_psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")), "Z_res",
#'   paste0("YZ_psi", c(c("00", "10", "20", "g0", "01", "11", "21", "g1",
#'                        "02", "12", "22", "g2", "0g", "1g", "2g", "gg"))),"YZ_res"
#'   )
#' ## Fit model
#' BLS_PLGCM_f <- getMGM(
#'   dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS", intrinsic = TRUE,
#'   records = list(1:9, 1:9), y_model = "LGCM", res_scale = c(0.1, 0.1), res_cor = 0.3,
#'   paramOut = TRUE, names = paraBLS_PLGCM.f
#'   )
#' printTable(BLS_PLGCM_f)
#' }
#'
#' @importFrom OpenMx mxTryHard mxRun
#' @importFrom methods new
#'
getMGM <- function(dat, t_var, y_var, curveFun, intrinsic = TRUE, records, y_model, starts = NULL, res_scale = NULL,
                   res_cor = NULL, tries = NULL, OKStatus = 0, jitterD = "runif", loc = 1, scale = 0.25,
                   paramOut = FALSE, names = NULL){
  if (I(paramOut & is.null(names))){
    stop("Please enter the original parameters if want to obtain them!")
  }
  if (I(any(res_scale <= 0) | any(res_scale >= 1))){
    stop("Please enter a value between 0 and 1 (exclusive) for res_scale!")
  }
  if (length(t_var) != length(y_var) || length(t_var) != length(records)) {
    stop("Lengths of t_var, y_var, and records must be equal!")
  }
  if (I(intrinsic & curveFun %in% c("linear", "LIN", "quadratic", "QUAD", "nonparametric", "NonP"))){
    stop("An intrinsic nonlinear function should be negative exponential, Jenss-Bayley, or bilinear spline!")
  }
  # Extract the first letter of each element in `y_var`, and generate the corresponding uppercase
  y_var <- sapply(y_var, function(x) toupper(substr(x, 1, 1)))
  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    starts <- getMULTI.initial(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, records = records,
                               res_scale = res_scale, res_cor = res_cor)
  }
  ## Build up a multivariate latent growth curve model or latent change score model
  ### Obtain manifest and latent variables, paths of the longitudinal outcomes
  model_mx <- getMGM.mxModel(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                             records = records, y_model = y_model, starts = starts)
  if (!is.null(tries)){
    model0 <- mxTryHard(model_mx, extraTries = tries, OKstatuscodes = OKStatus, jitterDistrib = jitterD,
                        loc = loc, scale = scale)
    model <- mxRun(model0)
  }
  else{
    model <- mxRun(model_mx)
  }
  if(paramOut){
    MGM_output <- getMGM.output(model = model, y_var = y_var, records = records, curveFun = curveFun, y_model = y_model,
                                names = names)
    model <- new("myMxOutput", mxOutput = model, Estimates = MGM_output)
  }
  else{
    model <- new("myMxOutput", mxOutput = model)
  }
  return(model)
}


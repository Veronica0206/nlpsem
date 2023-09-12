#' @title Fit a Latent Change Score Model with a Time-invariant Covariate (If Any)
#'
#' @description This function fits a latent change score model with or without time-invariant covariates to the provided data.
#' It manages model setup, optimization, and if requested, outputs parameter estimates and standard errors.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions, and time-invariant covariates (TICs) if any.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for latent change score
#' models include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}), \code{"Jenss-Bayley"}
#' (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}).
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
#' @return An object of class \code{myMxOutput}. Depending on the \code{paramOut} argument, the object may contain the following slots:
#' \itemize{
#'   \item \code{mxOutput}: This slot contains the fitted latent change score model. A summary of this model can be obtained using the
#'   \code{ModelSummary()} function.
#'   \item \code{Estimates} (optional): If \code{paramOut = TRUE}, a data frame with parameter estimates and standard errors. The content
#'   of this slot can be printed using the \code{printTable()} method for S4 objects.
#' }
#'
#' @references
#' \itemize{
#'   \item {Liu, J., & Perera, R. A. (2023). Estimating Rate of Change for Nonlinear Trajectories in the Framework of Individual Measurement
#'   Occasions: A New Perspective on Growth Curves. Behavior Research Methods. \doi{10.3758/s13428-023-02097-2}}
#'   \item {Liu, J. (2022). "Jenss–Bayley Latent Change Score Model With Individual Ratio of the Growth Acceleration in the Framework
#'   of Individual Measurement Occasions." Journal of Educational and Behavioral Statistics, 47(5), 507–543.
#'   \doi{10.3102/10769986221099919}}
#'   \item {Grimm, K. J., Zhang, Z., Hamagami, F., & Mazzocco, M. (2013). "Modeling Nonlinear Change via Latent Change and Latent
#'   Acceleration Frameworks: Examining Velocity and Acceleration of Growth Trajectories." Multivariate Behavioral Research, 48(1),
#'   117-143. \doi{10.1080/00273171.2012.755111}}
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
#' RMS_dat0$T1 <- (RMS_dat0$T1 - baseT)/12
#' RMS_dat0$T2 <- (RMS_dat0$T2 - baseT)/12
#' RMS_dat0$T3 <- (RMS_dat0$T3 - baseT)/12
#' RMS_dat0$T4 <- (RMS_dat0$T4 - baseT)/12
#' RMS_dat0$T5 <- (RMS_dat0$T5 - baseT)/12
#' RMS_dat0$T6 <- (RMS_dat0$T6 - baseT)/12
#' RMS_dat0$T7 <- (RMS_dat0$T7 - baseT)/12
#' RMS_dat0$T8 <- (RMS_dat0$T8 - baseT)/12
#' RMS_dat0$T9 <- (RMS_dat0$T9 - baseT)/12
#' # Standardized time-invariant covariates
#' RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#' RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
#'
#' \donttest{
#' # Fit nonparametric change score model for reading development
#' ## Fit model
#' NonP_LCSM <- getLCSM(
#'   dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "nonparametric",
#'   intrinsic = FALSE, records = 1:9, res_scale = 0.1
#'   )
#' }
#'
#' @importFrom OpenMx mxTryHard mxRun
#' @importFrom methods new
#'
getLCSM <- function(dat, t_var, y_var, curveFun, intrinsic = TRUE, records, growth_TIC = NULL, starts = NULL,
                    res_scale = NULL, tries = NULL, OKStatus = 0, jitterD = "runif", loc = 1, scale = 0.25,
                    paramOut = FALSE, names = NULL){
  if (I(paramOut & is.null(names))){
    stop("Please enter the original parameters if want to obtain them!")
  }
  if (I(res_scale <= 0 | res_scale >= 1)){
    stop("Please enter a value between 0 and 1 (exclusive) for res_scale!")
  }
  if (I(intrinsic & curveFun %in% c("nonparametric", "NonP", "quadratic", "QUAD"))){
    stop("An intrinsic nonlinear function should be negative exponential or Jenss-Bayley for a LCSM!")
  }
  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    starts <- getUNI.initial(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, records = records,
                             growth_TIC = growth_TIC, res_scale = res_scale)
  }
  ## Build up a latent change score model, with or without TICs
  model_mx <- getLCSM.mxModel(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                              records = records, growth_TIC = growth_TIC, starts = starts)
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
    LCSM_output <- getLCSM.output(model = model, curveFun = curveFun, growth_TIC = growth_TIC, names = names)
    model <- new("myMxOutput", mxOutput = model, Estimates = LCSM_output)
  }
  else{
    model <- new("myMxOutput", mxOutput = model)
  }
  return(model)
}

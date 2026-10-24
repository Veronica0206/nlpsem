#' @title Fit a Longitudinal Multiple Group Model
#'
#' @description This function fits a longitudinal multiple group model based on the specified sub-model. Supported submodels
#' include:
#' \itemize{
#'   \item Latent growth curve models,
#'   \item Latent change score models,
#'   \item Latent growth curve models or latent change score models with a time-varying covariate,
#'   \item Multivariate latent growth curve models or multivariate latent change score models,
#'   \item Longitudinal mediation models.
#' } For the first three submodels, time-invariant covariates are allowed.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions for each longitudinal process, time-invariant covariates (TICs) if any, and a variable
#' that indicates manifested group.
#' @param grp_var A string specifying the column that indicates manifested classes. Any coding (numeric, character,
#' or factor) is accepted: the distinct non-missing values are sorted and relabelled \code{1, 2, ...} in that order for
#' estimation and in the output, and a message reports the mapping when the original codes differ. Rows with a missing
#' group value are removed with a warning.
#' @param sub_Model A string that specifies the sub-model for manifested classes. Supported sub-models include \code{"LGCM"}
#' (for latent growth curve models), \code{"LCSM"} (for latent change score models), \code{"TVC"} (for latent growth curve
#' models or latent change score models with a time-varying covariate), \code{"MGM"} (for multivariate latent growth curve
#' models or latent change score models), and \code{"MED"} (for longitudinal mediation models).
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
#' \code{NULL}, as it is unnecessary when \code{sub_Model} is \code{"MED"}.
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
#' \code{NULL}, in which case the residual correlations of \code{"TVC"} and \code{"MGM"} models are estimated from the
#' per-person curve residuals (with \code{0.3} as fallback) and those of \code{"MED"} models start at \code{0.3}.
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
#'   \item \code{mxOutput}: This slot contains the fitted longitudinal multiple group model. A summary of this model can be obtained using
#'   the \code{ModelSummary()} function.
#'   \item \code{Estimates} (optional): If \code{paramOut = TRUE}, a data frame with parameter estimates and standard errors. The content
#'   of this slot can be printed using the \code{printTable()} method for S4 objects.
#' }
#'
#' @references
#' \itemize{
#'   \item {Liu, J., & Perera, R. A. (2023). Estimating Rate of Change for Nonlinear Trajectories in the Framework of Individual Measurement
#'   Occasions: A New Perspective on Growth Curves. Behavior Research Methods, 56(3), 1349-1375. \doi{10.3758/s13428-023-02097-2}}
#'   }
#'
#' @seealso \code{\link{getLGCM}}, \code{\link{getLCSM}}, \code{\link{getMIX}}, \code{\link{getFigure}}
#'
#' @details New fits retain their resolved model specification and fitting provenance.
#' Use \code{\link{getSpecification}} to inspect it; standard \code{summary()}, \code{coef()},
#' \code{logLik()} and \code{plot()} methods are available.
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
#' # Multiple-group bilinear spline model with a fixed knot, groups defined by sex
#' MGroup_BLS_LGCM_r <- getMGroup(
#'   dat = RMS_dat0, grp_var = "SEX", sub_Model = "LGCM", y_var = "M", t_var = "T",
#'   records = 1:9, curveFun = "BLS", intrinsic = FALSE, paramOut = TRUE
#' )
#' printTable(MGroup_BLS_LGCM_r)
#' }
#' \dontrun{
#' # Random knot with time-invariant covariates and 10 extra optimization attempts;
#' # this takes several minutes
#' set.seed(20191029)
#' MGroup_BLS_LGCM.TIC_f <- getMGroup(
#'   dat = RMS_dat0, grp_var = "SEX", sub_Model = "LGCM", y_var = "M", t_var = "T",
#'   records = 1:9, curveFun = "BLS", intrinsic = TRUE,
#'   growth_TIC = c("ex1", "ex2"), tries = 10, paramOut = TRUE
#' )
#' printTable(MGroup_BLS_LGCM.TIC_f)
#' }
#'
#' @importFrom OpenMx mxTryHard mxRun
#' @importFrom methods new
#'
getMGroup <- function(dat, grp_var, sub_Model, t_var, records, y_var, curveFun, intrinsic = NULL, y_model = NULL,
                      m_var = NULL, x_type = NULL, x_var = NULL, TVC = NULL, decompose = NULL, growth_TIC = NULL,
                      starts = NULL, res_scale = NULL, res_cor = NULL, tries = NULL, OKStatus = 0, jitterD = "runif",
                      loc = 1, scale = 0.25, paramOut = FALSE, names = NULL){
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
                   x_var = x_var, x_type = x_type, grp_var = grp_var)
  metadata <- .capture_fit_metadata(match.call(), "getMGroup", environment(), dat)
  ## Downstream helpers index the groups 1..K: drop rows without a group and recode the grouping
  ## variable by the sorted order of its distinct values (a message reports the mapping).
  grp_recode <- .recode_group_codes(dat[[grp_var]], grp_var = grp_var)
  dat <- dat[!is.na(grp_recode$codes), , drop = FALSE]
  dat[[grp_var]] <- grp_recode$codes[!is.na(grp_recode$codes)]
  nClass <- grp_recode$nClass
  ## Derive initial values for the parameters of interest if not specified by users
  if (is.null(starts)){
    metadata$fitting$starts_source <- "automatic"
    starts <- getMGroup.initial(dat = dat, grp_var = grp_var, nClass = nClass, sub_Model = sub_Model, t_var = t_var,
                                y_var = y_var, curveFun = curveFun, intrinsic = intrinsic, records = records, m_var = m_var, x_var = x_var,
                                x_type = x_type, TVC = TVC, decompose = decompose, growth_TIC = growth_TIC,
                                res_scale = res_scale, res_cor = res_cor)
  }
  ## Build up a multiple group model
  model_mx <- getMGroup.mxModel(dat = dat, grp_var = grp_var, nClass = nClass, sub_Model = sub_Model, t_var = t_var, y_var = y_var,
                                curveFun = curveFun, intrinsic = intrinsic, records = records, y_model = y_model,
                                m_var = m_var, x_var = x_var, x_type = x_type, TVC = TVC, decompose = decompose,
                                growth_TIC = growth_TIC, starts = starts, res_cor = res_cor)
  ## Optimize the constructed multiple group model
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
    MGroup_output <- getMGroup.output(model = model, nClass = nClass, sub_Model = sub_Model, y_var = y_var,
                                      curveFun = curveFun, x_type = x_type, records = records, growth_TIC = growth_TIC,
                                      y_model = y_model, decompose = decompose, names = names)
    model <- new("myMxOutput", mxOutput = model, Estimates = MGroup_output)
  }
  else{
    model <- new("myMxOutput", mxOutput = model)
  }
  model <- .attach_fit_metadata(model, metadata, starts, nClass = nClass, group_recode = grp_recode, frame = environment())
  return(model)
}


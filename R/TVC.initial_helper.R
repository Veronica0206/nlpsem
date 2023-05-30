#' @title Compute Initial Values for Parameters of Latent Growth Curve Models or Latent Change Score Models with a Time-varying
#' Covariate and Time-invariant Covariates (if any)
#'
#' @description This function computes the initial values of the parameters for a latent growth curve model or a latent change
#' score model with a time-varying covariate and time-invariant covariates (if any).
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements (for the longitudinal outcome and time-varying covariates), occasions, and time-invariant covariates
#' (TICs) if any. It takes the value passed from \code{getTVCmodel()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave. It takes
#' the value passed from \code{getTVCmodel()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave. It
#' takes the value passed from \code{getTVCmodel()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the value
#' passed from \code{getTVCmodel()}.
#' @param records A numeric vector specifying the indices of the observed study waves. It takes the value passed from \code{getTVCmodel()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) that account for the
#' variability of growth factors, if any. It takes the value passed from \code{getTVCmodel()}.
#' @param TVC A string specifying the prefix of the column names corresponding to the time-varying covariate at each study wave.
#' It takes the value passed from \code{getTVCmodel()}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). It takes the value
#' passed from \code{getTVCmodel()}.
#' @param res_scale A numeric value or numeric vector. For a model with \code{decompose = 0}, it is a numeric value representing
#' the scaling factor used to calculate the initial value for the residual variance of the longitudinal outcome. In cases where
#' \code{decompose != 0}, it is a numeric vector of user-specified scaling factors used to calculate the initial values for the
#' residual variance of both the longitudinal outcome and the time-varying covariate. It takes the value passed from \code{getTVCmodel()}.
#' @param res_cor A numeric value. When \code{decompose != 0}, this represents the user-specified residual correlation between the
#' longitudinal outcome and the time-varying covariate, which is used to calculate the corresponding initial value. If \code{decompose = 0},
#' this should be \code{NULL}. It takes the value passed from \code{getTVCmodel()}.
#'
#' @return A list containing the initial values for parameters related to growth factors, TVC, TICs (if any), and path
#' coefficients (if any) for a latent growth curve model or a latent change score model with a time-varying covariate and
#' time-invariant covariates (if any).
#'
#' @importFrom stats var cov lm na.exclude
#'
getTVC.initial <- function(dat, t_var, y_var, curveFun, records, growth_TIC, TVC, decompose, res_scale, res_cor){
  # Calculate the number of time points
  nT <- length(records)
  # Extract the measurement values and occasions for longitudinal outcome and TVC
  dat_traj <- dat[, paste0(y_var, records)]
  dat_TVC <- dat[, paste0(TVC, records)]
  dat_time <- dat[, paste0(t_var, records)]
  # Initial values of parameters related to baseline TVC and TICs if any
  if (!is.null(growth_TIC)){
    # For a decomposed TVC
    if (decompose != 0){
      dat_covariate <- dat[, c(growth_TIC, paste0(TVC, records[1]))]
      BL_mean <- apply(dat_covariate, 2, mean, na.rm = TRUE)
      BL_var <- var(dat_covariate, na.rm = TRUE)
      BL_starts <- list(BL_mean, BL_var)
    }
    # For a TVC in its original form
    else{
      dat_covariate <- dat[, growth_TIC]
      BL_mean <- if (length(growth_TIC) == 1){
        mean(dat_covariate, na.rm = TRUE)
      } else{
        apply(dat_covariate, 2, mean, na.rm = TRUE)
      }
      BL_var <- var(dat_covariate, na.rm = TRUE)
      BL_starts <- list(BL_mean, BL_var[lower.tri(BL_var, diag = TRUE)])
    }
  }
  else if (is.null(growth_TIC)){
    if (decompose != 0){
      dat_covariate <- dat[, paste0(TVC, records[1])]
      BL_mean <- mean(dat_covariate, na.rm = TRUE)
      BL_var <- var(dat_covariate, na.rm = TRUE)
      BL_starts <- list(BL_mean, BL_var)
    }
  }
  # Calculate initial values for TVC-related parameters
  if (decompose == 0){
    TVC_starts <- list(apply(dat_TVC, 2, mean), apply(dat_TVC, 2, var))
  }
  else if (decompose != 0){
    TVC_starts <- getUNI.initial(dat = dat, y_var = TVC, t_var = "T", records = records, growth_TIC = NULL,
                                 curveFun = "NonP", res_scale = res_scale[2])[[1]]
  }
  # Calculate initial values for growth factor-related parameters
  growth_factor <- getUNI.GF(dat_traj = dat_traj, dat_time = dat_time, nT = nT, curveFun = curveFun)
  # Calculate the number of growth factors of the specified functional form
  nGF <- ncol(growth_factor)
  # Calculate the initial value of residual variance
  Y_residuals <- var(growth_factor[, 1]) * res_scale[1]
  # If `dat_covariate` exists (i.e., `growth_TIC` is provided or the TVC is added in the model in its decomposed form),
  # calculate initial values for alpha0 (growth factor intercepts), psi_r (unexplained var-cov matrix of growth
  # factors), and beta0 (path coefficients)
  if (exists("dat_covariate")){
    Y_alpha0 <- sapply(1:nGF, function(p){
      as.numeric(lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p], dat_covariate)),
                    na.action = na.exclude)$coefficients[1])
    })
    beta0 <- matrix(sapply(1:nGF, function(p){
      as.numeric(lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p], dat_covariate)),
                    na.action = na.exclude)$coefficients[-1])
    }), nrow = nGF, byrow = TRUE)
    Y_psi_r <- cov(growth_factor) - beta0 %*% BL_var %*% t(beta0)
    Y_starts <- list(Y_alpha0, Y_psi_r, Y_residuals)
  }
  # If `dat_covariate` dose not exists (i.e., `growth_TIC` is not provided and TVC is added in the model in its original
  # form), calculate initial values for mean0 (growth factor means) and psi0 (var-cov matrix of growth factors)
  else{
    Y_mean0 <- sapply(1:nGF, function(p){
      as.numeric(lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p])),
                    na.action = na.exclude)$coefficients[1])
    })
    Y_psi0 <- cov(growth_factor)
    Y_starts <- list(Y_mean0, Y_psi0, Y_residuals)
  }
  if (curveFun %in% c("nonparametric", "NonP")){
    abs_rate <- rep(0, nT - 1)
    for (j in 1:(nT - 1)){
      abs_rate[j] <- as.numeric(lm(traj_delta ~ .,
                                   data = data.frame(traj_delta = dat[, paste0(y_var, records[j + 1])] - dat[, paste0(y_var, records[j])],
                                                     time_delta = dat[, paste0(t_var, records[j + 1])] - dat[, paste0(t_var, records[j])]),
                                   na.action = na.exclude)$coefficients[-1])
    }
    Y_starts <- c(Y_starts, abs_rate/Y_starts[[1]][2])
  }
  # Calculate initial value for state effect
  if (decompose == 0){
    kappa <- as.numeric(lm(Y_unexp ~ .,
                           data = data.frame(Y_unexp = dat_traj[, 1] - Y_starts[[1]][1], TVC = dat_TVC[, 1]),
                           na.action = na.exclude)$coefficients[-1])
  }
  else if (decompose == 1){
    kappa <- as.numeric(lm(Y_unexp ~ .,
                           data = data.frame(Y_unexp = dat_traj[, 2] - dat_traj[, 1] -
                                               TVC_starts[[1]][2] * dat_time[, 2],
                                             TVC_slp = (dat_TVC[, 2] - dat_TVC[, 1])/(dat_time[, 2] - dat_time[, 1])),
                           na.action = na.exclude)$coefficients[-1])
  }
  else if (I(decompose == 2 | decompose == 3)){
    kappa <- as.numeric(lm(Y_unexp ~ .,
                           data = data.frame(Y_unexp = dat_traj[, 2] - dat_traj[, 1] -
                                               TVC_starts[[1]][2] * dat_time[, 2],
                                             TVC_chg = dat_TVC[, 2] - dat_TVC[, 1]),
                           na.action = na.exclude)$coefficients[-1])
  }
  # Define list `starts` to include initial values for all parameters
  if (decompose == 0){
    starts <- if (exists("dat_covariate")){
      list(Y_starts, TVC_starts, BL_starts, beta0, kappa)
    }
    else{
      list(Y_starts, TVC_starts, kappa)
    }
  }
  else if (decompose != 0){
    # Calculate initial value of residual covariance between TVC and longitudinal outcome if TVC is added in the model
    # in a decomposed form
    XY_cov <- res_cor * sqrt(Y_residuals * TVC_starts[[3]])
    starts <- list(Y_starts, TVC_starts, BL_starts, beta0, kappa, XY_cov)
  }
  if (curveFun %in% c("bilinear spline", "BLS")){
    # Transformed matrices obtained by multivariate Delta method
    ## For mean vector
    func0 <- matrix(c(1, starts[[1]][[1]][4], 0, 0,
                      0, 0.5, 0.5, 0,
                      0, -0.5, 0.5, 0,
                      0, 0, 0, 1), nrow = 4, byrow = T)
    ## For var-cov matrix
    grad0 <- matrix(c(1, starts[[1]][[1]][4], 0, starts[[1]][[1]][2],
                      0, 0.5, 0.5, 0,
                      0, -0.5, 0.5, 0,
                      0, 0, 0, 1), nrow = 4, byrow = T)
    starts[[1]][[1]] <- func0 %*% starts[[1]][[1]]
    starts[[1]][[2]] <- grad0 %*% starts[[1]][[2]] %*% t(grad0)
    if (!is.null(growth_TIC)){
      starts[[4]] <- grad0 %*% starts[[4]]
    }
  }
  starts[[1]][[2]] <- starts[[1]][[2]][lower.tri(starts[[1]][[2]], diag = TRUE)]
  return(starts)
}

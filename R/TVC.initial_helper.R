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
#' @param res_scale An optional numeric value or numeric vector. For a model with \code{decompose = 0}, it is a numeric value
#' representing the scaling factor used to calculate the initial value for the residual variance of the longitudinal outcome. In
#' cases where \code{decompose != 0}, it is a numeric vector of scaling factors used to calculate the initial values for the
#' residual variance of both the longitudinal outcome and the time-varying covariate. Default is \code{NULL}, in which case
#' data-driven estimation is used with a heuristic of \code{0.1} as fallback. It takes the value passed from \code{getTVCmodel()}.
#' @param res_cor An optional numeric value. When \code{decompose != 0}, this represents the user-specified residual correlation
#' between the longitudinal outcome and the time-varying covariate, which is used to calculate the corresponding initial value. If
#' \code{decompose = 0}, this should be \code{NULL}. When \code{NULL} and \code{decompose != 0}, data-driven estimation is used.
#' It takes the value passed from \code{getTVCmodel()}.
#'
#' @return A list containing the initial values for parameters related to growth factors, TVC, TICs (if any), and path
#' coefficients (if any) for a latent growth curve model or a latent change score model with a time-varying covariate and
#' time-invariant covariates (if any).
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats var cov cor lm na.exclude
#'
.stabilize_kappa_start <- function(kappa, intrinsic, decompose){
  if (is.null(kappa) || !isTRUE(intrinsic) || decompose == 0) return(kappa)
  sign(kappa) * min(abs(kappa) * 0.5, 10)
}

.predict_curve_from_gf <- function(gf, time, curveFun){
  if (any(!is.finite(gf)) || !is.finite(time)) return(NA_real_)

  if (curveFun %in% c("linear", "LIN")){
    gf[1] + gf[2] * time
  } else if (curveFun %in% c("quadratic", "QUAD")){
    gf[1] + gf[2] * time + gf[3] * time^2
  } else if (curveFun %in% c("negative exponential", "EXP")){
    gf[1] + gf[2] * (1 - exp(-gf[3] * time))
  } else if (curveFun %in% c("Jenss-Bayley", "JB")){
    gf[1] + gf[2] * time + gf[3] * (exp(gf[4] * time) - 1)
  } else if (curveFun %in% c("bilinear spline", "BLS")){
    gf[1] + gf[2] * pmin(time, gf[4]) + gf[3] * pmax(time - gf[4], 0)
  } else{
    NA_real_
  }
}

.predict_curve_matrix <- function(growth_factor, dat_time, curveFun, rel_rate = NULL){
  n_person <- nrow(dat_time)
  n_wave <- ncol(dat_time)
  fitted <- matrix(NA_real_, nrow = n_person, ncol = n_wave)
  if (!n_person || !n_wave){
    return(fitted)
  }

  if (curveFun %in% c("nonparametric", "NonP")){
    if (is.null(rel_rate) || length(rel_rate) < (n_wave - 1)){
      return(fitted)
    }
    fitted[, 1] <- growth_factor[, 1]
    for (j in 2:n_wave){
      time_delta <- dat_time[, j] - dat_time[, j - 1]
      fitted[, j] <- fitted[, j - 1] + growth_factor[, 2] * rel_rate[j - 1] * time_delta
    }
    return(fitted)
  }

  for (i in seq_len(n_person)){
    for (j in seq_len(n_wave)){
      fitted[i, j] <- .predict_curve_from_gf(as.numeric(growth_factor[i, ]), dat_time[i, j], curveFun)
    }
  }
  fitted
}

.estimate_decompose3_kappa <- function(dat_traj, dat_TVC, dat_time, growth_factor, curveFun,
                                       outcome_rel_rate = NULL, tvc_rel_rate = NULL){
  y_fit <- .predict_curve_matrix(
    growth_factor = growth_factor,
    dat_time = dat_time,
    curveFun = curveFun,
    rel_rate = outcome_rel_rate
  )
  tvc_growth_factor <- getUNI.GF(dat_traj = dat_TVC, dat_time = dat_time, nT = ncol(dat_time), curveFun = "NonP")
  tvc_fit <- .predict_curve_matrix(
    growth_factor = as.matrix(tvc_growth_factor),
    dat_time = dat_time,
    curveFun = "NonP",
    rel_rate = tvc_rel_rate
  )

  if (all(!is.finite(y_fit[, -1, drop = FALSE]))){
    y_unexp <- c(sweep(as.matrix(dat_traj[, -1, drop = FALSE]), 1, dat_traj[, 1], "-"))
  } else{
    y_unexp <- c(as.matrix(dat_traj[, -1, drop = FALSE] - y_fit[, -1, drop = FALSE]))
  }

  if (all(!is.finite(tvc_fit[, -1, drop = FALSE]))){
    tvc_state <- c(sweep(as.matrix(dat_TVC[, -1, drop = FALSE]), 1, dat_TVC[, 1], "-"))
  } else{
    tvc_state <- c(sweep(as.matrix(tvc_fit[, -1, drop = FALSE]), 1, tvc_fit[, 1], "-"))
  }

  kappa <- suppressWarnings(as.numeric(lm(y_unexp ~ tvc_state, na.action = na.exclude)$coefficients[-1]))
  if (!is.finite(kappa)) {
    kappa <- suppressWarnings(as.numeric(lm((dat_traj[, 2] - dat_traj[, 1]) ~ I(dat_TVC[, 2] - dat_TVC[, 1]),
                                           na.action = na.exclude)$coefficients[-1]))
  }
  kappa
}

getTVC.initial <- function(dat, t_var, y_var, curveFun, records, growth_TIC, TVC, decompose, res_scale, res_cor = NULL,
                           intrinsic = FALSE){
  # Calculate the number of time points
  nT <- length(records)
  # Extract the measurement values and occasions for longitudinal outcome and TVC
  dat_traj <- dat[, paste0(y_var, records)]
  dat_TVC <- dat[, paste0(TVC, records)]
  dat_time <- dat[, paste0(t_var, records)]
  # Initialize dat_covariate to NULL; set below if growth_TIC is provided or decompose != 0
  dat_covariate <- NULL
  BL_starts <- NULL
  # Initial values of parameters related to baseline TVC and TICs if any
  if (!is.null(growth_TIC)){
    # For a decomposed TVC
    if (decompose != 0){
      dat_covariate <- dat[, c(growth_TIC, paste0(TVC, records[1]))]
      BL_mean <- apply(dat_covariate, 2, mean, na.rm = TRUE)
      BL_var <- var(dat_covariate, use = "pairwise.complete.obs")
      BL_starts <- list(means = BL_mean, covMatrix = BL_var)
    }
    # For a TVC in its original form
    else{
      dat_covariate <- dat[, growth_TIC]
      BL_mean <- if (length(growth_TIC) == 1){
        mean(dat_covariate, na.rm = TRUE)
      } else{
        apply(dat_covariate, 2, mean, na.rm = TRUE)
      }
      BL_var <- var(dat_covariate, use = "pairwise.complete.obs")
      BL_starts <- list(means = BL_mean, covMatrix = BL_var[lower.tri(BL_var, diag = TRUE)])
    }
  }
  else if (is.null(growth_TIC)){
    if (decompose != 0){
      dat_covariate <- dat[, paste0(TVC, records[1])]
      BL_mean <- mean(dat_covariate, na.rm = TRUE)
      BL_var <- var(dat_covariate, use = "pairwise.complete.obs")
      BL_starts <- list(means = BL_mean, covMatrix = BL_var)
    }
  }
  # Calculate initial values for TVC-related parameters
  if (decompose == 0){
    TVC_starts <- list(means = apply(dat_TVC, 2, mean, na.rm = TRUE), covMatrix = apply(dat_TVC, 2, var, na.rm = TRUE))
  }
  else if (decompose != 0){
    TVC_starts <- getUNI.initial(dat = dat, y_var = TVC, t_var = t_var, records = records, growth_TIC = NULL,
                                 curveFun = "NonP", res_scale = if (!is.null(res_scale)) res_scale[2] else NULL,
                                 intrinsic = FALSE)$Y_starts
  }
  # Calculate initial values for growth factor-related parameters
  growth_factor <- getUNI.GF(dat_traj = dat_traj, dat_time = dat_time, nT = nT, curveFun = curveFun)
  # Calculate the number of growth factors of the specified functional form
  nGF <- ncol(growth_factor)
  gf_cov_total <- cov(growth_factor, use = "pairwise.complete.obs")
  # Calculate the initial value of Y residual variance.
  # Priority: (1) user-specified res_scale  (2) data-driven mean_res_var  (3) 0.1 heuristic
  if (!is.null(res_scale)){
    Y_residuals <- var(growth_factor[, 1], na.rm = TRUE) * res_scale[1]
  } else {
    data_res_var <- attr(growth_factor, "mean_res_var")
    if (!is.null(data_res_var) && !is.na(data_res_var) && data_res_var > 0){
      Y_residuals <- data_res_var
    } else{
      warning("Data-driven residual variance estimation produced invalid result for outcome ", y_var,
              "; falling back to heuristic (0.1 * var(eta0)).", call. = FALSE)
      Y_residuals <- var(growth_factor[, 1], na.rm = TRUE) * 0.1
    }
  }
  # If `dat_covariate` exists (i.e., `growth_TIC` is provided or the TVC is added in the model in its decomposed form),
  # calculate initial values for alpha0 (growth factor intercepts), psi_r (unexplained var-cov matrix of growth
  # factors), and beta0 (path coefficients)
  if (!is.null(dat_covariate)){
    gf_means <- sapply(1:nGF, function(p){
      as.numeric(lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p], dat_covariate)),
                    na.action = na.exclude)$coefficients[1])
    })
    gf_beta <- matrix(sapply(1:nGF, function(p){
      as.numeric(lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p], dat_covariate)),
                    na.action = na.exclude)$coefficients[-1])
    }), nrow = nGF, byrow = TRUE)
  }
  # If `dat_covariate` does not exist (i.e., `growth_TIC` is not provided and TVC is added in the model in its original
  # form), calculate initial values for mean0 (growth factor means) and psi0 (var-cov matrix of growth factors)
  else{
    gf_means <- sapply(1:nGF, function(p){
      as.numeric(lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p])),
                    na.action = na.exclude)$coefficients[1])
    })
  }
  rel_rate <- NULL
  if (curveFun %in% c("nonparametric", "NonP")){
    abs_rate <- rep(0, nT - 1)
    for (j in 1:(nT - 1)){
      traj_delta <- dat[, paste0(y_var, records[j + 1])] - dat[, paste0(y_var, records[j])]
      time_delta <- dat[, paste0(t_var, records[j + 1])] - dat[, paste0(t_var, records[j])]
      # NonP-6: If time intervals have near-zero variance (fixed-occasion design),
      # skip regression and use mean(abs_rate) directly as the slope
      if (var(time_delta, na.rm = TRUE) < .Machine$double.eps){
        abs_rate[j] <- mean(traj_delta / time_delta, na.rm = TRUE)
      } else{
        abs_rate[j] <- as.numeric(lm(traj_delta ~ .,
                                     data = data.frame(traj_delta = traj_delta,
                                                       time_delta = time_delta),
                                     na.action = na.exclude)$coefficients[-1])
      }
    }
    # NonP-3: Safeguard against division by near-zero mean slope (eta1)
    mean_eta1 <- gf_means[2]
    if (abs(mean_eta1) < 1e-10){
      warning("Mean slope (eta1) is near zero; using abs_rate directly as rel_rate.", call. = FALSE)
      rel_rate <- abs_rate
    } else{
      rel_rate <- abs_rate / mean_eta1
    }
  }
  # Calculate initial value for state effect
  if (decompose == 0){
    kappa <- as.numeric(lm(Y_unexp ~ .,
                           data = data.frame(Y_unexp = dat_traj[, 1] - gf_means[1], TVC = dat_TVC[, 1]),
                           na.action = na.exclude)$coefficients[-1])
  }
  else if (decompose == 1){
    kappa <- as.numeric(lm(Y_unexp ~ .,
                           data = data.frame(Y_unexp = dat_traj[, 2] - dat_traj[, 1] -
                                               TVC_starts$means[2] * dat_time[, 2],
                                             TVC_slp = (dat_TVC[, 2] - dat_TVC[, 1])/(dat_time[, 2] - dat_time[, 1])),
                           na.action = na.exclude)$coefficients[-1])
  }
  else if (decompose == 2){
    kappa <- as.numeric(lm(Y_unexp ~ .,
                           data = data.frame(Y_unexp = dat_traj[, 2] - dat_traj[, 1] -
                                               TVC_starts$means[2] * dat_time[, 2],
                                             TVC_chg = dat_TVC[, 2] - dat_TVC[, 1]),
                           na.action = na.exclude)$coefficients[-1])
  } else if (decompose == 3){
    kappa <- .estimate_decompose3_kappa(
      dat_traj = as.matrix(dat_traj),
      dat_TVC = as.matrix(dat_TVC),
      dat_time = as.matrix(dat_time),
      growth_factor = as.matrix(growth_factor),
      curveFun = curveFun,
      outcome_rel_rate = rel_rate,
      tvc_rel_rate = TVC_starts$rel_rate
    )
  }
  kappa <- .stabilize_kappa_start(kappa, intrinsic = intrinsic, decompose = decompose)
  if (curveFun %in% c("bilinear spline", "BLS")){
    # Transformed matrices obtained by multivariate Delta method
    ## For mean vector
    func0 <- matrix(c(1, gf_means[4], 0, 0,
                      0, 0.5, 0.5, 0,
                      0, -0.5, 0.5, 0,
                      0, 0, 0, 1), nrow = 4, byrow = TRUE)
    ## For var-cov matrix
    grad0 <- matrix(c(1, gf_means[4], 0, gf_means[2],
                      0, 0.5, 0.5, 0,
                      0, -0.5, 0.5, 0,
                      0, 0, 0, 1), nrow = 4, byrow = TRUE)
    gf_means <- drop(func0 %*% gf_means)
    gf_cov_total <- grad0 %*% gf_cov_total %*% t(grad0)
    if (!is.null(dat_covariate)){
      gf_beta <- grad0 %*% gf_beta
    }
  }
  if (!is.null(dat_covariate)){
    baseline_tvc_col <- if (decompose != 0) ncol(gf_beta) else integer()
    gf_beta <- .stabilize_covariate_effect_starts(
      gf_beta, curveFun = curveFun, intrinsic = intrinsic, baseline_tvc_col = baseline_tvc_col
    )
    Y_cov <- gf_cov_total - gf_beta %*% BL_var %*% t(gf_beta)
    Y_cov <- .ensure_pd(Y_cov, label = "Y_psi_r")
    Y_starts <- list(means = gf_means, covMatrix = Y_cov, residuals = Y_residuals)
  } else{
    Y_starts <- list(means = gf_means, covMatrix = gf_cov_total, residuals = Y_residuals)
  }
  if (!is.null(rel_rate)){
    Y_starts$rel_rate <- rel_rate
  }
  # Define list `starts` to include initial values for all parameters
  if (decompose == 0){
    starts <- if (!is.null(dat_covariate)){
      list(Y_starts = Y_starts, TVC_starts = TVC_starts, BL_starts = BL_starts, beta = gf_beta, kappa = kappa)
    } else{
      list(Y_starts = Y_starts, TVC_starts = TVC_starts, kappa = kappa)
    }
  } else if (decompose != 0){
    # Calculate initial value of residual covariance between TVC and longitudinal outcome if TVC is added in the model
    # in a decomposed form. If res_cor not provided, estimate from data using within-person residuals.
    if (is.null(res_cor)){
      dat_y <- dat[, paste0(y_var, records), drop = FALSE]
      dat_x <- dat[, paste0(TVC, records), drop = FALSE]
      resid_y <- dat_y - rowMeans(dat_y, na.rm = TRUE)
      resid_x <- dat_x - rowMeans(dat_x, na.rm = TRUE)
      res_cor <- cor(unlist(resid_y), unlist(resid_x), use = "complete.obs")
      if (is.na(res_cor)) {
        warning("Data-driven residual correlation estimation failed for outcome ", y_var,
                " and TVC ", TVC, "; falling back to heuristic (res_cor = 0.3).", call. = FALSE)
        res_cor <- 0.3
      }
      res_cor <- max(min(res_cor, 0.99), -0.99)
    }
    XY_cov <- res_cor * sqrt(Y_residuals * TVC_starts$residuals)
    starts <- list(Y_starts = Y_starts, TVC_starts = TVC_starts, BL_starts = BL_starts, beta = gf_beta, kappa = kappa, XY_cov = XY_cov)
  }
  starts <- .stabilize_intrinsic_nonlinear_starts(
    starts = starts, curveFun = curveFun, intrinsic = intrinsic,
    time_values = as.numeric(as.matrix(dat_time))
  )
  starts$Y_starts$covMatrix <- starts$Y_starts$covMatrix[lower.tri(starts$Y_starts$covMatrix, diag = TRUE)]
  return(starts)
}

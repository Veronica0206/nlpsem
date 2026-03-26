#' @title Compute Initial Values for Parameters of Latent Growth Curve Models or Latent Change Score Models with Time-invariant
#' Covariates (If Any)
#'
#' @description This function computes the initial values of the parameters for a latent growth curve model or a latent change
#' score model with time-invariant covariates (if any).
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions, and time-invariant covariates (TICs) if any. It takes the value passed from \code{getLGCM()}
#' or \code{getLCSM()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' It takes the value passed from \code{getLGCM()} or \code{getLCSM()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getLGCM()} or \code{getLCSM()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for latent growth curve models include:
#' \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported options for
#' latent change score models include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the value passed from \code{getLGCM()}
#' or \code{getLCSM()}.
#' @param records A numeric vector specifying indices of the study waves. It takes the value passed from \code{getLGCM()} or \code{getLCSM()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing to the
#' variability of growth factors if any. It takes the value passed from \code{getLGCM()} or \code{getLCSM()}.
#' @param res_scale An optional numeric value representing the scaling factor for the initial calculation of the residual
#' variance. This value should be between \code{0} and \code{1}, exclusive, or \code{NULL} (default) to use data-driven estimation
#' with a heuristic of \code{0.1} as fallback. It takes the value passed from \code{getLGCM()} or \code{getLCSM()}.
#'
#' @return A list containing the initial values for parameters related to growth factors, TICs (if any), and path
#' coefficients (if any) of a latent growth curve model or a latent change score model. The returned list has the
#' following structure:
#' \describe{
#'   \item{Y_starts:}{A list containing three elements:
#'   \describe{
#'     \item{alpha0 or mean0:}{Depends on whether \code{growth_TIC} is provided,}
#'     \item{psi_r or psi0:}{Depends on whether \code{growth_TIC} is provided,}
#'     \item{residuals.}{}
#'   }}
#'   \item{TIC_starts:}{Only provided when \code{growth_TIC} is not NULL.}
#'   \item{beta0:}{Only provided when \code{growth_TIC} is not NULL.}
#' }
#' Each of these elements is a numeric vector or matrix containing the initial parameter estimates.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats var cov lm na.exclude
#'
.stabilize_covariate_effect_starts <- function(beta, curveFun, intrinsic, baseline_tvc_col = integer()){
  if (is.null(beta) || !isTRUE(intrinsic)) return(beta)

  beta <- as.matrix(beta)
  if (!curveFun %in% c("negative exponential", "EXP", "Jenss-Bayley", "JB", "bilinear spline", "BLS")){
    return(beta)
  }

  beta <- beta * 0.7

  if (length(baseline_tvc_col) > 0){
    beta[, baseline_tvc_col] <- beta[, baseline_tvc_col, drop = FALSE] * 0.5
  }

  shape_idx <- nrow(beta)
  if (curveFun %in% c("bilinear spline", "BLS")){
    beta[shape_idx, ] <- beta[shape_idx, ] * (0.25 / 0.7)
  } else if (curveFun %in% c("Jenss-Bayley", "JB")){
    beta[shape_idx, ] <- beta[shape_idx, ] * (0.15 / 0.7)
  }

  beta
}

.stabilize_intrinsic_nonlinear_starts <- function(starts, curveFun, intrinsic, time_values = NULL){
  if (!isTRUE(intrinsic) || !curveFun %in% c("Jenss-Bayley", "JB", "bilinear spline", "BLS")){
    return(starts)
  }

  shape_idx <- length(starts$Y_starts$means)
  if (shape_idx < 3 || is.null(dim(starts$Y_starts$covMatrix))) return(starts)

  span <- .curve_span(as.numeric(time_values))
  cov_mat <- starts$Y_starts$covMatrix

  if (curveFun %in% c("bilinear spline", "BLS")){
    time_range <- range(as.numeric(time_values), na.rm = TRUE)
    lower <- time_range[1] + 0.15 * span
    upper <- time_range[2] - 0.15 * span
    if (is.finite(lower) && is.finite(upper) && lower < upper){
      starts$Y_starts$means[shape_idx] <- min(max(starts$Y_starts$means[shape_idx], lower), upper)
    }
    cov_mat[-shape_idx, shape_idx] <- 0
    cov_mat[shape_idx, -shape_idx] <- 0
    cov_mat[shape_idx, shape_idx] <- min(cov_mat[shape_idx, shape_idx] * 0.1, (span / 8)^2)
  } else {
    cov_mat[shape_idx, ] <- cov_mat[shape_idx, ] * 0.15
    cov_mat[, shape_idx] <- cov_mat[, shape_idx] * 0.15
    if (curveFun %in% c("Jenss-Bayley", "JB")){
      min_abs <- max(1e-4 / span, max(0.05, 3 / span) / 100)
      if (abs(starts$Y_starts$means[shape_idx]) < min_abs){
        starts$Y_starts$means[shape_idx] <- if (starts$Y_starts$means[shape_idx] < 0) -min_abs else min_abs
      }
    }
  }

  starts$Y_starts$covMatrix <- .ensure_pd(cov_mat, label = "intrinsic nonlinear starts")
  starts
}

getUNI.initial <- function(dat, t_var, y_var, curveFun, records, growth_TIC, res_scale, intrinsic = FALSE){
  # Extract the measurement values and occasions
  dat_traj <- dat[, paste0(y_var, records)]
  dat_time <- dat[, paste0(t_var, records)]
  # Calculate the number of time points
  nT <- length(records)
  # If growth TICs are provided, calculate initial values for TIC-related parameters
  if (!is.null(growth_TIC)){
    dat_covariate <- dat[, growth_TIC]
    # Compute mean(s) and variance (variance-covariance matrix) of TICs
    TIC_mean <- if (length(growth_TIC) == 1){
      mean(dat_covariate, na.rm = TRUE)
    } else{
      apply(dat_covariate, 2, mean, na.rm = TRUE)
    }
    TIC_var <- var(dat_covariate, use = "pairwise.complete.obs")
    TIC_starts <- list(means = TIC_mean, covMatrix = if (length(growth_TIC) == 1){
      TIC_var
    } else{
      TIC_var[lower.tri(TIC_var, diag = TRUE)]
    })
  }
  # Calculate initial values for growth factor-related parameters
  growth_factor <- getUNI.GF(dat_traj = dat_traj, dat_time = dat_time, nT = nT, curveFun = curveFun)
  # Calculate the number of growth factors of the specified functional form
  nGF <- ncol(growth_factor)
  gf_cov_total <- cov(growth_factor, use = "pairwise.complete.obs")
  # Calculate the initial value of residual variance.
  # Priority: (1) user-specified res_scale  (2) data-driven mean_res_var  (3) 0.1 heuristic
  if (!is.null(res_scale)){
    residuals <- var(growth_factor[, 1], na.rm = TRUE) * res_scale
  } else {
    data_res_var <- attr(growth_factor, "mean_res_var")
    if (!is.null(data_res_var) && !is.na(data_res_var) && data_res_var > 0){
      residuals <- data_res_var
    } else{
      warning("Data-driven residual variance estimation produced invalid result; falling back to heuristic (0.1 * var(eta0)).", call. = FALSE)
      residuals <- var(growth_factor[, 1], na.rm = TRUE) * 0.1
    }
  }
  # If growth TICs are provided, calculate initial values for alpha0 (growth factor intercepts), psi_r (unexplained
  # var-cov matrix of growth factors), and beta0 (path coefficients)
  if (!is.null(growth_TIC)){
    gf_means <- sapply(1:nGF, function(p){
      lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p], dat_covariate)),
         na.action = na.exclude)$coefficients[1]
    })
    gf_beta <- matrix(sapply(1:nGF, function(p){
      lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p], dat_covariate)),
         na.action = na.exclude)$coefficients[-1]
    }), nrow = nGF, byrow = TRUE)
  }
  # If no growth TICs are provided, calculate initial values for mean0 (growth factor means) and psi0 (var-cov matrix
  # of growth factors)
  else if (is.null(growth_TIC)){
    gf_means <- sapply(1:nGF, function(p){
      lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p])),
         na.action = na.exclude)$coefficients[1]
    })
  }
  # If curveFun is nonparametric, calculate initial values for interval-specific absolute rate
  rel_rate <- NULL
  if (curveFun %in% c("nonparametric", "NonP")){
    abs_rate <- rep(0, nT - 1)
    for (j in 1:(nT - 1)){
      traj_delta <- dat_traj[, j + 1] - dat_traj[, j]
      time_delta <- dat_time[, j + 1] - dat_time[, j]
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
      rel_rate <- unlist(abs_rate)
    } else{
      rel_rate <- unlist(abs_rate / mean_eta1)
    }
  }
  # If curveFun is bilinear spline, calculate initial values for reparameterized functional form
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
    if (!is.null(growth_TIC)){
      gf_beta <- grad0 %*% gf_beta
    }
  }
  if (!is.null(growth_TIC)){
    gf_beta <- .stabilize_covariate_effect_starts(gf_beta, curveFun = curveFun, intrinsic = intrinsic)
    psi_r <- gf_cov_total - gf_beta %*% TIC_var %*% t(gf_beta)
    psi_r <- .ensure_pd(psi_r, label = "psi_r")
    Y_starts <- list(means = gf_means, covMatrix = psi_r, residuals = residuals)
    starts <- list(Y_starts = Y_starts, TIC_starts = TIC_starts, beta = gf_beta)
  } else{
    Y_starts <- list(means = gf_means, covMatrix = gf_cov_total, residuals = residuals)
    starts <- list(Y_starts = Y_starts)
  }
  if (!is.null(rel_rate)){
    starts$Y_starts$rel_rate <- rel_rate
  }
  starts <- .stabilize_intrinsic_nonlinear_starts(
    starts = starts, curveFun = curveFun, intrinsic = intrinsic,
    time_values = as.numeric(as.matrix(dat_time))
  )
  starts$Y_starts$covMatrix <- starts$Y_starts$covMatrix[lower.tri(starts$Y_starts$covMatrix, diag = TRUE)]
  return(starts)
}

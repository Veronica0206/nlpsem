#' @title Compute Initial Values for Parameters of Multivariate Latent Growth Curve Models or Latent Change Score Models
#'
#' @description This function computes the initial values for the parameters for a multivariate latent growth curve model
#' or a latent change score model.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables
#' with repeated measurements and occasions for multiple longitudinal outcomes. It takes the value passed from \code{getMGM()}.
#' @param t_var A vector of strings, with each element representing the prefix for column names related to the time
#' variable for the corresponding outcome variable at each study wave. It takes the value passed from \code{getMGM()}.
#' @param y_var A vector of strings, with each element representing the prefix for column names corresponding to a
#' particular outcome variable at each study wave. It takes the value passed from \code{getMGM()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the
#' value passed from \code{getMGM()}.
#' @param records A list of numeric vectors, with each vector specifying the indices of the observed study waves for
#' the corresponding outcome variable. It takes the value passed from \code{getMGM()}.
#' @param res_scale An optional numeric vector with each element representing the scaling factor for the initial calculation of
#' the residual variance. These values should be between \code{0} and \code{1}, exclusive, or \code{NULL} (default) to use
#' data-driven estimation with a heuristic of \code{0.1} as fallback. It takes the value passed from \code{getMGM()}.
#' @param res_cor An optional numeric value or vector for user-specified residual correlation between any two longitudinal outcomes
#' to calculate the corresponding initial value. When \code{NULL}, data-driven estimation is used. It takes the value passed from
#' \code{getMGM()}.
#'
#' @return A list containing the initial values for the parameters in the multivariate latent growth curve
#' model or a latent change score model growth curve model.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats var cov cor lm na.exclude
#'
.stabilize_multi_intrinsic_starts <- function(starts, curveFun, intrinsic, dat, t_var, y_var, records){
  if (!isTRUE(intrinsic) || !curveFun %in% c("Jenss-Bayley", "JB", "bilinear spline", "BLS")){
    return(starts)
  }

  cov_mat <- starts$covMatrix
  n_traj <- length(y_var)
  n_gf <- if (curveFun %in% c("bilinear spline", "BLS", "Jenss-Bayley", "JB")) 4 else 3
  shape_idx <- seq(n_gf, n_traj * n_gf, by = n_gf)

  for (traj in seq_len(n_traj)){
    time_values <- as.numeric(as.matrix(dat[, paste0(t_var[traj], records[[traj]])]))
    span <- .curve_span(time_values)
    idx <- shape_idx[traj]

    if (curveFun %in% c("bilinear spline", "BLS")){
      time_range <- range(time_values, na.rm = TRUE)
      lower <- time_range[1] + 0.15 * span
      upper <- time_range[2] - 0.15 * span
      if (is.finite(lower) && is.finite(upper) && lower < upper){
        starts$means[[traj]][4] <- min(max(starts$means[[traj]][4], lower), upper)
      }
      non_shape_idx <- setdiff(seq_len(nrow(cov_mat)), shape_idx)
      cov_mat[idx, non_shape_idx] <- 0
      cov_mat[non_shape_idx, idx] <- 0
      cov_mat[idx, idx] <- min(cov_mat[idx, idx] * 0.1, (span / 8)^2)
    } else {
      cov_mat[idx, ] <- cov_mat[idx, ] * 0.15
      cov_mat[, idx] <- cov_mat[, idx] * 0.15
      min_abs <- max(1e-4 / span, max(0.05, 3 / span) / 100)
      if (abs(starts$means[[traj]][4]) < min_abs){
        starts$means[[traj]][4] <- if (starts$means[[traj]][4] < 0) -min_abs else min_abs
      }
    }
  }

  if (length(shape_idx) > 1){
    shape_block <- cov_mat[shape_idx, shape_idx, drop = FALSE]
    off_diag <- row(shape_block) != col(shape_block)
    if (curveFun %in% c("bilinear spline", "BLS")){
      shape_block[off_diag] <- 0
    } else{
      shape_block[off_diag] <- shape_block[off_diag] * 0.25
    }
    cov_mat[shape_idx, shape_idx] <- shape_block
  }

  starts$covMatrix <- .ensure_pd(cov_mat, label = "intrinsic multivariate starts")
  starts
}

getMULTI.initial <- function(dat, t_var, y_var, curveFun, records, res_scale, res_cor, intrinsic = FALSE){
  # Growth factors from raw data
  GF_list <- uni_mean0 <- list(); uni_residual <- rep(0, length(y_var))
  # Extract the measurement values and occasions for each y_var
  for (traj in seq_along(y_var)){
    nT <- length(records[[traj]])
    dat_traj <- dat[, paste0(y_var[traj], records[[traj]])]
    dat_time <- dat[, paste0(t_var[traj], records[[traj]])]
    # Calculate individual growth factors for each longitudinal outcome using `getUNI.GF()` function
    uni_GF <- getUNI.GF(dat_traj = dat_traj, dat_time = dat_time, nT = nT, curveFun = curveFun)
    colnames(uni_GF) <- paste0(colnames(uni_GF), y_var[traj])
    GF_list[[length(GF_list) + 1]] <- uni_GF
    uni_mean0[[length(uni_mean0) + 1]] <- apply(uni_GF, 2, mean)
    # Calculate the residual variance for each longitudinal outcome.
    # Priority: (1) user-specified res_scale  (2) data-driven mean_res_var  (3) 0.1 heuristic
    if (!is.null(res_scale)){
      uni_residual[traj] <- var(uni_GF[, 1], na.rm = TRUE) * res_scale[traj]
    } else {
      data_res_var <- attr(uni_GF, "mean_res_var")
      if (!is.null(data_res_var) && !is.na(data_res_var) && data_res_var > 0){
        uni_residual[traj] <- data_res_var
      } else{
        warning("Data-driven residual variance estimation produced invalid result for outcome ", y_var[traj],
                "; falling back to heuristic (0.1 * var(eta0)).", call. = FALSE)
        uni_residual[traj] <- var(uni_GF[, 1], na.rm = TRUE) * 0.1
      }
    }
  }
  # If res_cor not provided, estimate from data using within-person residuals
  if (is.null(res_cor)){
    res_cor_est <- numeric()
    for (traj_i in 1:(length(y_var) - 1)){
      for (traj_j in (traj_i + 1):length(y_var)){
        nT_common <- min(length(records[[traj_i]]), length(records[[traj_j]]))
        dat_i <- dat[, paste0(y_var[traj_i], records[[traj_i]][1:nT_common]), drop = FALSE]
        dat_j <- dat[, paste0(y_var[traj_j], records[[traj_j]][1:nT_common]), drop = FALSE]
        # Remove person means to get within-person deviations
        resid_i <- dat_i - rowMeans(dat_i, na.rm = TRUE)
        resid_j <- dat_j - rowMeans(dat_j, na.rm = TRUE)
        # Correlate all residuals across person-timepoint observations
        r <- cor(unlist(resid_i), unlist(resid_j), use = "complete.obs")
        # Clamp to valid range and use 0.3 fallback if estimation fails
        if (is.na(r)) {
          warning("Data-driven residual correlation estimation failed for outcomes ",
                  y_var[traj_i], " and ", y_var[traj_j],
                  "; falling back to heuristic (res_cor = 0.3).", call. = FALSE)
          r <- 0.3
        }
        r <- max(min(r, 0.99), -0.99)
        res_cor_est <- c(res_cor_est, r)
      }
    }
    res_cor <- res_cor_est
  }
  # Combine the univariate growth factors into a single matrix
  growth_factor <- do.call("cbind", GF_list)
  # Calculate the initial values for var-cov matrix for the growth factors of all longitudinal outcome
  multi_var0 <- cov(growth_factor, use = "pairwise.complete.obs")
  multi_var0 <- .ensure_pd(multi_var0, label = "multi_var0")
  # Calculate the initial values for residual var-cov matrix
  multi_residuals <- diag(uni_residual)
  rc_idx <- 1
  for (traj_i in 1:(length(y_var) - 1)){
    for (traj_j in (traj_i + 1):length(y_var)){
      multi_residuals[traj_i, traj_j] <- res_cor[rc_idx] *
        sqrt(uni_residual[traj_i] * uni_residual[traj_j])
      rc_idx <- rc_idx + 1
    }
  }
  # Initial values of parameters related to growth factors
  starts <- list(means = uni_mean0, covMatrix = multi_var0, residuals = multi_residuals)
  # Transformed matrices if bilinear spline (BLS) functional form specified
  if (curveFun %in% c("bilinear spline", "BLS")){
    uni_func0 <- uni_grad0 <- list()
    for (traj in seq_along(y_var)){
      # Transformed matrices obtained by multivariate Delta method
      ## For mean vector
      uni_func0[[length(uni_func0) + 1]] <- matrix(c(1, uni_mean0[[traj]][4], 0, 0,
                                                     0, 0.5, 0.5, 0,
                                                     0, -0.5, 0.5, 0,
                                                     0, 0, 0, 1), nrow = 4, byrow = TRUE)
      ## For var-cov matrix
      uni_grad0[[length(uni_grad0) + 1]] <- matrix(c(1, uni_mean0[[traj]][4], 0, uni_mean0[[traj]][2],
                                                     0, 0.5, 0.5, 0,
                                                     0, -0.5, 0.5, 0,
                                                     0, 0, 0, 1), nrow = 4, byrow = TRUE)
      starts$means[[traj]] <- drop(uni_func0[[traj]] %*% starts$means[[traj]])
    }
    multi_grad0 <- as.matrix(Matrix::bdiag(uni_grad0))
    starts$covMatrix <- multi_grad0 %*% starts$covMatrix %*% t(multi_grad0)
  }
  # Absolute rate if nonparametric functional form specified
  if (curveFun %in% c("nonparametric", "NonP")){
    uni_rate <- list()
    for (traj in seq_along(y_var)){
      nT <- length(records[[traj]])
      dat_traj <- dat[, paste0(y_var[traj], records[[traj]])]
      dat_time <- dat[, paste0(t_var[traj], records[[traj]])]
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
      mean_eta1 <- starts$means[[traj]][2]
      if (abs(mean_eta1) < 1e-10){
        warning("Mean slope (eta1) is near zero for outcome ", y_var[traj],
                "; using abs_rate directly as rel_rate.", call. = FALSE)
        uni_rate[[length(uni_rate) + 1]] <- abs_rate
      } else{
        uni_rate[[length(uni_rate) + 1]] <- abs_rate / mean_eta1
      }
    }
    starts$rel_rate <- uni_rate
  }
  starts <- .stabilize_multi_intrinsic_starts(
    starts = starts, curveFun = curveFun, intrinsic = intrinsic,
    dat = dat, t_var = t_var, y_var = y_var, records = records
  )
  return(starts)
}

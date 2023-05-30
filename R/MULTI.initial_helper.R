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
#' @param res_scale A numeric vector with each element representing the scaling factor for the initial calculation of the residual
#' variance. These values should be between \code{0} and \code{1}, exclusive. It takes the value passed from \code{getMGM()}.
#' @param res_cor A numeric value or vector for user-specified residual correlation between any two longitudinal outcomes to calculate
#' the corresponding initial value. It takes the value passed from \code{getMGM()}.
#'
#' @return A list containing the initial values for the parameters in the multivariate latent growth curve
#' model or a latent change score model growth curve model.
#'
#' @importFrom stats var cov lm na.exclude
#'
getMULTI.initial <- function(dat, t_var, y_var, curveFun, records, res_scale, res_cor){
  # Growth factors from raw data
  GF_list <- uni_mean0 <- list(); uni_residual <- rep(0, length(y_var))
  # Extract the measurement values and occasions for each y_var
  for (traj in 1:length(y_var)){
    nT <- length(records[[traj]])
    dat_traj <- dat[, paste0(y_var[traj], records[[traj]])]
    dat_time <- dat[, paste0(t_var[traj], records[[traj]])]
    # Calculate individual growth factors for each longitudinal outcome using `getUNI.GF()` function
    uni_GF <- getUNI.GF(dat_traj = dat_traj, dat_time = dat_time, nT = nT, curveFun = curveFun)
    colnames(uni_GF) <- paste0(colnames(uni_GF), y_var[traj])
    GF_list[[length(GF_list) + 1]] <- uni_GF
    uni_mean0[[length(uni_mean0) + 1]] <- apply(uni_GF, 2, mean)
    # Calculate the residual variance for each longitudinal outcome
    uni_residual[traj] <- var(uni_GF[, 1]) * res_scale[traj]
  }
  # Combine the univariate growth factors into a single matrix
  growth_factor <- do.call("cbind", GF_list)
  # Calculate the initial values for var-cov matrix for the growth factors of all longitudinal outcome
  multi_var0 <- cov(growth_factor)
  # Calculate the initial values for residual var-cov matrix
  multi_residuals <- diag(uni_residual)
  for (traj_i in 1:(length(y_var) - 1)){
    for (traj_j in traj_i:(length(y_var) - 1)){
      multi_residuals[traj_i, (traj_j + 1)] <- res_cor[traj_i + traj_j - 1] *
        sqrt(uni_residual[[traj_i]] * uni_residual[[traj_j + 1]])
    }
  }
  # Initial values of parameters related to growth factors
  starts <- list(uni_mean0, multi_var0, multi_residuals)
  # Transformed matrices if bilinear spline (BLS) functional form specified
  if (curveFun %in% c("bilinear spline", "BLS")){
    uni_func0 <- uni_grad0 <- list()
    for (traj in 1:length(y_var)){
      # Transformed matrices obtained by multivariate Delta method
      ## For mean vector
      uni_func0[[length(uni_func0) + 1]] <- matrix(c(1, uni_mean0[[traj]][4], 0, 0,
                                                     0, 0.5, 0.5, 0,
                                                     0, -0.5, 0.5, 0,
                                                     0, 0, 0, 1), nrow = 4, byrow = T)
      ## For var-cov matrix
      uni_grad0[[length(uni_grad0) + 1]] <- matrix(c(1, uni_mean0[[traj]][4], 0, uni_mean0[[traj]][2],
                                                     0, 0.5, 0.5, 0,
                                                     0, -0.5, 0.5, 0,
                                                     0, 0, 0, 1), nrow = 4, byrow = T)
      starts[[1]][[traj]] <- uni_func0[[traj]] %*% starts[[1]][[traj]]
    }
    multi_grad0 <- as.matrix(Matrix::bdiag(uni_grad0))
    starts[[2]] <- multi_grad0 %*% starts[[2]] %*% t(multi_grad0)
  }
  # Absolute rate if nonparametric functional form specified
  if (curveFun %in% c("nonparametric", "NonP")){
    uni_rate <- list()
    for (traj in 1:length(y_var)){
      nT <- length(records[[traj]])
      dat_traj <- dat[, paste0(y_var[traj], records[[traj]])]
      dat_time <- dat[, paste0(t_var[traj], records[[traj]])]
      abs_rate <- rep(0, nT - 1)
      for (j in 1:(nT - 1)){
        abs_rate[j] <- as.numeric(lm(traj_delta ~ .,
                                     data = data.frame(traj_delta = dat_traj[, j + 1] - dat_traj[, j],
                                                       time_delta = dat_time[, j + 1] - dat_time[, j]),
                                     na.action = na.exclude)$coefficients[-1])
      }
      uni_rate[[length(uni_rate) + 1]] <- abs_rate/starts[[1]][[traj]][2]
    }
    starts <- append(starts, list(uni_rate))
  }
  return(starts)
}


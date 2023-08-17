#' @title Derive Individual Growth Factors for Latent Growth Curve Models or Latent Change Score Models with Time-Invariant Covariates (If Any)
#'
#' @description This function derives individual growth factors for the specified latent growth curve model or
#' latent change score model from raw data. These individual growth factors help further compute initial values
#' for parameters related to growth factors, time-invariant covariates (if any), and path coefficients (if any).
#'
#' @param dat_traj A data frame containing the records for the repeated measurements.
#' @param dat_time A data frame containing the records for measurement occasions associated with the repeated
#' measurements.
#' @param nT An integer representing the number of repeated measurements.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the
#' value passed from \code{getUNI.initial()} or \code{getMULTI.initial()} or \code{getTVC.initial()} or \code{getMIX.initial()}.
#'
#' @return A data frame containing the derived individual growth factors from the raw data.
#'
#' @keywords internal
#'
#' @importFrom stats coef nls
#'
getUNI.GF <- function(dat_traj, dat_time, nT, curveFun){
  # Calculate individual growth factors eta0 (initial values) and eta1 (slopes) for linear trajectories.
  if (curveFun %in% c("linear", "LIN")){
    # Apply the linear model to each row
    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- ind[1:nT]
      times <- ind[-(1:nT)]
      ind_model <- lm(measures ~ times)
      c(eta0 = coef(ind_model)[1], eta1 = coef(ind_model)[2])
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2])
  }
  # Calculate individual growth factors eta0 (initial values), eta1 (linear slopes), and eta2 (quadratic slopes)
  # for quadratic trajectories.
  else if (curveFun %in% c("quadratic", "QUAD")){
    # Apply the quadratic model to each row
    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- ind[1:nT]
      times <- ind[-(1:nT)]
      times_sq <- times^2
      ind_model <- lm(measures ~ times + times_sq)
      c(eta0 = coef(ind_model)[1], eta1 = coef(ind_model)[2], eta2 = coef(ind_model)[3])
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2], eta2 = coef_mat[, 3])
  }
  # Calculate individual growth factors eta0 (initial values), eta1 (asymptotic change between final and initial
  # values), and gamma (growth rate parameter) for negative exponential trajectories.
  else if (curveFun %in% c("negative exponential", "EXP")){
    # Apply the negative exponential model to each row
    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- as.numeric(ind[1:nT])
      times <- as.numeric(ind[-(1:nT)])
      a_init <- measures[1]
      b_init <- measures[nT] - measures[1]
      delta_t <- diff(times)
      delta_y <- diff(measures)
      log_slope <- log(abs(delta_y/delta_t) + 1e-9)
      c_init <- abs(coef(lm(log_slope ~ times[-nT]))[2])
      ind_model <- try(nls(measures ~ a + b * (1 - exp(-c * times)),
                           start = list(a = a_init, b = b_init, c = c_init)), silent = TRUE)
      if ("try-error" %in% class(ind_model)){
        c(eta0 = a_init, eta1 = b_init, gamma = c_init)
      }
      else{
        c(eta0 = coef(ind_model)[1], eta1 = coef(ind_model)[2], gamma = coef(ind_model)[3])
      }
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2], gamma = coef_mat[, 3])
  }
  # Calculate individual growth factors eta0 (initial values), eta1 (slope of linear asymptote), eta2 (change
  # between intercept of linear asymptote and initial status), and gamma (growth acceleration parameter) for
  # Jenss-Bayley trajectories.
  else if (curveFun %in% c("Jenss-Bayley", "JB")){
    # Apply the Jenss-Bayley model to each row
    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- as.numeric(ind[1:nT])
      times <- as.numeric(ind[-(1:nT)])
      a_init <- measures[1]
      b_init <- (measures[nT] - measures[1])/(times[nT] - times[1])
      c_init <- (measures[nT] - measures[1])
      # Compute the second difference of the response values
      delta2_y <- diff(measures, differences = 2)
      delta_t2 <- diff(times)^2
      curvature <- delta2_y/delta_t2[-1]
      lin_model <- lm(curvature ~ times[-c(1, length(times))])
      d_init <- -coef(lin_model)[2]
      ind_model <- try(nls(measures ~ a + b * times + c * (exp(d * times) - 1),
                           start = list(a = a_init, b = b_init, c = c_init, d = d_init)), silent = TRUE)
      if ("try-error" %in% class(ind_model)){
        c(eta0 = a_init, eta1 = b_init, eta2 = c_init, gamma = d_init)
      }
      else{
        c(eta0 = coef(ind_model)[1], eta1 = coef(ind_model)[2], eta2 = coef(ind_model)[3],
          gamma = coef(ind_model)[4])
      }
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2], eta2 = coef_mat[, 3],
                                gamma = coef_mat[, 4])
  }
  # Calculate individual growth factors eta0 (initial values), eta1 (slopes in the 1st stage), eta2 (slopes in
  # the 2nd stage), and gamma (knot locations) for bilinear spline trajectories.
  else if (curveFun %in% c("bilinear spline", "BLS")){
    # Apply the bilinear spline model to each row
    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- as.numeric(ind[1:nT])
      times <- as.numeric(ind[-(1:nT)])
      d_init <- (times[1] + times[nT])/2
      a_init <- measures[1]
      b_init <- (measures[which.max(times[times <= d_init])] - measures[1])/(times[which.max(times[times <= d_init])] - times[1])
      c_init <- (measures[nT] - measures[which.min(times[times >= d_init])])/(times[nT] - times[which.min(times[times >= d_init])])
      ind_model <- try(nls(measures ~ a + b * pmin(times, d) + c * pmax(times - d, 0),
                           start = list(a = a_init, b = b_init, c = c_init, d = d_init)), silent = TRUE)
      if ("try-error" %in% class(ind_model)){
        c(eta0 = a_init, eta1 = b_init, eta2 = c_init, gamma = d_init)
      }
      else{
        c(eta0 = coef(ind_model)[1], eta1 = coef(ind_model)[2], eta2 = coef(ind_model)[3],
          gamma = coef(ind_model)[4])
      }
    }))
    # Check if the rejection rate is exceeded and stop the function if necessary
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2], eta2 = coef_mat[, 3],
                                gamma = coef_mat[, 4])
  }
  # Calculate individual growth factors eta0 (initial values) and eta1 (slopes in the 1st time interval) for
  # nonparametric trajectories.
  else if (curveFun %in% c("nonparametric", "NonP")){
    slp <- rep(0, nT - 1)
    for (j in 1:(nT - 1)){
      slp[j] <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = (dat_traj[, j + 1] - dat_traj[, j])/
                                                                  (time_delta = dat_time[, j + 1] - dat_time[, j])),
                              na.action = na.exclude)$coefficients)
    }
    eta0 <- dat_traj[, 1]
    eta1 <- (dat_traj[, 2] - dat_traj[, 1])/(time_delta = dat_time[, 2] - dat_time[, 1])
    growth_factor <- data.frame(eta0 = eta0, eta1 = eta1)
  }
  return(growth_factor)
}

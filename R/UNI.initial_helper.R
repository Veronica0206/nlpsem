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
#' @param res_scale A numeric value representing the scaling factor for the initial calculation of the residual variance. This
#' value should be between \code{0} and \code{1}, exclusive. It takes the value passed from \code{getLGCM()} or \code{getLCSM()}.
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
#'
#' @importFrom stats var cov lm na.exclude
#'
getUNI.initial <- function(dat, t_var, y_var, curveFun, records, growth_TIC, res_scale){
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
    TIC_var <- var(dat_covariate, na.rm = TRUE)
    TIC_starts <- list(TIC_mean, if (length(growth_TIC) == 1){
      TIC_var
    } else{
      TIC_var[lower.tri(TIC_var, diag = TRUE)]
    })
  }
  # Calculate initial values for growth factor-related parameters
  growth_factor <- getUNI.GF(dat_traj = dat_traj, dat_time = dat_time, nT = nT, curveFun = curveFun)
  # Calculate the number of growth factors of the specified functional form
  nGF <- ncol(growth_factor)
  # Calculate the initial value of residual variance
  residuals <- var(growth_factor[, 1]) * res_scale
  # If growth TICs are provided, calculate initial values for alpha0 (growth factor intercepts), psi_r (unexplained
  # var-cov matrix of growth factors), and beta0 (path coefficients)
  if (!is.null(growth_TIC)){
    alpha0 <- sapply(1:nGF, function(p){
      lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p], dat_covariate)),
         na.action = na.exclude)$coefficients[1]
    })
    beta0 <- matrix(sapply(1:nGF, function(p){
      lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p], dat_covariate)),
         na.action = na.exclude)$coefficients[-1]
    }), nrow = nGF, byrow = TRUE)
    psi_r <- cov(growth_factor) - beta0 %*% TIC_var %*% t(beta0)
    Y_starts <- list(alpha0, psi_r, residuals)
    starts <- list(Y_starts, TIC_starts, beta0)
  }
  # If no growth TICs are provided, calculate initial values for mean0 (growth factor means) and psi0 (var-cov matrix
  # of growth factors)
  else if (is.null(growth_TIC)){
    mean0 <- sapply(1:nGF, function(p){
      lm(GF ~ ., data = data.frame(cbind(GF = growth_factor[, p])),
         na.action = na.exclude)$coefficients[1]
    })
    psi0 <- cov(growth_factor)
    Y_starts <- list(mean0, psi0, residuals)
    starts <- list(Y_starts)
  }
  # If curveFun is nonparametric, calculate initial values for interval-specific absolute rate
  if (curveFun %in% c("nonparametric", "NonP")){
    abs_rate <- rep(0, nT - 1)
    for (j in 1:(nT - 1)){
      abs_rate[j] <- as.numeric(lm(traj_delta ~ .,
                                   data = data.frame(traj_delta = dat_traj[, j + 1] - dat_traj[, j],
                                                     time_delta = dat_time[, j + 1] - dat_time[, j]),
                                   na.action = na.exclude)$coefficients[-1])
    }
    starts[[1]][[length(starts[[1]]) + 1]] <- unlist(abs_rate/starts[[1]][[1]][2])
  }
  # If curveFun is bilinear spline, calculate initial values for reparameterized functional form
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
      starts[[3]] <- grad0 %*% starts[[3]]
    }
  }
  starts[[1]][[2]] <- starts[[1]][[2]][lower.tri(starts[[1]][[2]], diag = TRUE)]
  return(starts)
}

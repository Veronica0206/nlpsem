#' @title Compute Initial Values for Parameters of Longitudinal Mediation Models
#'
#' @description This function computes the initial values of the parameters for a longitudinal mediation model.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables
#' with repeated measurements and occasions for multiple longitudinal processes and a baseline predictor when applicable.
#' It takes the value passed from \code{getMediation()}.
#' @param t_var A vector of strings, with each element representing the prefix for column names related to the time
#' variable for the corresponding longitudinal variable at each study wave. It takes the value passed from \code{getMediation()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getMediation()}.
#' @param m_var A string specifying the prefix of the column names corresponding to the mediator variable at each study wave.
#' It takes the value passed from \code{getMediation()}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. It takes the value passed from \code{getMediation()}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. It takes the value passed from
#' \code{getMediation()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options include: \code{"linear"}
#' (or \code{"LIN"}), and \code{"bilinear spline"} (or \code{"BLS"}). It takes the value passed from \code{getMediation()}.
#' @param records A list of numeric vectors, with each vector specifying the indices of the observed study waves for
#' the corresponding longitudinal variable. It takes the value passed from \code{getMediation()}.
#' @param res_scale A numeric vector with each element representing the scaling factor for the initial calculation of the residual
#' variance. These values should be between \code{0} and \code{1}, exclusive. It takes the value passed from \code{getMediation()}.
#' @param res_cor A numeric value or vector for user-specified residual correlation between any two longitudinal processes to calculate
#' the corresponding initial value. It takes the value passed from \code{getMediation()}.
#'
#' @return A list containing the initial values for parameters related to growth factors and path coefficients of a
#' longitudinal mediation model.
#'
getMED.initial <- function(dat, t_var, y_var, m_var, x_type, x_var, curveFun, records, res_scale, res_cor){
  # Define the number of growth factors for each functional form
  nGF <- ifelse(curveFun %in% c("linear", "LIN"), 2, 3)
  # Calculate individual growth factors for longitudinal outcome using `getUNI.GF()` function
  Y_growth_factor <- getUNI.GF(dat_traj = dat[, paste0(y_var, records[[1]])], dat_time = dat[, paste0(t_var[1], records[[1]])],
                               nT = length(records[[1]]), curveFun = curveFun)
  # Calculate individual growth factors for longitudinal mediator using `getUNI.GF()` function
  M_growth_factor <- getUNI.GF(dat_traj = dat[, paste0(m_var, records[[2]])], dat_time = dat[, paste0(t_var[2], records[[2]])],
                               nT = length(records[[2]]), curveFun = curveFun)
  # Reparameterize growth factors to first slope, knot measurement, second slope, and knot if bilinear spline
  # functional form specified
  if (curveFun %in% c("bilinear spline", "BLS")){
    Y_tmp_etar <- Y_growth_factor[, 1] + Y_growth_factor[, 2] * Y_growth_factor[, 4]
    Y_tmp_eta1 <- Y_growth_factor[, 2]
    Y_growth_factor[, 1] <- Y_tmp_eta1
    Y_growth_factor[, 2] <- Y_tmp_etar
    names(Y_growth_factor)[1:2] <- c("eta1", "etar")
    M_tmp_etar <- M_growth_factor[, 1] + M_growth_factor[, 2] * M_growth_factor[, 4]
    M_tmp_eta1 <- M_growth_factor[, 2]
    M_growth_factor[, 1] <- M_tmp_eta1
    M_growth_factor[, 2] <- M_tmp_etar
    names(M_growth_factor)[1:2] <- c("eta1", "etar")
  }
  # Calculate initial values for all parameters for a model with baseline predictor
  if (x_type == "baseline"){
    # Calculate initial values for parameters related to baseline predictor
    X_mean <- mean(dat[, x_var], na.rm = T)
    X_var <- var(dat[, x_var], na.rm = T)
    starts.X <- list(X_mean, X_var)
    # Calculate initial values for parameters related to longitudinal mediator
    M_alpha <- sapply(1:nGF, function(p){
      lm(GF ~ ., data = data.frame(cbind(GF = M_growth_factor[, p], dat[, x_var])),
         na.action = na.exclude)$coefficients[1]
    })
    M_beta <- matrix(sapply(1:nGF, function(p){
      lm(GF ~ ., data = data.frame(cbind(GF = M_growth_factor[, p], dat[, x_var])),
         na.action = na.exclude)$coefficients[-1]
    }), nrow = nGF, byrow = TRUE)
    M_alpha <- if (curveFun %in% c("bilinear spline", "BLS")){
      c(M_alpha, mean(M_growth_factor[, 4]))
    } else{
      M_alpha
    }
    M_psi_r <- if (curveFun %in% c("bilinear spline", "BLS")){
      var(M_growth_factor[, 1:3]) - M_beta %*% X_var %*% t(M_beta)
    } else{
      var(M_growth_factor[, 1:2]) - M_beta %*% X_var %*% t(M_beta)
    }
    M_residuals <- ifelse(curveFun %in% c("linear", "LIN"), var(M_growth_factor[, 1]) * res_scale[2],
                          var(M_growth_factor[, 2]) * res_scale[2])
    starts.M <- list(M_alpha, M_beta, M_psi_r, M_residuals)
    # Calculate initial values for parameters related to longitudinal outcome
    if (curveFun %in% c("linear", "LIN")){
      Y_reg0 <- as.numeric(lm(Yeta0 ~ ., data = data.frame(cbind(Yeta0 = Y_growth_factor[, 1], dat[, x_var],
                                                                 M_growth_factor[, 1])),
                              na.action = na.exclude)$coefficients)
      Y_reg1 <- as.numeric(lm(Yeta1 ~ ., data = data.frame(cbind(Yeta1 = Y_growth_factor[, 2], dat[, x_var],
                                                                 M_growth_factor[, 1:2])),
                              na.action = na.exclude)$coefficients)
      Y_alpha <- c(Y_reg0[1], Y_reg1[1])
      Y_beta <- matrix(c(Y_reg0[2], Y_reg0[3], 0, Y_reg1[-1]), byrow = T, nrow = 2, ncol = 3)
      Y_psi_r <- var(Y_growth_factor) - t(t(Y_beta[, 1])) %*% X_var %*% t(Y_beta[, 1]) -
        Y_beta[, 2:3] %*% M_psi_r %*% t(Y_beta[, 2:3])
      Y_residuals <- var(Y_growth_factor[, 1]) * res_scale[1]
      starts.Y <- list(Y_alpha, Y_beta[, 1], Y_beta[, 2:3], Y_psi_r, Y_residuals)
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      Y_reg_1 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = Y_growth_factor[, 1],
                                                                 dat[, x_var], M_growth_factor[, 1]),
                               na.action = na.exclude)$coefficients)
      Y_reg_r <- as.numeric(lm(Y_gammaV ~ ., data = data.frame(Y_gammaV = Y_growth_factor[, 2],
                                                               dat[, x_var], M_growth_factor[, 1:2]),
                               na.action = na.exclude)$coefficients)
      Y_reg_2 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = Y_growth_factor[, 3],
                                                                 dat[, x_var], M_growth_factor[, 1:3]),
                               na.action = na.exclude)$coefficients)
      Y_alpha <- c(Y_reg_1[1], Y_reg_r[1], Y_reg_2[1], mean(Y_growth_factor[, 4]))
      Y_beta <- matrix(c(Y_reg_1[2:3], 0, 0,
                         Y_reg_r[2:4], 0,
                         Y_reg_2[-1]), byrow = T, nrow = 3, ncol = 4)
      Y_psi_r <- cov(Y_growth_factor[, 1:3]) - t(t(Y_beta[, 1])) %*% X_var %*% t(Y_beta[, 1]) - Y_beta[, 2:4] %*% M_psi_r %*% t(Y_beta[, 2:4])
      Y_residuals <- var(Y_growth_factor[, 2]) * res_scale[1]
      starts.Y <- list(Y_alpha, Y_beta[, 1], Y_beta[, 2:4], Y_psi_r, Y_residuals)
    }
  }
  # Calculate initial values for all parameters for a model with longitudinal predictor
  else if (x_type == "longitudinal"){
    # Calculate individual growth factors for longitudinal predictor using `getUNI.GF()` function if the predictor
    # is longitudinal
    X_growth_factor <- getUNI.GF(dat_traj = dat[, paste0(x_var, records[[3]])], dat_time = dat[, paste0(t_var[3], records[[3]])],
                                 nT = length(records[[3]]), curveFun = curveFun)
    # Reparameterize growth factors to first slope, knot measurement, second slope, and knot if bilinear spline
    # functional form specified
    if (curveFun %in% c("bilinear spline", "BLS")){
      X_tmp_etar <- X_growth_factor[, 1] + X_growth_factor[, 2] * X_growth_factor[, 4]
      X_tmp_eta1 <- X_growth_factor[, 2]
      X_growth_factor[, 1] <- X_tmp_eta1
      X_growth_factor[, 2] <- X_tmp_etar
      names(X_growth_factor)[1:2] <- c("eta1", "etar")
    }
    # Calculate initial values for parameters related to baseline predictor
    X_mean <- sapply(1:nGF, function(p){
      lm(GF ~ ., data = data.frame(cbind(GF = X_growth_factor[, p])),
         na.action = na.exclude)$coefficients[1]
    })
    X_mean <- if (curveFun %in% c("bilinear spline", "BLS")){
      c(X_mean, mean(X_growth_factor[, 4]))
    } else{
      X_mean
    }
    X_psi0 <- if (curveFun %in% c("bilinear spline", "BLS")){
      var(X_growth_factor[, 1:3])
    } else{
      var(X_growth_factor[, 1:2])
    }
    X_residuals <- ifelse(curveFun %in% c("linear", "LIN"), var(X_growth_factor[, 1]) * res_scale[3],
                          var(X_growth_factor[, 2]) * res_scale[3])
    starts.X <- list(X_mean, X_psi0, X_residuals)
    if (curveFun %in% c("linear", "LIN")){
      # Calculate initial values for parameters related to longitudinal mediator
      M_reg0 <- as.numeric(lm(Meta0 ~ ., data = data.frame(cbind(Meta0 = M_growth_factor[, 1],
                                                                 X_growth_factor[, 1])),
                              na.action = na.exclude)$coefficients)
      M_reg1 <- as.numeric(lm(Meta1 ~ ., data = data.frame(cbind(Meta1 = M_growth_factor[, 2],
                                                                 X_growth_factor[, 1:2])),
                              na.action = na.exclude)$coefficients)
      M_alpha <- c(M_reg0[1], M_reg1[1])
      M_beta <- matrix(c(M_reg0[-1], 0, M_reg1[-1]), byrow = T, nrow = 2, ncol = 2)
      M_psi_r <- var(M_growth_factor) - M_beta %*% X_psi0 %*% t(M_beta)
      M_residuals <- var(M_growth_factor[, 1]) * res_scale[2]
      starts.M <- list(M_alpha, M_beta, M_psi_r, M_residuals)
      # Calculate initial values for parameters related to longitudinal outcome
      Y_reg0 <- as.numeric(lm(Yeta0 ~ ., data = data.frame(cbind(Yeta0 = Y_growth_factor[, 1],
                                                                 M_growth_factor[, 1],
                                                                 X_growth_factor[, 1])),
                              na.action = na.exclude)$coefficients)
      Y_reg1 <- as.numeric(lm(Yeta1 ~ ., data = data.frame(cbind(Yeta1 = Y_growth_factor[, 2],
                                                                 M_growth_factor[, 1:2],
                                                                 X_growth_factor[, 1:2])),
                              na.action = na.exclude)$coefficients)
      Y_alpha <- c(Y_reg0[1], Y_reg1[1])
      Y_beta <- matrix(c(Y_reg0[2], 0, Y_reg0[3], 0,
                         Y_reg1[-1]), byrow = T, nrow = 2, ncol = 4)
      Y_psi_r <- var(Y_growth_factor) - Y_beta[, 1:2] %*% X_psi0 %*% t(Y_beta[, 1:2]) -
        Y_beta[, 3:4] %*% M_psi_r %*% t(Y_beta[, 3:4])
      Y_residuals <- var(Y_growth_factor[, 1]) * res_scale[1]
      starts.Y <- list(Y_alpha, Y_beta[, 1:2], Y_beta[, 3:4], Y_psi_r, Y_residuals)
      starts <- list(starts.X, starts.M, starts.Y)
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      # Calculate initial values for parameters related to longitudinal mediator
      M_reg_1 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = M_growth_factor[, 1],
                                                                 X_growth_factor[, 1]),
                               na.action = na.exclude)$coefficients)
      M_reg_r <- as.numeric(lm(M_gammaV ~ ., data = data.frame(M_gammaV = M_growth_factor[, 2],
                                                               X_growth_factor[, 1:2]),
                               na.action = na.exclude)$coefficients)
      M_reg_2 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = M_growth_factor[, 3],
                                                                 X_growth_factor[, 1:3]),
                               na.action = na.exclude)$coefficients)
      M_alpha <- c(M_reg_1[1], M_reg_r[1], M_reg_2[1], mean(M_growth_factor[, 4]))
      M_beta <- matrix(c(M_reg_1[-1], 0, 0, M_reg_r[-1], 0, M_reg_2[-1]), byrow = T, nrow = 3, ncol = 3)
      M_psi_r <- var(M_growth_factor[, 1:3]) - M_beta %*% X_psi0 %*% t(M_beta)
      M_residuals <- var(M_growth_factor[, 2]) * res_scale[2]
      starts.M <- list(M_alpha, M_beta, M_psi_r, M_residuals)
      # Calculate initial values for parameters related to longitudinal outcome
      Y_reg_1 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = Y_growth_factor[, 1],
                                                                 M_growth_factor[, 1],
                                                                 X_growth_factor[, 1]),
                               na.action = na.exclude)$coefficients)
      Y_reg_r <- as.numeric(lm(Y_gammaV ~ ., data = data.frame(Y_gammaV = Y_growth_factor[, 2],
                                                               M_growth_factor[, 1:2],
                                                               X_growth_factor[, 1:2]),
                               na.action = na.exclude)$coefficients)
      Y_reg_2 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = Y_growth_factor[, 3],
                                                                 M_growth_factor[, 1:3],
                                                                 X_growth_factor[, 1:3]),
                               na.action = na.exclude)$coefficients)
      Y_alpha <- c(Y_reg_1[1], Y_reg_r[1], Y_reg_2[1], mean(Y_growth_factor[, 4]))
      Y_beta <- matrix(c(Y_reg_1[2], rep(0, 2), Y_reg_1[3], rep(0, 2),
                         Y_reg_r[2:3], 0, Y_reg_r[4:5], 0, Y_reg_2[-1]),
                       byrow = T, nrow = 3, ncol = 6)
      Y_psi_r <- var(Y_growth_factor[, 1:3]) - Y_beta[, 1:3] %*% M_psi_r %*% t(Y_beta[, 1:3]) -
        Y_beta[, 4:6] %*% X_psi0 %*% t(Y_beta[, 4:6])
      Y_residuals <- var(Y_growth_factor[, 2]) * res_scale[1]
      starts.Y <- list(Y_alpha, Y_beta[, 1:3], Y_beta[, 4:6], Y_psi_r, Y_residuals)
    }
  }
  starts <- list(starts.X, starts.M, starts.Y)
  return(starts)
}


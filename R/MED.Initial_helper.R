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
#' @param res_scale An optional numeric vector with each element representing the scaling factor for the initial calculation of
#' the residual variance. These values should be between \code{0} and \code{1}, exclusive, or \code{NULL} (default) to use
#' data-driven estimation with a heuristic of \code{0.1} as fallback. It takes the value passed from \code{getMediation()}.
#' @param res_cor An optional numeric value or vector for user-specified residual correlation between any two longitudinal processes
#' to calculate the corresponding initial value. When \code{NULL}, data-driven estimation is used. It takes the value passed from
#' \code{getMediation()}.
#'
#' @return A list containing the initial values for parameters related to growth factors and path coefficients of a
#' longitudinal mediation model.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats cor
#'
.stabilize_med_path_starts <- function(beta, curveFun){
  if (is.null(beta) || !curveFun %in% c("bilinear spline", "BLS")) return(beta)

  beta <- as.matrix(beta) * 0.7
  if (nrow(beta) >= 3){
    beta[3, ] <- beta[3, ] * (0.5 / 0.7)
  }
  beta
}

getMED.initial <- function(dat, t_var, y_var, m_var, x_type, x_var, curveFun, records, res_scale, res_cor){
  # If res_cor not provided, estimate from data using within-person residuals
  if (is.null(res_cor)){
    if (x_type == "baseline"){
      traj_var <- c(y_var, m_var)
      traj_records <- records[1:2]
    } else {
      traj_var <- c(y_var, m_var, x_var)
      traj_records <- records[1:3]
    }
    res_cor_est <- numeric()
    for (traj_i in 1:(length(traj_var) - 1)){
      for (traj_j in (traj_i + 1):length(traj_var)){
        nT_common <- min(length(traj_records[[traj_i]]), length(traj_records[[traj_j]]))
        dat_i <- dat[, paste0(traj_var[traj_i], traj_records[[traj_i]][1:nT_common]), drop = FALSE]
        dat_j <- dat[, paste0(traj_var[traj_j], traj_records[[traj_j]][1:nT_common]), drop = FALSE]
        resid_i <- dat_i - rowMeans(dat_i, na.rm = TRUE)
        resid_j <- dat_j - rowMeans(dat_j, na.rm = TRUE)
        r <- cor(unlist(resid_i), unlist(resid_j), use = "complete.obs")
        if (is.na(r)) r <- 0.3
        r <- max(min(r, 0.99), -0.99)
        res_cor_est <- c(res_cor_est, r)
      }
    }
    res_cor <- res_cor_est
  }
  # Helper: use data-driven residual variance if available, otherwise fall back to heuristic
  .get_res_init <- function(gf, col_idx, scale){
    data_res_var <- attr(gf, "mean_res_var")
    if (!is.null(data_res_var) && !is.na(data_res_var) && data_res_var > 0){
      data_res_var
    } else if (!is.null(scale)){
      var(gf[, col_idx]) * scale
    } else{
      var(gf[, col_idx]) * 0.1
    }
  }
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
    X_mean <- mean(dat[, x_var], na.rm = TRUE)
    X_var <- var(dat[, x_var], na.rm = TRUE)
    starts.X <- list(means = X_mean, covMatrix = X_var)
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
      M_beta <- .stabilize_med_path_starts(M_beta, curveFun = curveFun)
      .ensure_pd(var(M_growth_factor[, 1:3]) - M_beta %*% X_var %*% t(M_beta), label = "MED M_psi_r")
    } else{
      var(M_growth_factor[, 1:2]) - M_beta %*% X_var %*% t(M_beta)
    }
    M_residuals <- .get_res_init(M_growth_factor,
                                  ifelse(curveFun %in% c("linear", "LIN"), 1, 2), res_scale[2])
    starts.M <- list(means = M_alpha, beta = M_beta, covMatrix = M_psi_r, residuals = M_residuals)
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
      Y_residuals <- .get_res_init(Y_growth_factor, 1, res_scale[1])
      starts.Y <- list(means = Y_alpha, beta_x = Y_beta[, 1], beta_m = Y_beta[, 2:3], covMatrix = Y_psi_r, residuals = Y_residuals)
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
      Y_beta <- .stabilize_med_path_starts(Y_beta, curveFun = curveFun)
      Y_psi_r <- cov(Y_growth_factor[, 1:3]) - t(t(Y_beta[, 1])) %*% X_var %*% t(Y_beta[, 1]) - Y_beta[, 2:4] %*% M_psi_r %*% t(Y_beta[, 2:4])
      Y_psi_r <- .ensure_pd(Y_psi_r, label = "MED Y_psi_r")
      Y_residuals <- .get_res_init(Y_growth_factor, 2, res_scale[1])
      starts.Y <- list(means = Y_alpha, beta_x = Y_beta[, 1], beta_m = Y_beta[, 2:4], covMatrix = Y_psi_r, residuals = Y_residuals)
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
      .ensure_pd(var(X_growth_factor[, 1:3]), label = "MED X_psi0")
    } else{
      var(X_growth_factor[, 1:2])
    }
    X_residuals <- .get_res_init(X_growth_factor,
                                  ifelse(curveFun %in% c("linear", "LIN"), 1, 2), res_scale[3])
    starts.X <- list(means = X_mean, covMatrix = X_psi0, residuals = X_residuals)
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
      M_residuals <- .get_res_init(M_growth_factor, 1, res_scale[2])
      starts.M <- list(means = M_alpha, beta = M_beta, covMatrix = M_psi_r, residuals = M_residuals)
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
      Y_residuals <- .get_res_init(Y_growth_factor, 1, res_scale[1])
      starts.Y <- list(means = Y_alpha, beta_x = Y_beta[, 3:4], beta_m = Y_beta[, 1:2], covMatrix = Y_psi_r, residuals = Y_residuals)
      starts <- list(X = starts.X, M = starts.M, Y = starts.Y)
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
      M_beta <- .stabilize_med_path_starts(M_beta, curveFun = curveFun)
      M_psi_r <- var(M_growth_factor[, 1:3]) - M_beta %*% X_psi0 %*% t(M_beta)
      M_psi_r <- .ensure_pd(M_psi_r, label = "MED M_psi_r")
      M_residuals <- .get_res_init(M_growth_factor, 2, res_scale[2])
      starts.M <- list(means = M_alpha, beta = M_beta, covMatrix = M_psi_r, residuals = M_residuals)
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
      Y_beta <- .stabilize_med_path_starts(Y_beta, curveFun = curveFun)
      Y_psi_r <- var(Y_growth_factor[, 1:3]) - Y_beta[, 1:3] %*% M_psi_r %*% t(Y_beta[, 1:3]) -
        Y_beta[, 4:6] %*% X_psi0 %*% t(Y_beta[, 4:6])
      Y_psi_r <- .ensure_pd(Y_psi_r, label = "MED Y_psi_r")
      Y_residuals <- .get_res_init(Y_growth_factor, 2, res_scale[1])
      starts.Y <- list(means = Y_alpha, beta_x = Y_beta[, 4:6], beta_m = Y_beta[, 1:3], covMatrix = Y_psi_r, residuals = Y_residuals)
    }
  }
  starts <- list(X = starts.X, M = starts.M, Y = starts.Y)
  return(starts)
}

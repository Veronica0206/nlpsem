#' @title Get Additional Parameters Related to Interval-specific Slopes, Interval-specific Changes and Values of Change-from-
#' baseline for Mixture Model with Multivariate Latent Change Score Models as Submodels
#'
#' @description This function derives additional parameters for mixture model with multivariate latent change score models as submodels.
#' It specifies the means and variances of interval-specific slopes, interval-specific changes, and values of change-from-
#' baseline.
#'
#' @param dat Data frame. Contains the observed variables, including repeated measurements and occasions, and covariates
#' if any. It takes the value passed from `getMIX()`.
#' @param nClass Numeric. Indicates the number of latent classes. It takes the value passed from `getMIX()`.
#' @param t_var String. Prefix of the column names corresponding to the time variable at each time point. It takes the value
#' passed from `getMIX()`.
#' @param records Numeric vector. Indices of the observed time points. It takes the value passed from `getMIX()`.
#' @param y_var String. Prefix of the column names corresponding to the outcome variable at each time point. It takes the
#' value passed from `getMIX()`.
#' @param curveFun String. The functional form of the growth curve. Supported options include: "linear" (or "LIN"),
#' "quadratic" (or "QUAD"), "negative exponential" (or "EXP"), "Jenss-Bayley" (or "JB"), "bilinear spline" (or "BLS"),
#' and "nonparametric" (or "NonP"). It takes the value passed from `getMIX()`.
#' @param intrinsic Logical. Whether an intrinsically nonlinear longitudinal model is built up. It takes the value passed
#' from `getMIX()`.
#' @param starts A list of initial values for free parameters, either takes the value passed from `getMIX()` or derived by
#' the helper function `getMIX.initial()`.
#'
#' @return A list containing the specification of the means and variances of interval-specific slopes, interval-specific
#' changes, and values of change-from-baseline for a mixture model with multivariate latent change score models as submodels.
#'
#' @keywords internal
#'
#' @importFrom OpenMx mxMatrix mxAlgebra mxAlgebraFromString diag2vec
#'
getMIX_MULTI.addpara <- function(dat, nClass, t_var, records, y_var, curveFun, intrinsic = NULL, starts){
  AddPara_CL <- AddPara_L <- list()
  for (k in 1:nClass){
    for (traj in 1:length(y_var)){
      # Define mxMatrix to include the mean vector of each measurement occasion and the mean vector of
      # each middle time point, which are useful to derive the means and variances of interval-specific
      # slopes, interval-specific changes and change from baseline.
      m_time <- sapply(paste0(t_var[traj], records[[traj]]), function(x) mean(dat[, x]))
      m_mid_time <- c(0, (m_time[-1] + m_time[-length(records[[traj]])])/2)
      M_TIME <- mxMatrix("Full", length(records[[traj]]), 1, free = FALSE, values = m_time,
                         name = paste0(y_var[traj], "M_TIME"))
      M_MID_TIME <- mxMatrix("Full", length(records[[traj]]), 1, free = FALSE, values = m_mid_time,
                             name = paste0(y_var[traj], "M_MID_TIME"))
      # Define mxMatrix to include the mean value time lag matrix, which is useful to derive the means and
      # variances of interval-specific slopes, interval-specific changes and change from baseline.
      m_lag <- diff(m_time)
      m_lag_sq <- m_lag^2
      M_LAG <- mxMatrix("Full", length(records[[traj]]), 1, free = FALSE, values = c(0, m_lag),
                        name = paste0(y_var[traj], "M_LAG"))
      M_LAG_SQ <- mxMatrix("Full", length(records[[traj]]), 1, free = FALSE, values = c(0, m_lag_sq),
                           name = paste0(y_var[traj], "M_LAG_SQ"))
      m_Omega_val <- matrix(0, nrow = length(m_time), ncol = length(m_lag))
      for (j in records[[traj]][-1]){
        m_Omega_val[j, ] <- c(m_lag[1:(j - 1)], rep(0, length(m_time) - j))
      }
      # Define Omega matrix using mxMatrix
      m_Omega <- mxMatrix(type = "Full", nrow = length(m_time), ncol = length(m_lag),
                          values = m_Omega_val, free = FALSE, name = paste0(y_var[traj], "Omega"))
      if (curveFun %in% c("nonparametric", "NonP")){
        rate_loads <- mxMatrix("Full", nrow = length(m_lag), ncol = 1, c(F, rep(T, length(m_lag) - 1)),
                               values = c(1, starts[[k]][[4]][[traj]][-1]),
                               labels = paste0("c", k, y_var[traj], "_rel_rate", 1:length(m_lag)),
                               byrow = T, name = paste0("c", k, y_var[traj], "r_loads"))
        status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "c", k, y_var[traj], "r_loads"),
                                            name = paste0("c", k, y_var[traj], "s_loads"))
        slp_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                            "c", k, y_var[traj], "_mean0[2, ]"),
                                     name = paste0("c", k, y_var[traj], "slp_m"))
        slp_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                            "c", k, y_var[traj], "_psi0[2, 2] %*% t(",
                                            "c", k, y_var[traj], "r_loads)"),
                                     name = paste0("c", k, y_var[traj], "slp_v"))
        chg_inv_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                         name = paste0("c", k, y_var[traj], "chg_inv_m"))
        chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(c",
                                                k, y_var[traj], "slp_v))"),
                                         name = paste0("c", k, y_var[traj], "chg_inv_v"))
        chg_bl_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                               "c", k, y_var[traj], "_mean0[2, ]"),
                                        name = paste0("c", k, y_var[traj], "chg_bl_m"))
        chg_bl_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                               "c", k, y_var[traj], "_psi0[2, 2] %*% t(c", k,
                                               y_var[traj], "s_loads[-1, ])"),
                                        name = paste0("c", k, y_var[traj], "chg_bl_v"))
        AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                  M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
      }
      else if (curveFun %in% c("quadratic", "QUAD")){
        ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                         free = FALSE, name = "ADD1")
        rate_loads <- mxAlgebraFromString(paste0("cbind(ADD1, ", y_var[traj], "M_MID_TIME[-1, ] * 2)"),
                                          name = paste0("c", k, y_var[traj], "r_loads"))
        status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "c", k, y_var[traj], "r_loads"),
                                            name = paste0("c", k, y_var[traj], "s_loads"))
        slp_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                            "c", k, y_var[traj], "_mean0[2:3, ]"),
                                     name = paste0("c", k, y_var[traj], "slp_m"))
        slp_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                            "c", k, y_var[traj], "_psi0[2:3, 2:3] %*% t(",
                                            "c", k, y_var[traj], "r_loads)"),
                                     name = paste0("c", k, y_var[traj], "slp_v"))
        chg_inv_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                         name = paste0("c", k, y_var[traj], "chg_inv_m"))
        chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(c",
                                                k, y_var[traj], "slp_v))"),
                                         name = paste0("c", k, y_var[traj], "chg_inv_v"))
        chg_bl_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                               "c", k, y_var[traj], "_mean0[2:3, ]"),
                                        name = paste0("c", k, y_var[traj], "chg_bl_m"))
        chg_bl_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                               "c", k, y_var[traj], "_psi0[2:3, 2:3] %*% t(c", k,
                                               y_var[traj], "s_loads[-1, ])"),
                                        name = paste0("c", k, y_var[traj], "chg_bl_v"))
        AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                  M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
      }
      else if (curveFun %in% c("negative exponential", "EXP")){
        if (intrinsic){
          rate_loads <- mxAlgebraFromString(paste0("cbind(c", k, y_var[traj], "_mug * exp(-c", k, y_var[traj], "_mug * ",
                                                   y_var[traj], "M_MID_TIME[-1, ]), c", k, y_var[traj], "_mueta1 * exp(-c", k,
                                                   y_var[traj], "_mug * ", y_var[traj], "M_MID_TIME[-1, ]) * (1 - c", k,
                                                   y_var[traj], "_mug * ", y_var[traj], "M_MID_TIME[-1, ]))"),
                                            name = paste0("c", k, y_var[traj], "r_loads"))
          status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "c", k, y_var[traj], "r_loads"),
                                              name = paste0("c", k, y_var[traj], "s_loads"))
          slp_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                              "rbind(c", k, y_var[traj], "_mean0[2, ], 0)"),
                                       name = paste0("c", k, y_var[traj], "slp_m"))
          slp_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                              "c", k, y_var[traj], "_psi0[2:3, 2:3] %*% t(",
                                              "c", k, y_var[traj], "r_loads)"),
                                       name = paste0("c", k, y_var[traj], "slp_v"))
          chg_inv_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                           name = paste0("c", k, y_var[traj], "chg_inv_m"))
          chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(c",
                                                  k, y_var[traj], "slp_v))"),
                                           name = paste0("c", k, y_var[traj], "chg_inv_v"))
          chg_bl_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                                 "rbind(c", k, y_var[traj], "_mean0[2, ], 0)"),
                                          name = paste0("c", k, y_var[traj], "chg_bl_m"))
          chg_bl_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                                 "c", k, y_var[traj], "_psi0[2:3, 2:3] %*% t(c", k,
                                                 y_var[traj], "s_loads[-1, ])"),
                                          name = paste0("c", k, y_var[traj], "chg_bl_v"))
          AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                    M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
        }
        else if (!intrinsic){
          rate_loads <- mxAlgebraFromString(paste0("c", k, y_var[traj], "_mug * exp(-c", k, y_var[traj], "_mug * ",
                                                   y_var[traj], "M_MID_TIME[-1, ])"),
                                            name = paste0("c", k, y_var[traj], "r_loads"))
          status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "c", k, y_var[traj], "r_loads"),
                                              name = paste0("c", k, y_var[traj], "s_loads"))
          slp_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                              "c", k, y_var[traj], "_mean0[2, ]"),
                                       name = paste0("c", k, y_var[traj], "slp_m"))
          slp_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                              "c", k, y_var[traj], "_psi0[2, 2] %*% t(",
                                              "c", k, y_var[traj], "r_loads)"),
                                       name = paste0("c", k, y_var[traj], "slp_v"))
          chg_inv_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                           name = paste0("c", k, y_var[traj], "chg_inv_m"))
          chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(c",
                                                  k, y_var[traj], "slp_v))"),
                                           name = paste0("c", k, y_var[traj], "chg_inv_v"))
          chg_bl_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                                 "c", k, y_var[traj], "_mean0[2, ]"),
                                          name = paste0("c", k, y_var[traj], "chg_bl_m"))
          chg_bl_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                                 "c", k, y_var[traj], "_psi0[2, 2] %*% t(c", k,
                                                 y_var[traj], "s_loads[-1, ])"),
                                          name = paste0("c", k, y_var[traj], "chg_bl_v"))
          AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                    M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
        }
      }
      else if (curveFun %in% c("Jenss-Bayley", "JB")){
        if (intrinsic){
          ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                           free = FALSE, name = "ADD1")
          rate_loads <- mxAlgebraFromString(paste0("cbind(ADD1, c", k, y_var[traj], "_mug * exp(c", k, y_var[traj], "_mug * ",
                                                   y_var[traj], "M_MID_TIME[-1, ]), c", k, y_var[traj], "_mueta2 * exp(c", k,
                                                   y_var[traj], "_mug * ", y_var[traj], "M_MID_TIME[-1, ]) * (1 + c", k,
                                                   y_var[traj], "_mug * ", y_var[traj], "M_MID_TIME[-1, ]))"),
                                            name = paste0("c", k, y_var[traj], "r_loads"))
          status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "c", k, y_var[traj], "r_loads"),
                                              name = paste0("c", k, y_var[traj], "s_loads"))
          slp_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                              "rbind(c", k, y_var[traj], "_mean0[2:3, ], 0)"),
                                       name = paste0("c", k, y_var[traj], "slp_m"))
          slp_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                              "c", k, y_var[traj], "_psi0[2:4, 2:4] %*% t(",
                                              "c", k, y_var[traj], "r_loads)"),
                                       name = paste0("c", k, y_var[traj], "slp_v"))
          chg_inv_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                           name = paste0("c", k, y_var[traj], "chg_inv_m"))
          chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(c",
                                                  k, y_var[traj], "slp_v))"),
                                           name = paste0("c", k, y_var[traj], "chg_inv_v"))
          chg_bl_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                                 "rbind(c", k, y_var[traj], "_mean0[2:3, ], 0)"),
                                          name = paste0("c", k, y_var[traj], "chg_bl_m"))
          chg_bl_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                                 "c", k, y_var[traj], "_psi0[2:4, 2:4] %*% t(c", k,
                                                 y_var[traj], "s_loads[-1, ])"),
                                          name = paste0("c", k, y_var[traj], "chg_bl_v"))
          AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                    M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
        }
        else if (!intrinsic){
          ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                           free = FALSE, name = "ADD1")
          rate_loads <- mxAlgebraFromString(paste0("cbind(ADD1, c", k, y_var[traj], "_mug * exp(c", k, y_var[traj], "_mug * ",
                                                   y_var[traj], "M_MID_TIME[-1, ]))"),
                                            name = paste0("c", k, y_var[traj], "r_loads"))
          status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "c", k, y_var[traj], "r_loads"),
                                              name = paste0("c", k, y_var[traj], "s_loads"))
          slp_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                              "c", k, y_var[traj], "_mean0[2:3, ]"),
                                       name = paste0("c", k, y_var[traj], "slp_m"))
          slp_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "r_loads %*% ",
                                              "c", k, y_var[traj], "_psi0[2:3, 2:3] %*% t(",
                                              "c", k, y_var[traj], "r_loads)"),
                                       name = paste0("c", k, y_var[traj], "slp_v"))
          chg_inv_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                           name = paste0("c", k, y_var[traj], "chg_inv_m"))
          chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(c",
                                                  k, y_var[traj], "slp_v))"),
                                           name = paste0("c", k, y_var[traj], "chg_inv_v"))
          chg_bl_m <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                                 "c", k, y_var[traj], "_mean0[2:3, ]"),
                                          name = paste0("c", k, y_var[traj], "chg_bl_m"))
          chg_bl_v <- mxAlgebraFromString(paste0("c", k, y_var[traj], "s_loads[-1, ] %*% ",
                                                 "c", k, y_var[traj], "_psi0[2:3, 2:3] %*% t(c", k,
                                                 y_var[traj], "s_loads[-1, ])"),
                                          name = paste0("c", k, y_var[traj], "chg_bl_v"))
          AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                    M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
        }
      }
    }
    AddPara_CL[[k]] <- AddPara_L
  }
  return(AddPara_CL)
}

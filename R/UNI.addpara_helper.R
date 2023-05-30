#' @title Get Additional Parameters Related to Interval-specific Slopes, Interval-specific Changes and Values of
#' Change-from-baseline for Latent Change Score Models for Longitudinal Outcome
#'
#' @description This function derives additional parameters for latent change score models. In particular, it
#' specifies the means and variances of interval-specific slopes, interval-specific changes, and values of change-
#' from-baseline.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables
#' with repeated measurements and occasions, and time-invariant covariates (TICs) if any. It takes the value passed from
#' \code{getLCSM()} or \code{getTVCmodel()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for latent change score
#' models include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}), \code{"Jenss-Bayley"}
#' (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the value passed from \code{getLCSM()} or
#' \code{getTVCmodel()}.
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. It takes the
#' value passed from \code{getLCSM()} or \code{getTVCmodel()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' It takes the value passed from \code{getLCSM()} or \code{getTVCmodel()}.
#' @param records A numeric vector specifying indices of the study waves. It takes the value passed from \code{getLCSM()}
#' or \code{getTVCmodel()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing
#' to the variability of growth factors if any. It takes the value passed from \code{getLCSM()} or \code{getTVCmodel()}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with
#' interval-specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states).
#' It takes the value passed from \code{getTVCmodel()}. The argument is set as \code{NULL} when fitting a latent change score
#' model.
#' @param starts A list containing initial values for the parameters, either takes the value passed from \code{getLCSM()} or
#' \code{getTVCmodel()}, or derived by the helper function \code{getUNI.initial()} or \code{getTVC.initial()}.

#' @return A list containing the specification of the means and variances of interval-specific slopes, interval-specific
#' changes, and values of change-from-baseline for latent change score models.
#'
#' @importFrom OpenMx mxMatrix mxAlgebra diag2vec
#'
getUNI.addpara <- function(dat, curveFun, intrinsic = NULL, t_var, records, growth_TIC, decompose, starts){
  # Define mxMatrix to include the mean vector of each measurement occasion and the mean vector of
  # each middle time point, which are useful to derive the means and variances of interval-specific
  # slopes, interval-specific changes and change from baseline.
  m_time <- sapply(paste0(t_var, records), function(x) mean(dat[, x]))
  m_mid_time <- c(0, (m_time[-1] + m_time[-length(records)])/2)
  M_TIME <- mxMatrix("Full", length(records), 1, free = FALSE, values = m_time, name = "M_TIME")
  M_MID_TIME <- mxMatrix("Full", length(records), 1, free = FALSE, values = m_mid_time,
                         name = "M_MID_TIME")
  # Define mxMatrix to include the mean value time lag matrix, which is useful to derive the means and
  # variances of interval-specific slopes, interval-specific changes and change from baseline.
  m_lag <- diff(m_time)
  m_lag_sq <- m_lag^2
  M_LAG <- mxMatrix("Full", length(records), 1, free = FALSE, values = c(0, m_lag), name = "M_LAG")
  M_LAG_SQ <- mxMatrix("Full", length(records), 1, free = FALSE, values = c(0, m_lag_sq),
                       name = "M_LAG_SQ")
  m_Omega_val <- matrix(0, nrow = length(m_time), ncol = length(m_lag))
  for (j in records[-1]){
    m_Omega_val[j, ] <- c(m_lag[1:(j - 1)], rep(0, length(m_time) - j))
  }
  # Define Omega matrix using mxMatrix
  m_Omega <- mxMatrix(type = "Full", nrow = length(m_time), ncol = length(m_lag),
                      values = m_Omega_val, free = FALSE, name = "Omega")
  if (curveFun %in% c("nonparametric", "NonP")){
    rate_loads <- mxMatrix("Full", nrow = length(m_lag), ncol = 1, c(F, rep(T, length(m_lag) - 1)),
                           values = c(1, starts[[1]][[4]][-1]),
                           labels = paste0("Y_rel_rate", 1:length(m_lag)),
                           byrow = T, name = "r_loads")
    status_loads <- mxAlgebra(Omega %*% r_loads, name = "s_loads")
    slp_m <- mxAlgebra(r_loads %*% Y_mean0[2, ], name = "Yslp_m")
    slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
      mxAlgebra(r_loads %*% Y_psi0[2, 2] %*% t(r_loads), name = "Yslp_v")
    } else {
      mxAlgebra(r_loads %*% Y_psi_r[2, 2] %*% t(r_loads), name = "Yslp_v")
    }
    chg_inv_m <- mxAlgebra(Yslp_m * M_LAG[-1, ], name = "Ychg_inv_m")
    chg_inv_v <- mxAlgebra(t(M_LAG_SQ[-1, ] * diag2vec(Yslp_v)), name = "Ychg_inv_v")
    chg_bl_m <- mxAlgebra(s_loads[-1, ] %*% Y_mean0[2, ], name = "Ychg_bl_m")
    chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
      mxAlgebra(s_loads[-1, ] %*% Y_psi0[2, 2] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
    } else {
      mxAlgebra(s_loads[-1, ] %*% Y_psi_r[2, 2] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
    }
    AddPara <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                    M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
  }
  else if (curveFun %in% c("quadratic", "QUAD")){
    ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                     free = FALSE, name = "ADD1")
    rate_loads <- mxAlgebra(cbind(ADD1, M_MID_TIME[-1, ] * 2), name = "r_loads")
    status_loads <- mxAlgebra(Omega %*% r_loads, name = "s_loads")
    slp_m <- mxAlgebra(r_loads %*% Y_mean0[2:3, ], name = "Yslp_m")
    slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
      mxAlgebra(r_loads %*% Y_psi0[2:3, 2:3] %*% t(r_loads), name = "Yslp_v")
    } else {
      mxAlgebra(r_loads %*% Y_psi_r[2:3, 2:3] %*% t(r_loads), name = "Yslp_v")
    }
    chg_inv_m <- mxAlgebra(Yslp_m * M_LAG[-1, ], name = "Ychg_inv_m")
    chg_inv_v <- mxAlgebra(t(M_LAG_SQ[-1, ] * diag2vec(Yslp_v)), name = "Ychg_inv_v")
    chg_bl_m <- mxAlgebra(s_loads[-1, ] %*% Y_mean0[2:3, ], name = "Ychg_bl_m")
    chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
      mxAlgebra(s_loads[-1, ] %*% Y_psi0[2:3, 2:3] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
    } else {
      mxAlgebra(s_loads[-1, ] %*% Y_psi_r[2:3, 2:3] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
    }
    AddPara <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                    M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
  }
  else if (curveFun %in% c("negative exponential", "EXP")){
    if (intrinsic){
      rate_loads <- mxAlgebra(cbind(Y_mug * exp(-Y_mug * M_MID_TIME[-1, ]),
                                    Y_mueta1 * exp(-Y_mug * M_MID_TIME[-1, ]) * (1 - Y_mug * M_MID_TIME[-1, ])),
                              name = "r_loads")
      status_loads <- mxAlgebra(Omega %*% r_loads, name = "s_loads")
      slp_m <- mxAlgebra(r_loads %*% rbind(Y_mean0[2, ], 0), name = "Yslp_m")
      slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebra(r_loads %*% Y_psi0[2:3, 2:3] %*% t(r_loads), name = "Yslp_v")
      } else {
        mxAlgebra(r_loads %*% Y_psi_r[2:3, 2:3] %*% t(r_loads), name = "Yslp_v")
      }
      chg_inv_m <- mxAlgebra(Yslp_m * M_LAG[-1, ], name = "Ychg_inv_m")
      chg_inv_v <- mxAlgebra(t(M_LAG_SQ[-1, ] * diag2vec(Yslp_v)), name = "Ychg_inv_v")
      chg_bl_m <- mxAlgebra(s_loads[-1, ] %*% rbind(Y_mean0[2, ], 0), name = "Ychg_bl_m")
      chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebra(s_loads[-1, ] %*% Y_psi0[2:3, 2:3] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
      } else {
        mxAlgebra(s_loads[-1, ] %*% Y_psi_r[2:3, 2:3] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
      }
      AddPara <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                      M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
    }
    else if (!intrinsic){
      rate_loads <- mxAlgebra(Y_mug * exp(-Y_mug * M_MID_TIME[-1, ]), name = "r_loads")
      status_loads <- mxAlgebra(Omega %*% r_loads, name = "s_loads")
      slp_m <- mxAlgebra(r_loads %*% Y_mean0[2, ], name = "Yslp_m")
      slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebra(r_loads %*% Y_psi0[2, 2] %*% t(r_loads), name = "Yslp_v")
      } else {
        mxAlgebra(r_loads %*% Y_psi_r[2, 2] %*% t(r_loads), name = "Yslp_v")
      }
      chg_inv_m <- mxAlgebra(Yslp_m * M_LAG[-1, ], name = "Ychg_inv_m")
      chg_inv_v <- mxAlgebra(t(M_LAG_SQ[-1, ] * diag2vec(Yslp_v)), name = "Ychg_inv_v")
      chg_bl_m <- mxAlgebra(s_loads[-1, ] %*% Y_mean0[2, ], name = "Ychg_bl_m")
      chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebra(s_loads[-1, ] %*% Y_psi0[2, 2] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
      } else {
        mxAlgebra(s_loads[-1, ] %*% Y_psi_r[2, 2] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
      }
      AddPara <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                      M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
    }
  }
  else if (curveFun %in% c("Jenss-Bayley", "JB")){
    if (intrinsic){
      ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                       free = FALSE, name = "ADD1")
      rate_loads <- mxAlgebra(cbind(ADD1, Y_mug * exp(Y_mug * M_MID_TIME[-1, ]),
                                    Y_mueta2 * exp(Y_mug * M_MID_TIME[-1, ]) * (1 + Y_mug * M_MID_TIME[-1, ])),
                              name = "r_loads")
      status_loads <- mxAlgebra(Omega %*% r_loads, name = "s_loads")
      slp_m <- mxAlgebra(r_loads %*% rbind(Y_mean0[2:3, ], 0), name = "Yslp_m")
      slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebra(r_loads %*% Y_psi0[2:4, 2:4] %*% t(r_loads), name = "Yslp_v")
      } else {
        mxAlgebra(r_loads %*% Y_psi_r[2:4, 2:4] %*% t(r_loads), name = "Yslp_v")
      }
      chg_inv_m <- mxAlgebra(Yslp_m * M_LAG[-1, ], name = "Ychg_inv_m")
      chg_inv_v <- mxAlgebra(t(M_LAG_SQ[-1, ] * diag2vec(Yslp_v)), name = "Ychg_inv_v")
      chg_bl_m <- mxAlgebra(s_loads[-1, ] %*% rbind(Y_mean0[2:3, ], 0), name = "Ychg_bl_m")
      chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebra(s_loads[-1, ] %*% Y_psi0[2:4, 2:4] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
      } else {
        mxAlgebra(s_loads[-1, ] %*% Y_psi_r[2:4, 2:4] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
      }
      AddPara <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                      M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
    }
    else if (!intrinsic){
      ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                       free = FALSE, name = "ADD1")
      rate_loads <- mxAlgebra(cbind(ADD1, Y_mug * exp(Y_mug * M_MID_TIME[-1, ])), name = "r_loads")
      status_loads <- mxAlgebra(Omega %*% r_loads, name = "s_loads")
      slp_m <- mxAlgebra(r_loads %*% Y_mean0[2:3, ], name = "Yslp_m")
      slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebra(r_loads %*% Y_psi0[2:3, 2:3] %*% t(r_loads), name = "Yslp_v")
      } else {
        mxAlgebra(r_loads %*% Y_psi_r[2:3, 2:3] %*% t(r_loads), name = "Yslp_v")
      }
      chg_inv_m <- mxAlgebra(Yslp_m * M_LAG[-1, ], name = "Ychg_inv_m")
      chg_inv_v <- mxAlgebra(t(M_LAG_SQ[-1, ] * diag2vec(Yslp_v)), name = "Ychg_inv_v")
      chg_bl_m <- mxAlgebra(s_loads[-1, ] %*% Y_mean0[2:3, ], name = "Ychg_bl_m")
      chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebra(s_loads[-1, ] %*% Y_psi0[2:3, 2:3] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
      } else {
        mxAlgebra(s_loads[-1, ] %*% Y_psi_r[2:3, 2:3] %*% t(s_loads[-1, ]), name = "Ychg_bl_v")
      }
      AddPara <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                      M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
    }
  }
  return(AddPara)
}

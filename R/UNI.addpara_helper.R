#' @title Derived change-score quantities of a univariate outcome, one set per prefix
#' @description \code{pfx = ""} gives the single-group objects; \code{pfx = c("c1", "c2", ...)} the class-specific
#' ones used by the multiple-group and mixture models. The names of matrices and algebras carry the prefix.
#' @keywords internal
#' @noRd
#'
#' @importFrom OpenMx mxMatrix mxAlgebra diag2vec
#' @importFrom OpenMx mxMatrix mxAlgebra mxAlgebraFromString diag2vec
#'
.addpara_UNI_classes <- function(pfx, dat, curveFun, intrinsic, t_var, records, growth_TIC, decompose, starts,
                                  grp_var = NULL){
  output_CL <- list()
  for (k in seq_along(pfx)){
    # Manifest groups use their own observed occasions. Latent mixture classes share
    # the pooled reference grid because their group membership is not observed.
    time_dat <- if (is.null(grp_var)) dat else dat[dat[[grp_var]] == k, , drop = FALSE]
    # Define mxMatrix to include the mean vector of each measurement occasion and the mean vector of
    # each middle time point, which are useful to derive the means and variances of interval-specific
    # slopes, interval-specific changes and change from baseline.
    m_time <- sapply(paste0(t_var, records), function(x) mean(time_dat[[x]]))
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
    for (j in seq_along(records)[-1]){
      m_Omega_val[j, seq_len(j - 1)] <- m_lag[seq_len(j - 1)]
    }
    # Define Omega matrix using mxMatrix
    m_Omega <- mxMatrix(type = "Full", nrow = length(m_time), ncol = length(m_lag),
                        values = m_Omega_val, free = FALSE, name = "Omega")
    if (curveFun %in% c("nonparametric", "NonP")){
      rate_loads <- mxMatrix("Full", nrow = length(m_lag), ncol = 1, c(FALSE, rep(TRUE, length(m_lag) - 1)),
                             values = c(1, starts[[k]]$Y_starts$rel_rate[-1]),
                             labels = paste0(pfx[k], "Y_rel_rate", seq_along(m_lag)),
                             byrow = TRUE, name = paste0(pfx[k], "r_loads"))
      status_loads <- mxAlgebraFromString(paste0("Omega %*% ", pfx[k], "r_loads"),
                                          name = paste0(pfx[k], "s_loads"))
      slp_m <- mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_mean0[2, ]"),
                                   name = paste0(pfx[k], "Yslp_m"))
      slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi0[2, 2] %*% t(", pfx[k], "r_loads)"),
                            name = paste0(pfx[k], "Yslp_v"))
      } else {
        mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi_r[2, 2] %*% t(", pfx[k], "r_loads)"),
                            name = paste0(pfx[k], "Yslp_v"))
      }
      chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], "Yslp_m * M_LAG[-1, ]"),
                                       name = paste0(pfx[k], "Ychg_inv_m"))
      chg_inv_v <- mxAlgebraFromString(paste0("t(M_LAG_SQ[-1, ] * diag2vec(", pfx[k], "Yslp_v))"),
                                       name = paste0(pfx[k], "Ychg_inv_v"))
      chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_mean0[2, ]"),
                                      name = paste0(pfx[k], "Ychg_bl_m"))
      chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi0[2, 2] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
      } else {
        mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi_r[2, 2] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
      }
      output_CL[[k]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                             M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
    }
    else if (curveFun %in% c("quadratic", "QUAD")){
      ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                       free = FALSE, name = "ADD1")
      rate_loads <- mxAlgebra(cbind(ADD1, M_MID_TIME[-1, ] * 2), name = paste0(pfx[k], "r_loads"))
      status_loads <- mxAlgebraFromString(paste0("Omega %*% ", pfx[k], "r_loads"),
                                          name = paste0(pfx[k], "s_loads"))
      slp_m <- mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_mean0[2:3, ]"),
                                   name = paste0(pfx[k], "Yslp_m"))
      slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi0[2:3, 2:3] %*% t(", pfx[k], "r_loads)"),
                            name = paste0(pfx[k], "Yslp_v"))
      } else {
        mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi_r[2:3, 2:3] %*% t(", pfx[k], "r_loads)"),
                            name = paste0(pfx[k], "Yslp_v"))
      }
      chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], "Yslp_m * M_LAG[-1, ]"),
                                       name = paste0(pfx[k], "Ychg_inv_m"))
      chg_inv_v <- mxAlgebraFromString(paste0("t(M_LAG_SQ[-1, ] * diag2vec(", pfx[k], "Yslp_v))"),
                                       name = paste0(pfx[k], "Ychg_inv_v"))
      chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_mean0[2:3, ]"),
                                      name = paste0(pfx[k], "Ychg_bl_m"))
      chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
        mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi0[2:3, 2:3] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
      } else {
        mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi_r[2:3, 2:3] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
      }
      output_CL[[k]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                             M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
    }
    else if (curveFun %in% c("negative exponential", "EXP")){
      if (intrinsic){
        rate_loads <- mxAlgebraFromString(paste0("cbind(", pfx[k], "Y_mug * exp(-", pfx[k], "Y_mug * M_MID_TIME[-1, ]), ",
                                                 "", pfx[k], "Y_mueta1 * exp(-", pfx[k], "Y_mug * M_MID_TIME[-1, ]) * ",
                                                 "(1 - ", pfx[k], "Y_mug * M_MID_TIME[-1, ]))"),
                                          name = paste0(pfx[k], "r_loads"))
        status_loads <- mxAlgebraFromString(paste0("Omega %*% ", pfx[k], "r_loads"),
                                            name = paste0(pfx[k], "s_loads"))
        slp_m <- mxAlgebraFromString(paste0(pfx[k], "r_loads %*% rbind(", pfx[k], "Y_mean0[2, ], 0)"),
                                     name = paste0(pfx[k], "Yslp_m"))
        slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
          mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi0[2:3, 2:3] %*% t(", pfx[k], "r_loads)"),
                              name = paste0(pfx[k], "Yslp_v"))
        } else {
          mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi_r[2:3, 2:3] %*% t(", pfx[k], "r_loads)"),
                              name = paste0(pfx[k], "Yslp_v"))
        }
        chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], "Yslp_m * M_LAG[-1, ]"),
                                         name = paste0(pfx[k], "Ychg_inv_m"))
        chg_inv_v <- mxAlgebraFromString(paste0("t(M_LAG_SQ[-1, ] * diag2vec(", pfx[k], "Yslp_v))"),
                                         name = paste0(pfx[k], "Ychg_inv_v"))
        chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% rbind(", pfx[k], "Y_mean0[2, ], 0)"),
                                        name = paste0(pfx[k], "Ychg_bl_m"))
        chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
          mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi0[2:3, 2:3] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
        } else {
          mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi_r[2:3, 2:3] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
        }
        output_CL[[k]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                               M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
      }
      else if (!intrinsic){
        rate_loads <- mxAlgebraFromString(paste0(pfx[k], "Y_mug * exp(-", pfx[k], "Y_mug * M_MID_TIME[-1, ])"),
                                          name = paste0(pfx[k], "r_loads"))
        status_loads <- mxAlgebraFromString(paste0("Omega %*% ", pfx[k], "r_loads"),
                                            name = paste0(pfx[k], "s_loads"))
        slp_m <- mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_mean0[2, ]"),
                                     name = paste0(pfx[k], "Yslp_m"))
        slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
          mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi0[2, 2] %*% t(", pfx[k], "r_loads)"),
                              name = paste0(pfx[k], "Yslp_v"))
        } else {
          mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi_r[2, 2] %*% t(", pfx[k], "r_loads)"),
                              name = paste0(pfx[k], "Yslp_v"))
        }
        chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], "Yslp_m * M_LAG[-1, ]"),
                                         name = paste0(pfx[k], "Ychg_inv_m"))
        chg_inv_v <- mxAlgebraFromString(paste0("t(M_LAG_SQ[-1, ] * diag2vec(", pfx[k], "Yslp_v))"),
                                         name = paste0(pfx[k], "Ychg_inv_v"))
        chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_mean0[2, ]"),
                                        name = paste0(pfx[k], "Ychg_bl_m"))
        chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
          mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi0[2, 2] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
        } else {
          mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi_r[2, 2] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
        }
        output_CL[[k]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                               M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
      }
    }
    else if (curveFun %in% c("Jenss-Bayley", "JB")){
      if (intrinsic){
        ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                         free = FALSE, name = "ADD1")
        rate_loads <- mxAlgebraFromString(paste0("cbind(ADD1, ", pfx[k], "Y_mug * exp(", pfx[k], "Y_mug * M_MID_TIME[-1, ]), ",
                                                 "", pfx[k], "Y_mueta2 * exp(", pfx[k], "Y_mug * M_MID_TIME[-1, ]) * ",
                                                 "(1 + ", pfx[k], "Y_mug * M_MID_TIME[-1, ]))"),
                                          name = paste0(pfx[k], "r_loads"))
        status_loads <- mxAlgebraFromString(paste0("Omega %*% ", pfx[k], "r_loads"),
                                            name = paste0(pfx[k], "s_loads"))
        slp_m <- mxAlgebraFromString(paste0(pfx[k], "r_loads %*% rbind(", pfx[k], "Y_mean0[2:3, ], 0)"),
                                     name = paste0(pfx[k], "Yslp_m"))
        slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
          mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi0[2:4, 2:4] %*% t(", pfx[k], "r_loads)"),
                              name = paste0(pfx[k], "Yslp_v"))
        } else {
          mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi_r[2:4, 2:4] %*% t(", pfx[k], "r_loads)"),
                              name = paste0(pfx[k], "Yslp_v"))
        }
        chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], "Yslp_m * M_LAG[-1, ]"),
                                         name = paste0(pfx[k], "Ychg_inv_m"))
        chg_inv_v <- mxAlgebraFromString(paste0("t(M_LAG_SQ[-1, ] * diag2vec(", pfx[k], "Yslp_v))"),
                                         name = paste0(pfx[k], "Ychg_inv_v"))
        chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% rbind(", pfx[k], "Y_mean0[2:3, ], 0)"),
                                        name = paste0(pfx[k], "Ychg_bl_m"))
        chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
          mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi0[2:4, 2:4] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
        } else {
          mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi_r[2:4, 2:4] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
        }
        output_CL[[k]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                               M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
      }
      else if (!intrinsic){
        ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                         free = FALSE, name = "ADD1")
        rate_loads <- mxAlgebraFromString(paste0("cbind(ADD1, ", pfx[k], "Y_mug * exp(", pfx[k], "Y_mug * M_MID_TIME[-1, ]))"),
                                          name = paste0(pfx[k], "r_loads"))
        status_loads <- mxAlgebraFromString(paste0("Omega %*% ", pfx[k], "r_loads"),
                                            name = paste0(pfx[k], "s_loads"))
        slp_m <- mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_mean0[2:3, ]"),
                                     name = paste0(pfx[k], "Yslp_m"))
        slp_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
          mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi0[2:3, 2:3] %*% t(", pfx[k], "r_loads)"),
                              name = paste0(pfx[k], "Yslp_v"))
        } else {
          mxAlgebraFromString(paste0(pfx[k], "r_loads %*% ", pfx[k], "Y_psi_r[2:3, 2:3] %*% t(", pfx[k], "r_loads)"),
                              name = paste0(pfx[k], "Yslp_v"))
        }
        chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], "Yslp_m * M_LAG[-1, ]"),
                                         name = paste0(pfx[k], "Ychg_inv_m"))
        chg_inv_v <- mxAlgebraFromString(paste0("t(M_LAG_SQ[-1, ] * diag2vec(", pfx[k], "Yslp_v))"),
                                         name = paste0(pfx[k], "Ychg_inv_v"))
        chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_mean0[2:3, ]"),
                                        name = paste0(pfx[k], "Ychg_bl_m"))
        chg_bl_v <- if ((is.null(decompose)||decompose == 0) && is.null(growth_TIC)) {
          mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi0[2:3, 2:3] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
        } else {
          mxAlgebraFromString(paste0(pfx[k], "s_loads[-1, ] %*% ", pfx[k], "Y_psi_r[2:3, 2:3] %*% t(", pfx[k], "s_loads[-1, ])"), name = paste0(pfx[k], "Ychg_bl_v"))
        }
        output_CL[[k]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                               M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
      }
    }
  }
  return(output_CL)
}

#' @keywords internal
#' @noRd
getUNI.addpara <- function(dat, curveFun, intrinsic = NULL, t_var, records, growth_TIC, decompose, starts){
  .addpara_UNI_classes(pfx = "", dat = dat, curveFun = curveFun, intrinsic = intrinsic, t_var = t_var,
                       records = records, growth_TIC = growth_TIC, decompose = decompose, starts = list(starts))[[1]]
}
getMIX_UNI.addpara <- function(dat, nClass, curveFun, intrinsic, t_var, records, growth_TIC, decompose, starts,
                                grp_var = NULL){
  .addpara_UNI_classes(pfx = paste0("c", seq_len(nClass)), dat = dat, curveFun = curveFun, intrinsic = intrinsic,
                       t_var = t_var, records = records, growth_TIC = growth_TIC, decompose = decompose, starts = starts,
                       grp_var = grp_var)
}

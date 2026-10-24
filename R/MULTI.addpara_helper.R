#' @title Derived change-score quantities of multivariate outcomes, one set per prefix
#' @description \code{pfx = ""} gives the single-group objects; \code{pfx = c("c1", "c2", ...)} the class-specific
#' ones used by the multiple-group and mixture models. The names of matrices and algebras carry the prefix.
#' @keywords internal
#' @noRd
#'
#' @importFrom OpenMx mxMatrix mxAlgebra mxAlgebraFromString diag2vec
#'
.addpara_MULTI_classes <- function(pfx, dat, t_var, records, y_var, curveFun, intrinsic = NULL, starts,
                                    grp_var = NULL){
  AddPara_CL <- AddPara_L <- list()
  for (k in seq_along(pfx)){
    # Group-specific summaries use observed group occasions; latent mixtures retain
    # a common pooled time grid for comparisons between classes.
    time_dat <- if (is.null(grp_var)) dat else dat[dat[[grp_var]] == k, , drop = FALSE]
    for (traj in seq_along(y_var)){
      # Define mxMatrix to include the mean vector of each measurement occasion and the mean vector of
      # each middle time point, which are useful to derive the means and variances of interval-specific
      # slopes, interval-specific changes and change from baseline.
      m_time <- sapply(paste0(t_var[traj], records[[traj]]), function(x) mean(time_dat[[x]]))
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
      for (j in seq_along(records[[traj]])[-1]){
        m_Omega_val[j, seq_len(j - 1)] <- m_lag[seq_len(j - 1)]
      }
      # Define Omega matrix using mxMatrix
      m_Omega <- mxMatrix(type = "Full", nrow = length(m_time), ncol = length(m_lag),
                          values = m_Omega_val, free = FALSE, name = paste0(y_var[traj], "Omega"))
      if (curveFun %in% c("nonparametric", "NonP")){
        rate_loads <- mxMatrix("Full", nrow = length(m_lag), ncol = 1, c(FALSE, rep(TRUE, length(m_lag) - 1)),
                               values = c(1, starts[[k]]$rel_rate[[traj]][-1]),
                               labels = paste0(pfx[k], y_var[traj], "_rel_rate", seq_along(m_lag)),
                               byrow = TRUE, name = paste0(pfx[k], y_var[traj], "r_loads"))
        status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "", pfx[k], y_var[traj], "r_loads"),
                                            name = paste0(pfx[k], y_var[traj], "s_loads"))
        slp_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                            "", pfx[k], y_var[traj], "_mean0[2, ]"),
                                     name = paste0(pfx[k], y_var[traj], "slp_m"))
        slp_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                            "", pfx[k], y_var[traj], "_psi0[2, 2] %*% t(",
                                            "", pfx[k], y_var[traj], "r_loads)"),
                                     name = paste0(pfx[k], y_var[traj], "slp_v"))
        chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                         name = paste0(pfx[k], y_var[traj], "chg_inv_m"))
        chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(", pfx[k], y_var[traj], "slp_v))"),
                                         name = paste0(pfx[k], y_var[traj], "chg_inv_v"))
        chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                               "", pfx[k], y_var[traj], "_mean0[2, ]"),
                                        name = paste0(pfx[k], y_var[traj], "chg_bl_m"))
        chg_bl_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                               "", pfx[k], y_var[traj], "_psi0[2, 2] %*% t(", pfx[k], y_var[traj], "s_loads[-1, ])"),
                                        name = paste0(pfx[k], y_var[traj], "chg_bl_v"))
        AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                  M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
      }
      else if (curveFun %in% c("quadratic", "QUAD")){
        ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                         free = FALSE, name = paste0(pfx[k], y_var[traj], "ADD1"))
        rate_loads <- mxAlgebraFromString(paste0("cbind(", pfx[k], y_var[traj], "ADD1, ", y_var[traj], "M_MID_TIME[-1, ] * 2)"),
                                          name = paste0(pfx[k], y_var[traj], "r_loads"))
        status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "", pfx[k], y_var[traj], "r_loads"),
                                            name = paste0(pfx[k], y_var[traj], "s_loads"))
        slp_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                            "", pfx[k], y_var[traj], "_mean0[2:3, ]"),
                                     name = paste0(pfx[k], y_var[traj], "slp_m"))
        slp_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                            "", pfx[k], y_var[traj], "_psi0[2:3, 2:3] %*% t(",
                                            "", pfx[k], y_var[traj], "r_loads)"),
                                     name = paste0(pfx[k], y_var[traj], "slp_v"))
        chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                         name = paste0(pfx[k], y_var[traj], "chg_inv_m"))
        chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(", pfx[k], y_var[traj], "slp_v))"),
                                         name = paste0(pfx[k], y_var[traj], "chg_inv_v"))
        chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                               "", pfx[k], y_var[traj], "_mean0[2:3, ]"),
                                        name = paste0(pfx[k], y_var[traj], "chg_bl_m"))
        chg_bl_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                               "", pfx[k], y_var[traj], "_psi0[2:3, 2:3] %*% t(", pfx[k], y_var[traj], "s_loads[-1, ])"),
                                        name = paste0(pfx[k], y_var[traj], "chg_bl_v"))
        AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                  M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
      }
      else if (curveFun %in% c("negative exponential", "EXP")){
        if (intrinsic){
          rate_loads <- mxAlgebraFromString(paste0("cbind(", pfx[k], y_var[traj], "_mug * exp(-", pfx[k], y_var[traj], "_mug * ",
                                                   y_var[traj], "M_MID_TIME[-1, ]), ", pfx[k], y_var[traj], "_mueta1 * exp(-", pfx[k], y_var[traj], "_mug * ", y_var[traj], "M_MID_TIME[-1, ]) * (1 - ", pfx[k], y_var[traj], "_mug * ", y_var[traj], "M_MID_TIME[-1, ]))"),
                                            name = paste0(pfx[k], y_var[traj], "r_loads"))
          status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "", pfx[k], y_var[traj], "r_loads"),
                                              name = paste0(pfx[k], y_var[traj], "s_loads"))
          slp_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                              "rbind(", pfx[k], y_var[traj], "_mean0[2, ], 0)"),
                                       name = paste0(pfx[k], y_var[traj], "slp_m"))
          slp_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                              "", pfx[k], y_var[traj], "_psi0[2:3, 2:3] %*% t(",
                                              "", pfx[k], y_var[traj], "r_loads)"),
                                       name = paste0(pfx[k], y_var[traj], "slp_v"))
          chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                           name = paste0(pfx[k], y_var[traj], "chg_inv_m"))
          chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(", pfx[k], y_var[traj], "slp_v))"),
                                           name = paste0(pfx[k], y_var[traj], "chg_inv_v"))
          chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                                 "rbind(", pfx[k], y_var[traj], "_mean0[2, ], 0)"),
                                          name = paste0(pfx[k], y_var[traj], "chg_bl_m"))
          chg_bl_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                                 "", pfx[k], y_var[traj], "_psi0[2:3, 2:3] %*% t(", pfx[k], y_var[traj], "s_loads[-1, ])"),
                                          name = paste0(pfx[k], y_var[traj], "chg_bl_v"))
          AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                    M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
        }
        else if (!intrinsic){
          rate_loads <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "_mug * exp(-", pfx[k], y_var[traj], "_mug * ",
                                                   y_var[traj], "M_MID_TIME[-1, ])"),
                                            name = paste0(pfx[k], y_var[traj], "r_loads"))
          status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "", pfx[k], y_var[traj], "r_loads"),
                                              name = paste0(pfx[k], y_var[traj], "s_loads"))
          slp_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                              "", pfx[k], y_var[traj], "_mean0[2, ]"),
                                       name = paste0(pfx[k], y_var[traj], "slp_m"))
          slp_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                              "", pfx[k], y_var[traj], "_psi0[2, 2] %*% t(",
                                              "", pfx[k], y_var[traj], "r_loads)"),
                                       name = paste0(pfx[k], y_var[traj], "slp_v"))
          chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                           name = paste0(pfx[k], y_var[traj], "chg_inv_m"))
          chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(", pfx[k], y_var[traj], "slp_v))"),
                                           name = paste0(pfx[k], y_var[traj], "chg_inv_v"))
          chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                                 "", pfx[k], y_var[traj], "_mean0[2, ]"),
                                          name = paste0(pfx[k], y_var[traj], "chg_bl_m"))
          chg_bl_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                                 "", pfx[k], y_var[traj], "_psi0[2, 2] %*% t(", pfx[k], y_var[traj], "s_loads[-1, ])"),
                                          name = paste0(pfx[k], y_var[traj], "chg_bl_v"))
          AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                    M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads)
        }
      }
      else if (curveFun %in% c("Jenss-Bayley", "JB")){
        if (intrinsic){
          ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                           free = FALSE, name = paste0(pfx[k], y_var[traj], "ADD1"))
          rate_loads <- mxAlgebraFromString(paste0("cbind(", pfx[k], y_var[traj], "ADD1, ", pfx[k], y_var[traj], "_mug * exp(", pfx[k], y_var[traj], "_mug * ",
                                                   y_var[traj], "M_MID_TIME[-1, ]), ", pfx[k], y_var[traj], "_mueta2 * exp(", pfx[k], y_var[traj], "_mug * ", y_var[traj], "M_MID_TIME[-1, ]) * (1 + ", pfx[k], y_var[traj], "_mug * ", y_var[traj], "M_MID_TIME[-1, ]))"),
                                            name = paste0(pfx[k], y_var[traj], "r_loads"))
          status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "", pfx[k], y_var[traj], "r_loads"),
                                              name = paste0(pfx[k], y_var[traj], "s_loads"))
          slp_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                              "rbind(", pfx[k], y_var[traj], "_mean0[2:3, ], 0)"),
                                       name = paste0(pfx[k], y_var[traj], "slp_m"))
          slp_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                              "", pfx[k], y_var[traj], "_psi0[2:4, 2:4] %*% t(",
                                              "", pfx[k], y_var[traj], "r_loads)"),
                                       name = paste0(pfx[k], y_var[traj], "slp_v"))
          chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                           name = paste0(pfx[k], y_var[traj], "chg_inv_m"))
          chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(", pfx[k], y_var[traj], "slp_v))"),
                                           name = paste0(pfx[k], y_var[traj], "chg_inv_v"))
          chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                                 "rbind(", pfx[k], y_var[traj], "_mean0[2:3, ], 0)"),
                                          name = paste0(pfx[k], y_var[traj], "chg_bl_m"))
          chg_bl_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                                 "", pfx[k], y_var[traj], "_psi0[2:4, 2:4] %*% t(", pfx[k], y_var[traj], "s_loads[-1, ])"),
                                          name = paste0(pfx[k], y_var[traj], "chg_bl_v"))
          AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                    M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
        }
        else if (!intrinsic){
          ADD1 <- mxMatrix(type = "Full", nrow = length(m_lag), ncol = 1, values = rep(1, length(m_lag)),
                           free = FALSE, name = paste0(pfx[k], y_var[traj], "ADD1"))
          rate_loads <- mxAlgebraFromString(paste0("cbind(", pfx[k], y_var[traj], "ADD1, ", pfx[k], y_var[traj], "_mug * exp(", pfx[k], y_var[traj], "_mug * ",
                                                   y_var[traj], "M_MID_TIME[-1, ]))"),
                                            name = paste0(pfx[k], y_var[traj], "r_loads"))
          status_loads <- mxAlgebraFromString(paste0(y_var[traj], "Omega %*% ", "", pfx[k], y_var[traj], "r_loads"),
                                              name = paste0(pfx[k], y_var[traj], "s_loads"))
          slp_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                              "", pfx[k], y_var[traj], "_mean0[2:3, ]"),
                                       name = paste0(pfx[k], y_var[traj], "slp_m"))
          slp_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "r_loads %*% ",
                                              "", pfx[k], y_var[traj], "_psi0[2:3, 2:3] %*% t(",
                                              "", pfx[k], y_var[traj], "r_loads)"),
                                       name = paste0(pfx[k], y_var[traj], "slp_v"))
          chg_inv_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "slp_m * ", y_var[traj], "M_LAG[-1, ]"),
                                           name = paste0(pfx[k], y_var[traj], "chg_inv_m"))
          chg_inv_v <- mxAlgebraFromString(paste0("t(", y_var[traj], "M_LAG_SQ[-1, ] * diag2vec(", pfx[k], y_var[traj], "slp_v))"),
                                           name = paste0(pfx[k], y_var[traj], "chg_inv_v"))
          chg_bl_m <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                                 "", pfx[k], y_var[traj], "_mean0[2:3, ]"),
                                          name = paste0(pfx[k], y_var[traj], "chg_bl_m"))
          chg_bl_v <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "s_loads[-1, ] %*% ",
                                                 "", pfx[k], y_var[traj], "_psi0[2:3, 2:3] %*% t(", pfx[k], y_var[traj], "s_loads[-1, ])"),
                                          name = paste0(pfx[k], y_var[traj], "chg_bl_v"))
          AddPara_L[[traj]] <- list(slp_m, slp_v, chg_inv_m, chg_inv_v, chg_bl_m, chg_bl_v,
                                    M_TIME, M_MID_TIME, M_LAG, M_LAG_SQ, m_Omega, rate_loads, status_loads, ADD1)
        }
      }
    }
    AddPara_CL[[k]] <- AddPara_L
  }
  return(AddPara_CL)
}

#' @keywords internal
#' @noRd
getMULTI.addpara <- function(dat, t_var, y_var, curveFun, intrinsic, records, starts){
  .addpara_MULTI_classes(pfx = "", dat = dat, t_var = t_var, records = records, y_var = y_var, curveFun = curveFun,
                         intrinsic = intrinsic, starts = list(starts))[[1]]
}
getMIX_MULTI.addpara <- function(dat, nClass, t_var, records, y_var, curveFun, intrinsic = NULL, starts,
                                  grp_var = NULL){
  .addpara_MULTI_classes(pfx = paste0("c", seq_len(nClass)), dat = dat, t_var = t_var, records = records, y_var = y_var,
                         curveFun = curveFun, intrinsic = intrinsic, starts = starts, grp_var = grp_var)
}

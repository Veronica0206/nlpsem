#' @title Define Longitudinal Mediation Models as Class-specific Models (Submodels) for a Longitudinal Mixture Model
#'
#' @description This function defines longitudinal mediation models as class-specific models (submodels) for a longitudinal mixture model.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions for each longitudinal process, and time-invariant covariates (TICs) if any. It takes
#' the value passed from \code{getMIX()}.
#' @param nClass An integer specifying the number of latent classes for the mixture model. It takes the value passed
#' from \code{getMIX()}.
#' @param t_var A vector of strings, with each element representing the prefix for column names related to the time variable for the
#' corresponding longitudinal variable at each study wave. It takes the value passed from \code{getMIX()}.
#' @param records A list of numeric vectors, with each vector specifying the indices of the observed study waves for
#' the corresponding longitudinal variable. It takes the value passed from \code{getMIX()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave. It takes the value
#' passed from \code{getMIX()}.
#' @param m_var A string specifying the prefix of the column names corresponding to the mediator variable at each study wave.
#' It takes the value passed from \code{getMIX()}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. It takes the value passed from \code{getMIX()}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. It takes the value passed from \code{getMIX()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options include: "linear" (or "LIN"),
#' and "bilinear spline" (or "BLS"). It takes the value passed from \code{getMIX()}.
#' @param starts A list of initial values for the parameters, either takes the value passed from \code{getMIX()} or derived by the
#' helper function \code{getMIX.initial()}.
#' @param res_cor A numeric value or vector for user-specified residual correlation between any two longitudinal processes to calculate
#' the corresponding initial value. It takes the value passed from \code{getMIX()}.
#'
#' @return A list of manifest and latent variables and paths for an mxModel object.
#'
#' @keywords internal
#'
#' @importFrom OpenMx mxPath mxModel mxAlgebraFromString mxMatrix mxFitFunctionML
#'
getsub.MED_l <- function(dat, nClass, t_var, records, y_var, curveFun, m_var, x_var, x_type, starts, res_cor){
  RES_L <- COV_L <- list()
  for (k in 1:nClass){
    RES <- COV <- list()
    if (x_type == "baseline"){
      y_records <- records[[1]]
      m_records <- records[[2]]
      traj_var <- c(y_var, m_var)
      traj_list <- list()
      for (traj in 1:length(traj_var)){
        traj_list[[length(traj_list) + 1]] <- paste0(traj_var[traj], records[[traj]])
      }
      manifests <- c(unlist(traj_list), x_var)
      var0 <- c(starts[[k]][[3]][[5]], starts[[k]][[2]][[4]])
      var_1 <- c(starts[[k]][[3]][[5]])
      var_2 <- c(starts[[k]][[2]][[4]])
    }
    else if (x_type == "longitudinal"){
      y_records <- records[[1]]
      m_records <- records[[2]]
      x_records <- records[[3]]
      traj_var <- c(y_var, m_var, x_var)
      traj_list <- list()
      for (traj in 1:length(traj_var)){
        traj_list[[length(traj_list) + 1]] <- paste0(traj_var[traj], records[[traj]])
      }
      manifests <- unlist(traj_list)
      var0 <- c(starts[[k]][[3]][[5]], starts[[k]][[2]][[4]], starts[[k]][[1]][[3]])
      var_1 <- c(starts[[k]][[3]][[5]], starts[[k]][[3]][[5]], starts[[k]][[2]][[4]])
      var_2 <- c(starts[[k]][[2]][[4]], starts[[k]][[1]][[3]], starts[[k]][[1]][[3]])
    }
    for (traj in 1:length(traj_var)){
      RES[[length(RES) + 1]] <- mxPath(from = traj_list[[traj]], to = traj_list[[traj]], arrows = 2, free = TRUE, values = var0[traj],
                                       labels = paste0("c", k, traj_var[traj], "_residuals"))
    }
    for (traj_i in 1:(length(traj_var) - 1)){
      for (traj_j in traj_i:(length(traj_var) - 1)){
        if (setequal(readr::parse_number(traj_list[[traj_i]]), readr::parse_number(traj_list[[traj_j + 1]]))){
          COV[[length(COV) + 1]] <- mxPath(from = traj_list[[traj_i]], to = traj_list[[traj_j + 1]],
                                           arrows = 2, free = TRUE, values = res_cor[[k]][traj_i + traj_j - 1] * sqrt(var_1[traj_i] * var_2[traj_j]),
                                           labels = paste0("c", k, traj_var[traj_i], traj_var[traj_j + 1], "_RES"))
        }
        else{
          T_common <- Reduce(intersect, list(readr::parse_number(traj_list[[traj_i]]), readr::parse_number(traj_list[[traj_j + 1]])))
          COV[[length(COV) + 1]] <- mxPath(from = paste0(traj_var[traj_i], T_common),
                                           to = paste0(traj_var[traj_j + 1], T_common),
                                           arrows = 2, free = TRUE, values = res_cor[[k]][traj_i + traj_j - 1] * sqrt(var_1[traj_i] * var_2[traj_j]),
                                           labels = paste0("c", k, traj_var[traj_i], traj_var[traj_j + 1], "_RES"))
        }
      }
    }
    RES_L[[k]] <- RES
    COV_L[[k]] <- COV
    RES <- COV <- list()
  }
  class.list <- list()
  if (x_type == "baseline"){
    GF_loadings <- getMIX_MED.loadings(nClass = nClass, t_var = t_var, y_var = y_var, m_var = m_var,
                                       x_type = x_type, x_var = x_var, curveFun = curveFun,
                                       y_records = y_records, m_records = m_records)
    if (curveFun %in% c("linear", "LIN")){
      latents <- c("eta0Y", "eta1Y", "eta0M", "eta1M")
      for (k in 1:nClass){
        GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE, values = c(starts[[k]][[3]][[1]][1:2], starts[[k]][[2]][[1]][1:2]),
                          labels = paste0("c", k, c("Y_alpha0", "Y_alpha1", "M_alpha0", "M_alpha1")))
        GF_VAR <- list(mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]][[3]][[4]][row(starts[[k]][[3]][[4]]) >= col(starts[[k]][[3]][[4]])],
                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                       mxPath(from = latents[3:4], to = latents[3:4], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]][[2]][[3]][row(starts[[k]][[2]][[3]]) >= col(starts[[k]][[2]][[3]])],
                              labels = paste0("c", k, c("M_psi00", "M_psi01", "M_psi11"))))
        GF_LOADINGS <- list(mxPath(from = "eta0Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L1", y_records, "Y[1,1]")),
                            mxPath(from = "eta0M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L1", m_records, "M[1,1]")))
        X_BS <- list(mxPath(from = "one", to = "X", arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]], labels = paste0("c", k, "muX")),
                     mxPath(from = "X", to = "X", connect = "unique.pairs", arrows = 2, free = TRUE, values = c(starts[[k]][[1]][[2]]),
                            labels = paste0("c", k, "phi11")))
        BETA <- list(mxPath(from = "X", to = latents, arrows = 1, free = TRUE, values = c(starts[[k]][[3]][[2]], starts[[k]][[2]][[2]]),
                            labels = paste0("c", k, "beta", rep(c("Y", "M"), each = 2), rep(c(0, 1), 2))),
                     mxPath(from = latents[3], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][1, 1],
                            labels = paste0("c", k, "betaM0Y0")),
                     mxPath(from = latents[3], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][2, 1],
                            labels = paste0("c", k, "betaM0Y1")),
                     mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][2, 2],
                            labels = paste0("c", k, "betaM1Y1")))
        M_ALPHA <- mxAlgebraFromString(paste0("rbind(c", k, "M_alpha0, c", k, "M_alpha1", ")"),
                                       name = paste0("c", k, "M_alpha"))
        M_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "M_psi00, c", k, "M_psi01),",
                                              "cbind(c", k, "M_psi01, c", k, "M_psi11))"),
                                       name = paste0("c", k, "M_psi_r"))
        Y_ALPHA <- mxAlgebraFromString(paste0("rbind(c", k, "Y_alpha0, c", k, "Y_alpha1", ")"),
                                       name = paste0("c", k, "Y_alpha"))
        Y_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01),",
                                              "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                       name = paste0("c", k, "Y_psi_r"))
        BETA_XM <- mxAlgebraFromString(paste0("rbind(c", k, "betaM0, c", k, "betaM1)"),
                                       name = paste0("c", k, "beta_xm"))
        BETA_XY <- mxAlgebraFromString(paste0("rbind(c", k, "betaY0, c", k, "betaY1)"),
                                       name = paste0("c", k, "beta_xy"))
        BETA_MY <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "betaM0Y0,", 0, "),",
                                              "cbind(c", k, "betaM0Y1, c", k, "betaM1Y1))"),
                                       name = paste0("c", k, "beta_my"))
        M_MEAN0 <- mxAlgebraFromString(paste0("c", k, "M_alpha + c", k, "beta_xm %*% c", k, "muX"),
                                       name = paste0("c", k, "M_mean"))
        Y_MEAN0 <- mxAlgebraFromString(paste0("c", k, "Y_alpha + c", k, "beta_my %*% c", k,
                                              "M_mean + c", k, "beta_xy %*% c", k, "muX"),
                                       name = paste0("c", k, "Y_mean"))
        MED00 <- mxAlgebraFromString(paste0("c", k, "betaM0Y0 * c", k, "betaM0"),
                                     name = paste0("c", k, "mediator_00"))
        MED01 <- mxAlgebraFromString(paste0("c", k, "betaM0Y1 * c", k, "betaM0"),
                                     name = paste0("c", k, "mediator_01"))
        MED11 <- mxAlgebraFromString(paste0("c", k, "betaM1Y1 * c", k, "betaM1"),
                                     name = paste0("c", k, "mediator_11"))
        MED <- mxAlgebraFromString(paste0("cbind(c", k, "mediator_00, c", k, "mediator_01, c", k, "mediator_11)"),
                                   name = paste0("c", k, "mediator"))
        TOTAL <- mxAlgebraFromString(paste0("c", k, "beta_my %*% c", k, "beta_xm + c", k, "beta_xy"),
                                     name = paste0("c", k, "total"))
        class.list[[k]] <- mxModel(paste0("Class", k), type = "RAM", manifestVars = manifests, latentVars = latents,
                                   GF_loadings[[k]], GF_MEAN, GF_VAR, GF_LOADINGS, X_BS, BETA,
                                   RES_L[[k]], COV_L[[k]], M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r,
                                   BETA_XM, BETA_XY, BETA_MY, M_MEAN0, Y_MEAN0,
                                   MED00, MED01, MED11, MED, TOTAL, mxFitFunctionML(vector = T))
      }
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      latents <- c("eta1Y", "etaYr", "etaY2", "eta1M", "etaMr", "etaM2")
      for (k in 1:nClass){
        GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE,
                          values = c(starts[[k]][[3]][[1]][1:3], starts[[k]][[2]][[1]][1:3]),
                          labels = paste0("c", k, c("Y_alpha1", "Y_alphar", "Y_alpha2",
                                                    "M_alpha1", "M_alphar", "M_alpha2")))
        GF_VAR <- list(mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]][[3]][[4]][row(starts[[k]][[3]][[4]]) >= col(starts[[k]][[3]][[4]])],
                              labels = paste0("c", k, c("Y_psi11", "Y_psi1r", "Y_psi12",
                                                        "Y_psirr", "Y_psir2", "Y_psi22"))),
                       mxPath(from = latents[4:6], to = latents[4:6], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]][[2]][[3]][row(starts[[k]][[2]][[3]]) >= col(starts[[k]][[2]][[3]])],
                              labels = paste0("c", k, c("M_psi11", "M_psi1r", "M_psi12",
                                                        "M_psirr", "M_psir2", "M_psi22"))))
        GF_LOADINGS <- list(mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L1", y_records, "Y[1,1]")),
                            mxPath(from = "etaYr", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "etaY2", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L2", y_records, "Y[1,1]")),
                            mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L1", m_records, "M[1,1]")),
                            mxPath(from = "etaMr", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "etaM2", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L2", m_records, "M[1,1]")))
        GAMMA <- list(mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[3]][[1]][4],
                               labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                      mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[2]][[1]][4],
                               labels = paste0("c", k, "M_knot"), name = paste0("c", k, "M_mug")))
        X_BS <- list(mxPath(from = "one", to = "X", arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                            labels = paste0("c", k, "muX")),
                     mxPath(from = "X", to = "X", connect = "unique.pairs", arrows = 2, free = TRUE, values = c(starts[[k]][[1]][[2]]),
                            labels = paste0("c", k, "phi11")))
        BETA <- list(mxPath(from = "X", to = latents, arrows = 1, free = TRUE, values = c(starts[[k]][[3]][[2]], starts[[k]][[2]][[2]]),
                            labels = paste0("c", k, "beta", rep(c("Y", "M"), each = 3), rep(c(1, "r", 2), 2))),
                     mxPath(from = latents[4], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][1, 1],
                            labels = paste0("c", k, "betaM1Y1")),
                     mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][2, 1],
                            labels = paste0("c", k, "betaM1Yr")),
                     mxPath(from = latents[4], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][3, 1],
                            labels = paste0("c", k, "betaM1Y2")),
                     mxPath(from = latents[5], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][2, 2],
                            labels = paste0("c", k, "betaMrYr")),
                     mxPath(from = latents[5], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][3, 2],
                            labels = paste0("c", k, "betaMrY2")),
                     mxPath(from = latents[6], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][3, 3],
                            labels = paste0("c", k, "betaM2Y2")))
        M_ALPHA <- mxAlgebraFromString(paste0("rbind(c", k, "M_alpha1, c", k, "M_alphar, c", k, "M_alpha2, c", k, "M_mug)"),
                                       name = paste0("c", k, "M_alpha"))
        M_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "M_psi11, c", k, "M_psi1r, c", k, "M_psi12),",
                                                    "cbind(c", k, "M_psi1r, c", k, "M_psirr, c", k, "M_psir2),",
                                                    "cbind(c", k, "M_psi12, c", k, "M_psir2, c", k, "M_psi22))"),
                                       name = paste0("c", k, "M_psi_r"))
        Y_ALPHA <- mxAlgebraFromString(paste0("rbind(c", k, "Y_alpha1, c", k, "Y_alphar, c", k, "Y_alpha2, c", k, "Y_mug)"),
                                       name = paste0("c", k, "Y_alpha"))
        Y_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi11, c", k, "Y_psi1r, c", k, "Y_psi12),",
                                              "cbind(c", k, "Y_psi1r, c", k, "Y_psirr, c", k, "Y_psir2),",
                                              "cbind(c", k, "Y_psi12, c", k, "Y_psir2, c", k, "Y_psi22))"),
                                       name = paste0("c", k, "Y_psi_r"))
        BETA_XM <- mxAlgebraFromString(paste0("rbind(c", k, "betaM1, c", k, "betaMr, c", k, "betaM2)"),
                                       name = paste0("c", k, "beta_xm"))
        BETA_XY <- mxAlgebraFromString(paste0("rbind(c", k, "betaY1, c", k, "betaYr, c", k, "betaY2)"),
                                       name = paste0("c", k, "beta_xy"))
        BETA_MY <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "betaM1Y1, ", "0,", "0), ",
                                                    "cbind(c", k, "betaM1Yr, c", k, "betaMrYr,", "0),",
                                                    "cbind(c", k, "betaM1Y2, c", k, "betaMrY2, c", k, "betaM2Y2))"),
                                       name = paste0("c", k, "beta_my"))
        M_MEAN0 <- mxAlgebraFromString(paste0("c", k, "M_alpha[1:3, ] + c", k, "beta_xm %*% c", k, "muX"),
                                       name = paste0("c", k, "M_mean"))
        Y_MEAN0 <- mxAlgebraFromString(paste0("c", k, "Y_alpha[1:3, ] + c", k, "beta_my %*% c", k, "M_mean + c", k,
                                              "beta_xy %*% c", k, "muX"), name = paste0("c", k, "Y_mean"))
        MED11 <- mxAlgebraFromString(paste0("c", k, "betaM1Y1 * c", k, "betaM1"),
                                     name = paste0("c", k, "mediator_11"))
        MED1r <- mxAlgebraFromString(paste0("c", k, "betaM1Yr * c", k, "betaM1"),
                                     name = paste0("c", k, "mediator_1r"))
        MED12 <- mxAlgebraFromString(paste0("c", k, "betaM1Y2 * c", k, "betaM1"),
                                     name = paste0("c", k, "mediator_12"))
        MEDrr <- mxAlgebraFromString(paste0("c", k, "betaMrYr * c", k, "betaMr"),
                                     name = paste0("c", k, "mediator_rr"))
        MEDr2 <- mxAlgebraFromString(paste0("c", k, "betaMrY2 * c", k, "betaMr"),
                                     name = paste0("c", k, "mediator_r2"))
        MED22 <- mxAlgebraFromString(paste0("c", k, "betaM2Y2 * c", k, "betaM2"),
                                     name = paste0("c", k, "mediator_22"))
        MED <- mxAlgebraFromString(paste0("cbind(c", k, "mediator_11, c", k, "mediator_1r, c", k, "mediator_12, c", k,
                                          "mediator_rr, c", k, "mediator_r2, c", k, "mediator_22)"),
                                   name = paste0("c", k, "mediator"))
        TOTAL <- mxAlgebraFromString(paste0("c", k, "beta_my %*% c", k, "beta_xm + c", k, "beta_xy"),
                                     name = paste0("c", k, "total"))
        class.list[[k]] <- mxModel(paste0("Class", k), type = "RAM", manifestVars = manifests, latentVars = latents,
                                   GF_loadings[[k]], GF_MEAN, GF_VAR, GF_LOADINGS, GAMMA, X_BS, BETA,
                                   RES_L[[k]], COV_L[[k]], M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r,
                                   BETA_XM, BETA_XY, BETA_MY, M_MEAN0, Y_MEAN0,
                                   MED11, MED1r, MED12, MEDrr, MEDr2, MED22, MED, TOTAL, mxFitFunctionML(vector = T))
      }
    }
  }
  else if (x_type == "longitudinal"){
    GF_loadings <- getMIX_MED.loadings(nClass = nClass, t_var = t_var, y_var = y_var, m_var = m_var,
                                       x_type = x_type, x_var = x_var, curveFun = curveFun,
                                       y_records = y_records, m_records = m_records, x_records = x_records)
    if (curveFun %in% c("linear", "LIN")){
      latents <- c("eta0Y", "eta1Y", "eta0M", "eta1M", "eta0X", "eta1X")
      for (k in 1:nClass){
        GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE,
                          values = c(starts[[k]][[3]][[1]][1:2], starts[[k]][[2]][[1]][1:2], starts[[k]][[1]][[1]][1:2]),
                          labels = paste0("c", k, c("Y_alpha0", "Y_alpha1",
                                                    "M_alpha0", "M_alpha1",
                                                    "X_mean0", "X_mean1")))
        GF_VAR <- list(mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]][[3]][[4]][row(starts[[k]][[3]][[4]]) >= col(starts[[k]][[3]][[4]])],
                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                       mxPath(from = latents[3:4], to = latents[3:4], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]][[2]][[3]][row(starts[[k]][[2]][[3]]) >= col(starts[[k]][[2]][[3]])],
                              labels = paste0("c", k, c("M_psi00", "M_psi01", "M_psi11"))),
                       mxPath(from = latents[5:6], to = latents[5:6], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]][[1]][[2]][row(starts[[k]][[1]][[2]]) >= col(starts[[k]][[1]][[2]])],
                              labels = paste0("c", k, c("X_psi00", "X_psi01", "X_psi11"))))
        GF_LOADINGS <- list(mxPath(from = "eta0Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L1", y_records, "Y[1,1]")),
                            mxPath(from = "eta0M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L1", m_records, "M[1,1]")),
                            mxPath(from = "eta0X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L1", x_records, "X[1,1]")))
        BETA <- list(mxPath(from = latents[5], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]][[3]][[2]][1, 1],
                            labels = paste0("c", k, "betaX0Y0")),
                     mxPath(from = latents[5], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[2]][2, 1],
                            labels = paste0("c", k, "betaX0Y1")),
                     mxPath(from = latents[6], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[2]][2, 2],
                            labels = paste0("c", k, "betaX1Y1")),
                     mxPath(from = latents[5], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]][[2]][[2]][1, 1],
                            labels = paste0("c", k, "betaX0M0")),
                     mxPath(from = latents[5], to = latents[4], arrows = 1, free = TRUE, values = starts[[k]][[2]][[2]][2, 1],
                            labels = paste0("c", k, "betaX0M1")),
                     mxPath(from = latents[6], to = latents[4], arrows = 1, free = TRUE, values = starts[[k]][[2]][[2]][2, 2],
                            labels = paste0("c", k, "betaX1M1")),
                     mxPath(from = latents[3], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][1, 1],
                            labels = paste0("c", k, "betaM0Y0")),
                     mxPath(from = latents[3], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][2, 1],
                            labels = paste0("c", k, "betaM0Y1")),
                     mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][2, 2],
                            labels = paste0("c", k, "betaM1Y1")))
        X_MEAN <- mxAlgebraFromString(paste0("rbind(c", k, "X_mean0, c", k, "X_mean1)"),
                                      name = paste0("c", k, "X_mean"))
        X_PSI <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01), ",
                                                  "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                     name = paste0("c", k, "X_psi0"))
        M_ALPHA <- mxAlgebraFromString(paste0("rbind(c", k, "M_alpha0, c", k, "M_alpha1)"),
                                       name = paste0("c", k, "M_alpha"))
        M_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "M_psi00, c", k, "M_psi01), ",
                                              "cbind(c", k, "M_psi01, c", k, "M_psi11))"),
                                       name = paste0("c", k, "M_psi_r"))
        Y_ALPHA <- mxAlgebraFromString(paste0("rbind(c", k, "Y_alpha0, c", k, "Y_alpha1)"),
                                       name = paste0("c", k, "Y_alpha"))
        Y_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                              "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                       name = paste0("c", k, "Y_psi_r"))
        BETA_XM <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "betaX0M0, 0), ",
                                                    "cbind(c", k, "betaX0M1, c", k, "betaX1M1))"),
                                       name = paste0("c", k, "beta_xm"))
        BETA_XY <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "betaX0Y0, 0), ",
                                              "cbind(c", k, "betaX0Y1, c", k, "betaX1Y1))"),
                                       name = paste0("c", k, "beta_xy"))
        BETA_MY <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "betaM0Y0, 0), ",
                                              "cbind(c", k, "betaM0Y1, c", k, "betaM1Y1))"),
                                       name = paste0("c", k, "beta_my"))
        M_MEAN0 <- mxAlgebraFromString(paste0("c", k, "M_alpha + c", k, "beta_xm %*% c", k, "X_mean"),
                                       name = paste0("c", k, "M_mean"))
        Y_MEAN0 <- mxAlgebraFromString(paste0("c", k, "Y_alpha + c", k, "beta_my %*% c", k, "M_mean + c", k,
                                              "beta_xy %*% c", k, "X_mean"), name = paste0("c", k, "Y_mean"))
        MED000 <- mxAlgebraFromString(paste0("c", k, "betaX0M0 * c", k, "betaM0Y0"),
                                      name = paste0("c", k, "mediator_000"))
        MED001 <- mxAlgebraFromString(paste0("c", k, "betaX0M0 * c", k, "betaM0Y1"),
                                      name = paste0("c", k, "mediator_001"))
        MED011 <- mxAlgebraFromString(paste0("c", k, "betaX0M1 * c", k, "betaM1Y1"),
                                      name = paste0("c", k, "mediator_011"))
        MED111 <- mxAlgebraFromString(paste0("c", k, "betaX1M1 * c", k, "betaM1Y1"),
                                      name = paste0("c", k, "mediator_111"))
        MED <- mxAlgebraFromString(paste0("cbind(c", k, "mediator_000, c", k, "mediator_001, c", k,
                                          "mediator_011, c", k, "mediator_111)"),
                                   name = paste0("c", k, "mediator"))
        TOTAL <- mxAlgebraFromString(paste0("c", k, "beta_my %*% c", k, "beta_xm + c", k, "beta_xy"),
                                     name = paste0("c", k, "total"))
        class.list[[k]] <- mxModel(paste0("Class", k), type = "RAM", manifestVars = manifests, latentVars = latents,
                                   GF_loadings[[k]], GF_MEAN, GF_VAR, GF_LOADINGS, BETA,
                                   RES_L[[k]], COV_L[[k]], X_MEAN, X_PSI, M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r,
                                   BETA_XM, BETA_XY, BETA_MY, M_MEAN0, Y_MEAN0,
                                   MED000, MED001, MED011, MED111, MED, TOTAL, mxFitFunctionML(vector = T))
      }
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      latents <- c("eta1Y", "etaYr", "etaY2", "eta1M", "etaMr", "etaM2", "eta1X", "etaXr", "etaX2")
      for (k in 1:nClass){
        GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE,
                          values = c(starts[[k]][[3]][[1]][1:3], starts[[k]][[2]][[1]][1:3], starts[[k]][[1]][[1]][1:3]),
                          labels = paste0("c", k, c("Y_alpha1", "Y_alphar", "Y_alpha2",
                                                    "M_alpha1", "M_alphar", "M_alpha2",
                                                    "X_mean1", "X_meanr", "X_mean2")))
        GF_VAR <- list(mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]][[3]][[4]][row(starts[[k]][[3]][[4]]) >= col(starts[[k]][[3]][[4]])],
                              labels = paste0("c", k, c("Y_psi11", "Y_psi1r", "Y_psi12",
                                                        "Y_psirr", "Y_psir2",
                                                        "Y_psi22"))),
                       mxPath(from = latents[4:6], to = latents[4:6], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]][[2]][[3]][row(starts[[k]][[2]][[3]]) >= col(starts[[k]][[2]][[3]])],
                              labels = paste0("c", k, c("M_psi11", "M_psi1r", "M_psi12",
                                                        "M_psirr", "M_psir2",
                                                        "M_psi22"))),
                       mxPath(from = latents[7:9], to = latents[7:9], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]][[1]][[2]][row(starts[[k]][[1]][[2]]) >= col(starts[[k]][[1]][[2]])],
                              labels = paste0("c", k, c("X_psi11", "X_psi1r", "X_psi12",
                                                        "X_psirr", "X_psir2",
                                                        "X_psi22"))))
        GF_LOADINGS <- list(mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L1", y_records, "Y[1,1]")),
                            mxPath(from = "etaYr", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "etaY2", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L2", y_records, "Y[1,1]")),
                            mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L1", m_records, "M[1,1]")),
                            mxPath(from = "etaMr", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "etaM2", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L2", m_records, "M[1,1]")),
                            mxPath(from = "eta1X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L1", x_records, "X[1,1]")),
                            mxPath(from = "etaXr", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "etaX2", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("c", k, "L2", x_records, "X[1,1]")))
        GAMMA <- list(mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[3]][[1]][4],
                               labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                      mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[2]][[1]][4],
                               labels = paste0("c", k, "M_knot"), name = paste0("c", k, "M_mug")),
                      mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                               labels = paste0("c", k, "X_knot"), name = paste0("c", k, "X_mug")))
        BETA <- list(mxPath(from = latents[7], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]][[3]][[2]][1, 1],
                            labels = paste0("c", k, "betaX1Y1")),
                     ##### Slope 1 of X to intercept of Y
                     mxPath(from = latents[7], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[2]][2, 1],
                            labels = paste0("c", k, "betaX1Yr")),
                     ##### Slope 1 of X to slope 2 of Y
                     mxPath(from = latents[7], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]][[3]][[2]][3, 1],
                            labels = paste0("c", k, "betaX1Y2")),
                     ##### Intercept of X to intercept of Y
                     mxPath(from = latents[8], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[2]][2, 2],
                            labels = paste0("c", k, "betaXrYr")),
                     ##### Intercept of X to slope 2 of Y
                     mxPath(from = latents[8], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]][[3]][[2]][3, 2],
                            labels = paste0("c", k, "betaXrY2")),
                     ##### Slope 2 of X to slope 2 of Y
                     mxPath(from = latents[9], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]][[3]][[2]][3, 3],
                            labels = paste0("c", k, "betaX2Y2")),
                     ##### Slope 1 of X to slope 1 of M
                     mxPath(from = latents[7], to = latents[4], arrows = 1, free = TRUE, values = starts[[k]][[2]][[2]][1, 1],
                            labels = paste0("c", k, "betaX1M1")),
                     ##### Slope 1 of X to intercept of M
                     mxPath(from = latents[7], to = latents[5], arrows = 1, free = TRUE, values = starts[[k]][[2]][[2]][2, 1],
                            labels = paste0("c", k, "betaX1Mr")),
                     ##### Slope 1 of X to slope 2 of M
                     mxPath(from = latents[7], to = latents[6], arrows = 1, free = TRUE, values = starts[[k]][[2]][[2]][3, 1],
                            labels = paste0("c", k, "betaX1M2")),
                     ##### Intercept of X to intercept of M
                     mxPath(from = latents[8], to = latents[5], arrows = 1, free = TRUE, values = starts[[k]][[2]][[2]][2, 2],
                            labels = paste0("c", k, "betaXrMr")),
                     ##### Intercept of X to slope 2 of M
                     mxPath(from = latents[8], to = latents[6], arrows = 1, free = TRUE, values = starts[[k]][[2]][[2]][3, 2],
                            labels = paste0("c", k, "betaXrM2")),
                     ##### Slope 2 of X to slope 2 of M
                     mxPath(from = latents[9], to = latents[6], arrows = 1, free = TRUE, values = starts[[k]][[2]][[2]][3, 3],
                            labels = paste0("c", k, "betaX2M2")),

                     ##### Slope 1 of M to slope 1 of Y
                     mxPath(from = latents[4], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][1, 1],
                            labels = paste0("c", k, "betaM1Y1")),
                     ##### Slope 1 of M to intercept of Y
                     mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][2, 1],
                            labels = paste0("c", k, "betaM1Yr")),
                     ##### Slope 1 of M to slope 2 of Y
                     mxPath(from = latents[4], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][3, 1],
                            labels = paste0("c", k, "betaM1Y2")),
                     ##### Intercept of M to intercept of Y
                     mxPath(from = latents[5], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][2, 2],
                            labels = paste0("c", k, "betaMrYr")),
                     ##### Intercept of M to slope 2 of Y
                     mxPath(from = latents[5], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][3, 2],
                            labels = paste0("c", k, "betaMrY2")),
                     ##### Slope 2 of M to slope 2 of Y
                     mxPath(from = latents[6], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]][[3]][[3]][3, 3],
                            labels = paste0("c", k, "betaM2Y2")))
        X_MEAN <- mxAlgebraFromString(paste0("rbind(c", k, "X_mean1, c", k, "X_meanr, c", k, "X_mean2, c", k, "X_mug)"),
                                      name = paste0("c", k, "X_mean"))
        X_PSI <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi11, c", k, "X_psi1r, c", k, "X_psi12), ",
                                                  "cbind(c", k, "X_psi1r, c", k, "X_psirr, c", k, "X_psir2), ",
                                                  "cbind(c", k, "X_psi12, c", k, "X_psir2, c", k, "X_psi22))"),
                                     name = paste0("c", k, "X_psi0"))
        M_ALPHA <- mxAlgebraFromString(paste0("rbind(c", k, "M_alpha1, c", k, "M_alphar, c", k, "M_alpha2, c", k, "M_mug)"),
                                       name = paste0("c", k, "M_alpha"))
        M_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "M_psi11, c", k, "M_psi1r, c", k, "M_psi12), ",
                                              "cbind(c", k, "M_psi1r, c", k, "M_psirr, c", k, "M_psir2), ",
                                              "cbind(c", k, "M_psi12, c", k, "M_psir2, c", k, "M_psi22))"),
                                       name = paste0("c", k, "M_psi_r"))
        Y_ALPHA <- mxAlgebraFromString(paste0("rbind(c", k, "Y_alpha1, c", k, "Y_alphar, c", k, "Y_alpha2, c", k, "Y_mug)"),
                                       name = paste0("c", k, "Y_alpha"))
        Y_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi11, c", k, "Y_psi1r, c", k, "Y_psi12), ",
                                              "cbind(c", k, "Y_psi1r, c", k, "Y_psirr, c", k, "Y_psir2), ",
                                              "cbind(c", k, "Y_psi12, c", k, "Y_psir2, c", k, "Y_psi22))"),
                                       name = paste0("c", k, "Y_psi_r"))
        BETA_XM <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "betaX1M1, 0, 0), ",
                                                    "cbind(c", k, "betaX1Mr, c", k, "betaXrMr, 0), ",
                                                    "cbind(c", k, "betaX1M2, c", k, "betaXrM2, c", k, "betaX2M2))"),
                                       name = paste0("c", k, "beta_xm"))
        BETA_XY <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "betaX1Y1, 0, 0), ",
                                              "cbind(c", k, "betaX1Yr, c", k, "betaXrYr, 0), ",
                                              "cbind(c", k, "betaX1Y2, c", k, "betaXrY2, c", k, "betaX2Y2))"),
                                       name = paste0("c", k, "beta_xy"))
        BETA_MY <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "betaM1Y1, 0, 0), ",
                                              "cbind(c", k, "betaM1Yr, c", k, "betaMrYr, 0), ",
                                              "cbind(c", k, "betaM1Y2, c", k, "betaMrY2, c", k, "betaM2Y2))"),
                                       name = paste0("c", k, "beta_my"))
        M_MEAN0 <- mxAlgebraFromString(paste0("c", k, "M_alpha[1:3, ] + c", k, "beta_xm %*% c", k, "X_mean[1:3, ]"),
                                              name = paste0("c", k, "M_mean"))
        Y_MEAN0 <- mxAlgebraFromString(paste0("c", k, "Y_alpha[1:3, ] + c", k, "beta_my %*% c", k, "M_mean + c", k,
                                              "beta_xy %*% c", k, "X_mean[1:3, ]"),name = paste0("c", k, "Y_mean"))
        MED111 <- mxAlgebraFromString(paste0("c", k, "betaX1M1 * c", k, "betaM1Y1"), name = paste0("c", k, "mediator_111"))
        MED11r <- mxAlgebraFromString(paste0("c", k, "betaX1M1 * c", k, "betaM1Yr"), name = paste0("c", k, "mediator_11r"))
        MED112 <- mxAlgebraFromString(paste0("c", k, "betaX1M1 * c", k, "betaM1Y2"), name = paste0("c", k, "mediator_112"))
        MED1rr <- mxAlgebraFromString(paste0("c", k, "betaX1Mr * c", k, "betaMrYr"), name = paste0("c", k, "mediator_1rr"))
        MED1r2 <- mxAlgebraFromString(paste0("c", k, "betaX1Mr * c", k, "betaMrY2"), name = paste0("c", k, "mediator_1r2"))
        MED122 <- mxAlgebraFromString(paste0("c", k, "betaX1M2 * c", k, "betaM2Y2"), name = paste0("c", k, "mediator_122"))
        MEDrrr <- mxAlgebraFromString(paste0("c", k, "betaXrMr * c", k, "betaMrYr"), name = paste0("c", k, "mediator_rrr"))
        MEDrr2 <- mxAlgebraFromString(paste0("c", k, "betaXrMr * c", k, "betaMrY2"), name = paste0("c", k, "mediator_rr2"))
        MEDr22 <- mxAlgebraFromString(paste0("c", k, "betaXrM2 * c", k, "betaM2Y2"), name = paste0("c", k, "mediator_r22"))
        MED222 <- mxAlgebraFromString(paste0("c", k, "betaX2M2 * c", k, "betaM2Y2"), name = paste0("c", k, "mediator_222"))
        MED <- mxAlgebraFromString(paste0("cbind(c", k, "mediator_111, c", k, "mediator_11r, c", k, "mediator_112, c", k,
                                                "mediator_1rr, c", k, "mediator_1r2, c", k, "mediator_122, c", k,
                                                "mediator_rrr, c", k, "mediator_rr2, c", k, "mediator_r22, c", k,
                                                "mediator_222)"), name = paste0("c", k, "mediator"))
        TOTAL <- mxAlgebraFromString(paste0("c", k, "beta_my %*% c", k, "beta_xm + c", k, "beta_xy"),
                                     name = paste0("c", k, "total"))
        class.list[[k]] <- mxModel(paste0("Class", k), type = "RAM", manifestVars = manifests, latentVars = latents,
                                   GF_loadings[[k]], GF_MEAN, GF_VAR, GF_LOADINGS, GAMMA, BETA,
                                   RES_L[[k]], COV_L[[k]], X_MEAN, X_PSI, M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r,
                                   BETA_XM, BETA_XY, BETA_MY, M_MEAN0, Y_MEAN0,
                                   MED111, MED11r, MED112, MED1rr, MED1r2, MED122, MEDrrr, MEDrr2,
                                   MEDr22, MED222, MED, TOTAL, mxFitFunctionML(vector = T))
      }
    }
  }
  return(class.list)
}






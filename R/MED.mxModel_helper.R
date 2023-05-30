#' @title Construct An Object of mxModel for Longitudinal Mediation Models To Be Evaluated
#'
#' @description This function builds up an object of mxModel for a longitudinal mediation model with user-specified functional form.
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
#' @param starts A list of initial values for the parameters, either takes the value passed from \code{getMediation()} or
#' derived by the helper function \code{getMED.initial()}.
#' @param res_cor A numeric value or vector for user-specified residual correlation between any two longitudinal processes to calculate
#' the corresponding initial value. It takes the value passed from \code{getMediation()}.
#'
#' @return A pre-optimized mxModel for a longitudinal mediation model.
#'
#' @importFrom OpenMx mxPath mxModel mxAlgebra mxMatrix mxData
#'
getMED.mxModel <- function(dat, t_var, y_var, m_var, x_type, x_var, curveFun, records, starts, res_cor){
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
    var0 <- c(starts[[3]][[5]], starts[[2]][[4]])
    var_1 <- c(starts[[3]][[5]])
    var_2 <- c(starts[[2]][[4]])
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
    var0 <- c(starts[[3]][[5]], starts[[2]][[4]], starts[[1]][[3]])
    var_1 <- c(starts[[3]][[5]], starts[[3]][[5]], starts[[2]][[4]])
    var_2 <- c(starts[[2]][[4]], starts[[1]][[3]], starts[[1]][[3]])
  }
  for (traj in 1:length(traj_var)){
    RES[[length(RES) + 1]] <- mxPath(from = traj_list[[traj]], to = traj_list[[traj]], arrows = 2, free = TRUE, values = var0[traj],
                                     labels = paste0(traj_var[traj], "_residuals"))
  }
  for (traj_i in 1:(length(traj_var) - 1)){
    for (traj_j in traj_i:(length(traj_var) - 1)){
      if (setequal(readr::parse_number(traj_list[[traj_i]]), readr::parse_number(traj_list[[traj_j + 1]]))){
        COV[[length(COV) + 1]] <- mxPath(from = traj_list[[traj_i]], to = traj_list[[traj_j + 1]],
                                         arrows = 2, free = TRUE, values = res_cor[traj_i + traj_j - 1] * sqrt(var_1[traj_i] * var_2[traj_j]),
                                         labels = paste0(traj_var[traj_i], traj_var[traj_j + 1], "_RES"))
      }
      else{
        T_common <- Reduce(intersect, list(readr::parse_number(traj_list[[traj_i]]), readr::parse_number(traj_list[[traj_j + 1]])))
        COV[[length(COV) + 1]] <- mxPath(from = paste0(traj_var[traj_i], T_common),
                                         to = paste0(traj_var[traj_j + 1], T_common),
                                         arrows = 2, free = TRUE, values = res_cor[traj_i + traj_j - 1] * sqrt(var_1[traj_i] * var_2[traj_j]),
                                         labels = paste0(traj_var[traj_i], traj_var[traj_j + 1], "_RES"))
      }
    }
  }
  if (x_type == "baseline"){
    GF_loadings <- getMED.loadings(t_var = t_var, y_var = y_var, m_var = m_var, x_type = x_type, x_var = x_var,
                                   curveFun = curveFun, y_records = y_records, m_records = m_records)
    if (curveFun %in% c("linear", "LIN")){
      latents <- c("eta0Y", "eta1Y", "eta0M", "eta1M")
      GF_MEAN <- GF_VAR <- GF_LOADINGS <- X_BS <- BETA <- list()
      GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE, values = c(starts[[3]][[1]][1:2], starts[[2]][[1]][1:2]),
                        labels = c("Y_alpha0", "Y_alpha1", "M_alpha0", "M_alpha1"))
      GF_VAR <- list(mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                            free = TRUE, values = starts[[3]][[4]][row(starts[[3]][[4]]) >= col(starts[[3]][[4]])],
                            labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                     mxPath(from = latents[3:4], to = latents[3:4], arrows = 2, connect = "unique.pairs",
                            free = TRUE, values = starts[[2]][[3]][row(starts[[2]][[3]]) >= col(starts[[2]][[3]])],
                            labels = c("M_psi00", "M_psi01", "M_psi11")))
      GF_LOADINGS <- list(mxPath(from = "eta0Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", y_records, "Y[1,1]")),
                          mxPath(from = "eta0M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", m_records, "M[1,1]")))
      X_BS <- list(mxPath(from = "one", to = x_var, arrows = 1, free = TRUE, values = starts[[1]][[1]], labels = c("muX")),
                   mxPath(from = x_var, to = x_var, connect = "unique.pairs", arrows = 2, free = TRUE, values = c(starts[[1]][[2]]),
                          labels = "phi11"))
      BETA <- list(mxPath(from = x_var, to = latents, arrows = 1, free = TRUE, values = c(starts[[3]][[2]], starts[[2]][[2]]),
                          labels = paste0("beta", rep(c("Y", "M"), each = 2), rep(c(0, 1), 2))),
                   mxPath(from = latents[3], to = latents[1], arrows = 1, free = TRUE, values = starts[[3]][[3]][1, 1],
                          labels = "betaM0Y0"),
                   mxPath(from = latents[3], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[3]][2, 1],
                          labels = "betaM0Y1"),
                   mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[3]][2, 2],
                          labels = "betaM1Y1"))
      M_ALPHA <- mxAlgebra(rbind(M_alpha0, M_alpha1), name = "M_alpha")
      M_PSI_r <- mxAlgebra(rbind(cbind(M_psi00, M_psi01),
                                 cbind(M_psi01, M_psi11)), name = "M_psi_r")
      Y_ALPHA <- mxAlgebra(rbind(Y_alpha0, Y_alpha1), name = "Y_alpha")
      Y_PSI_r <- mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                 cbind(Y_psi01, Y_psi11)), name = "Y_psi_r")
      BETA_XM <- mxAlgebra(rbind(betaM0, betaM1), name = "beta_xm")
      BETA_XY <- mxAlgebra(rbind(betaY0, betaY1), name = "beta_xy")
      BETA_MY <- mxAlgebra(rbind(cbind(betaM0Y0, 0),
                                 cbind(betaM0Y1, betaM1Y1)), name = "beta_my")
      M_MEAN0 <- mxAlgebra(M_alpha + beta_xm %*% muX, name = "M_mean")
      Y_MEAN0 <- mxAlgebra(Y_alpha + beta_my %*% M_mean + beta_xy %*% muX, name = "Y_mean")
      MED00 <- mxAlgebra(betaM0Y0 * betaM0, name = "mediator_00")
      MED01 <- mxAlgebra(betaM0Y1 * betaM0, name = "mediator_01")
      MED11 <- mxAlgebra(betaM1Y1 * betaM1, name = "mediator_11")
      MED <- mxAlgebra(cbind(mediator_00, mediator_01, mediator_11), name = "mediator")
      TOTAL <- mxAlgebra(beta_my %*% beta_xm + beta_xy, name = "total")
      model_mx <- mxModel("Longitudinal Mediation Model with Baseline Predictor, Linear", type = "RAM",
                          manifestVars = manifests, latentVars = latents, mxData(observed = dat, type = "raw"),
                          GF_loadings, GF_MEAN, GF_VAR, GF_LOADINGS, X_BS, BETA, RES, COV,
                          M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r, BETA_XM, BETA_XY, BETA_MY, M_MEAN0, Y_MEAN0,
                          MED00, MED01, MED11, MED, TOTAL)
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      latents <- c("eta1Y", "etarY", "eta2Y", "eta1M", "etarM", "eta2M")
      GF_MEAN <- GF_VAR <- GF_LOADINGS <- GAMMA <- X_BS <- BETA <- list()
      GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE,
                        values = c(starts[[3]][[1]][1:3], starts[[2]][[1]][1:3]),
                        labels = c("Y_alpha1", "Y_alphar", "Y_alpha2", "M_alpha1", "M_alphar", "M_alpha2"))
      GF_VAR <- list(mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                            free = TRUE, values = starts[[3]][[4]][row(starts[[3]][[4]]) >= col(starts[[3]][[4]])],
                            labels = c("Y_psi11", "Y_psi1r", "Y_psi12", "Y_psirr", "Y_psir2", "Y_psi22")),
                     mxPath(from = latents[4:6], to = latents[4:6], arrows = 2, connect = "unique.pairs",
                            free = TRUE, values = starts[[2]][[3]][row(starts[[2]][[3]]) >= col(starts[[2]][[3]])],
                            labels = c("M_psi11", "M_psi1r", "M_psi12", "M_psirr", "M_psir2", "M_psi22")))
      GF_LOADINGS <- list(mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", y_records, "Y[1,1]")),
                          mxPath(from = "etarY", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta2Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L2", y_records, "Y[1,1]")),
                          mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", m_records, "M[1,1]")),
                          mxPath(from = "etarM", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta2M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L2", m_records, "M[1,1]")))
      GAMMA <- list(mxMatrix("Full", 1, 1, free = TRUE, values = starts[[3]][[1]][4], labels = "Y_knot", name = "Y_mug"),
                    mxMatrix("Full", 1, 1, free = TRUE, values = starts[[2]][[1]][4], labels = "M_knot", name = "M_mug"))
      X_BS <- list(mxPath(from = "one", to = x_var, arrows = 1, free = TRUE, values = starts[[1]][[1]], labels = c("muX")),
                   mxPath(from = x_var, to = x_var, connect = "unique.pairs", arrows = 2, free = TRUE, values = c(starts[[1]][[2]]),
                          labels = "phi11"))
      BETA <- list(mxPath(from = x_var, to = latents, arrows = 1, free = TRUE, values = c(starts[[3]][[2]], starts[[2]][[2]]),
                          labels = paste0("beta", rep(c("Y", "M"), each = 3), rep(c(1, "r", 2), 2))),
                   mxPath(from = latents[4], to = latents[1], arrows = 1, free = TRUE, values = starts[[3]][[3]][1, 1],
                          labels = "betaM1Y1"),
                   mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[3]][2, 1],
                          labels = "betaM1Yr"),
                   mxPath(from = latents[4], to = latents[3], arrows = 1, free = TRUE, values = starts[[3]][[3]][3, 1],
                          labels = "betaM1Y2"),
                   mxPath(from = latents[5], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[3]][2, 2],
                          labels = "betaMrYr"),
                   mxPath(from = latents[5], to = latents[3], arrows = 1, free = TRUE, values = starts[[3]][[3]][3, 2],
                          labels = "betaMrY2"),
                   mxPath(from = latents[6], to = latents[3], arrows = 1, free = TRUE, values = starts[[3]][[3]][3, 3],
                          labels = "betaM2Y2"))
      M_ALPHA <- mxAlgebra(rbind(M_alpha1, M_alphar, M_alpha2, M_mug), name = "M_alpha")
      M_PSI_r <- mxAlgebra(rbind(cbind("M_psi11", "M_psi1r", "M_psi12"),
                                 cbind("M_psi1r", "M_psirr", "M_psir2"),
                                 cbind("M_psi12", "M_psir2", "M_psi22")), name = "M_psi_r")
      Y_ALPHA <- mxAlgebra(rbind(Y_alpha1, Y_alphar, Y_alpha2, Y_mug), name = "Y_alpha")
      Y_PSI_r <- mxAlgebra(rbind(cbind("Y_psi11", "Y_psi1r", "Y_psi12"),
                                 cbind("Y_psi1r", "Y_psirr", "Y_psir2"),
                                 cbind("Y_psi12", "Y_psir2", "Y_psi22")), name = "Y_psi_r")
      BETA_XM <- mxAlgebra(rbind(betaM1, betaMr, betaM2), name = "beta_xm")
      BETA_XY <- mxAlgebra(rbind(betaY1, betaYr, betaY2), name = "beta_xy")
      BETA_MY <- mxAlgebra(rbind(cbind(betaM1Y1, 0, 0),
                                 cbind(betaM1Yr, betaMrYr, 0),
                                 cbind(betaM1Y2, betaMrY2, betaM2Y2)), name = "beta_my")
      M_MEAN0 <- mxAlgebra(M_alpha[1:3, ] + beta_xm %*% muX, name = "M_mean")
      Y_MEAN0 <- mxAlgebra(Y_alpha[1:3, ] + beta_my %*% M_mean + beta_xy %*% muX, name = "Y_mean")
      MED11 <- mxAlgebra(betaM1Y1 * betaM1, name = "mediator_11")
      MED1r <- mxAlgebra(betaM1Yr * betaM1, name = "mediator_1r")
      MED12 <- mxAlgebra(betaM1Y2 * betaM1, name = "mediator_12")
      MEDrr <- mxAlgebra(betaMrYr * betaMr, name = "mediator_rr")
      MEDr2 <- mxAlgebra(betaMrY2 * betaMr, name = "mediator_r2")
      MED22 <- mxAlgebra(betaM2Y2 * betaM2, name = "mediator_22")
      MED <- mxAlgebra(cbind(mediator_11, mediator_1r, mediator_12,
                             mediator_rr, mediator_r2, mediator_22), name = "mediator")
      TOTAL <- mxAlgebra(beta_my %*% beta_xm + beta_xy, name = "total")
      model_mx <- mxModel("Longitudinal Mediation Model with Baseline Predictor, Bilinear Spline", type = "RAM",
                          manifestVars = manifests, latentVars = latents, mxData(observed = dat, type = "raw"),
                          GF_loadings, GF_MEAN, GF_VAR, GF_LOADINGS, GAMMA, X_BS, BETA, RES, COV,
                          M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r, BETA_XM, BETA_XY, BETA_MY, M_MEAN0, Y_MEAN0,
                          MED11, MED1r, MED12, MEDrr, MEDr2, MED22, MED, TOTAL)
    }
  }
  else if (x_type == "longitudinal"){
    GF_loadings <- getMED.loadings(t_var = t_var, y_var = y_var, m_var = m_var, x_type = x_type, x_var = x_var,
                                   curveFun = curveFun, y_records = y_records, m_records = m_records,
                                   x_records = x_records)
    if (curveFun %in% c("linear", "LIN")){
      latents <- c("eta0Y", "eta1Y", "eta0M", "eta1M", "eta0X", "eta1X")
      GF_MEAN <- GF_VAR <- GF_LOADINGS <- BETA <- list()
      GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE,
                        values = c(starts[[3]][[1]][1:2], starts[[2]][[1]][1:2], starts[[1]][[1]][1:2]),
                        labels = c("Y_alpha0", "Y_alpha1", "M_alpha0", "M_alpha1", "X_mean0", "X_mean1"))
      GF_VAR <- list(mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                            free = TRUE, values = starts[[3]][[4]][row(starts[[3]][[4]]) >= col(starts[[3]][[4]])],
                            labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                     mxPath(from = latents[3:4], to = latents[3:4], arrows = 2, connect = "unique.pairs",
                            free = TRUE, values = starts[[2]][[3]][row(starts[[2]][[3]]) >= col(starts[[2]][[3]])],
                            labels = c("M_psi00", "M_psi01", "M_psi11")),
                     mxPath(from = latents[5:6], to = latents[5:6], arrows = 2, connect = "unique.pairs",
                            free = TRUE, values = starts[[1]][[2]][row(starts[[1]][[2]]) >= col(starts[[1]][[2]])],
                            labels = c("X_psi00", "X_psi01", "X_psi11")))
      GF_LOADINGS <- list(mxPath(from = "eta0Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", y_records, "Y[1,1]")),
                          mxPath(from = "eta0M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", m_records, "M[1,1]")),
                          mxPath(from = "eta0X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta1X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", x_records, "X[1,1]")))
      BETA <- list(mxPath(from = latents[5], to = latents[1], arrows = 1, free = TRUE, values = starts[[3]][[2]][1, 1],
                          labels = "betaX0Y0"),
                   mxPath(from = latents[5], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[2]][2, 1],
                          labels = "betaX0Y1"),
                   mxPath(from = latents[6], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[2]][2, 2],
                          labels = "betaX1Y1"),
                   mxPath(from = latents[5], to = latents[3], arrows = 1, free = TRUE, values = starts[[2]][[2]][1, 1],
                          labels = "betaX0M0"),
                   mxPath(from = latents[5], to = latents[4], arrows = 1, free = TRUE, values = starts[[2]][[2]][2, 1],
                          labels = "betaX0M1"),
                   mxPath(from = latents[6], to = latents[4], arrows = 1, free = TRUE, values = starts[[2]][[2]][2, 2],
                          labels = "betaX1M1"),
                   mxPath(from = latents[3], to = latents[1], arrows = 1, free = TRUE, values = starts[[3]][[3]][1, 1],
                          labels = "betaM0Y0"),
                   mxPath(from = latents[3], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[3]][2, 1],
                          labels = "betaM0Y1"),
                   mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[3]][2, 2],
                          labels = "betaM1Y1"))
      X_MEAN <- mxAlgebra(rbind(X_mean0, X_mean1), name = "X_mean")
      X_PSI <- mxAlgebra(rbind(cbind("X_psi00", "X_psi01"),
                               cbind("X_psi01", "X_psi11")), name = "X_psi0")
      M_ALPHA <- mxAlgebra(rbind(M_alpha0, M_alpha1), name = "M_alpha")
      M_PSI_r <- mxAlgebra(rbind(cbind("M_psi00", "M_psi01"),
                                 cbind("M_psi01", "M_psi11")), name = "M_psi_r")
      Y_ALPHA <- mxAlgebra(rbind(Y_alpha0, Y_alpha1), name = "Y_alpha")
      Y_PSI_r <- mxAlgebra(rbind(cbind("Y_psi00", "Y_psi01"),
                                 cbind("Y_psi01", "Y_psi11")), name = "Y_psi_r")
      BETA_XM <- mxAlgebra(rbind(cbind(betaX0M0, 0),
                                 cbind(betaX0M1, betaX1M1)), name = "beta_xm")
      BETA_XY <- mxAlgebra(rbind(cbind(betaX0Y0, 0),
                                 cbind(betaX0Y1, betaX1Y1)), name = "beta_xy")
      BETA_MY <- mxAlgebra(rbind(cbind(betaM0Y0, 0),
                                 cbind(betaM0Y1, betaM1Y1)), name = "beta_my")
      M_MEAN0 <- mxAlgebra(M_alpha + beta_xm %*% X_mean, name = "M_mean")
      Y_MEAN0 <- mxAlgebra(Y_alpha + beta_my %*% M_mean + beta_xy %*% X_mean, name = "Y_mean")
      MED000 <- mxAlgebra(betaX0M0 * betaM0Y0, name = "mediator_000")
      MED001 <- mxAlgebra(betaX0M0 * betaM0Y1, name = "mediator_001")
      MED011 <- mxAlgebra(betaX0M1 * betaM1Y1, name = "mediator_011")
      MED111 <- mxAlgebra(betaX1M1 * betaM1Y1, name = "mediator_111")
      MED <- mxAlgebra(cbind(mediator_000, mediator_001, mediator_011, mediator_111), name = "mediator")
      TOTAL <- mxAlgebra(beta_my %*% beta_xm + beta_xy, name = "total")
      model_mx <- mxModel("Longitudinal Mediation Model with Longitudinal Predictor, Linear", type = "RAM",
                          manifestVars = manifests, latentVars = latents, mxData(observed = dat, type = "raw"),
                          GF_loadings, GF_MEAN, GF_VAR, GF_LOADINGS, BETA, RES, COV,
                          X_MEAN, X_PSI, M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r, BETA_XM, BETA_XY, BETA_MY,
                          M_MEAN0, Y_MEAN0, MED000, MED001, MED011, MED111, MED, TOTAL)
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      latents <- c("eta1Y", "etarY", "eta2Y", "eta1M", "etarM", "eta2M", "eta1X", "etarX", "eta2X")
      GF_MEAN <- GF_VAR <- GF_LOADINGS <- GAMMA <- BETA <- list()
      GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE,
                        values = c(starts[[3]][[1]][1:3], starts[[2]][[1]][1:3], starts[[1]][[1]][1:3]),
                        labels = c("Y_alpha1", "Y_alphar", "Y_alpha2", "M_alpha1", "M_alphar", "M_alpha2",
                                   "X_mean1", "X_meanr", "X_mean2"))
      GF_VAR <- list(mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                            free = TRUE, values = starts[[3]][[4]][row(starts[[3]][[4]]) >= col(starts[[3]][[4]])],
                            labels = c("Y_psi11", "Y_psi1r", "Y_psi12",
                                       "Y_psirr", "Y_psir2",
                                       "Y_psi22")),
                     mxPath(from = latents[4:6], to = latents[4:6], arrows = 2, connect = "unique.pairs",
                            free = TRUE, values = starts[[2]][[3]][row(starts[[2]][[3]]) >= col(starts[[2]][[3]])],
                            labels = c("M_psi11", "M_psi1r", "M_psi12",
                                       "M_psirr", "M_psir2",
                                       "M_psi22")),
                     mxPath(from = latents[7:9], to = latents[7:9], arrows = 2, connect = "unique.pairs",
                            free = TRUE, values = starts[[1]][[2]][row(starts[[1]][[2]]) >= col(starts[[1]][[2]])],
                            labels = c("X_psi11", "X_psi1r", "X_psi12",
                                       "X_psirr", "X_psir2",
                                       "X_psi22")))
      GF_LOADINGS <- list(mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", y_records, "Y[1,1]")),
                          mxPath(from = "etarY", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta2Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L2", y_records, "Y[1,1]")),
                          mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", m_records, "M[1,1]")),
                          mxPath(from = "etarM", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta2M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L2", m_records, "M[1,1]")),
                          mxPath(from = "eta1X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", x_records, "X[1,1]")),
                          mxPath(from = "etarX", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta2X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L2", x_records, "X[1,1]")))
      GAMMA <- list(mxMatrix("Full", 1, 1, free = TRUE, values = starts[[3]][[1]][4],
                             labels = "Y_knot", name = "Y_mug"),
                    mxMatrix("Full", 1, 1, free = TRUE, values = starts[[2]][[1]][4],
                             labels = "M_knot", name = "M_mug"),
                    mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4],
                             labels = "X_knot", name = "X_mug"))
      BETA <- list(mxPath(from = latents[7], to = latents[1], arrows = 1, free = TRUE, values = starts[[3]][[2]][1, 1],
                          labels = "betaX1Y1"),
                   ##### Slope 1 of X to intercept of Y
                   mxPath(from = latents[7], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[2]][2, 1],
                          labels = "betaX1Yr"),
                   ##### Slope 1 of X to slope 2 of Y
                   mxPath(from = latents[7], to = latents[3], arrows = 1, free = TRUE, values = starts[[3]][[2]][3, 1],
                          labels = "betaX1Y2"),
                   ##### Intercept of X to intercept of Y
                   mxPath(from = latents[8], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[2]][2, 2],
                          labels = "betaXrYr"),
                   ##### Intercept of X to slope 2 of Y
                   mxPath(from = latents[8], to = latents[3], arrows = 1, free = TRUE, values = starts[[3]][[2]][3, 2],
                          labels = "betaXrY2"),
                   ##### Slope 2 of X to slope 2 of Y
                   mxPath(from = latents[9], to = latents[3], arrows = 1, free = TRUE, values = starts[[3]][[2]][3, 3],
                          labels = "betaX2Y2"),

                   ##### Slope 1 of X to slope 1 of M
                   mxPath(from = latents[7], to = latents[4], arrows = 1, free = TRUE, values = starts[[2]][[2]][1, 1],
                          labels = "betaX1M1"),
                   ##### Slope 1 of X to intercept of M
                   mxPath(from = latents[7], to = latents[5], arrows = 1, free = TRUE, values = starts[[2]][[2]][2, 1],
                          labels = "betaX1Mr"),
                   ##### Slope 1 of X to slope 2 of M
                   mxPath(from = latents[7], to = latents[6], arrows = 1, free = TRUE, values = starts[[2]][[2]][3, 1],
                          labels = "betaX1M2"),
                   ##### Intercept of X to intercept of M
                   mxPath(from = latents[8], to = latents[5], arrows = 1, free = TRUE, values = starts[[2]][[2]][2, 2],
                          labels = "betaXrMr"),
                   ##### Intercept of X to slope 2 of M
                   mxPath(from = latents[8], to = latents[6], arrows = 1, free = TRUE, values = starts[[2]][[2]][3, 2],
                          labels = "betaXrM2"),
                   ##### Slope 2 of X to slope 2 of M
                   mxPath(from = latents[9], to = latents[6], arrows = 1, free = TRUE, values = starts[[2]][[2]][3, 3],
                          labels = "betaX2M2"),

                   ##### Slope 1 of M to slope 1 of Y
                   mxPath(from = latents[4], to = latents[1], arrows = 1, free = TRUE, values = starts[[3]][[3]][1, 1],
                          labels = "betaM1Y1"),
                   ##### Slope 1 of M to intercept of Y
                   mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[3]][2, 1],
                          labels = "betaM1Yr"),
                   ##### Slope 1 of M to slope 2 of Y
                   mxPath(from = latents[4], to = latents[3], arrows = 1, free = TRUE, values = starts[[3]][[3]][3, 1],
                          labels = "betaM1Y2"),
                   ##### Intercept of M to intercept of Y
                   mxPath(from = latents[5], to = latents[2], arrows = 1, free = TRUE, values = starts[[3]][[3]][2, 2],
                          labels = "betaMrYr"),
                   ##### Intercept of M to slope 2 of Y
                   mxPath(from = latents[5], to = latents[3], arrows = 1, free = TRUE, values = starts[[3]][[3]][3, 2],
                          labels = "betaMrY2"),
                   ##### Slope 2 of M to slope 2 of Y
                   mxPath(from = latents[6], to = latents[3], arrows = 1, free = TRUE, values = starts[[3]][[3]][3, 3],
                          labels = "betaM2Y2"))
      X_MEAN <- mxAlgebra(rbind(X_mean1, X_meanr, X_mean2, X_mug), name = "X_mean")
      X_PSI <- mxAlgebra(rbind(cbind("X_psi11", "X_psi1r", "X_psi12"),
                               cbind("X_psi1r", "X_psirr", "X_psir2"),
                               cbind("X_psi12", "X_psir2", "X_psi22")), name = "X_psi0")
      M_ALPHA <- mxAlgebra(rbind(M_alpha1, M_alphar, M_alpha2, M_mug), name = "M_alpha")
      M_PSI_r <- mxAlgebra(rbind(cbind("M_psi11", "M_psi1r", "M_psi12"),
                                 cbind("M_psi1r", "M_psirr", "M_psir2"),
                                 cbind("M_psi12", "M_psir2", "M_psi22")), name = "M_psi_r")
      Y_ALPHA <- mxAlgebra(rbind(Y_alpha1, Y_alphar, Y_alpha2, Y_mug), name = "Y_alpha")
      Y_PSI_r <- mxAlgebra(rbind(cbind("Y_psi11", "Y_psi1r", "Y_psi12"),
                                 cbind("Y_psi1r", "Y_psirr", "Y_psir2"),
                                 cbind("Y_psi12", "Y_psir2", "Y_psi22")), name = "Y_psi_r")
      BETA_XM <- mxAlgebra(rbind(cbind(betaX1M1, 0, 0),
                                 cbind(betaX1Mr, betaXrMr, 0),
                                 cbind(betaX1M2, betaXrM2, betaX2M2)), name = "beta_xm")
      BETA_XY <- mxAlgebra(rbind(cbind(betaX1Y1, 0, 0),
                                 cbind(betaX1Yr, betaXrYr, 0),
                                 cbind(betaX1Y2, betaXrY2, betaX2Y2)), name = "beta_xy")
      BETA_MY <- mxAlgebra(rbind(cbind(betaM1Y1, 0, 0),
                                 cbind(betaM1Yr, betaMrYr, 0),
                                 cbind(betaM1Y2, betaMrY2, betaM2Y2)), name = "beta_my")
      M_MEAN0 <- mxAlgebra(M_alpha[1:3, ] + beta_xm %*% X_mean[1:3, ], name = "M_mean")
      Y_MEAN0 <- mxAlgebra(Y_alpha[1:3, ] + beta_my %*% M_mean + beta_xy %*% X_mean[1:3, ], name = "Y_mean")
      MED111 <- mxAlgebra(betaX1M1 * betaM1Y1, name = "mediator_111")
      MED11r <- mxAlgebra(betaX1M1 * betaM1Yr, name = "mediator_11r")
      MED112 <- mxAlgebra(betaX1M1 * betaM1Y2, name = "mediator_112")
      MED1rr <- mxAlgebra(betaX1Mr * betaMrYr, name = "mediator_1rr")
      MED1r2 <- mxAlgebra(betaX1Mr * betaMrY2, name = "mediator_1r2")
      MED122 <- mxAlgebra(betaX1M2 * betaM2Y2, name = "mediator_122")
      MEDrrr <- mxAlgebra(betaXrMr * betaMrYr, name = "mediator_rrr")
      MEDrr2 <- mxAlgebra(betaXrMr * betaMrY2, name = "mediator_rr2")
      MEDr22 <- mxAlgebra(betaXrM2 * betaM2Y2, name = "mediator_r22")
      MED222 <- mxAlgebra(betaX2M2 * betaM2Y2, name = "mediator_222")
      MED <- mxAlgebra(cbind(mediator_111, mediator_11r, mediator_112,
                             mediator_1rr, mediator_1r2, mediator_122,
                             mediator_rrr, mediator_rr2, mediator_r22,
                             mediator_222), name = "mediator")
      TOTAL <- mxAlgebra(beta_my %*% beta_xm + beta_xy, name = "total")
      model_mx <- mxModel("Longitudinal Mediation Model with Longitudinal Predictor, Bilinear Spline", type = "RAM",
                          manifestVars = manifests, latentVars = latents,
                          mxData(observed = dat, type = "raw"),
                          GF_loadings, GF_MEAN, GF_VAR, GF_LOADINGS, GAMMA, BETA, RES, COV,
                          X_MEAN, X_PSI, M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r, BETA_XM, BETA_XY, BETA_MY,
                          M_MEAN0, Y_MEAN0, MED111, MED11r, MED112, MED1rr, MED1r2, MED122, MEDrrr, MEDrr2,
                          MEDr22, MED222, MED, TOTAL)
    }
  }
  return(model_mx)
}






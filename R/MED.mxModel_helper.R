#' @title One builder for the MED sub-model, single-group or class-specific
#' @description Builds one mxModel per element of \code{pfx}: with \code{pfx = ""} the single-group model,
#' with \code{pfx = c("c1", "c2", ...)} the class-specific models of a multiple-group or mixture model.
#' Parameter labels, matrix and algebra names carry the prefix. \code{starts}, \code{model_names},
#' \code{extra} and the loading/derived-parameter objects are per-class lists; \code{extra[[k]]} holds
#' the components the container adds to class k (its \code{mxData} for a multiple-group model, or
#' \code{mxFitFunctionML(vector = TRUE)} for a mixture model).
#' @keywords internal
#' @noRd
#'
#' @importFrom OpenMx mxPath mxModel mxAlgebraFromString mxMatrix mxData mxFitFunctionML
#'
.build_MED <- function(dat, t_var, y_var, curveFun, records, m_var, x_var, x_type, starts, res_cor, pfx, model_names, extra, GF_loadings){
  nClass <- length(pfx)
  # If res_cor not provided, default to 0.3 for mxPath starting values (freely estimated parameter)
  if (is.null(res_cor)){
    n_traj <- if (x_type == "baseline") 2 else 3
    n_pairs <- n_traj * (n_traj - 1) / 2
    res_cor <- lapply(1:nClass, function(k) rep(0.3, n_pairs))
  }
  RES_L <- COV_L <- list()
  for (k in 1:nClass){
    RES <- COV <- list()
    if (x_type == "baseline"){
      y_records <- records[[1]]
      m_records <- records[[2]]
      traj_var <- c(y_var, m_var)
      traj_list <- list()
      for (traj in seq_along(traj_var)){
        traj_list[[length(traj_list) + 1]] <- paste0(traj_var[traj], records[[traj]])
      }
      manifests <- c(unlist(traj_list), x_var)
      var0 <- c(starts[[k]]$Y$residuals, starts[[k]]$M$residuals)
    }
    else if (x_type == "longitudinal"){
      y_records <- records[[1]]
      m_records <- records[[2]]
      x_records <- records[[3]]
      traj_var <- c(y_var, m_var, x_var)
      traj_list <- list()
      for (traj in seq_along(traj_var)){
        traj_list[[length(traj_list) + 1]] <- paste0(traj_var[traj], records[[traj]])
      }
      manifests <- unlist(traj_list)
      var0 <- c(starts[[k]]$Y$residuals, starts[[k]]$M$residuals, starts[[k]]$X$residuals)
    }
    for (traj in seq_along(traj_var)){
      RES[[length(RES) + 1]] <- mxPath(from = traj_list[[traj]], to = traj_list[[traj]], arrows = 2, free = TRUE, values = var0[traj],
                                       labels = paste0(pfx[k], traj_var[traj], "_residuals"))
    }
    # Pair contemporaneous residuals by the supplied record IDs, independently
    # of column ordering or digits in variable prefixes. Every process pair
    # has one freely estimated residual covariance and therefore needs a shared
    # record; silently omitting it would change the output/model contract.
    for (traj_i in seq_len(length(traj_var) - 1L)){
      for (traj_j in seq.int(traj_i + 1L, length(traj_var))){
        common_records <- intersect(records[[traj_i]], records[[traj_j]])
        if (!length(common_records)){
          stop("Mediation processes '", traj_var[traj_i], "' and '", traj_var[traj_j],
               "' must share at least one record ID to estimate their residual covariance.",
               call. = FALSE)
        }
        pair <- .outcome_pair_index(traj_i, traj_j, length(traj_var))
        COV[[pair]] <- mxPath(from = paste0(traj_var[traj_i], common_records),
                              to = paste0(traj_var[traj_j], common_records),
                              arrows = 2, free = TRUE,
                              values = res_cor[[k]][pair] * sqrt(var0[traj_i] * var0[traj_j]),
                              labels = paste0(pfx[k], traj_var[traj_i], traj_var[traj_j], "_RES"))
      }
    }
    RES_L[[k]] <- RES
    COV_L[[k]] <- COV
    RES <- COV <- list()
  }
  class.list <- list()
  if (x_type == "baseline"){
    if (curveFun %in% c("linear", "LIN")){
      latents <- c("eta0Y", "eta1Y", "eta0M", "eta1M")
      for (k in 1:nClass){
        GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE, values = c(starts[[k]]$Y$means[1:2], starts[[k]]$M$means[1:2]),
                          labels = paste0(pfx[k], c("Y_alpha0", "Y_alpha1", "M_alpha0", "M_alpha1")))
        GF_VAR <- list(mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]]$Y$covMatrix[row(starts[[k]]$Y$covMatrix) >= col(starts[[k]]$Y$covMatrix)],
                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                       mxPath(from = latents[3:4], to = latents[3:4], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]]$M$covMatrix[row(starts[[k]]$M$covMatrix) >= col(starts[[k]]$M$covMatrix)],
                              labels = paste0(pfx[k], c("M_psi00", "M_psi01", "M_psi11"))))
        GF_LOADINGS <- list(mxPath(from = "eta0Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L1", y_records, "Y[1,1]")),
                            mxPath(from = "eta0M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L1", m_records, "M[1,1]")))
        X_BS <- list(mxPath(from = "one", to = x_var, arrows = 1, free = TRUE, values = starts[[k]]$X$means, labels = paste0(pfx[k], "muX")),
                     mxPath(from = x_var, to = x_var, connect = "unique.pairs", arrows = 2, free = TRUE, values = c(starts[[k]]$X$covMatrix),
                            labels = paste0(pfx[k], "phi11")))
        BETA <- list(mxPath(from = x_var, to = latents, arrows = 1, free = TRUE, values = c(starts[[k]]$Y$beta_x, starts[[k]]$M$beta),
                            labels = paste0(pfx[k], "beta", rep(c("Y", "M"), each = 2), rep(c(0, 1), 2))),
                     mxPath(from = latents[3], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[1, 1],
                            labels = paste0(pfx[k], "betaM0Y0")),
                     mxPath(from = latents[3], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[2, 1],
                            labels = paste0(pfx[k], "betaM0Y1")),
                     mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[2, 2],
                            labels = paste0(pfx[k], "betaM1Y1")))
        M_ALPHA <- mxAlgebraFromString(paste0("rbind(", pfx[k], "M_alpha0, ", pfx[k], "M_alpha1", ")"),
                                       name = paste0(pfx[k], "M_alpha"))
        M_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "M_psi00, ", pfx[k], "M_psi01),",
                                              "cbind(", pfx[k], "M_psi01, ", pfx[k], "M_psi11))"),
                                       name = paste0(pfx[k], "M_psi_r"))
        Y_ALPHA <- mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_alpha0, ", pfx[k], "Y_alpha1", ")"),
                                       name = paste0(pfx[k], "Y_alpha"))
        Y_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01),",
                                              "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                       name = paste0(pfx[k], "Y_psi_r"))
        BETA_XM <- mxAlgebraFromString(paste0("rbind(", pfx[k], "betaM0, ", pfx[k], "betaM1)"),
                                       name = paste0(pfx[k], "beta_xm"))
        BETA_XY <- mxAlgebraFromString(paste0("rbind(", pfx[k], "betaY0, ", pfx[k], "betaY1)"),
                                       name = paste0(pfx[k], "beta_xy"))
        BETA_MY <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "betaM0Y0,", 0, "),",
                                              "cbind(", pfx[k], "betaM0Y1, ", pfx[k], "betaM1Y1))"),
                                       name = paste0(pfx[k], "beta_my"))
        M_MEAN0 <- mxAlgebraFromString(paste0(pfx[k], "M_alpha + ", pfx[k], "beta_xm %*% ", pfx[k], "muX"),
                                       name = paste0(pfx[k], "M_mean"))
        Y_MEAN0 <- mxAlgebraFromString(paste0(pfx[k], "Y_alpha + ", pfx[k], "beta_my %*% ", pfx[k], "M_mean + ", pfx[k], "beta_xy %*% ", pfx[k], "muX"),
                                       name = paste0(pfx[k], "Y_mean"))
        MED00 <- mxAlgebraFromString(paste0(pfx[k], "betaM0Y0 * ", pfx[k], "betaM0"),
                                     name = paste0(pfx[k], "mediator_00"))
        MED01 <- mxAlgebraFromString(paste0(pfx[k], "betaM0Y1 * ", pfx[k], "betaM0"),
                                     name = paste0(pfx[k], "mediator_01"))
        MED11 <- mxAlgebraFromString(paste0(pfx[k], "betaM1Y1 * ", pfx[k], "betaM1"),
                                     name = paste0(pfx[k], "mediator_11"))
        MED <- mxAlgebraFromString(paste0("cbind(", pfx[k], "mediator_00, ", pfx[k], "mediator_01, ", pfx[k], "mediator_11)"),
                                   name = paste0(pfx[k], "mediator"))
        TOTAL <- mxAlgebraFromString(paste0(pfx[k], "beta_my %*% ", pfx[k], "beta_xm + ", pfx[k], "beta_xy"),
                                     name = paste0(pfx[k], "total"))
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   GF_loadings[[k]], GF_MEAN, GF_VAR, GF_LOADINGS, X_BS, BETA,
                                   RES_L[[k]], COV_L[[k]], M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r,
                                   BETA_XM, BETA_XY, BETA_MY, M_MEAN0, Y_MEAN0,
                                   MED00, MED01, MED11, MED, TOTAL)
      }
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      latents <- c("eta1Y", "etarY", "eta2Y", "eta1M", "etarM", "eta2M")
      for (k in 1:nClass){
        GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE,
                          values = c(starts[[k]]$Y$means[1:3], starts[[k]]$M$means[1:3]),
                          labels = paste0(pfx[k], c("Y_alpha1", "Y_alphar", "Y_alpha2",
                                                    "M_alpha1", "M_alphar", "M_alpha2")))
        GF_VAR <- list(mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]]$Y$covMatrix[row(starts[[k]]$Y$covMatrix) >= col(starts[[k]]$Y$covMatrix)],
                              labels = paste0(pfx[k], c("Y_psi11", "Y_psi1r", "Y_psi12",
                                                        "Y_psirr", "Y_psir2", "Y_psi22"))),
                       mxPath(from = latents[4:6], to = latents[4:6], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]]$M$covMatrix[row(starts[[k]]$M$covMatrix) >= col(starts[[k]]$M$covMatrix)],
                              labels = paste0(pfx[k], c("M_psi11", "M_psi1r", "M_psi12",
                                                        "M_psirr", "M_psir2", "M_psi22"))))
        GF_LOADINGS <- list(mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L1", y_records, "Y[1,1]")),
                            mxPath(from = "etarY", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta2Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L2", y_records, "Y[1,1]")),
                            mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L1", m_records, "M[1,1]")),
                            mxPath(from = "etarM", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta2M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L2", m_records, "M[1,1]")))
        GAMMA <- list(mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y$means[4],
                               labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                      mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$M$means[4],
                               labels = paste0(pfx[k], "M_knot"), name = paste0(pfx[k], "M_mug")))
        X_BS <- list(mxPath(from = "one", to = x_var, arrows = 1, free = TRUE, values = starts[[k]]$X$means,
                            labels = paste0(pfx[k], "muX")),
                     mxPath(from = x_var, to = x_var, connect = "unique.pairs", arrows = 2, free = TRUE, values = c(starts[[k]]$X$covMatrix),
                            labels = paste0(pfx[k], "phi11")))
        BETA <- list(mxPath(from = x_var, to = latents, arrows = 1, free = TRUE, values = c(starts[[k]]$Y$beta_x, starts[[k]]$M$beta),
                            labels = paste0(pfx[k], "beta", rep(c("Y", "M"), each = 3), rep(c(1, "r", 2), 2))),
                     mxPath(from = latents[4], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[1, 1],
                            labels = paste0(pfx[k], "betaM1Y1")),
                     mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[2, 1],
                            labels = paste0(pfx[k], "betaM1Yr")),
                     mxPath(from = latents[4], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[3, 1],
                            labels = paste0(pfx[k], "betaM1Y2")),
                     mxPath(from = latents[5], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[2, 2],
                            labels = paste0(pfx[k], "betaMrYr")),
                     mxPath(from = latents[5], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[3, 2],
                            labels = paste0(pfx[k], "betaMrY2")),
                     mxPath(from = latents[6], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[3, 3],
                            labels = paste0(pfx[k], "betaM2Y2")))
        M_ALPHA <- mxAlgebraFromString(paste0("rbind(", pfx[k], "M_alpha1, ", pfx[k], "M_alphar, ", pfx[k], "M_alpha2, ", pfx[k], "M_mug)"),
                                       name = paste0(pfx[k], "M_alpha"))
        M_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "M_psi11, ", pfx[k], "M_psi1r, ", pfx[k], "M_psi12),",
                                              "cbind(", pfx[k], "M_psi1r, ", pfx[k], "M_psirr, ", pfx[k], "M_psir2),",
                                              "cbind(", pfx[k], "M_psi12, ", pfx[k], "M_psir2, ", pfx[k], "M_psi22))"),
                                       name = paste0(pfx[k], "M_psi_r"))
        Y_ALPHA <- mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_alpha1, ", pfx[k], "Y_alphar, ", pfx[k], "Y_alpha2, ", pfx[k], "Y_mug)"),
                                       name = paste0(pfx[k], "Y_alpha"))
        Y_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1r, ", pfx[k], "Y_psi12),",
                                              "cbind(", pfx[k], "Y_psi1r, ", pfx[k], "Y_psirr, ", pfx[k], "Y_psir2),",
                                              "cbind(", pfx[k], "Y_psi12, ", pfx[k], "Y_psir2, ", pfx[k], "Y_psi22))"),
                                       name = paste0(pfx[k], "Y_psi_r"))
        BETA_XM <- mxAlgebraFromString(paste0("rbind(", pfx[k], "betaM1, ", pfx[k], "betaMr, ", pfx[k], "betaM2)"),
                                       name = paste0(pfx[k], "beta_xm"))
        BETA_XY <- mxAlgebraFromString(paste0("rbind(", pfx[k], "betaY1, ", pfx[k], "betaYr, ", pfx[k], "betaY2)"),
                                       name = paste0(pfx[k], "beta_xy"))
        BETA_MY <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "betaM1Y1, ", "0,", "0), ",
                                              "cbind(", pfx[k], "betaM1Yr, ", pfx[k], "betaMrYr,", "0),",
                                              "cbind(", pfx[k], "betaM1Y2, ", pfx[k], "betaMrY2, ", pfx[k], "betaM2Y2))"),
                                       name = paste0(pfx[k], "beta_my"))
        M_MEAN0 <- mxAlgebraFromString(paste0(pfx[k], "M_alpha[1:3, ] + ", pfx[k], "beta_xm %*% ", pfx[k], "muX"),
                                       name = paste0(pfx[k], "M_mean"))
        Y_MEAN0 <- mxAlgebraFromString(paste0(pfx[k], "Y_alpha[1:3, ] + ", pfx[k], "beta_my %*% ", pfx[k], "M_mean + ", pfx[k], "beta_xy %*% ", pfx[k], "muX"), name = paste0(pfx[k], "Y_mean"))
        MED11 <- mxAlgebraFromString(paste0(pfx[k], "betaM1Y1 * ", pfx[k], "betaM1"),
                                     name = paste0(pfx[k], "mediator_11"))
        MED1r <- mxAlgebraFromString(paste0(pfx[k], "betaM1Yr * ", pfx[k], "betaM1"),
                                     name = paste0(pfx[k], "mediator_1r"))
        MED12 <- mxAlgebraFromString(paste0(pfx[k], "betaM1Y2 * ", pfx[k], "betaM1"),
                                     name = paste0(pfx[k], "mediator_12"))
        MEDrr <- mxAlgebraFromString(paste0(pfx[k], "betaMrYr * ", pfx[k], "betaMr"),
                                     name = paste0(pfx[k], "mediator_rr"))
        MEDr2 <- mxAlgebraFromString(paste0(pfx[k], "betaMrY2 * ", pfx[k], "betaMr"),
                                     name = paste0(pfx[k], "mediator_r2"))
        MED22 <- mxAlgebraFromString(paste0(pfx[k], "betaM2Y2 * ", pfx[k], "betaM2"),
                                     name = paste0(pfx[k], "mediator_22"))
        MED <- mxAlgebraFromString(paste0("cbind(", pfx[k], "mediator_11, ", pfx[k], "mediator_1r, ", pfx[k], "mediator_12, ", pfx[k], "mediator_rr, ", pfx[k], "mediator_r2, ", pfx[k], "mediator_22)"),
                                   name = paste0(pfx[k], "mediator"))
        TOTAL <- mxAlgebraFromString(paste0(pfx[k], "beta_my %*% ", pfx[k], "beta_xm + ", pfx[k], "beta_xy"),
                                     name = paste0(pfx[k], "total"))
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   GF_loadings[[k]], GF_MEAN, GF_VAR, GF_LOADINGS, GAMMA, X_BS, BETA,
                                   RES_L[[k]], COV_L[[k]], M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r,
                                   BETA_XM, BETA_XY, BETA_MY, M_MEAN0, Y_MEAN0,
                                   MED11, MED1r, MED12, MEDrr, MEDr2, MED22, MED, TOTAL)
      }
    }
  }
  else if (x_type == "longitudinal"){
    if (curveFun %in% c("linear", "LIN")){
      latents <- c("eta0Y", "eta1Y", "eta0M", "eta1M", "eta0X", "eta1X")
      for (k in 1:nClass){
        GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE,
                          values = c(starts[[k]]$Y$means[1:2], starts[[k]]$M$means[1:2], starts[[k]]$X$means[1:2]),
                          labels = paste0(pfx[k], c("Y_alpha0", "Y_alpha1",
                                                    "M_alpha0", "M_alpha1",
                                                    "X_mean0", "X_mean1")))
        GF_VAR <- list(mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]]$Y$covMatrix[row(starts[[k]]$Y$covMatrix) >= col(starts[[k]]$Y$covMatrix)],
                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                       mxPath(from = latents[3:4], to = latents[3:4], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]]$M$covMatrix[row(starts[[k]]$M$covMatrix) >= col(starts[[k]]$M$covMatrix)],
                              labels = paste0(pfx[k], c("M_psi00", "M_psi01", "M_psi11"))),
                       mxPath(from = latents[5:6], to = latents[5:6], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]]$X$covMatrix[row(starts[[k]]$X$covMatrix) >= col(starts[[k]]$X$covMatrix)],
                              labels = paste0(pfx[k], c("X_psi00", "X_psi01", "X_psi11"))))
        GF_LOADINGS <- list(mxPath(from = "eta0Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L1", y_records, "Y[1,1]")),
                            mxPath(from = "eta0M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L1", m_records, "M[1,1]")),
                            mxPath(from = "eta0X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L1", x_records, "X[1,1]")))
        BETA <- list(mxPath(from = latents[5], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_x[1, 1],
                            labels = paste0(pfx[k], "betaX0Y0")),
                     mxPath(from = latents[5], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_x[2, 1],
                            labels = paste0(pfx[k], "betaX0Y1")),
                     mxPath(from = latents[6], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_x[2, 2],
                            labels = paste0(pfx[k], "betaX1Y1")),
                     mxPath(from = latents[5], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]]$M$beta[1, 1],
                            labels = paste0(pfx[k], "betaX0M0")),
                     mxPath(from = latents[5], to = latents[4], arrows = 1, free = TRUE, values = starts[[k]]$M$beta[2, 1],
                            labels = paste0(pfx[k], "betaX0M1")),
                     mxPath(from = latents[6], to = latents[4], arrows = 1, free = TRUE, values = starts[[k]]$M$beta[2, 2],
                            labels = paste0(pfx[k], "betaX1M1")),
                     mxPath(from = latents[3], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[1, 1],
                            labels = paste0(pfx[k], "betaM0Y0")),
                     mxPath(from = latents[3], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[2, 1],
                            labels = paste0(pfx[k], "betaM0Y1")),
                     mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[2, 2],
                            labels = paste0(pfx[k], "betaM1Y1")))
        X_MEAN <- mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mean0, ", pfx[k], "X_mean1)"),
                                      name = paste0(pfx[k], "X_mean"))
        X_PSI <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01), ",
                                            "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                     name = paste0(pfx[k], "X_psi0"))
        M_ALPHA <- mxAlgebraFromString(paste0("rbind(", pfx[k], "M_alpha0, ", pfx[k], "M_alpha1)"),
                                       name = paste0(pfx[k], "M_alpha"))
        M_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "M_psi00, ", pfx[k], "M_psi01), ",
                                              "cbind(", pfx[k], "M_psi01, ", pfx[k], "M_psi11))"),
                                       name = paste0(pfx[k], "M_psi_r"))
        Y_ALPHA <- mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_alpha0, ", pfx[k], "Y_alpha1)"),
                                       name = paste0(pfx[k], "Y_alpha"))
        Y_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                              "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                       name = paste0(pfx[k], "Y_psi_r"))
        BETA_XM <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "betaX0M0, 0), ",
                                              "cbind(", pfx[k], "betaX0M1, ", pfx[k], "betaX1M1))"),
                                       name = paste0(pfx[k], "beta_xm"))
        BETA_XY <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "betaX0Y0, 0), ",
                                              "cbind(", pfx[k], "betaX0Y1, ", pfx[k], "betaX1Y1))"),
                                       name = paste0(pfx[k], "beta_xy"))
        BETA_MY <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "betaM0Y0, 0), ",
                                              "cbind(", pfx[k], "betaM0Y1, ", pfx[k], "betaM1Y1))"),
                                       name = paste0(pfx[k], "beta_my"))
        M_MEAN0 <- mxAlgebraFromString(paste0(pfx[k], "M_alpha + ", pfx[k], "beta_xm %*% ", pfx[k], "X_mean"),
                                       name = paste0(pfx[k], "M_mean"))
        Y_MEAN0 <- mxAlgebraFromString(paste0(pfx[k], "Y_alpha + ", pfx[k], "beta_my %*% ", pfx[k], "M_mean + ", pfx[k], "beta_xy %*% ", pfx[k], "X_mean"), name = paste0(pfx[k], "Y_mean"))
        MED000 <- mxAlgebraFromString(paste0(pfx[k], "betaX0M0 * ", pfx[k], "betaM0Y0"),
                                      name = paste0(pfx[k], "mediator_000"))
        MED001 <- mxAlgebraFromString(paste0(pfx[k], "betaX0M0 * ", pfx[k], "betaM0Y1"),
                                      name = paste0(pfx[k], "mediator_001"))
        MED011 <- mxAlgebraFromString(paste0(pfx[k], "betaX0M1 * ", pfx[k], "betaM1Y1"),
                                      name = paste0(pfx[k], "mediator_011"))
        MED111 <- mxAlgebraFromString(paste0(pfx[k], "betaX1M1 * ", pfx[k], "betaM1Y1"),
                                      name = paste0(pfx[k], "mediator_111"))
        MED <- mxAlgebraFromString(paste0("cbind(", pfx[k], "mediator_000, ", pfx[k], "mediator_001, ", pfx[k], "mediator_011, ", pfx[k], "mediator_111)"),
                                   name = paste0(pfx[k], "mediator"))
        TOTAL <- mxAlgebraFromString(paste0(pfx[k], "beta_my %*% ", pfx[k], "beta_xm + ", pfx[k], "beta_xy"),
                                     name = paste0(pfx[k], "total"))
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   GF_loadings[[k]], GF_MEAN, GF_VAR, GF_LOADINGS, BETA,
                                   RES_L[[k]], COV_L[[k]], X_MEAN, X_PSI, M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r,
                                   BETA_XM, BETA_XY, BETA_MY, M_MEAN0, Y_MEAN0,
                                   MED000, MED001, MED011, MED111, MED, TOTAL)
      }
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      latents <- c("eta1Y", "etarY", "eta2Y", "eta1M", "etarM", "eta2M", "eta1X", "etarX", "eta2X")
      for (k in 1:nClass){
        GF_MEAN <- mxPath(from = "one", to = latents, arrows = 1, free = TRUE,
                          values = c(starts[[k]]$Y$means[1:3], starts[[k]]$M$means[1:3], starts[[k]]$X$means[1:3]),
                          labels = paste0(pfx[k], c("Y_alpha1", "Y_alphar", "Y_alpha2",
                                                    "M_alpha1", "M_alphar", "M_alpha2",
                                                    "X_mean1", "X_meanr", "X_mean2")))
        GF_VAR <- list(mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]]$Y$covMatrix[row(starts[[k]]$Y$covMatrix) >= col(starts[[k]]$Y$covMatrix)],
                              labels = paste0(pfx[k], c("Y_psi11", "Y_psi1r", "Y_psi12",
                                                        "Y_psirr", "Y_psir2",
                                                        "Y_psi22"))),
                       mxPath(from = latents[4:6], to = latents[4:6], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]]$M$covMatrix[row(starts[[k]]$M$covMatrix) >= col(starts[[k]]$M$covMatrix)],
                              labels = paste0(pfx[k], c("M_psi11", "M_psi1r", "M_psi12",
                                                        "M_psirr", "M_psir2",
                                                        "M_psi22"))),
                       mxPath(from = latents[7:9], to = latents[7:9], arrows = 2, connect = "unique.pairs",
                              free = TRUE, values = starts[[k]]$X$covMatrix[row(starts[[k]]$X$covMatrix) >= col(starts[[k]]$X$covMatrix)],
                              labels = paste0(pfx[k], c("X_psi11", "X_psi1r", "X_psi12",
                                                        "X_psirr", "X_psir2",
                                                        "X_psi22"))))
        GF_LOADINGS <- list(mxPath(from = "eta1Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L1", y_records, "Y[1,1]")),
                            mxPath(from = "etarY", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta2Y", to = paste0(y_var, y_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L2", y_records, "Y[1,1]")),
                            mxPath(from = "eta1M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L1", m_records, "M[1,1]")),
                            mxPath(from = "etarM", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta2M", to = paste0(m_var, m_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L2", m_records, "M[1,1]")),
                            mxPath(from = "eta1X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L1", x_records, "X[1,1]")),
                            mxPath(from = "etarX", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta2X", to = paste0(x_var, x_records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(pfx[k], "L2", x_records, "X[1,1]")))
        GAMMA <- list(mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y$means[4],
                               labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                      mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$M$means[4],
                               labels = paste0(pfx[k], "M_knot"), name = paste0(pfx[k], "M_mug")),
                      mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$X$means[4],
                               labels = paste0(pfx[k], "X_knot"), name = paste0(pfx[k], "X_mug")))
        BETA <- list(mxPath(from = latents[7], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_x[1, 1],
                            labels = paste0(pfx[k], "betaX1Y1")),
                     ##### Slope 1 of X to intercept of Y
                     mxPath(from = latents[7], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_x[2, 1],
                            labels = paste0(pfx[k], "betaX1Yr")),
                     ##### Slope 1 of X to slope 2 of Y
                     mxPath(from = latents[7], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_x[3, 1],
                            labels = paste0(pfx[k], "betaX1Y2")),
                     ##### Intercept of X to intercept of Y
                     mxPath(from = latents[8], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_x[2, 2],
                            labels = paste0(pfx[k], "betaXrYr")),
                     ##### Intercept of X to slope 2 of Y
                     mxPath(from = latents[8], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_x[3, 2],
                            labels = paste0(pfx[k], "betaXrY2")),
                     ##### Slope 2 of X to slope 2 of Y
                     mxPath(from = latents[9], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_x[3, 3],
                            labels = paste0(pfx[k], "betaX2Y2")),
                     ##### Slope 1 of X to slope 1 of M
                     mxPath(from = latents[7], to = latents[4], arrows = 1, free = TRUE, values = starts[[k]]$M$beta[1, 1],
                            labels = paste0(pfx[k], "betaX1M1")),
                     ##### Slope 1 of X to intercept of M
                     mxPath(from = latents[7], to = latents[5], arrows = 1, free = TRUE, values = starts[[k]]$M$beta[2, 1],
                            labels = paste0(pfx[k], "betaX1Mr")),
                     ##### Slope 1 of X to slope 2 of M
                     mxPath(from = latents[7], to = latents[6], arrows = 1, free = TRUE, values = starts[[k]]$M$beta[3, 1],
                            labels = paste0(pfx[k], "betaX1M2")),
                     ##### Intercept of X to intercept of M
                     mxPath(from = latents[8], to = latents[5], arrows = 1, free = TRUE, values = starts[[k]]$M$beta[2, 2],
                            labels = paste0(pfx[k], "betaXrMr")),
                     ##### Intercept of X to slope 2 of M
                     mxPath(from = latents[8], to = latents[6], arrows = 1, free = TRUE, values = starts[[k]]$M$beta[3, 2],
                            labels = paste0(pfx[k], "betaXrM2")),
                     ##### Slope 2 of X to slope 2 of M
                     mxPath(from = latents[9], to = latents[6], arrows = 1, free = TRUE, values = starts[[k]]$M$beta[3, 3],
                            labels = paste0(pfx[k], "betaX2M2")),

                     ##### Slope 1 of M to slope 1 of Y
                     mxPath(from = latents[4], to = latents[1], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[1, 1],
                            labels = paste0(pfx[k], "betaM1Y1")),
                     ##### Slope 1 of M to intercept of Y
                     mxPath(from = latents[4], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[2, 1],
                            labels = paste0(pfx[k], "betaM1Yr")),
                     ##### Slope 1 of M to slope 2 of Y
                     mxPath(from = latents[4], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[3, 1],
                            labels = paste0(pfx[k], "betaM1Y2")),
                     ##### Intercept of M to intercept of Y
                     mxPath(from = latents[5], to = latents[2], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[2, 2],
                            labels = paste0(pfx[k], "betaMrYr")),
                     ##### Intercept of M to slope 2 of Y
                     mxPath(from = latents[5], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[3, 2],
                            labels = paste0(pfx[k], "betaMrY2")),
                     ##### Slope 2 of M to slope 2 of Y
                     mxPath(from = latents[6], to = latents[3], arrows = 1, free = TRUE, values = starts[[k]]$Y$beta_m[3, 3],
                            labels = paste0(pfx[k], "betaM2Y2")))
        X_MEAN <- mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mean1, ", pfx[k], "X_meanr, ", pfx[k], "X_mean2, ", pfx[k], "X_mug)"),
                                      name = paste0(pfx[k], "X_mean"))
        X_PSI <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi11, ", pfx[k], "X_psi1r, ", pfx[k], "X_psi12), ",
                                            "cbind(", pfx[k], "X_psi1r, ", pfx[k], "X_psirr, ", pfx[k], "X_psir2), ",
                                            "cbind(", pfx[k], "X_psi12, ", pfx[k], "X_psir2, ", pfx[k], "X_psi22))"),
                                     name = paste0(pfx[k], "X_psi0"))
        M_ALPHA <- mxAlgebraFromString(paste0("rbind(", pfx[k], "M_alpha1, ", pfx[k], "M_alphar, ", pfx[k], "M_alpha2, ", pfx[k], "M_mug)"),
                                       name = paste0(pfx[k], "M_alpha"))
        M_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "M_psi11, ", pfx[k], "M_psi1r, ", pfx[k], "M_psi12), ",
                                              "cbind(", pfx[k], "M_psi1r, ", pfx[k], "M_psirr, ", pfx[k], "M_psir2), ",
                                              "cbind(", pfx[k], "M_psi12, ", pfx[k], "M_psir2, ", pfx[k], "M_psi22))"),
                                       name = paste0(pfx[k], "M_psi_r"))
        Y_ALPHA <- mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_alpha1, ", pfx[k], "Y_alphar, ", pfx[k], "Y_alpha2, ", pfx[k], "Y_mug)"),
                                       name = paste0(pfx[k], "Y_alpha"))
        Y_PSI_r <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1r, ", pfx[k], "Y_psi12), ",
                                              "cbind(", pfx[k], "Y_psi1r, ", pfx[k], "Y_psirr, ", pfx[k], "Y_psir2), ",
                                              "cbind(", pfx[k], "Y_psi12, ", pfx[k], "Y_psir2, ", pfx[k], "Y_psi22))"),
                                       name = paste0(pfx[k], "Y_psi_r"))
        BETA_XM <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "betaX1M1, 0, 0), ",
                                              "cbind(", pfx[k], "betaX1Mr, ", pfx[k], "betaXrMr, 0), ",
                                              "cbind(", pfx[k], "betaX1M2, ", pfx[k], "betaXrM2, ", pfx[k], "betaX2M2))"),
                                       name = paste0(pfx[k], "beta_xm"))
        BETA_XY <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "betaX1Y1, 0, 0), ",
                                              "cbind(", pfx[k], "betaX1Yr, ", pfx[k], "betaXrYr, 0), ",
                                              "cbind(", pfx[k], "betaX1Y2, ", pfx[k], "betaXrY2, ", pfx[k], "betaX2Y2))"),
                                       name = paste0(pfx[k], "beta_xy"))
        BETA_MY <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "betaM1Y1, 0, 0), ",
                                              "cbind(", pfx[k], "betaM1Yr, ", pfx[k], "betaMrYr, 0), ",
                                              "cbind(", pfx[k], "betaM1Y2, ", pfx[k], "betaMrY2, ", pfx[k], "betaM2Y2))"),
                                       name = paste0(pfx[k], "beta_my"))
        M_MEAN0 <- mxAlgebraFromString(paste0(pfx[k], "M_alpha[1:3, ] + ", pfx[k], "beta_xm %*% ", pfx[k], "X_mean[1:3, ]"),
                                       name = paste0(pfx[k], "M_mean"))
        Y_MEAN0 <- mxAlgebraFromString(paste0(pfx[k], "Y_alpha[1:3, ] + ", pfx[k], "beta_my %*% ", pfx[k], "M_mean + ", pfx[k], "beta_xy %*% ", pfx[k], "X_mean[1:3, ]"),name = paste0(pfx[k], "Y_mean"))
        MED111 <- mxAlgebraFromString(paste0(pfx[k], "betaX1M1 * ", pfx[k], "betaM1Y1"), name = paste0(pfx[k], "mediator_111"))
        MED11r <- mxAlgebraFromString(paste0(pfx[k], "betaX1M1 * ", pfx[k], "betaM1Yr"), name = paste0(pfx[k], "mediator_11r"))
        MED112 <- mxAlgebraFromString(paste0(pfx[k], "betaX1M1 * ", pfx[k], "betaM1Y2"), name = paste0(pfx[k], "mediator_112"))
        MED1rr <- mxAlgebraFromString(paste0(pfx[k], "betaX1Mr * ", pfx[k], "betaMrYr"), name = paste0(pfx[k], "mediator_1rr"))
        MED1r2 <- mxAlgebraFromString(paste0(pfx[k], "betaX1Mr * ", pfx[k], "betaMrY2"), name = paste0(pfx[k], "mediator_1r2"))
        MED122 <- mxAlgebraFromString(paste0(pfx[k], "betaX1M2 * ", pfx[k], "betaM2Y2"), name = paste0(pfx[k], "mediator_122"))
        MEDrrr <- mxAlgebraFromString(paste0(pfx[k], "betaXrMr * ", pfx[k], "betaMrYr"), name = paste0(pfx[k], "mediator_rrr"))
        MEDrr2 <- mxAlgebraFromString(paste0(pfx[k], "betaXrMr * ", pfx[k], "betaMrY2"), name = paste0(pfx[k], "mediator_rr2"))
        MEDr22 <- mxAlgebraFromString(paste0(pfx[k], "betaXrM2 * ", pfx[k], "betaM2Y2"), name = paste0(pfx[k], "mediator_r22"))
        MED222 <- mxAlgebraFromString(paste0(pfx[k], "betaX2M2 * ", pfx[k], "betaM2Y2"), name = paste0(pfx[k], "mediator_222"))
        MED <- mxAlgebraFromString(paste0("cbind(", pfx[k], "mediator_111, ", pfx[k], "mediator_11r, ", pfx[k], "mediator_112, ", pfx[k], "mediator_1rr, ", pfx[k], "mediator_1r2, ", pfx[k], "mediator_122, ", pfx[k], "mediator_rrr, ", pfx[k], "mediator_rr2, ", pfx[k], "mediator_r22, ", pfx[k], "mediator_222)"), name = paste0(pfx[k], "mediator"))
        TOTAL <- mxAlgebraFromString(paste0(pfx[k], "beta_my %*% ", pfx[k], "beta_xm + ", pfx[k], "beta_xy"),
                                     name = paste0(pfx[k], "total"))
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                    manifestVars = manifests, latentVars = latents,
                                   GF_loadings[[k]], GF_MEAN, GF_VAR, GF_LOADINGS, GAMMA, BETA,
                                   RES_L[[k]], COV_L[[k]], X_MEAN, X_PSI, M_ALPHA, M_PSI_r, Y_ALPHA, Y_PSI_r,
                                   BETA_XM, BETA_XY, BETA_MY, M_MEAN0, Y_MEAN0,
                                   MED111, MED11r, MED112, MED1rr, MED1r2, MED122, MEDrrr, MEDrr2,
                                   MEDr22, MED222, MED, TOTAL)
      }
    }
  }
  for (k in seq_len(nClass)) class.list[[k]] <- do.call(mxModel, c(list(class.list[[k]]), extra[[k]]))
  class.list
}






#' @keywords internal
#' @noRd
#' @importFrom OpenMx mxData
getMED.mxModel <- function(dat, t_var, y_var, m_var, x_type, x_var, curveFun, records, starts, res_cor){
  loadings <- if (x_type == "baseline"){
    getMED.loadings(t_var = t_var, y_var = y_var, m_var = m_var, x_type = x_type, x_var = x_var,
                    curveFun = curveFun, y_records = records[[1]], m_records = records[[2]])
  } else{
    getMED.loadings(t_var = t_var, y_var = y_var, m_var = m_var, x_type = x_type, x_var = x_var,
                    curveFun = curveFun, y_records = records[[1]], m_records = records[[2]], x_records = records[[3]])
  }
  .build_MED(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, records = records, m_var = m_var,
             x_var = x_var, x_type = x_type, starts = list(starts),
             res_cor = if (is.null(res_cor)) NULL else list(res_cor), pfx = "",
             model_names = .single_model_name_med(x_type, curveFun),
             extra = list(list(mxData(observed = dat, type = "raw"))), GF_loadings = list(loadings))[[1]]
}

#' @keywords internal
#' @noRd
getsub.MED <- function(dat, nClass, t_var, y_var, curveFun, records, m_var, x_var, x_type, starts, res_cor, container,
                       grp_var = NULL){
  loadings <- if (x_type == "baseline"){
    getMIX_MED.loadings(nClass = nClass, t_var = t_var, y_var = y_var, m_var = m_var, x_type = x_type, x_var = x_var,
                        curveFun = curveFun, y_records = records[[1]], m_records = records[[2]])
  } else{
    getMIX_MED.loadings(nClass = nClass, t_var = t_var, y_var = y_var, m_var = m_var, x_type = x_type, x_var = x_var,
                        curveFun = curveFun, y_records = records[[1]], m_records = records[[2]], x_records = records[[3]])
  }
  .build_MED(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, records = records, m_var = m_var,
             x_var = x_var, x_type = x_type, starts = starts, res_cor = res_cor, pfx = paste0("c", seq_len(nClass)),
             model_names = paste0("Class", seq_len(nClass)), extra = .class_extras(dat, nClass, container, grp_var),
             GF_loadings = loadings)
}

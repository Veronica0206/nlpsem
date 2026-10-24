#' @title One builder for the MGM sub-model, single-group or class-specific
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
.build_MGM <- function(dat, t_var, y_var, curveFun, intrinsic, records, y_model, starts, pfx, model_names, extra, GF_loadings, AddPara){
  nClass <- length(pfx)
  y_var_L <- tolower(y_var)
  # Define manifest variables
  traj_list <- list()
  for (traj in seq_along(y_var)){
    traj_list[[length(traj_list) + 1]] <- paste0(y_var[traj], records[[traj]])
  }
  manifests <- unlist(traj_list)
  # Define variances and covariance(s) of residuals
  RES_L <- COV_L <- list()
  for (k in 1:nClass){
    RES <- COV <- list()
    for (traj in seq_along(y_var)){
      RES[[traj]] <- mxPath(from = traj_list[[traj]], to = traj_list[[traj]], arrows = 2, free = TRUE,
                            values = starts[[k]]$residuals[traj, traj],
                            labels = paste0(pfx[k], y_var[traj], "_residuals"))
    }
    for (traj_i in 1:(length(y_var) - 1)){
      for (traj_j in traj_i:(length(y_var) - 1)){
        common_records <- Reduce(intersect, list(records[[traj_i]], records[[traj_j + 1]]))
        if (setequal(records[[traj_i]], records[[traj_j + 1]])){
          COV[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxPath(from = traj_list[[traj_i]], to = traj_list[[traj_j + 1]],
                                               arrows = 2, free = TRUE,
                                               values = starts[[k]]$residuals[traj_i, traj_j + 1],
                                               labels = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_RES"))
        }
        else{
          COV[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxPath(from = paste0(y_var[traj_i], common_records),
                                               to = paste0(y_var[traj_j + 1], common_records),
                                               arrows = 2, free = TRUE, values = starts[[k]]$residuals[traj_i, traj_j + 1],
                                               labels = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_RES"))
        }
      }
    }
    RES_L[[k]] <- RES
    COV_L[[k]] <- COV
  }
  # Obtain factor loadings for the specified trajectory model (either "LGCM" or "LCSM") and functional form
  class.list <- list()
  if (y_model == "LGCM"){
    if (curveFun %in% c("linear", "LIN")){
      latents <- paste0(rep(c("eta0", "eta1"), length(y_var)), rep(y_var, each = 2))
      for (k in 1:nClass){
        multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
        GF_MEAN <- GF_LOADINGS <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in seq_along(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]), arrows = 1, free = TRUE,
                                    values = starts[[k]]$means[[traj]][1:2],
                                    labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1")))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]], arrows = 1,
                                             free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]], arrows = 1,
                                             free = FALSE, labels = paste0(pfx[k], "L1", records[[traj]], y_var[traj], "[1,1]")))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj], "_mueta1)"),
                                                name = paste0(pfx[k], y_var[traj], "_mean0"))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj], "_psi01), ",
                                                      "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi11))"),
                                               name = paste0(pfx[k], y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "01", "11")),
                                         byrow = TRUE, nrow = 2, ncol = 2)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01), ",
                                                                           "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11))"),
                                                                    name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi", c("00", "01", "10", "11")),
                                                          byrow = TRUE, nrow = 2, ncol = 2)
          }
        }
        for (traj in seq_along(y_var)){
          multi_label[((traj - 1) * 2 + 1):((traj - 1) * 2 + 2), ((traj - 1) * 2 + 1):((traj - 1) * 2 + 2)] <- gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 2 + 1):((traj_i - 1) * 2 + 2), (traj_j * 2 + 1):(traj_j * 2 + 2)] <- gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
          }
        }
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = t(starts[[k]]$covMatrix)[row(t(starts[[k]]$covMatrix)) >= col(t(starts[[k]]$covMatrix))],
                                          labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                   GF_loadings[[k]], GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]], mean_L, psi_L,
                                   psi_btw_L)
      }
    }
    else if (curveFun %in% c("quadratic", "QUAD")){
      latents <- paste0(rep(c("eta0", "eta1", "eta2"), length(y_var)), rep(y_var, each = 3))
      for (k in 1:nClass){
        multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
        GF_MEAN <- GF_LOADINGS <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in seq_along(y_var)) {
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:3],
                                    labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE,
                                             labels = paste0(pfx[k], "L1", records[[traj]], y_var[traj], "[1,1]")),
                                      mxPath(from = paste0("eta2", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE,
                                             labels = paste0(pfx[k], "L2", records[[traj]], y_var[traj], "[1,1]")))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj],
                                                       "_mueta1, ", pfx[k], y_var[traj], "_mueta2)"),
                                                name = paste0(pfx[k], y_var[traj], "_mean0"))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi02), ",
                                                      "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi11, ", pfx[k], y_var[traj], "_psi12), ",
                                                      "cbind(", pfx[k], y_var[traj], "_psi02, ", pfx[k], y_var[traj], "_psi12, ", pfx[k], y_var[traj], "_psi22))"),
                                               name = paste0(pfx[k], y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "02", "01", "11", "12", "02", "12", "22")),
                                         byrow = TRUE, nrow = 3, ncol = 3)
        }
        for (traj_i in 1:(length(y_var) - 1)) {
          for (traj_j in traj_i:(length(y_var) - 1)) {
            psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi02), ",
                                                                           "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi12), ",
                                                                           "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi20, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi21, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi22))"),
                                                                    name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi", c("00", "01", "02", "10", "11", "12", "20", "21", "22")),
                                                          byrow = TRUE, nrow = 3, ncol = 3)
          }
        }
        for (traj in seq_along(y_var)) {
          multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <- gf_var_label[[traj]]
        }

        for (traj_i in 1:(length(y_var) - 1)) {
          for (traj_j in traj_i:(length(y_var) - 1)) {
            multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <- gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
          }
        }
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = t(starts[[k]]$covMatrix)[row(t(starts[[k]]$covMatrix)) >= col(t(starts[[k]]$covMatrix))],
                                          labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                   GF_loadings[[k]], GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]], mean_L, psi_L,
                                   psi_btw_L)
      }
    }
    else if (curveFun %in% c("negative exponential", "EXP")){
      if (intrinsic){
        latents <- paste0(rep(c("eta0", "eta1", "deltag"), length(y_var)), rep(y_var, each = 3))
        for (k in 1:nClass){
          multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
          GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
          mean_L <- psi_L <- list()
          gf_var_label <- gf_cov_label <- list()
          psi_btw_L <- list()
          for (traj in seq_along(y_var)){
            GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                      arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:2],
                                      labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1")))
            GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$means[[traj]][3],
                                      labels = paste0(pfx[k], y_var[traj], "_slp_ratio"),
                                      name = paste0(pfx[k], y_var[traj], "_mug"))
            mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj], "_mueta1, ",
                                                         pfx[k], y_var[traj], "_slp_ratio)"),
                                                  name = paste0(pfx[k], y_var[traj], "_mean0"))
            GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 1),
                                        mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE,
                                               labels = paste0(pfx[k], "L1", records[[traj]], y_var[traj], "[1,1]")),
                                        mxPath(from = paste0("deltag", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L2", records[[traj]], y_var[traj], "[1,1]")))
            psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj],
                                                        "_psi01, ", pfx[k], y_var[traj], "_psi0g), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj],
                                                        "_psi11, ", pfx[k], y_var[traj], "_psi1g), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi0g, ", pfx[k], y_var[traj],
                                                        "_psi1g, ", pfx[k], y_var[traj], "_psigg))"),
                                                 name = paste0(pfx[k], y_var[traj], "_psi0"))
            gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "0g", "01", "11", "1g", "0g", "1g", "gg")),
                                           byrow = TRUE, nrow = 3, ncol = 3)
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1],
                                                                             "_psi01, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi0g),",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi1g), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig0, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig1, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psigg))"),
                                                                      name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
              gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                   c("00", "01", "0g", "10", "11", "1g", "g0", "g1", "gg")),
                                                            byrow = TRUE, nrow = 3, ncol = 3)
            }
          }
          for (traj in seq_along(y_var)){
            multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
              gf_var_label[[traj]]
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
                gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
            }
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = t(starts[[k]]$covMatrix)[row(t(starts[[k]]$covMatrix)) >= col(t(starts[[k]]$covMatrix))],
                                            labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                     GF_loadings[[k]], GAMMA, GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]],
                                     mean_L, psi_L, psi_btw_L)
        }
      }
      else if (!intrinsic){
        latents <- paste0(rep(c("eta0", "eta1"), length(y_var)), rep(y_var, each = 2))
        for (k in 1:nClass){
          multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
          GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
          mean_L <- psi_L <- list()
          gf_var_label <- gf_cov_label <- list()
          psi_btw_L <- list()
          for (traj in seq_along(y_var)){
            GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                      arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:2],
                                      labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1")))
            GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$means[[traj]][3],
                                      labels = paste0(pfx[k], y_var[traj], "_slp_ratio"),
                                      name = paste0(pfx[k], y_var[traj], "_mug"))
            mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj], "_mueta1, ",
                                                         pfx[k], y_var[traj], "_slp_ratio)"),
                                                  name = paste0(pfx[k], y_var[traj], "_mean0"))
            GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 1),
                                        mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE,
                                               labels = paste0(pfx[k], "L1", records[[traj]], y_var[traj], "[1,1]")))
            psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj], "_psi01), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi11))"),
                                                 name = paste0(pfx[k], y_var[traj], "_psi0"))
            gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "01", "11")),
                                           byrow = TRUE, nrow = 2, ncol = 2)
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01),",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11))"),
                                                                      name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
              gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                   c("00", "01", "10", "11")),
                                                            byrow = TRUE, nrow = 2, ncol = 2)
            }
          }
          for (traj in seq_along(y_var)){
            multi_label[((traj - 1) * 2 + 1):((traj - 1) * 2 + 2), ((traj - 1) * 2 + 1):((traj - 1) * 2 + 2)] <-
              gf_var_label[[traj]]
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              multi_label[((traj_i - 1) * 2 + 1):((traj_i - 1) * 2 + 2), (traj_j * 2 + 1):(traj_j * 2 + 2)] <-
                gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
            }
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs", free = TRUE,
                                            values = t(starts[[k]]$covMatrix[c(1:2, 4:5), c(1:2, 4:5)])[row(t(starts[[k]]$covMatrix[c(1:2, 4:5), c(1:2, 4:5)])) >=
                                                                                                     col(t(starts[[k]]$covMatrix[c(1:2, 4:5), c(1:2, 4:5)]))],
                                            labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                     GF_loadings[[k]], GAMMA, GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]],
                                     mean_L, psi_L, psi_btw_L)
        }
      }
    }
    else if (curveFun %in% c("Jenss-Bayley", "JB")){
      if (intrinsic){
        latents <- paste0(rep(c("eta0", "eta1", "eta2", "deltag"), length(y_var)), rep(y_var, each = 4))
        for (k in 1:nClass){
          multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
          GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
          mean_L <- psi_L <- list()
          gf_var_label <- gf_cov_label <- list()
          psi_btw_L <- list()
          for (traj in seq_along(y_var)){
            GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                      arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:3],
                                      labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
            GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$means[[traj]][4],
                                      labels = paste0(pfx[k], y_var[traj], "_acc_ratio"),
                                      name = paste0(pfx[k], y_var[traj], "_mug"))
            mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj], "_mueta1, ", pfx[k], y_var[traj], "_mueta2, ", pfx[k], y_var[traj], "_acc_ratio)"),
                                                  name = paste0(pfx[k], y_var[traj], "_mean0"))
            GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 1),
                                        mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L1", records[[traj]], y_var[traj], "[1,1]")),
                                        mxPath(from = paste0("eta2", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L2", records[[traj]], y_var[traj], "[1,1]")),
                                        mxPath(from = paste0("deltag", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L3", records[[traj]], y_var[traj], "[1,1]")))
            psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi02, ", pfx[k], y_var[traj], "_psi0g), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi11, ", pfx[k], y_var[traj], "_psi12, ", pfx[k], y_var[traj], "_psi1g), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi02, ", pfx[k], y_var[traj], "_psi12, ", pfx[k], y_var[traj], "_psi22, ", pfx[k], y_var[traj], "_psi2g), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi0g, ", pfx[k], y_var[traj], "_psi1g, ", pfx[k], y_var[traj], "_psi2g, ", pfx[k], y_var[traj], "_psigg))"),
                                                 name = paste0(pfx[k], y_var[traj], "_psi0"))
            gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "02", "0g", "01", "11", "12", "1g",
                                                                                 "02", "12", "22", "2g", "0g", "1g", "2g", "gg")),
                                           byrow = TRUE, nrow = 4, ncol = 4)
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi02, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi0g), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi12, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi1g), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi20, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi21, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi22, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi2g), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig0, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig1, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig2, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psigg))"),
                                                                      name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
              gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                   c("00", "01", "02", "0g", "10", "11", "12", "1g",
                                                                     "20", "21", "22", "2g", "g0", "g1", "g2", "gg")),
                                                            byrow = TRUE, nrow = 4, ncol = 4)
            }
          }
          for (traj in seq_along(y_var)){
            multi_label[((traj - 1) * 4 + 1):((traj - 1) * 4 + 4), ((traj - 1) * 4 + 1):((traj - 1) * 4 + 4)] <-
              gf_var_label[[traj]]
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              multi_label[((traj_i - 1) * 4 + 1):((traj_i - 1) * 4 + 4), (traj_j * 4 + 1):(traj_j * 4 + 4)] <-
                gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
            }
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = t(starts[[k]]$covMatrix)[row(t(starts[[k]]$covMatrix)) >= col(t(starts[[k]]$covMatrix))],
                                            labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                     GF_loadings[[k]], GAMMA, GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]],
                                     mean_L, psi_L, psi_btw_L)
        }
      }
      else if (!intrinsic){
        latents <- paste0(rep(c("eta0", "eta1", "eta2"), length(y_var)), rep(y_var, each = 3))
        for (k in 1:nClass){
          multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
          GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
          mean_L <- psi_L <- list()
          gf_var_label <- gf_cov_label <- list()
          psi_btw_L <- list()
          for (traj in seq_along(y_var)){
            GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                      arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:3],
                                      labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
            GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$means[[traj]][4],
                                      labels = paste0(pfx[k], y_var[traj], "_acc_ratio"),
                                      name = paste0(pfx[k], y_var[traj], "_mug"))
            mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj], "_mueta1, ",
                                                         pfx[k], y_var[traj], "_mueta2, ", pfx[k], y_var[traj], "_acc_ratio)"),
                                                  name = paste0(pfx[k], y_var[traj], "_mean0"))
            GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 1),
                                        mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L1", records[[traj]], y_var[traj], "[1,1]")),
                                        mxPath(from = paste0("eta2", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L2", records[[traj]], y_var[traj], "[1,1]")))
            psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj], "_psi01, ",
                                                        pfx[k], y_var[traj], "_psi02), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi11, ",
                                                        pfx[k], y_var[traj], "_psi12), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi02, ", pfx[k], y_var[traj], "_psi12, ",
                                                        pfx[k], y_var[traj], "_psi22))"),
                                                 name = paste0(pfx[k], y_var[traj], "_psi0"))
            gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "02", "01", "11", "12",
                                                                                 "02", "12", "22")),
                                           byrow = TRUE, nrow = 3, ncol = 3)
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi02), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi12), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi20, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi21, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi22))"),
                                                                      name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
              gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                   c("00", "01", "02", "10", "11", "12",
                                                                     "20", "21", "22")),
                                                            byrow = TRUE, nrow = 3, ncol = 3)
            }
          }
          for (traj in seq_along(y_var)){
            multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
              gf_var_label[[traj]]
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
                gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
            }
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs", free = TRUE,
                                            values = .multi_reduced_cov_starts(starts[[k]]$covMatrix, n_traj = length(y_var)),
                                            labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                     GF_loadings[[k]], GAMMA, GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]],
                                     mean_L, psi_L, psi_btw_L)
        }
      }
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      if (intrinsic){
        latents <- paste0(rep(c("eta0s", "eta1s", "eta2s", "deltag"), length(y_var)), rep(y_var, each = 4))
        for (k in 1:nClass){
          multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
          GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
          func_L <- grad_L <- mean_s_L <- psi_s_L <- mean_L <- psi_L <- list()
          gf_var_label <- gf_cov_label <- list()
          psi_btw_s_L <- psi_btw_L <- list()
          for (traj in seq_along(y_var)){
            GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0s", "eta1s", "eta2s"), y_var[traj]),
                                      arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:3],
                                      labels = paste0(pfx[k], y_var[traj], c("_mueta0s", "_mueta1s", "_mueta2s")))
            GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$means[[traj]][4],
                                      labels = paste0(pfx[k], y_var[traj], "_knot"),
                                      name = paste0(pfx[k], y_var[traj], "_mug"))
            mean_s_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0s, ", pfx[k], y_var[traj], "_mueta1s, ", pfx[k], y_var[traj], "_mueta2s, ", pfx[k], y_var[traj], "_knot)"),
                                                    name = paste0(pfx[k], y_var[traj], "_mean_s"))
            GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0s", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 1),
                                        mxPath(from = paste0("eta1s", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L1", records[[traj]], y_var[traj], "[1,1]")),
                                        mxPath(from = paste0("eta2s", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L2", records[[traj]], y_var[traj], "[1,1]")),
                                        mxPath(from = paste0("deltag", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L3", records[[traj]], y_var[traj], "[1,1]")))
            func_L[[length(func_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], y_var[traj], "_mug, ", pfx[k], y_var[traj], "_mug, 0), ",
                                                                       "cbind(0, 1, -1, 0), ", "cbind(0, 1, 1, 0)", ", cbind(0, 0, 0, 1))"),
                                                                name = paste0(pfx[k], "func", y_var[traj]))
            grad_L[[length(grad_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], y_var[traj], "_mug, ", pfx[k], y_var[traj], "_mug, ", pfx[k], y_var[traj], "_mueta2s - ", pfx[k], y_var[traj], "_mueta1s), ",
                                                                       "cbind(0, 1, -1, 0), ", "cbind(0, 1, 1, 0)", ", cbind(0, 0, 0, 1))"),
                                                                name = paste0(pfx[k], "grad", y_var[traj]))
            psi_s_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00s, ", pfx[k], y_var[traj], "_psi01s, ",
                                                          pfx[k], y_var[traj], "_psi02s, ", pfx[k], y_var[traj], "_psi0gs), ",
                                                          "cbind(", pfx[k], y_var[traj], "_psi01s, ", pfx[k], y_var[traj], "_psi11s, ",
                                                          pfx[k], y_var[traj], "_psi12s, ", pfx[k], y_var[traj], "_psi1gs), ",
                                                          "cbind(", pfx[k], y_var[traj], "_psi02s, ", pfx[k], y_var[traj], "_psi12s, ",
                                                          pfx[k], y_var[traj], "_psi22s, ", pfx[k], y_var[traj], "_psi2gs), ",
                                                          "cbind(", pfx[k], y_var[traj], "_psi0gs, ", pfx[k], y_var[traj], "_psi1gs, ",
                                                          pfx[k], y_var[traj], "_psi2gs, ", pfx[k], y_var[traj], "_psiggs))"),
                                                   name = paste0(pfx[k], y_var[traj], "_psi_s"))
            mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], "func", y_var[traj], " %*% ", pfx[k], y_var[traj],
                                                         "_mean_s)"), name = paste0(pfx[k], y_var[traj], "_mean0"))
            psi_L[[traj]] <- mxAlgebraFromString(paste0(pfx[k], "grad", y_var[traj], " %*% ", pfx[k], y_var[traj], "_psi_s",
                                                        " %*% t(", pfx[k], "grad", y_var[traj], ")"),
                                                 name = paste0(pfx[k], y_var[traj], "_psi0"))
            gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00s", "01s", "02s", "0gs", "01s", "11s", "12s", "1gs",
                                                                                 "02s", "12s", "22s", "2gs", "0gs", "1gs", "2gs", "ggs")),
                                           byrow = TRUE, nrow = 4, ncol = 4)
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              psi_btw_s_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi02s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi0gs), ",
                                                                               "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi12s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi1gs), ",
                                                                               "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi20s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi21s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi22s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi2gs), ",
                                                                               "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig0s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig1s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig2s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psiggs))"),
                                                                        name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi_s"))
              psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0(pfx[k], "grad", y_var[traj_i], " %*% ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi_s",
                                                                             " %*% t(", pfx[k], "grad", y_var[traj_j + 1], ")"),
                                                                      name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
              gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                   c("00s", "01s", "02s", "0gs", "10s", "11s", "12s", "1gs",
                                                                     "20s", "21s", "22s", "2gs", "g0s", "g1s", "g2s", "ggs")),
                                                            byrow = TRUE, nrow = 4, ncol = 4)
            }
          }
          for (traj in seq_along(y_var)){
            multi_label[((traj - 1) * 4 + 1):((traj - 1) * 4 + 4), ((traj - 1) * 4 + 1):((traj - 1) * 4 + 4)] <-
              gf_var_label[[traj]]
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              multi_label[((traj_i - 1) * 4 + 1):((traj_i - 1) * 4 + 4), (traj_j * 4 + 1):(traj_j * 4 + 4)] <-
                gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
            }
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = t(starts[[k]]$covMatrix)[row(t(starts[[k]]$covMatrix)) >= col(t(starts[[k]]$covMatrix))],
                                            labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                     GF_loadings[[k]], GAMMA, GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]],
                                     func_L, grad_L, mean_s_L, psi_s_L, psi_btw_s_L,
                                     mean_L, psi_L, psi_btw_L)
        }
      }
      else if (!intrinsic){
        latents <- paste0(rep(c("eta0s", "eta1s", "eta2s"), length(y_var)), rep(y_var, each = 3))
        for (k in 1:nClass){
          multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
          GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
          func_L <- grad_L <- mean_s_L <- psi_s_L <- mean_L <- psi_L <- list()
          gf_var_label <- gf_cov_label <- list()
          psi_btw_s_L <- psi_btw_L <- list()
          for (traj in seq_along(y_var)){
            GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0s", "eta1s", "eta2s"), y_var[traj]),
                                      arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:3],
                                      labels = paste0(pfx[k], y_var[traj], c("_mueta0s", "_mueta1s", "_mueta2s")))
            GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$means[[traj]][4],
                                      labels = paste0(pfx[k], y_var[traj], "_knot"),
                                      name = paste0(pfx[k], y_var[traj], "_mug"))
            mean_s_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0s, ", pfx[k], y_var[traj], "_mueta1s, ", pfx[k], y_var[traj], "_mueta2s)"),
                                                    name = paste0(pfx[k], y_var[traj], "_mean_s"))
            GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0s", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 1),
                                        mxPath(from = paste0("eta1s", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L1", records[[traj]], y_var[traj], "[1,1]")),
                                        mxPath(from = paste0("eta2s", y_var[traj]), to = traj_list[[traj]],
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], "L2", records[[traj]], y_var[traj], "[1,1]")))
            func_L[[length(func_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], y_var[traj], "_mug, ", pfx[k], y_var[traj], "_mug), ",
                                                                       "cbind(0, 1, -1), ", "cbind(0, 1, 1))"),
                                                                name = paste0(pfx[k], "func", y_var[traj]))
            grad_L[[length(grad_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], y_var[traj], "_mug, ", pfx[k], y_var[traj], "_mug), ",
                                                                       "cbind(0, 1, -1), ", "cbind(0, 1, 1))"),
                                                                name = paste0(pfx[k], "grad", y_var[traj]))
            psi_s_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00s, ", pfx[k], y_var[traj], "_psi01s, ",
                                                          pfx[k], y_var[traj], "_psi02s), ",
                                                          "cbind(", pfx[k], y_var[traj], "_psi01s, ", pfx[k], y_var[traj], "_psi11s, ",
                                                          pfx[k], y_var[traj], "_psi12s), ",
                                                          "cbind(", pfx[k], y_var[traj], "_psi02s, ", pfx[k], y_var[traj], "_psi12s, ",
                                                          pfx[k], y_var[traj], "_psi22s))"),
                                                   name = paste0(pfx[k], y_var[traj], "_psi_s"))
            mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], "func", y_var[traj], " %*% ", pfx[k], y_var[traj],
                                                         "_mean_s, ", pfx[k], y_var[traj], "_knot)"), name = paste0(pfx[k], y_var[traj], "_mean0"))
            psi_L[[traj]] <- mxAlgebraFromString(paste0(pfx[k], "grad", y_var[traj], " %*% ", pfx[k], y_var[traj], "_psi_s",
                                                        " %*% t(", pfx[k], "grad", y_var[traj], ")"),
                                                 name = paste0(pfx[k], y_var[traj], "_psi0"))
            gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00s", "01s", "02s", "01s", "11s", "12s",
                                                                                 "02s", "12s", "22s")),
                                           byrow = TRUE, nrow = 3, ncol = 3)
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              psi_btw_s_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi02s), ",
                                                                               "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi12s), ",
                                                                               "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi20s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi21s, ",
                                                                               pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi22s))"),
                                                                        name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi_s"))
              psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0(pfx[k], "grad", y_var[traj_i], " %*% ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi_s",
                                                                             " %*% t(", pfx[k], "grad", y_var[traj_j + 1], ")"),
                                                                      name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
              gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                   c("00s", "01s", "02s", "10s", "11s", "12s",
                                                                     "20s", "21s", "22s")),
                                                            byrow = TRUE, nrow = 3, ncol = 3)
            }
          }
          for (traj in seq_along(y_var)){
            multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
              gf_var_label[[traj]]
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
                gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
            }
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs", free = TRUE,
                                            values = .multi_reduced_cov_starts(starts[[k]]$covMatrix, n_traj = length(y_var)),
                                            labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                     GF_loadings[[k]], GAMMA, GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]],
                                     func_L, grad_L, mean_s_L, psi_s_L, psi_btw_s_L,
                                     mean_L, psi_L, psi_btw_L)
        }
      }
    }
  }
  else if (y_model == "LCSM"){
    # Obtain additional parameters (derived) for the specified functional form
    # Define latent variables in addition to growth factors
    L_list <- D_list <- list()
    for (traj in seq_along(y_var)){
      L_list[[traj]] <- paste0("l", y_var_L[traj], records[[traj]])
      D_list[[traj]] <- paste0("d", y_var_L[traj], records[[traj]][-1])
    }
    PATH_L_L_CL <- PATH_SLP_L_CL <- PATH_AUTO_L_CL <- list()
    for (k in 1:nClass){
      PATH_L_L <- PATH_SLP_L <- PATH_AUTO_L <- list()
      for (traj in seq_along(y_var)){
        ## Define paths from latent true scores to observed scores
        PATH_L_L[[traj]] <- mxPath(from = paste0("l", tolower(y_var[traj]), records[[traj]]),
                                   to = paste0(y_var[[traj]], records[[traj]]), arrows = 1, free = FALSE, values = 1)
        ## Define paths from latent instantaneous rate of change at each measurement to true scores
        PATH_SLP_L[[traj]] <- mxPath(from = paste0("d", tolower(y_var[traj]), records[[traj]][-1]),
                                     to = paste0("l", tolower(y_var[traj]), records[[traj]][-1]),
                                     arrows = 1, free = FALSE, values = 0,
                                     labels = paste0(y_var[traj], "lag", records[[traj]][-1], "[1,1]"))
        #### Define autoregressive paths
        PATH_AUTO_L[[traj]] <- mxPath(from = paste0("l", tolower(y_var[traj]), records[[traj]][-length(records[[traj]])]),
                                      to = paste0("l", tolower(y_var[traj]), records[[traj]][-1]), arrows = 1,
                                      free = FALSE, values = 1)
      }
      PATH_L_L_CL[[k]] <- PATH_L_L
      PATH_SLP_L_CL[[k]] <- PATH_SLP_L
      PATH_AUTO_L_CL[[k]] <- PATH_AUTO_L
      PATH_L_L <- PATH_SLP_L <- PATH_AUTO_L <- list()
    }
    if (curveFun %in% c("nonparametric", "NonP")){
      latents <- c(paste0(rep(c("eta0", "eta1"), length(y_var)), rep(y_var, each = 2)),
                   unlist(L_list), unlist(D_list))
      for (k in 1:nClass){
        multi_label <- matrix(NA, nrow = 2 * length(y_var), ncol = 2 * length(y_var))
        GF_MEAN <- GF_LOADINGS <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in seq_along(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:2],
                                    labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1")))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                             to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]),
                                             to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                             arrows = 1, free = c(FALSE, rep(TRUE, length(records[[traj]]) - 2)),
                                             values = c(1, starts[[k]]$rel_rate[[traj]][-1]),
                                             labels = paste0(pfx[k], y_var[traj], "_rel_rate",
                                                             1:(length(records[[traj]]) - 1))))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj], "_mueta1)"),
                                                name = paste0(pfx[k], y_var[traj], "_mean0"))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj], "_psi01), ",
                                                      "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi11))"),
                                               name = paste0(pfx[k], y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "01", "11")),
                                         byrow = TRUE, nrow = 2, ncol = 2)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01), ",
                                                                           "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11))"),
                                                                    name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi", c("00", "01", "10", "11")),
                                                          byrow = TRUE, nrow = 2, ncol = 2)
          }
        }
        for (traj in seq_along(y_var)){
          multi_label[((traj - 1) * 2 + 1):((traj - 1) * 2 + 2), ((traj - 1) * 2 + 1):((traj - 1) * 2 + 2)] <- gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 2 + 1):((traj_i - 1) * 2 + 2), (traj_j * 2 + 1):(traj_j * 2 + 2)] <- gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
          }
        }
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   mxPath(from = latents[1:(length(y_var) * 2)], to = latents[1:(length(y_var) * 2)], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = t(starts[[k]]$covMatrix)[row(t(starts[[k]]$covMatrix)) >= col(t(starts[[k]]$covMatrix))],
                                          labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                   GF_loadings[[k]], AddPara[[k]], GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]], mean_L, psi_L,
                                   psi_btw_L, PATH_L_L_CL[[k]], PATH_SLP_L_CL[[k]], PATH_AUTO_L_CL[[k]])
      }
    }
    else if (curveFun %in% c("quadratic", "QUAD")){
      latents <- c(paste0(rep(c("eta0", "eta1", "eta2"), length(y_var)), rep(y_var, each = 3)),
                   unlist(L_list), unlist(D_list))
      for (k in 1:nClass){
        multi_label <- matrix(NA, nrow = 3 * length(y_var), ncol = 3 * length(y_var))
        GF_MEAN <- GF_LOADINGS <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in seq_along(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:3],
                                    labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                             to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]),
                                             to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta2", y_var[traj]),
                                             to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0(pfx[k], y_var[traj], "L2", records[[traj]][-1], "[1,1]")))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj],
                                                       "_mueta1, ", pfx[k], y_var[traj], "_mueta2)"),
                                                name = paste0(pfx[k], y_var[traj], "_mean0"))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi02), ",
                                                      "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi11, ", pfx[k], y_var[traj], "_psi12), ",
                                                      "cbind(", pfx[k], y_var[traj], "_psi02, ", pfx[k], y_var[traj], "_psi12, ", pfx[k], y_var[traj], "_psi22))"),
                                               name = paste0(pfx[k], y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "02", "01", "11", "12", "02", "12", "22")),
                                         byrow = TRUE, nrow = 3, ncol = 3)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)) {
            psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi02), ",
                                                                           "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi12), ",
                                                                           "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi20, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi21, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi22))"),
                                                                    name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi", c("00", "01", "02", "10", "11", "12", "20", "21", "22")),
                                                          byrow = TRUE, nrow = 3, ncol = 3)
          }
        }
        for (traj in seq_along(y_var)){
          multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <- gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <- gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
          }
        }
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   mxPath(from = latents[1:(length(y_var) * 3)], to = latents[1:(length(y_var) * 3)], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = t(starts[[k]]$covMatrix)[row(t(starts[[k]]$covMatrix)) >= col(t(starts[[k]]$covMatrix))],
                                          labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                   GF_loadings[[k]], AddPara[[k]], GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]], mean_L, psi_L,
                                   psi_btw_L, PATH_L_L_CL[[k]], PATH_SLP_L_CL[[k]], PATH_AUTO_L_CL[[k]])
      }
    }
    else if (curveFun %in% c("negative exponential", "EXP")){
      if (intrinsic){
        latents <- c(paste0(rep(c("eta0", "eta1", "deltag"), length(y_var)), rep(y_var, each = 3)),
                     unlist(L_list), unlist(D_list))
        for (k in 1:nClass){
          multi_label <- matrix(NA, nrow = 3 * length(y_var), ncol = 3 * length(y_var))
          GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
          mean_L <- psi_L <- list()
          gf_var_label <- gf_cov_label <- list()
          psi_btw_L <- list()
          for (traj in seq_along(y_var)){
            GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                      arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:2],
                                      labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1")))
            GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$means[[traj]][3],
                                      labels = paste0(pfx[k], y_var[traj], "_slp_ratio"),
                                      name = paste0(pfx[k], y_var[traj], "_mug"))
            mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj], "_mueta1, ",
                                                         pfx[k], y_var[traj], "_slp_ratio)"),
                                                  name = paste0(pfx[k], y_var[traj], "_mean0"))
            GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                               to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                               arrows = 1, free = FALSE, values = 1),
                                        mxPath(from = paste0("eta1", y_var[traj]),
                                               to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], y_var[traj], "L1", records[[traj]][-1], "[1,1]")),
                                        mxPath(from = paste0("deltag", y_var[traj]),
                                               to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], y_var[traj], "L2", records[[traj]][-1], "[1,1]")))
            psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj],
                                                        "_psi01, ", pfx[k], y_var[traj], "_psi0g), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj],
                                                        "_psi11, ", pfx[k], y_var[traj], "_psi1g), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi0g, ", pfx[k], y_var[traj],
                                                        "_psi1g, ", pfx[k], y_var[traj], "_psigg))"),
                                                 name = paste0(pfx[k], y_var[traj], "_psi0"))
            gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "0g", "01", "11", "1g", "0g", "1g", "gg")),
                                           byrow = TRUE, nrow = 3, ncol = 3)
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1],
                                                                             "_psi01, ", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi0g),",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi1g), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig0, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig1, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psigg))"),
                                                                      name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
              gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                   c("00", "01", "0g", "10", "11", "1g", "g0", "g1", "gg")),
                                                            byrow = TRUE, nrow = 3, ncol = 3)
            }
          }
          for (traj in seq_along(y_var)){
            multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
              gf_var_label[[traj]]
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
                gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
            }
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = latents[1:(length(y_var) * 3)], to = latents[1:(length(y_var) * 3)], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = t(starts[[k]]$covMatrix)[row(t(starts[[k]]$covMatrix)) >= col(t(starts[[k]]$covMatrix))],
                                            labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                     GF_loadings[[k]], AddPara[[k]], GAMMA, GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]],
                                     mean_L, psi_L, psi_btw_L, PATH_L_L_CL[[k]], PATH_SLP_L_CL[[k]], PATH_AUTO_L_CL[[k]])
        }
      }
      else if (!intrinsic){
        latents <- c(paste0(rep(c("eta0", "eta1"), length(y_var)), rep(y_var, each = 2)),
                     unlist(L_list), unlist(D_list))
        for (k in 1:nClass){
          multi_label <- matrix(NA, nrow = 2 * length(y_var), ncol = 2 * length(y_var))
          GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
          mean_L <- psi_L <- list()
          gf_var_label <- gf_cov_label <- list()
          psi_btw_L <- list()
          for (traj in seq_along(y_var)){
            GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                      arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:2],
                                      labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1")))
            GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$means[[traj]][3],
                                      labels = paste0(pfx[k], y_var[traj], "_slp_ratio"),
                                      name = paste0(pfx[k], y_var[traj], "_mug"))
            mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj], "_mueta1, ",
                                                         pfx[k], y_var[traj], "_slp_ratio)"),
                                                  name = paste0(pfx[k], y_var[traj], "_mean0"))
            GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                               to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                               arrows = 1, free = FALSE, values = 1),
                                        mxPath(from = paste0("eta1", y_var[traj]),
                                               to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], y_var[traj], "L1", records[[traj]][-1], "[1,1]")))
            psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj], "_psi01), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi11))"),
                                                 name = paste0(pfx[k], y_var[traj], "_psi0"))
            gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "01", "11")),
                                           byrow = TRUE, nrow = 2, ncol = 2)
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01),",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11))"),
                                                                      name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
              gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                   c("00", "01", "10", "11")),
                                                            byrow = TRUE, nrow = 2, ncol = 2)
            }
          }
          for (traj in seq_along(y_var)){
            multi_label[((traj - 1) * 2 + 1):((traj - 1) * 2 + 2), ((traj - 1) * 2 + 1):((traj - 1) * 2 + 2)] <-
              gf_var_label[[traj]]
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              multi_label[((traj_i - 1) * 2 + 1):((traj_i - 1) * 2 + 2), (traj_j * 2 + 1):(traj_j * 2 + 2)] <-
                gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
            }
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = latents[1:(length(y_var) * 2)], to = latents[1:(length(y_var) * 2)], arrows = 2, connect = "unique.pairs", free = TRUE,
                                            values = t(starts[[k]]$covMatrix[c(1:2, 4:5), c(1:2, 4:5)])[row(t(starts[[k]]$covMatrix[c(1:2, 4:5), c(1:2, 4:5)])) >=
                                                                                                     col(t(starts[[k]]$covMatrix[c(1:2, 4:5), c(1:2, 4:5)]))],
                                            labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                     GF_loadings[[k]], AddPara[[k]], GAMMA, GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]], mean_L, psi_L,
                                     psi_btw_L, PATH_L_L_CL[[k]], PATH_SLP_L_CL[[k]], PATH_AUTO_L_CL[[k]])
        }
      }
    }
    else if (curveFun %in% c("Jenss-Bayley", "JB")){
      if (intrinsic){
        latents <- c(paste0(rep(c("eta0", "eta1", "eta2", "deltag"), length(y_var)), rep(y_var, each = 4)),
                     unlist(L_list), unlist(D_list))
        for (k in 1:nClass){
          multi_label <- matrix(NA, nrow = 4 * length(y_var), ncol = 4 * length(y_var))
          GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
          mean_L <- psi_L <- list()
          gf_var_label <- gf_cov_label <- list()
          psi_btw_L <- list()
          for (traj in seq_along(y_var)){
            GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                      arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:3],
                                      labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
            GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$means[[traj]][4],
                                      labels = paste0(pfx[k], y_var[traj], "_acc_ratio"),
                                      name = paste0(pfx[k], y_var[traj], "_mug"))
            mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj], "_mueta1, ",
                                                         pfx[k], y_var[traj], "_mueta2, ", pfx[k], y_var[traj], "_acc_ratio)"),
                                                  name = paste0(pfx[k], y_var[traj], "_mean0"))
            GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                               to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                               arrows = 1, free = FALSE, values = 1),
                                        mxPath(from = paste0("eta1", y_var[traj]),
                                               to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                               arrows = 1, free = FALSE, values = 1),
                                        mxPath(from = paste0("eta2", y_var[traj]),
                                               to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                               arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], y_var[traj], "L2", records[[traj]][-1], "[1,1]")),
                                        mxPath(from = paste0("deltag", y_var[traj]),
                                               to = paste0("d", y_var_L[traj], records[[traj]][-1]), arrows = 1, free = FALSE, values = 0,
                                               labels = paste0(pfx[k], y_var[traj], "L3", records[[traj]][-1], "[1,1]")))
            psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj], "_psi01, ",
                                                        pfx[k], y_var[traj], "_psi02, ", pfx[k], y_var[traj], "_psi0g), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi11, ",
                                                        pfx[k], y_var[traj], "_psi12, ", pfx[k], y_var[traj], "_psi1g), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi02, ", pfx[k], y_var[traj], "_psi12, ",
                                                        pfx[k], y_var[traj], "_psi22, ", pfx[k], y_var[traj], "_psi2g), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi0g, ", pfx[k], y_var[traj], "_psi1g, ",
                                                        pfx[k], y_var[traj], "_psi2g, ", pfx[k], y_var[traj], "_psigg))"),
                                                 name = paste0(pfx[k], y_var[traj], "_psi0"))
            gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "02", "0g", "01", "11", "12", "1g",
                                                                                 "02", "12", "22", "2g", "0g", "1g", "2g", "gg")),
                                           byrow = TRUE, nrow = 4, ncol = 4)
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi02, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi0g), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi12, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi1g), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi20, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi21, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi22, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi2g), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig0, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig1, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psig2, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psigg))"),
                                                                      name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
              gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                   c("00", "01", "02", "0g", "10", "11", "12", "1g",
                                                                     "20", "21", "22", "2g", "g0", "g1", "g2", "gg")),
                                                            byrow = TRUE, nrow = 4, ncol = 4)
            }
          }
          for (traj in seq_along(y_var)){
            multi_label[((traj - 1) * 4 + 1):((traj - 1) * 4 + 4), ((traj - 1) * 4 + 1):((traj - 1) * 4 + 4)] <-
              gf_var_label[[traj]]
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              multi_label[((traj_i - 1) * 4 + 1):((traj_i - 1) * 4 + 4), (traj_j * 4 + 1):(traj_j * 4 + 4)] <-
                gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
            }
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = latents[1:(length(y_var) * 4)], to = latents[1:(length(y_var) * 4)], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = t(starts[[k]]$covMatrix)[row(t(starts[[k]]$covMatrix)) >= col(t(starts[[k]]$covMatrix))],
                                            labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                     GF_loadings[[k]], AddPara[[k]], GAMMA, GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]],
                                     mean_L, psi_L, psi_btw_L, PATH_L_L_CL[[k]], PATH_SLP_L_CL[[k]], PATH_AUTO_L_CL[[k]])
        }
      }
      else if (!intrinsic){
        latents <- c(paste0(rep(c("eta0", "eta1", "eta2"), length(y_var)), rep(y_var, each = 3)),
                     unlist(L_list), unlist(D_list))
        for (k in 1:nClass){
          multi_label <- matrix(NA, nrow = 3 * length(y_var), ncol = 3 * length(y_var))
          GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
          mean_L <- psi_L <- list()
          gf_var_label <- gf_cov_label <- list()
          psi_btw_L <- list()
          for (traj in seq_along(y_var)){
            GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                      arrows = 1, free = TRUE, values = starts[[k]]$means[[traj]][1:3],
                                      labels = paste0(pfx[k], y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
            GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$means[[traj]][4],
                                      labels = paste0(pfx[k], y_var[traj], "_acc_ratio"),
                                      name = paste0(pfx[k], y_var[traj], "_mug"))
            mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", pfx[k], y_var[traj], "_mueta0, ", pfx[k], y_var[traj], "_mueta1, ",
                                                         pfx[k], y_var[traj], "_mueta2, ", pfx[k], y_var[traj], "_acc_ratio)"),
                                                  name = paste0(pfx[k], y_var[traj], "_mean0"))
            GF_LOADINGS[[length(GF_LOADINGS) + 1]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                                                  to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                                                  arrows = 1, free = FALSE, values = 1),
                                                           mxPath(from = paste0("eta1", y_var[traj]),
                                                                  to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                                                  arrows = 1, free = FALSE, values = 1),
                                                           mxPath(from = paste0("eta2", y_var[traj]),
                                                                  to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                                                  arrows = 1, free = FALSE, values = 0,
                                                                  labels = paste0(pfx[k], y_var[traj], "L2", records[[traj]][-1], "[1,1]")))
            psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj], "_psi00, ", pfx[k], y_var[traj], "_psi01, ",
                                                        pfx[k], y_var[traj], "_psi02), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi01, ", pfx[k], y_var[traj], "_psi11, ",
                                                        pfx[k], y_var[traj], "_psi12), ",
                                                        "cbind(", pfx[k], y_var[traj], "_psi02, ", pfx[k], y_var[traj], "_psi12, ",
                                                        pfx[k], y_var[traj], "_psi22))"),
                                                 name = paste0(pfx[k], y_var[traj], "_psi0"))
            gf_var_label[[traj]] <- matrix(paste0(pfx[k], y_var[traj], "_psi", c("00", "01", "02", "01", "11", "12",
                                                                                 "02", "12", "22")),
                                           byrow = TRUE, nrow = 3, ncol = 3)
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              psi_btw_L[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi01, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi02), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi12), ",
                                                                             "cbind(", pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi20, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi21, ",
                                                                             pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi22))"),
                                                                      name = paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi"))
              gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]] <- matrix(paste0(pfx[k], y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                   c("00", "01", "02", "10", "11", "12",
                                                                     "20", "21", "22")),
                                                            byrow = TRUE, nrow = 3, ncol = 3)
            }
          }
          for (traj in seq_along(y_var)){
            multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
              gf_var_label[[traj]]
          }
          for (traj_i in 1:(length(y_var) - 1)){
            for (traj_j in traj_i:(length(y_var) - 1)){
              multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
                gf_cov_label[[.outcome_pair_index(traj_i, traj_j + 1L, length(y_var))]]
            }
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = latents[1:(length(y_var) * 3)], to = latents[1:(length(y_var) * 3)], arrows = 2, connect = "unique.pairs", free = TRUE,
                                            values = .multi_reduced_cov_starts(starts[[k]]$covMatrix, n_traj = length(y_var)),
                                            labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                     GF_loadings[[k]], AddPara[[k]], GAMMA, GF_MEAN, GF_LOADINGS, RES_L[[k]], COV_L[[k]],
                                     mean_L, psi_L, psi_btw_L, PATH_L_L_CL[[k]], PATH_SLP_L_CL[[k]], PATH_AUTO_L_CL[[k]])
        }
      }
    }
  }
  for (k in seq_len(nClass)) class.list[[k]] <- do.call(mxModel, c(list(class.list[[k]]), extra[[k]]))
  class.list
}

#' @keywords internal
#' @noRd
#' @importFrom OpenMx mxData
getMGM.mxModel <- function(dat, t_var, y_var, curveFun, intrinsic, records, y_model, starts){
  loadings <- getMULTI.loadings(y_model = y_model, t_var = t_var, y_var = y_var, curveFun = curveFun,
                                intrinsic = intrinsic, records = records)
  addpara <- if (y_model == "LCSM"){
    getMULTI.addpara(dat = dat, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic, t_var = t_var,
                     records = records, starts = starts)
  }
  .build_MGM(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic, records = records,
             y_model = y_model, starts = list(starts), pfx = "",
             model_names = .single_model_name_mgm(y_model, curveFun, intrinsic),
             extra = list(list(mxData(observed = dat, type = "raw"))), GF_loadings = list(loadings),
             AddPara = list(addpara))[[1]]
}

#' @keywords internal
#' @noRd
getsub.MGM <- function(dat, nClass, t_var, y_var, curveFun, intrinsic, records, y_model, starts, container,
                       grp_var = NULL){
  loadings <- getMIX_MULTI.loadings(nClass = nClass, y_model = y_model, t_var = t_var, y_var = y_var,
                                    curveFun = curveFun, intrinsic = intrinsic, records = records)
  addpara <- if (y_model == "LCSM"){
    getMIX_MULTI.addpara(dat = dat, nClass = nClass, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                         t_var = t_var, records = records, starts = starts,
                         grp_var = if (container == "multigroup") grp_var else NULL)
  }
  .build_MGM(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic, records = records,
             y_model = y_model, starts = starts, pfx = paste0("c", seq_len(nClass)),
             model_names = paste0("Class", seq_len(nClass)), extra = .class_extras(dat, nClass, container, grp_var),
             GF_loadings = loadings, AddPara = addpara)
}

#' @title One builder for the TVC sub-model, single-group or class-specific
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
.build_TVC <- function(dat, t_var, y_var, curveFun, intrinsic, records, y_model, TVC, decompose, growth_TIC, starts, pfx, model_names, extra, GF_loadings, AddPara, TVC_info){
  nClass <- length(pfx)
  ## Define manifest variables
  manifests <- c(paste0(TVC, records), paste0(y_var, records))
  if (!is.null(growth_TIC)){
    manifests <- c(manifests, growth_TIC)
  }
  # Define latent variables for the TVC for each type of decomposition
  if (decompose != 0){
    X_latents <- c("eta0x", "eta1x", paste0("lx", records), paste0("dx", records[-1]))
    if (decompose ==  2){
      X_latents <- c(X_latents, paste0("deltax", records[-1]))
    }
    else if (decompose ==  3){
      X_latents <- c(X_latents, paste0("Deltax", records[-1]))
    }
  }
  # Obtain factor loadings for the specified functional form
  BETA <- class.list <- list()
  ## Define latent variables, growth factor loadings, paths of the longitudinal outcome
  if (y_model == "LGCM"){
    if (curveFun %in% c("linear", "LIN")){
      latents <- c("eta0Y", "eta1Y")
      Y_nGF <- length(latents)
      if (decompose == 0){
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                  labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 0,
                                              labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                           name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                       mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                values = starts[[k]]$beta[1:2, seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta1TIC", seq_along(growth_TIC))),
                                                byrow = TRUE, name = paste0(pfx[k], "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                byrow = FALSE, name = paste0(pfx[k], "mux")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA)
          }
          else if (is.null(growth_TIC)){
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 0,
                                              labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                           name = paste0(pfx[k], "Y_psi0")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                       TVC_info[[k]], GF_loadings[[k]])
          }
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                  labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 0,
                                              labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                           name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                           name = paste0(pfx[k], "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                  "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                           name = paste0(pfx[k], "X_psi0")),
                                       mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                values = starts[[k]]$beta[1:2, seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta1TIC", seq_along(growth_TIC))),
                                                byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                       mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                  "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC))"),
                                                           name = paste0(pfx[k], "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                byrow = FALSE, name = paste0(pfx[k], "mux")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                           name = paste0(pfx[k], "BL_mean")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                  labels = paste0(pfx[k], "beta", p - 1, "TVC"))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 0,
                                              labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                           name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                           name = paste0(pfx[k], "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                  "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                           name = paste0(pfx[k], "X_psi0")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC)"),
                                                           name = paste0(pfx[k], "beta")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA)
          }
        }
      }
    }
    else if (curveFun %in% c("quadratic", "QUAD")){
      latents <- c("eta0Y", "eta1Y", "eta2Y")
      Y_nGF <- length(latents)
      if (decompose == 0){
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                  labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                       mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                              arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                              labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"), name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                  "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                       mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta2TIC", seq_along(growth_TIC))),
                                                byrow = TRUE, name = paste0(pfx[k], "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                byrow = FALSE, name = paste0(pfx[k], "mux")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA)
          }
          else if (is.null(growth_TIC)){
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                       mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                              arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                              labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"), name = paste0(pfx[k], "Y_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                  "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                           name = paste0(pfx[k], "Y_psi0")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                       TVC_info[[k]], GF_loadings[[k]])
          }
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                  labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                       mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                              arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                              labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"), name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                  "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                           name = paste0(pfx[k], "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                  "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                           name = paste0(pfx[k], "X_psi0")),
                                       mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta2TIC", seq_along(growth_TIC))),
                                                byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                       mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                  "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC))"),
                                                           name = paste0(pfx[k], "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                byrow = FALSE, name = paste0(pfx[k], "mux")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                           name = paste0(pfx[k], "BL_mean")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                  labels = paste0(pfx[k], "beta", p - 1, "TVC"))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                       mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                              arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                              labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"), name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                  "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                           name = paste0(pfx[k], "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                  "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                           name = paste0(pfx[k], "X_psi0")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC)"), name = paste0(pfx[k], "beta")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA)
          }
        }
      }
    }
    else if (curveFun %in% c("negative exponential", "EXP")){
      if (intrinsic){
        latents <- c("eta0Y", "eta1Y", "deltag")
        Y_nGF <- length(latents)
        if (decompose == 0){
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
              }
              BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", c(paste0("TIC", seq_along(growth_TIC)))))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"),
                                                  name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "betagTIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"),
                                                  name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]])
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (decompose == 1){
              Y_nGF <- length(latents) - (length(records) * 2 - 1) - 2
            }
            else if (decompose == 2 | decompose == 3) {
              Y_nGF <- length(latents) - (length(records) * 2 - 1) - (length(records) - 1) - 2
            }
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              }
              BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"),
                                                  name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "betagTIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                    "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "betagTVC))"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                             name = paste0(pfx[k], "BL_mean")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                    labels = paste0(pfx[k], "beta", p - 1, "TVC"))
              }
              BETA[[Y_nGF]] <- mxPath(from = paste0("lx", records[1]), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF],
                                      labels = paste0(pfx[k], "beta", "g", "TVC"))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"),
                                                  name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "betagTVC)"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
          }
        }
      }
      else if (!intrinsic){
        latents <- c("eta0Y", "eta1Y")
        Y_nGF <- length(latents)
        if (decompose == 0){
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:2, 4)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:2, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:2, ] + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:2, 4)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                             name = paste0(pfx[k], "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]])
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:2, 4)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:2, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                    "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC))"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                             name = paste0(pfx[k], "BL_mean")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:2, ] + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                    labels = paste0(pfx[k], "beta", p - 1, "TVC"))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:2, 4)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC)"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:2, ] + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
          }
        }
      }
    }
    else if (curveFun %in% c("Jenss-Bayley", "JB")){
      if (intrinsic){
        latents <- c("eta0Y", "eta1Y", "eta2Y", "deltag")
        Y_nGF <- length(latents)
        if (decompose == 0){
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
              }
              BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", c(paste0("TIC", seq_along(growth_TIC)))))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22, ", pfx[k], "Y_psi2g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psi2g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:4, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "betagTIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22, ", pfx[k], "Y_psi2g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psi2g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]])
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              }
              BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22, ", pfx[k], "Y_psi2g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psi2g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:4, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "betagTIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                    "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC, ", pfx[k], "betagTVC))"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                             name = paste0(pfx[k], "BL_mean")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                    labels = paste0(pfx[k], "beta", p - 1, "TVC"))
              }
              BETA[[Y_nGF]] <- mxPath(from = paste0("lx", records[1]), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF],
                                      labels = paste0(pfx[k], "beta", "g", "TVC"))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22, ", pfx[k], "Y_psi2g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psi2g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC, ", pfx[k], "betagTVC)"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
          }
        }
      }
      else if (!intrinsic){
        latents <- c("eta0Y", "eta1Y", "eta2Y")
        Y_nGF <- length(latents)
        if (decompose == 0){
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:3, ] + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                             name = paste0(pfx[k], "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]])
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                    "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC))"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                             name = paste0(pfx[k], "BL_mean")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:3, ] + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                    labels = paste0(pfx[k], "beta", p - 1, "TVC"))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC)"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:3, ] + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
          }
        }
      }
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      if (intrinsic){
        latents <- c("eta0sY", "eta1sY", "eta2sY", "deltag")
        Y_nGF <- length(latents)
        if (decompose == 0){
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
              }
              BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", c(paste0("TIC", seq_along(growth_TIC)))))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                          "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, ", pfx[k], "Y_mueta2s - ", pfx[k], "Y_mueta1s),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "grad")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                             name = paste0(pfx[k], "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi0gs), ",
                                                                    "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi1gs), ",
                                                                    "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s, ", pfx[k], "Y_psi2gs), ",
                                                                    "cbind(", pfx[k], "Y_psi0gs, ", pfx[k], "Y_psi1gs, ", pfx[k], "Y_psi2gs, ", pfx[k], "Y_psiggs))"),
                                                             name = paste0(pfx[k], "Y_psi_s")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:4, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "betagTIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_s")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "func[1:3, 1:3] %*% ", pfx[k], "Y_alpha_s, ", pfx[k], "Y_knot)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "beta_s"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                          "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, ", pfx[k], "Y_mueta2s - ", pfx[k], "Y_mueta1s),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "grad")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                             name = paste0(pfx[k], "Y_mean_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi0gs), ",
                                                                    "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi1gs), ",
                                                                    "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s, ", pfx[k], "Y_psi2gs), ",
                                                                    "cbind(", pfx[k], "Y_psi0gs, ", pfx[k], "Y_psi1gs, ", pfx[k], "Y_psi2gs, ", pfx[k], "Y_psiggs))"),
                                                             name = paste0(pfx[k], "Y_psi_s")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "func[1:3, 1:3] %*% ", pfx[k], "Y_mean_s, ", pfx[k], "Y_knot)"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                             name = paste0(pfx[k], "Y_psi0")),
                                         TVC_info[[k]], GF_loadings[[k]])
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              }
              BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                          "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, ", pfx[k], "Y_mueta2s - ", pfx[k], "Y_mueta1s),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "grad")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                             name = paste0(pfx[k], "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi0gs), ",
                                                                    "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi1gs), ",
                                                                    "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s, ", pfx[k], "Y_psi2gs), ",
                                                                    "cbind(", pfx[k], "Y_psi0gs, ", pfx[k], "Y_psi1gs, ", pfx[k], "Y_psi2gs, ", pfx[k], "Y_psiggs))"),
                                                             name = paste0(pfx[k], "Y_psi_s")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:4, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "betagTIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                    "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC, ", pfx[k], "betagTVC))"),
                                                             name = paste0(pfx[k], "beta_s")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                             name = paste0(pfx[k], "BL_mean")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "func[1:3, 1:3] %*% ", pfx[k], "Y_alpha_s, ", pfx[k], "Y_knot)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "beta_s"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, "TVC"))
              }
              BETA[[Y_nGF]] <- mxPath(from = paste0("lx", records[1]), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", "TVC"))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                          "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, ", pfx[k], "Y_mueta2s - ", pfx[k], "Y_mueta1s),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "grad")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                             name = paste0(pfx[k], "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi0gs), ",
                                                                    "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi1gs), ",
                                                                    "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s, ", pfx[k], "Y_psi2gs), ",
                                                                    "cbind(", pfx[k], "Y_psi0gs, ", pfx[k], "Y_psi1gs, ", pfx[k], "Y_psi2gs, ", pfx[k], "Y_psiggs))"),
                                                             name = paste0(pfx[k], "Y_psi_s")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC, ", pfx[k], "betagTVC)"),
                                                             name = paste0(pfx[k], "beta_s")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "func[1:3, 1:3] %*% ", pfx[k], "Y_alpha_s, ", pfx[k], "Y_knot)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "beta_s"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
          }
        }
      }
      else if (!intrinsic){
        latents <- c("eta0sY", "eta1sY", "eta2sY")
        Y_nGF <- length(latents)
        if (decompose == 0){
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s",
                                                                          "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0(pfx[k], "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0(pfx[k], "grad")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                             name = paste0(pfx[k], "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s), ",
                                                                    "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s), ",
                                                                    "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s))"),
                                                             name = paste0(pfx[k], "Y_psi_s")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_s")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "func %*% ", pfx[k], "Y_alpha_s, ", pfx[k], "Y_knot)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "beta_s"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:3, 1] + ", pfx[k], "beta %*% ", pfx[k], "mux"), name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s",
                                                                          "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0(pfx[k], "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0(pfx[k], "grad")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                             name = paste0(pfx[k], "Y_mean_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s), ",
                                                                    "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s), ",
                                                                    "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s))"),
                                                             name = paste0(pfx[k], "Y_psi_s")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "func %*% ", pfx[k], "Y_mean_s, ", pfx[k], "Y_knot)"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                             name = paste0(pfx[k], "Y_psi0")),
                                         TVC_info[[k]], GF_loadings[[k]])
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s",
                                                                          "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0(pfx[k], "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0(pfx[k], "grad")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                             name = paste0(pfx[k], "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s), ",
                                                                    "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s), ",
                                                                    "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s))"),
                                                             name = paste0(pfx[k], "Y_psi_s")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                    "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC))"),
                                                             name = paste0(pfx[k], "beta_s")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                             name = paste0(pfx[k], "BL_mean")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "func %*% ", pfx[k], "Y_alpha_s, ", pfx[k], "Y_knot)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "beta_s"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:3, 1] + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"), name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
            else if (is.null(growth_TIC)){
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                    labels = paste0(pfx[k], "beta", p - 1, "TVC"))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s",
                                                                          "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0(pfx[k], "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0(pfx[k], "grad")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                             name = paste0(pfx[k], "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s), ",
                                                                    "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s), ",
                                                                    "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s))"),
                                                             name = paste0(pfx[k], "Y_psi_s")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC)"),
                                                             name = paste0(pfx[k], "beta_s")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "func %*% ", pfx[k], "Y_alpha_s, ", pfx[k], "Y_knot)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "beta_s"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 [1:3, 1]+ ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"), name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA)
            }
          }
        }
      }
    }
  }
  else if (y_model == "LCSM"){
    Y_PATH_L_L <- Y_PATH_SLP_L <- Y_PATH_AUTO_L <- list()
    for (k in 1:nClass){
      ## Define paths from latent true scores to observed scores
      Y_PATH_L <- mxPath(from = paste0("ly", records), to = paste0(y_var, records), arrows = 1, free = FALSE,
                         values = 1)
      ## Define paths from latent instantaneous rate of change at each measurement to true scores
      Y_PATH_SLP <- mxPath(from = paste0("dy", records[-1]), to = paste0("ly", records[-1]), arrows = 1,
                           free = FALSE, values = 0, labels = paste0("lag", records[-1], "[1,1]"))
      #### Define autoregressive paths
      Y_PATH_AUTO <- mxPath(from = paste0("ly", records[-length(records)]), to = paste0("ly", records[-1]),
                            arrows = 1, free = FALSE, values = 1)
      Y_PATH_L_L[[k]] <- Y_PATH_L
      Y_PATH_SLP_L[[k]] <- Y_PATH_SLP
      Y_PATH_AUTO_L[[k]] <- Y_PATH_AUTO
    }
    # Obtain additional parameters (derived) for the specified functional form
    if (curveFun %in% c("nonparametric", "NonP")){
      latents <- c("eta0Y", "eta1Y", paste0("dy", records[-1]), paste0("ly", records))
      Y_nGF <- length(latents) - (length(records) * 2 - 1)
      if (decompose == 0){
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                  labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                              free = c(FALSE, rep(TRUE, length(records) - 2)), values = c(1, starts[[k]]$Y_starts$rel_rate[-1]),
                                              labels = paste0(pfx[k], "Y_rel_rate", 1:(length(records) - 1))),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                              labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                           name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                       mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                values = starts[[k]]$beta[1:2, seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta1TIC", seq_along(growth_TIC))),
                                                byrow = TRUE, name = paste0(pfx[k], "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                byrow = FALSE, name = paste0(pfx[k], "mux")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
          }
          else if (is.null(growth_TIC)){
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                              free = c(FALSE, rep(TRUE, length(records) - 2)), values = c(1, starts[[k]]$Y_starts$rel_rate[-1]),
                                              labels = paste0(pfx[k], "Y_rel_rate", 1:(length(records) - 1))),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                              labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                           name = paste0(pfx[k], "Y_psi0")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
          }
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                  labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                              free = c(FALSE, rep(TRUE, length(records) - 2)), values = c(1, starts[[k]]$Y_starts$rel_rate[-1]),
                                              labels = paste0(pfx[k], "Y_rel_rate", 1:(length(records) - 1))),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                              labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                           name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                           name = paste0(pfx[k], "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                  "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                           name = paste0(pfx[k], "X_psi0")),
                                       mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                values = starts[[k]]$beta[1:2, seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta1TIC", seq_along(growth_TIC))),
                                                byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                       mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                  "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC))"),
                                                           name = paste0(pfx[k], "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                byrow = FALSE, name = paste0(pfx[k], "mux")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                           name = paste0(pfx[k], "BL_mean")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                  labels = paste0(pfx[k], "beta", p - 1, "TVC"))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                              free = c(FALSE, rep(TRUE, length(records) - 2)), values = c(1, starts[[k]]$Y_starts$rel_rate[-1]),
                                              labels = paste0(pfx[k], "Y_rel_rate", 1:(length(records) - 1))),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                              labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                           name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                           name = paste0(pfx[k], "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                  "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                           name = paste0(pfx[k], "X_psi0")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC)"),
                                                           name = paste0(pfx[k], "beta")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
          }
        }
      }
    }
    else if (curveFun %in% c("quadratic", "QUAD")){
      latents <- c("eta0Y", "eta1Y", "eta2Y", paste0("dy", records[-1]), paste0("ly", records))
      Y_nGF <- length(latents) - (length(records) * 2 - 1)
      if (decompose == 0){
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                  labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                              labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"), name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                  "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                       mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta2TIC", seq_along(growth_TIC))),
                                                byrow = TRUE, name = paste0(pfx[k], "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                byrow = FALSE, name = paste0(pfx[k], "mux")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
          }
          else if (is.null(growth_TIC)){
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                              labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"), name = paste0(pfx[k], "Y_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                  "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                           name = paste0(pfx[k], "Y_psi0")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
          }
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                  labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                              labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"), name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                  "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                           name = paste0(pfx[k], "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                  "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                           name = paste0(pfx[k], "X_psi0")),
                                       mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                           paste0(pfx[k], "beta2TIC", seq_along(growth_TIC))),
                                                byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                       mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                  "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC))"),
                                                           name = paste0(pfx[k], "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                byrow = FALSE, name = paste0(pfx[k], "mux")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                           name = paste0(pfx[k], "BL_mean")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                  labels = paste0(pfx[k], "beta", p - 1, "TVC"))
            }
            class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means,
                                              labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                              labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                              labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"), name = paste0(pfx[k], "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                  "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                  "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                           name = paste0(pfx[k], "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                           name = paste0(pfx[k], "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                  "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                           name = paste0(pfx[k], "X_psi0")),
                                       mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC)"),
                                                           name = paste0(pfx[k], "beta")),
                                       mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                           name = paste0(pfx[k], "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
          }
        }
      }
    }
    else if (curveFun %in% c("negative exponential", "EXP")){
      if (intrinsic){
        latents <- c("eta0Y", "eta1Y", "deltag", paste0("dy", records[-1]), paste0("ly", records))
        Y_nGF <- length(latents) - (length(records) * 2 - 1)
        if (decompose == 0){
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
              }
              BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", c(paste0("TIC", seq_along(growth_TIC)))))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"),
                                                  name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "betagTIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"),
                                                  name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              }
              BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"),
                                                  name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "betagTIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                    "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "betagTVC))"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                             name = paste0(pfx[k], "BL_mean")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
            else if (is.null(growth_TIC)){
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                    labels = paste0(pfx[k], "beta", p - 1, "TVC"))
              }
              BETA[[Y_nGF]] <- mxPath(from = paste0("lx", records[1]), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF],
                                      labels = paste0(pfx[k], "beta", "g", "TVC"))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"),
                                                  name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "betagTVC)"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
          }
        }
      }
      else if (!intrinsic){
        latents <- c("eta0Y", "eta1Y", paste0("dy", records[-1]), paste0("ly", records))
        Y_nGF <- length(latents) - (length(records) * 2 - 1)
        if (decompose == 0){
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:2, 4)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:2, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:2, ] + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:2, 4)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                             name = paste0(pfx[k], "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:2, 4)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:2, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                    "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC))"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                             name = paste0(pfx[k], "BL_mean")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:2, ] + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
            else if (is.null(growth_TIC)){
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                    labels = paste0(pfx[k], "beta", p - 1, "TVC"))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:2],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                                  labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:2, 4)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                                labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"), name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC)"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:2, ] + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
          }
        }
      }
    }
    else if (curveFun %in% c("Jenss-Bayley", "JB")){
      if (intrinsic){
        latents <- c("eta0Y", "eta1Y", "eta2Y", "deltag", paste0("dy", records[-1]), paste0("ly", records))
        Y_nGF <- length(latents) - (length(records) * 2 - 1)
        if (decompose == 0){
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
              }
              BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", c(paste0("TIC", seq_along(growth_TIC)))))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L3", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22, ", pfx[k], "Y_psi2g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psi2g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:4, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "betagTIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L3", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22, ", pfx[k], "Y_psi2g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psi2g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              }
              BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF, ],
                                      labels = paste0(pfx[k], "beta", "g", c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L3", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22, ", pfx[k], "Y_psi2g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psi2g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:4, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "betagTIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                    "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC, ", pfx[k], "betagTVC))"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                             name = paste0(pfx[k], "BL_mean")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
            else if (is.null(growth_TIC)){
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                    labels = paste0(pfx[k], "beta", p - 1, "TVC"))
              }
              BETA[[Y_nGF]] <- mxPath(from = paste0("lx", records[1]), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[Y_nGF],
                                      labels = paste0(pfx[k], "beta", "g", "TVC"))
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]]$Y_starts$covMatrix,
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L3", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02, ", pfx[k], "Y_psi0g), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi1g), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22, ", pfx[k], "Y_psi2g), ",
                                                                    "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psi2g, ", pfx[k], "Y_psigg))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC, ", pfx[k], "betagTVC)"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
          }
        }
      }
      else if (!intrinsic){
        latents <- c("eta0Y", "eta1Y", "eta2Y", paste0("dy", records[-1]), paste0("ly", records))
        Y_nGF <- length(latents) - (length(records) * 2 - 1)
        if (decompose == 0){
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)))))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:3, ] + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                             name = paste0(pfx[k], "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]]$TVC_starts$means, labels = paste0(pfx[k], "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]]$TVC_starts$covMatrix, labels = paste0(pfx[k], "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = c(growth_TIC, paste0("lx", records[1])), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                    labels = paste0(pfx[k], "beta", p - 1, c(paste0("TIC", seq_along(growth_TIC)), "TVC")))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]]$beta[1:3, seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "beta0TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta1TIC", seq_along(growth_TIC)),
                                                             paste0(pfx[k], "beta2TIC", seq_along(growth_TIC))),
                                                  byrow = TRUE, name = paste0(pfx[k], "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(", pfx[k], "beta_TIC,",
                                                                    "rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC))"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]]$BL_starts$means[seq_along(growth_TIC)],
                                                  labels = c(paste0(pfx[k], "mux", seq_along(growth_TIC))),
                                                  byrow = FALSE, name = paste0(pfx[k], "mux")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "mux, ", pfx[k], "X_mueta0)"),
                                                             name = paste0(pfx[k], "BL_mean")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:3, ] + ", pfx[k], "beta %*% ", pfx[k], "BL_mean"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
            else if (is.null(growth_TIC)){
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = paste0("lx", records[1]), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p],
                                    labels = paste0(pfx[k], "beta", p - 1, "TVC"))
              }
              class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]]$Y_starts$means[1:3],
                                                labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                                  labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                                labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]]$Y_starts$residuals, labels = paste0(pfx[k], "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                             name = paste0(pfx[k], "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                    "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                    "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                             name = paste0(pfx[k], "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "X_mueta0, ", pfx[k], "X_mueta1)"),
                                                             name = paste0(pfx[k], "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "X_psi00, ", pfx[k], "X_psi01),",
                                                                    "cbind(", pfx[k], "X_psi01, ", pfx[k], "X_psi11))"),
                                                             name = paste0(pfx[k], "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(", pfx[k], "beta0TVC, ", pfx[k], "beta1TVC, ", pfx[k], "beta2TVC)"),
                                                             name = paste0(pfx[k], "beta")),
                                         mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:3, ] + ", pfx[k], "beta %*% ", pfx[k], "X_mueta0"),
                                                             name = paste0(pfx[k], "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]])
            }
          }
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
getTVC.mxModel <- function(dat, t_var, y_var, curveFun, intrinsic, records, y_model, TVC, decompose, growth_TIC,
                           starts){
  info <- getTVC.info(y_var = y_var, records = records, growth_TIC = growth_TIC, TVC = TVC, decompose = decompose,
                      starts = starts)
  loadings <- getUNI.loadings(y_model = y_model, t_var = t_var, y_var = y_var, curveFun = curveFun,
                              intrinsic = intrinsic, records = records)
  addpara <- if (y_model == "LCSM"){
    getUNI.addpara(dat = dat, curveFun = curveFun, intrinsic = intrinsic, t_var = t_var, records = records,
                   growth_TIC = growth_TIC, decompose = decompose, starts = starts)
  }
  .build_TVC(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic, records = records,
             y_model = y_model, TVC = TVC, decompose = decompose, growth_TIC = growth_TIC, starts = list(starts),
             pfx = "", model_names = .single_model_name_tvc(y_model, curveFun, intrinsic, growth_TIC, decompose),
             extra = list(list(mxData(observed = dat, type = "raw"))), GF_loadings = list(loadings),
             AddPara = list(addpara), TVC_info = list(info))[[1]]
}

#' @keywords internal
#' @noRd
getsub.TVC <- function(dat, nClass, t_var, y_var, curveFun, intrinsic, records, y_model, TVC, decompose, growth_TIC,
                       starts, container, grp_var = NULL){
  info <- getMIX_TVC.info(nClass = nClass, y_var = y_var, records = records, growth_TIC = growth_TIC, TVC = TVC,
                          decompose = decompose, starts = starts)
  loadings <- getMIX_UNI.loadings(nClass = nClass, y_model = y_model, t_var = t_var, records = records, y_var = y_var,
                                  curveFun = curveFun, intrinsic = intrinsic)
  addpara <- if (y_model == "LCSM"){
    getMIX_UNI.addpara(dat = dat, nClass = nClass, curveFun = curveFun, intrinsic = intrinsic, t_var = t_var,
                       records = records, growth_TIC = growth_TIC, decompose = decompose, starts = starts,
                       grp_var = if (container == "multigroup") grp_var else NULL)
  }
  .build_TVC(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic, records = records,
             y_model = y_model, TVC = TVC, decompose = decompose, growth_TIC = growth_TIC, starts = starts,
             pfx = paste0("c", seq_len(nClass)), model_names = paste0("Class", seq_len(nClass)),
             extra = .class_extras(dat, nClass, container, grp_var), GF_loadings = loadings, AddPara = addpara,
             TVC_info = info)
}

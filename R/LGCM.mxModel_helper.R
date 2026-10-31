#' @title One builder for the LGCM sub-model, single-group or class-specific
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
.build_LGCM <- function(dat, t_var, y_var, curveFun, intrinsic, records, growth_TIC, starts, pfx, model_names, extra, GF_loadings){
  nClass <- length(pfx)
  ## Define manifest variables
  manifests <- paste0(y_var, records)
  ## Define paths for adding growth TICs if any
  TIC_mean <- TIC_VAR <- BETA <- list()
  if (!is.null(growth_TIC)){
    for (k in 1:nClass){
      nTICs <- length(growth_TIC)
      ### Y_mean values of TIC(s)
      TIC_mean[[k]] <- mxPath(from = "one", to = growth_TIC, arrows = 1, free = TRUE, values = starts[[k]]$TIC_starts$means,
                              labels = paste0(pfx[k], "mux", 1:nTICs))
      ### Var-cov of TIC(s)
      TIC_VAR[[k]] <- mxPath(from = growth_TIC, to = growth_TIC, connect = "unique.pairs", arrows = 2, free = TRUE,
                             values = starts[[k]]$TIC_starts$covMatrix,
                             labels = paste0(pfx[k], "phi", 1:(nTICs * (nTICs + 1)/2)))
    }
  }
  class.list <- list()
  ## Define latent variables, growth factor loadings, paths of the longitudinal outcome
  if (curveFun %in% c("linear", "LIN")){
    latents <- c("eta0", "eta1")
    for (k in 1:nClass){
      if (!is.null(growth_TIC)){
        nGF <- length(latents)
        nTICs <- length(growth_TIC)
        for (p in 1:nGF){
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                              labels = paste0(pfx[k], "beta", p - 1, 1:nTICs))
        }
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                   mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                          values = starts[[k]]$Y_starts$means,
                                          labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                   mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                          labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                   mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                          values = 1),
                                   mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                          values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                          labels = paste0(pfx[k], "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                       name = paste0(pfx[k], "Y_alpha0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                              "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                       name = paste0(pfx[k], "Y_psi_r")),
                                   mxMatrix("Full", 2, length(growth_TIC), free = TRUE, values = starts[[k]]$beta,
                                            labels = c(paste0(pfx[k], "beta0", seq_along(growth_TIC)),
                                                       paste0(pfx[k], "beta1", seq_along(growth_TIC))),
                                            byrow = TRUE, name = paste0(pfx[k], "beta")),
                                   mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]]$TIC_starts$means,
                                            labels = paste0(pfx[k], "mux", seq_along(growth_TIC)),
                                            byrow = FALSE, name = paste0(pfx[k], "mux")),
                                   mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                       name = paste0(pfx[k], "Y_mean0")),
                                   GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
      }
      else if (is.null(growth_TIC)){
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                          values = starts[[k]]$Y_starts$means,
                                          labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                   mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                          labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                   mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                          values = 1),
                                   mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                          values = 0, labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                          labels = paste0(pfx[k], "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                       name = paste0(pfx[k], "Y_mean0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                              "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                       name = paste0(pfx[k], "Y_psi0")),
                                   GF_loadings[[k]])
      }
    }
  }
  else if (curveFun %in% c("quadratic", "QUAD")){
    latents <- c("eta0", "eta1", "eta2")
    for (k in 1:nClass){
      if (!is.null(growth_TIC)){
        nGF <- length(latents)
        nTICs <- length(growth_TIC)
        for (p in 1:nGF){
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                              labels = paste0(pfx[k], "beta", p - 1, 1:nTICs))
        }
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                   mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                          values = starts[[k]]$Y_starts$means,
                                          labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                   mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                          labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                   mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                   mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                          labels = paste0(pfx[k], "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"),
                                                       name = paste0(pfx[k], "Y_alpha0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                              "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                              "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                       name = paste0(pfx[k], "Y_psi_r")),
                                   mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[k]]$beta,
                                            labels = c(paste0(pfx[k], "beta0", seq_along(growth_TIC)),
                                                       paste0(pfx[k], "beta1", seq_along(growth_TIC)),
                                                       paste0(pfx[k], "beta2", seq_along(growth_TIC))),
                                            byrow = TRUE, name = paste0(pfx[k], "beta")),
                                   mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]]$TIC_starts$means,
                                            labels = paste0(pfx[k], "mux", seq_along(growth_TIC)),
                                            byrow = FALSE, name = paste0(pfx[k], "mux")),
                                   mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                       name = paste0(pfx[k], "Y_mean0")),
                                   GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
      }
      else if (is.null(growth_TIC)){
        class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                          values = starts[[k]]$Y_starts$means,
                                          labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                   mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                          labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                   mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                   mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                          labels = paste0(pfx[k], "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"),
                                                       name = paste0(pfx[k], "Y_mean0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                              "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                              "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                       name = paste0(pfx[k], "Y_psi0")),
                                   GF_loadings[[k]])
      }
    }
  }
  else if (curveFun %in% c("negative exponential", "EXP")){
    if (intrinsic){
      latents <- c("eta0", "eta1", "deltag")
      for (k in 1:nClass){
        if (!is.null(growth_TIC)){
          nGF <- length(latents)
          nTICs <- length(growth_TIC)
          for (p in 1:(nGF - 1)){
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                labels = paste0(pfx[k], "beta", p - 1, 1:nTICs))
          }
          BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[nGF, ],
                                labels = paste0(pfx[k], "beta", "g", 1:nTICs))
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:2],
                                            labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                              labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                            labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"),
                                                         name = paste0(pfx[k], "Y_alpha0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                         name = paste0(pfx[k], "Y_psi_r")),
                                     mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[k]]$beta,
                                              labels = c(paste0(pfx[k], "beta0", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "beta1", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "betag", seq_along(growth_TIC))),
                                              byrow = TRUE, name = paste0(pfx[k], "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]]$TIC_starts$means,
                                              labels = paste0(pfx[k], "mux", seq_along(growth_TIC)),
                                              byrow = FALSE, name = paste0(pfx[k], "mux")),
                                     mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:2],
                                            labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                              labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                            labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                         name = paste0(pfx[k], "Y_psi0")),
                                     GF_loadings[[k]])
        }
      }
    }
    else if (!intrinsic){
      latents <- c("eta0", "eta1")
      for (k in 1:nClass){
        if (!is.null(growth_TIC)){
          nGF <- length(latents)
          nTICs <- length(growth_TIC)
          for (p in 1:nGF){
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                labels = paste0(pfx[k], "beta", p - 1, 1:nTICs))
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:2],
                                            labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                              labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix[c(1:2, 4)],
                                            labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"),
                                                         name = paste0(pfx[k], "Y_alpha0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                         name = paste0(pfx[k], "Y_psi_r")),
                                     mxMatrix("Full", 2, length(growth_TIC), free = TRUE, values = starts[[k]]$beta[1:2, ],
                                              labels = c(paste0(pfx[k], "beta0", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "beta1", seq_along(growth_TIC))),
                                              byrow = TRUE, name = paste0(pfx[k], "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]]$TIC_starts$means,
                                              labels = paste0(pfx[k], "mux", seq_along(growth_TIC)),
                                              byrow = FALSE, name = paste0(pfx[k], "mux")),
                                     mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:2, ] + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:2],
                                            labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[3],
                                              labels = paste0(pfx[k], "Y_slp_ratio"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix[c(1:2, 4)],
                                            labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                         name = paste0(pfx[k], "Y_psi0")),
                                     GF_loadings[[k]])
        }
      }
    }
  }
  else if (curveFun %in% c("Jenss-Bayley", "JB")){
    if (intrinsic){
      latents <- c("eta0", "eta1", "eta2", "deltag")
      for (k in 1:nClass){
        if (!is.null(growth_TIC)){
          nGF <- length(latents)
          nTICs <- length(growth_TIC)
          for (p in 1:(nGF - 1)){
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                labels = paste0(pfx[k], "beta", p - 1, 1:nTICs))
          }
          BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[nGF, ],
                                labels = paste0(pfx[k], "beta", "g", 1:nTICs))
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:3],
                                            labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                              labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                            labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                      "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                         name = paste0(pfx[k], "Y_alpha0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02, ", pfx[k], "Y_psi0g), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi1g), ",
                                                                "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22, ", pfx[k], "Y_psi2g), ",
                                                                "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psi2g, ", pfx[k], "Y_psigg))"),
                                                         name = paste0(pfx[k], "Y_psi_r")),
                                     mxMatrix("Full", 4, length(growth_TIC), free = TRUE, values = starts[[k]]$beta,
                                              labels = c(paste0(pfx[k], "beta0", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "beta1", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "beta2", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "betag", seq_along(growth_TIC))),
                                              byrow = TRUE, name = paste0(pfx[k], "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]]$TIC_starts$means,
                                              labels = paste0(pfx[k], "mux", seq_along(growth_TIC)),
                                              byrow = FALSE, name = paste0(pfx[k], "mux")),
                                     mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:3],
                                            labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                              labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                            labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                      "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02, ", pfx[k], "Y_psi0g), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi1g), ",
                                                                "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22, ", pfx[k], "Y_psi2g), ",
                                                                "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psi2g, ", pfx[k], "Y_psigg))"),
                                                         name = paste0(pfx[k], "Y_psi0")),
                                     GF_loadings[[k]])
        }
      }
    }
    else if (!intrinsic){
      latents <- c("eta0", "eta1", "eta2")
      for (k in 1:nClass){
        if (!is.null(growth_TIC)){
          nGF <- length(latents)
          nTICs <- length(growth_TIC)
          for (p in 1:nGF){
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                labels = paste0(pfx[k], "beta", p - 1, 1:nTICs))
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:3],
                                            labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                              labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                            labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                         name = paste0(pfx[k], "Y_alpha0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                         name = paste0(pfx[k], "Y_psi_r")),
                                     mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[k]]$beta[1:3, ],
                                              labels = c(paste0(pfx[k], "beta0", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "beta1", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "beta2", seq_along(growth_TIC))),
                                              byrow = TRUE, name = paste0(pfx[k], "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]]$TIC_starts$means,
                                              labels = paste0(pfx[k], "mux", seq_along(growth_TIC)),
                                              byrow = FALSE, name = paste0(pfx[k], "mux")),
                                     mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:3, ] + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:3],
                                            labels = paste0(pfx[k], c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                              labels = paste0(pfx[k], "Y_acc_ratio"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                            labels = paste0(pfx[k], c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                         name = paste0(pfx[k], "Y_psi0")),
                                     GF_loadings[[k]])
        }
      }
    }
  }
  else if (curveFun %in% c("bilinear spline", "BLS")){
    if (intrinsic){
      latents <- c("eta0s", "eta1s", "eta2s", "deltag")
      for (k in 1:nClass){
        if (!is.null(growth_TIC)){
          nGF <- length(latents)
          nTICs <- length(growth_TIC)
          for (p in 1:(nGF - 1)){
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                labels = paste0(pfx[k], "beta", p - 1, 1:nTICs))
          }
          BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[k]]$beta[nGF, ],
                                labels = paste0(pfx[k], "beta", "g", 1:nTICs))
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:3],
                                            labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                              labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                            labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                      "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                     mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                         name = paste0(pfx[k], "Y_alpha_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi0gs), ",
                                                                "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi1gs), ",
                                                                "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s, ", pfx[k], "Y_psi2gs), ",
                                                                "cbind(", pfx[k], "Y_psi0gs, ", pfx[k], "Y_psi1gs, ", pfx[k], "Y_psi2gs, ", pfx[k], "Y_psiggs))"),
                                                         name = paste0(pfx[k], "Y_psi_s")),
                                     mxMatrix("Full", 4, length(growth_TIC), free = TRUE, values = starts[[k]]$beta,
                                              labels = c(paste0(pfx[k], "beta0", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "beta1", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "beta2", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "betag", seq_along(growth_TIC))),
                                              byrow = TRUE, name = paste0(pfx[k], "beta_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, 0),",
                                                                "cbind(0, 1, -1, 0),",
                                                                "cbind(0, 1, 1, 0),",
                                                                "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "func")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, ", pfx[k], "Y_mueta2s - ", pfx[k], "Y_mueta1s),",
                                                                "cbind(0, 1, -1, 0),",
                                                                "cbind(0, 1, 1, 0),",
                                                                "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "grad")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "func[1:3, 1:3] %*% ", pfx[k], "Y_alpha_s, ", pfx[k], "Y_knot)"),
                                                         name = paste0(pfx[k], "Y_alpha0")),
                                     mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                         name = paste0(pfx[k], "Y_psi_r")),
                                     mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "beta_s"), name = paste0(pfx[k], "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]]$TIC_starts$means,
                                              labels = paste0(pfx[k], "mux", seq_along(growth_TIC)),
                                              byrow = FALSE, name = paste0(pfx[k], "mux")),
                                     mxAlgebraFromString(paste0(pfx[k], "Y_alpha0 + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:3],
                                            labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                              labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix,
                                            labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                      "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                     mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L3", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                         name = paste0(pfx[k], "Y_mean_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi0gs), ",
                                                                "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi1gs), ",
                                                                "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s, ", pfx[k], "Y_psi2gs), ",
                                                                "cbind(", pfx[k], "Y_psi0gs, ", pfx[k], "Y_psi1gs, ", pfx[k], "Y_psi2gs, ", pfx[k], "Y_psiggs))"),
                                                         name = paste0(pfx[k], "Y_psi_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, 0),",
                                                                "cbind(0, 1, -1, 0),",
                                                                "cbind(0, 1, 1, 0),",
                                                                "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "func")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot, ", pfx[k], "Y_mueta2s - ", pfx[k], "Y_mueta1s),",
                                                                "cbind(0, 1, -1, 0),",
                                                                "cbind(0, 1, 1, 0),",
                                                                "cbind(0, 0, 0, 1))"), name = paste0(pfx[k], "grad")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "func[1:3, 1:3] %*% ", pfx[k], "Y_mean_s, ", pfx[k], "Y_knot)"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                         name = paste0(pfx[k], "Y_psi0")),
                                     GF_loadings[[k]])
        }
      }
    }
    else if (!intrinsic){
      latents <- c("eta0s", "eta1s", "eta2s")
      for (k in 1:nClass){
        if (!is.null(growth_TIC)){
          nGF <- length(latents)
          nTICs <- length(growth_TIC)
          for (p in 1:nGF){
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]]$beta[p, ],
                                labels = paste0(pfx[k], "beta", p - 1, 1:nTICs))
          }
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:3],
                                            labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                              labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                            labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                     mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                         name = paste0(pfx[k], "Y_alpha_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s), ",
                                                                "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s), ",
                                                                "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s))"),
                                                         name = paste0(pfx[k], "Y_psi_s")),
                                     mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[k]]$beta[1:3, ],
                                              labels = c(paste0(pfx[k], "beta0", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "beta1", seq_along(growth_TIC)),
                                                         paste0(pfx[k], "beta2", seq_along(growth_TIC))),
                                              byrow = TRUE, name = paste0(pfx[k], "beta_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                "cbind(0, 1, -1),",
                                                                "cbind(0, 1, 1))"), name = paste0(pfx[k], "func")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                "cbind(0, 1, -1),",
                                                                "cbind(0, 1, 1))"), name = paste0(pfx[k], "grad")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "func %*% ", pfx[k], "Y_alpha_s, ", pfx[k], "Y_knot)"),
                                                         name = paste0(pfx[k], "Y_alpha0")),
                                     mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                         name = paste0(pfx[k], "Y_psi_r")),
                                     mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "beta_s"), name = paste0(pfx[k], "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]]$TIC_starts$means,
                                              labels = paste0(pfx[k], "mux", seq_along(growth_TIC)),
                                              byrow = FALSE, name = paste0(pfx[k], "mux")),
                                     mxAlgebraFromString(paste0(pfx[k], "Y_alpha0[1:3, 1] + ", pfx[k], "beta %*% ", pfx[k], "mux"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          class.list[[k]] <- mxModel(name = model_names[k], type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]]$Y_starts$means[1:3],
                                            labels = paste0(pfx[k], c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]]$Y_starts$means[4],
                                              labels = paste0(pfx[k], "Y_knot"), name = paste0(pfx[k], "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]]$Y_starts$covMatrix[c(1:3, 5:6, 8)],
                                            labels = paste0(pfx[k], c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                     mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records, "[1,1]")),
                                     mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0s, ", pfx[k], "Y_mueta1s, ", pfx[k], "Y_mueta2s)"),
                                                         name = paste0(pfx[k], "Y_mean_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00s, ", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi02s), ",
                                                                "cbind(", pfx[k], "Y_psi01s, ", pfx[k], "Y_psi11s, ", pfx[k], "Y_psi12s), ",
                                                                "cbind(", pfx[k], "Y_psi02s, ", pfx[k], "Y_psi12s, ", pfx[k], "Y_psi22s))"),
                                                         name = paste0(pfx[k], "Y_psi_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                "cbind(0, 1, -1),",
                                                                "cbind(0, 1, 1))"), name = paste0(pfx[k], "func")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -", pfx[k], "Y_knot, ", pfx[k], "Y_knot),",
                                                                "cbind(0, 1, -1),",
                                                                "cbind(0, 1, 1))"), name = paste0(pfx[k], "grad")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "func %*% ", pfx[k], "Y_mean_s, ", pfx[k], "Y_knot)"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     mxAlgebraFromString(paste0(pfx[k], "grad %*% ", pfx[k], "Y_psi_s %*% t(", pfx[k], "grad)"),
                                                         name = paste0(pfx[k], "Y_psi0")),
                                     GF_loadings[[k]])
        }
      }
    }
  }
  ## Define the output of the function
  for (k in seq_len(nClass)) class.list[[k]] <- do.call(mxModel, c(list(class.list[[k]]), extra[[k]]))
  class.list
}

#' @title Single-group latent growth curve model
#' @description Thin wrapper around \code{.build_LGCM()} with an empty prefix and the data attached.
#' @keywords internal
#' @noRd
#' @importFrom OpenMx mxData
getLGCM.mxModel <- function(dat, t_var, y_var, curveFun, intrinsic, records, growth_TIC, starts){
  loadings <- getUNI.loadings(y_model = "LGCM", t_var = t_var, y_var = y_var, curveFun = curveFun,
                              intrinsic = intrinsic, records = records)
  .build_LGCM(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic, records = records,
              growth_TIC = growth_TIC, starts = list(starts), pfx = "",
              model_names = .single_model_name("LGCM", curveFun, intrinsic, growth_TIC),
              extra = list(list(mxData(observed = dat, type = "raw"))), GF_loadings = list(loadings))[[1]]
}

#' @title Class-specific latent growth curve models for a multiple-group or mixture model
#' @description One model per class with prefixes c1, c2, ...; \code{container} decides what each class carries
#' (\code{"multigroup"}: its rows of the data; \code{"mixture"}: the row-wise fit function).
#' @keywords internal
#' @noRd
getsub.LGCM <- function(dat, nClass, t_var, y_var, curveFun, intrinsic, records, growth_TIC, starts, container,
                        grp_var = NULL){
  loadings <- getMIX_UNI.loadings(nClass = nClass, y_model = "LGCM", t_var = t_var, records = records, y_var = y_var,
                                  curveFun = curveFun, intrinsic = intrinsic)
  .build_LGCM(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic, records = records,
              growth_TIC = growth_TIC, starts = starts, pfx = paste0("c", seq_len(nClass)),
              model_names = paste0("Class", seq_len(nClass)),
              extra = .class_extras(dat, nClass, container, grp_var), GF_loadings = loadings)
}

#' @title One builder for the LCSM sub-model, single-group or class-specific
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
.build_LCSM <- function(dat, t_var, y_var, curveFun, intrinsic, records, growth_TIC, starts, pfx, model_names, extra, GF_loadings, AddPara){
  nClass <- length(pfx)
  ## Define manifest variables
  manifests <- paste0(y_var, records)
  ## Define paths for adding growth TICs if any
  PATH_L <- PATH_SLP <- PATH_AUTO <- TIC_MEAN <- TIC_VAR <- BETA <- list()
  for (k in 1:nClass){
    ## Define paths from latent true scores to observed scores
    PATH_L[[k]] <- mxPath(from = paste0("ly", records), to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1)
    ## Define paths from latent instantaneous rate of change at each measurement to true scores
    PATH_SLP[[k]] <- mxPath(from = paste0("dy", records[-1]), to = paste0("ly", records[-1]), arrows = 1, free = FALSE,
                            values = 0, labels = paste0("lag", records[-1], "[1,1]"))
    #### Define autoregressive paths
    PATH_AUTO[[k]] <- mxPath(from = paste0("ly", records[-length(records)]), to = paste0("ly", records[-1]), arrows = 1,
                             free = FALSE, values = 1)
    if (!is.null(growth_TIC)){
      nTICs <- length(growth_TIC)
      ### Mean values of TIC(s)
      TIC_MEAN[[k]] <- mxPath(from = "one", to = growth_TIC, arrows = 1, free = TRUE, values = starts[[k]]$TIC_starts$means,
                              labels = paste0(pfx[k], "mux", 1:nTICs))
      ### Var-cov of TIC(s)
      TIC_VAR[[k]] <- mxPath(from = growth_TIC, to = growth_TIC, connect = "unique.pairs", arrows = 2, free = TRUE,
                             values = starts[[k]]$TIC_starts$covMatrix,
                             labels = paste0(pfx[k], "phi", 1:(nTICs * (nTICs + 1)/2)))
    }
  }
  # Obtain factor loadings for the specified functional form
  # Obtain additional parameters (derived) for the specified functional form
  class.list <- list()
  if (curveFun %in% c("nonparametric", "NonP")){
    latents <- c("eta0", "eta1", paste0("dy", records[-1]), paste0("ly", records))
    for (k in 1:nClass){
      if (!is.null(growth_TIC)){
        nGF <- length(latents) - length(records) * 2 + 1
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
                                   mxPath(from = "eta0", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1,
                                          free = c(FALSE, rep(TRUE, length(records) - 2)),
                                          values = c(1, starts[[k]]$Y_starts$rel_rate[-1]),
                                          labels = paste0(pfx[k], "Y_rel_rate", 1:(length(records) - 1))),
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
                                   PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                   BETA, GF_loadings[[k]], AddPara[[k]])
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
                                   mxPath(from = "eta0", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1,
                                          free = c(FALSE, rep(TRUE, length(records) - 2)),
                                          values = c(1, starts[[k]]$Y_starts$rel_rate[-1]),
                                          labels = paste0(pfx[k], "Y_rel_rate", 1:(length(records) - 1))),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                          labels = paste0(pfx[k], "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1)"),
                                                       name = paste0(pfx[k], "Y_mean0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                              "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                       name = paste0(pfx[k], "Y_psi0")),
                                   PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]])
      }
    }
  }
  else if (curveFun %in% c("quadratic", "QUAD")){
    latents <- c("eta0", "eta1", "eta2", paste0("dy", records[-1]), paste0("ly", records))
    for (k in 1:nClass){
      if (!is.null(growth_TIC)){
        nGF <- length(latents) - length(records) * 2 + 1
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
                                   mxPath(from = "eta0", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
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
                                   PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                   BETA, GF_loadings[[k]], AddPara[[k]])
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
                                   mxPath(from = "eta0", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                          labels = paste0(pfx[k], "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2)"),
                                                       name = paste0(pfx[k], "Y_mean0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                              "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                              "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                       name = paste0(pfx[k], "Y_psi0")),
                                   PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]])
      }
    }
  }
  else if (curveFun %in% c("negative exponential", "EXP")){
    if (intrinsic){
      latents <- c("eta0", "eta1", "deltag", paste0("dy", records[-1]), paste0("ly", records))
      for (k in 1:nClass){
        if (!is.null(growth_TIC)){
          nGF <- length(latents) - length(records) * 2 + 1
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
                                     mxPath(from = "eta0", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                     mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
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
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                     BETA, GF_loadings[[k]], AddPara[[k]])
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
                                     mxPath(from = "eta0", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                     mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi0g), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi1g), ",
                                                                "cbind(", pfx[k], "Y_psi0g, ", pfx[k], "Y_psi1g, ", pfx[k], "Y_psigg))"),
                                                         name = paste0(pfx[k], "Y_psi0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]])
        }
      }
    }
    else if (!intrinsic){
      latents <- c("eta0", "eta1", paste0("dy", records[-1]), paste0("ly", records))
      for (k in 1:nClass){
        if (!is.null(growth_TIC)){
          nGF <- length(latents) - length(records) * 2 + 1
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
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
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
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                     BETA, GF_loadings[[k]], AddPara[[k]])
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
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L1", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_slp_ratio)"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11))"),
                                                         name = paste0(pfx[k], "Y_psi0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]])
        }
      }
    }
  }
  else if (curveFun %in% c("Jenss-Bayley", "JB")){
    if (intrinsic){
      latents <- c("eta0", "eta1", "eta2", "deltag", paste0("dy", records[-1]), paste0("ly", records))
      for (k in 1:nClass){
        if (!is.null(growth_TIC)){
          nGF <- length(latents) - length(records) * 2 + 1
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
                                     mxPath(from = "eta0", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                     mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L3", records[-1], "[1,1]")),
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
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                     BETA, GF_loadings[[k]], AddPara[[k]])
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
                                     mxPath(from = "eta0", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                     mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L3", records[-1], "[1,1]")),
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
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]])
        }
      }
    }
    else if (!intrinsic){
      latents <- c("eta0", "eta1", "eta2", paste0("dy", records[-1]), paste0("ly", records))
      for (k in 1:nClass){
        if (!is.null(growth_TIC)){
          nGF <- length(latents) - length(records) * 2 + 1
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
                                     mxPath(from = "eta0", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
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
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                     BETA, GF_loadings[[k]], AddPara[[k]])
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
                                     mxPath(from = "eta0", to = paste0("ly", records[1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0(pfx[k], "L2", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]]$Y_starts$residuals,
                                            labels = paste0(pfx[k], "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(", pfx[k], "Y_mueta0, ", pfx[k], "Y_mueta1, ", pfx[k], "Y_mueta2, ", pfx[k], "Y_acc_ratio)"),
                                                         name = paste0(pfx[k], "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(", pfx[k], "Y_psi00, ", pfx[k], "Y_psi01, ", pfx[k], "Y_psi02), ",
                                                                "cbind(", pfx[k], "Y_psi01, ", pfx[k], "Y_psi11, ", pfx[k], "Y_psi12), ",
                                                                "cbind(", pfx[k], "Y_psi02, ", pfx[k], "Y_psi12, ", pfx[k], "Y_psi22))"),
                                                         name = paste0(pfx[k], "Y_psi0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]])
        }
      }
    }
  }
  ## Define the output of the function
  for (k in seq_len(nClass)) class.list[[k]] <- do.call(mxModel, c(list(class.list[[k]]), extra[[k]]))
  class.list
}

#' @keywords internal
#' @noRd
#' @importFrom OpenMx mxData
getLCSM.mxModel <- function(dat, t_var, y_var, curveFun, intrinsic, records, growth_TIC, starts){
  loadings <- getUNI.loadings(y_model = "LCSM", t_var = t_var, y_var = y_var, curveFun = curveFun,
                              intrinsic = intrinsic, records = records)
  addpara <- getUNI.addpara(dat = dat, curveFun = curveFun, intrinsic = intrinsic, t_var = t_var,
                            records = records, growth_TIC = growth_TIC, decompose = NULL, starts = starts)
  .build_LCSM(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic, records = records,
              growth_TIC = growth_TIC, starts = list(starts), pfx = "",
              model_names = .single_model_name("LCSM", curveFun, intrinsic, growth_TIC),
              extra = list(list(mxData(observed = dat, type = "raw"))), GF_loadings = list(loadings),
              AddPara = list(addpara))[[1]]
}

#' @keywords internal
#' @noRd
getsub.LCSM <- function(dat, nClass, t_var, y_var, curveFun, intrinsic, records, growth_TIC, starts, container,
                        grp_var = NULL){
  loadings <- getMIX_UNI.loadings(nClass = nClass, y_model = "LCSM", t_var = t_var, records = records, y_var = y_var,
                                  curveFun = curveFun, intrinsic = intrinsic)
  addpara <- getMIX_UNI.addpara(dat = dat, nClass = nClass, curveFun = curveFun, intrinsic = intrinsic,
                                t_var = t_var, records = records, growth_TIC = growth_TIC, decompose = NULL,
                                starts = starts, grp_var = if (container == "multigroup") grp_var else NULL)
  .build_LCSM(dat = dat, t_var = t_var, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic, records = records,
              growth_TIC = growth_TIC, starts = starts, pfx = paste0("c", seq_len(nClass)),
              model_names = paste0("Class", seq_len(nClass)), extra = .class_extras(dat, nClass, container, grp_var),
              GF_loadings = loadings, AddPara = addpara)
}

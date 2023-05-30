#' @title Define Latent Change Score Models as Class-specific Models (Submodels) for a Longitudinal Mixture Model
#'
#' @description This function defines latent change score models as class-specific models (submodels) for a longitudinal mixture model.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions for each longitudinal process, and time-invariant covariates (TICs) if any.
#' It takes the value passed from \code{getMIX()}.
#' @param nClass An integer specifying the number of latent classes for the mixture model. It takes the value passed from \code{getMIX()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' It takes the value passed from \code{getMIX()}.
#' @param records A numeric vector specifying indices of the study waves. It takes the value passed from \code{getMIX()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getMIX()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for latent change score
#' models include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}), \code{"Jenss-Bayley"}
#' (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the value passed from \code{getMIX()}.
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. It takes the value
#' passed from \code{getMIX()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing to the
#' variability of growth factors if any. It takes the value passed from \code{getMIX()}.
#' @param starts A list of initial values for the parameters, either takes the value passed from \code{getMIX()} or
#' derived by the helper function \code{getMIX.initial()}.
#'
#' @return A list of manifest and latent variables and paths for an mxModel object.
#'
#' @importFrom OpenMx mxPath mxModel mxAlgebraFromString mxMatrix mxFitFunctionML
#'
getsub.LCSM_l <- function(dat, nClass, t_var, records, y_var, curveFun, intrinsic, growth_TIC, starts){
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
      TIC_MEAN[[k]] <- mxPath(from = "one", to = growth_TIC, arrows = 1, free = TRUE, values = starts[[k]][[2]][[1]],
                              labels = paste0("c", k, "mux", 1:nTICs))
      ### Var-cov of TIC(s)
      TIC_VAR[[k]] <- mxPath(from = growth_TIC, to = growth_TIC, connect = "unique.pairs", arrows = 2, free = TRUE,
                             values = starts[[k]][[2]][[2]],
                             labels = paste0("c", k, "phi", 1:(nTICs * (nTICs + 1)/2)))
    }
  }
  # Obtain factor loadigs for the specified functional form
  GF_loadings <- getMIX_UNI.loadings(nClass = nClass, y_model = "LCSM", t_var = t_var, y_var = y_var,
                                     curveFun = curveFun, intrinsic = intrinsic, records = records)
  # Obtain additional parameters (derived) for the specified functional form
  AddPara <- getMIX_UNI.addpara(dat = dat, nClass = nClass, curveFun = curveFun, intrinsic = intrinsic,
                                t_var = t_var, records = records, growth_TIC = growth_TIC, decompose = NULL,
                                starts = starts)
  class.list <- list()
  if (curveFun %in% c("nonparametric", "NonP")){
    latents <- c("eta0", "eta1", paste0("dy", records[-1]), paste0("ly", records))
    for (k in 1:nClass){
      if (!is.null(growth_TIC)){
        nGF <- length(latents) - length(records) * 2 + 1
        nTICs <- length(growth_TIC)
        for (p in 1:nGF){
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                              labels = paste0("c", k, "beta", p - 1, 1:nTICs))
        }
        class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                   manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                   mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                          values = starts[[k]][[1]][[1]],
                                          labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                   mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]][[1]][[2]],
                                          labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                   mxPath(from = "eta0", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1,
                                          free = c(F, rep(T, length(records) - 2)),
                                          values = c(1, starts[[k]][[1]][[4]][-1]),
                                          labels = paste0("c", k, "Y_rel_rate", 1:(length(records) - 1))),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                          labels = paste0("c", k, "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                       name = paste0("c", k, "Y_alpha0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                              "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                       name = paste0("c", k, "Y_psi_r")),
                                   mxMatrix("Full", 2, length(growth_TIC), free = TRUE, values = starts[[k]][[3]],
                                            labels = c(paste0("c", k, "beta0", 1:length(growth_TIC)),
                                                       paste0("c", k, "beta1", 1:length(growth_TIC))),
                                            byrow = T, name = paste0("c", k, "beta")),
                                   mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]][[2]][[1]],
                                            labels = paste0("c", k, "mux", 1:length(growth_TIC)),
                                            byrow = F, name = paste0("c", k, "mux")),
                                   mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                       name = paste0("c", k, "Y_mean0")),
                                   PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                   BETA, GF_loadings[[k]], AddPara[[k]], mxFitFunctionML(vector = T))
      }
      else if (is.null(growth_TIC)){
        class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                          values = starts[[k]][[1]][[1]],
                                          labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                   mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]][[1]][[2]],
                                          labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                   mxPath(from = "eta0", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1,
                                          free = c(F, rep(T, length(records) - 2)),
                                          values = c(1, starts[[k]][[1]][[4]][-1]),
                                          labels = paste0("c", k, "Y_rel_rate", 1:(length(records) - 1))),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                          labels = paste0("c", k, "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                       name = paste0("c", k, "Y_mean0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                              "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                       name = paste0("c", k, "Y_psi0")),
                                   PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]],
                                   mxFitFunctionML(vector = T))
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
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                              labels = paste0("c", k, "beta", p - 1, 1:nTICs))
        }
        class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                   manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                   mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                          values = starts[[k]][[1]][[1]],
                                          labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                   mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]][[1]][[2]],
                                          labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                   mxPath(from = "eta0", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                          labels = paste0("c", k, "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2)"),
                                                       name = paste0("c", k, "Y_alpha0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                              "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                              "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                       name = paste0("c", k, "Y_psi_r")),
                                   mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[k]][[3]],
                                            labels = c(paste0("c", k, "beta0", 1:length(growth_TIC)),
                                                       paste0("c", k, "beta1", 1:length(growth_TIC)),
                                                       paste0("c", k, "beta2", 1:length(growth_TIC))),
                                            byrow = T, name = paste0("c", k, "beta")),
                                   mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]][[2]][[1]],
                                            labels = paste0("c", k, "mux", 1:length(growth_TIC)),
                                            byrow = F, name = paste0("c", k, "mux")),
                                   mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                       name = paste0("c", k, "Y_mean0")),
                                   PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                   BETA, GF_loadings[[k]], AddPara[[k]], mxFitFunctionML(vector = T))
      }
      else if (is.null(growth_TIC)){
        class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                   manifestVars = manifests, latentVars = latents,
                                   mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                          values = starts[[k]][[1]][[1]],
                                          labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                   mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]][[1]][[2]],
                                          labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                   mxPath(from = "eta0", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                          labels = paste0("c", k, "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2)"),
                                                       name = paste0("c", k, "Y_mean0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                              "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                              "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                       name = paste0("c", k, "Y_psi0")),
                                   PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]],
                                   mxFitFunctionML(vector = T))
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
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                                labels = paste0("c", k, "beta", p - 1, 1:nTICs))
          }
          BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[k]][[3]][nGF, ],
                                labels = paste0("c", k, "beta", "g", 1:nTICs))
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:2],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                              labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                     mxPath(from = "eta0", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                     mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_slp_ratio)"),
                                                         name = paste0("c", k, "Y_alpha0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                         name = paste0("c", k, "Y_psi_r")),
                                     mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[k]][[3]],
                                              labels = c(paste0("c", k, "beta0", 1:length(growth_TIC)),
                                                         paste0("c", k, "beta1", 1:length(growth_TIC)),
                                                         paste0("c", k, "betag", 1:length(growth_TIC))),
                                              byrow = T, name = paste0("c", k, "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]][[2]][[1]],
                                              labels = paste0("c", k, "mux", 1:length(growth_TIC)),
                                              byrow = F, name = paste0("c", k, "mux")),
                                     mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                     BETA, GF_loadings[[k]], AddPara[[k]], mxFitFunctionML(vector = T))
        }
        else if (is.null(growth_TIC)){
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:2],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                              labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                     mxPath(from = "eta0", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                     mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_slp_ratio)"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                         name = paste0("c", k, "Y_psi0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]],
                                     mxFitFunctionML(vector = T))
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
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                                labels = paste0("c", k, "beta", p - 1, 1:nTICs))
          }
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:2],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                              labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]][c(1:2, 4)],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_slp_ratio)"),
                                                         name = paste0("c", k, "Y_alpha0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                         name = paste0("c", k, "Y_psi_r")),
                                     mxMatrix("Full", 2, length(growth_TIC), free = TRUE, values = starts[[k]][[3]][1:2, ],
                                              labels = c(paste0("c", k, "beta0", 1:length(growth_TIC)),
                                                         paste0("c", k, "beta1", 1:length(growth_TIC))),
                                              byrow = T, name = paste0("c", k, "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]][[2]][[1]],
                                              labels = paste0("c", k, "mux", 1:length(growth_TIC)),
                                              byrow = F, name = paste0("c", k, "mux")),
                                     mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:2, ] + c", k, "beta %*% c", k, "mux"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                     BETA, GF_loadings[[k]], AddPara[[k]], mxFitFunctionML(vector = T))
        }
        else if (is.null(growth_TIC)){
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:2],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                              labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]][c(1:2, 4)],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_slp_ratio)"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                         name = paste0("c", k, "Y_psi0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]],
                                     mxFitFunctionML(vector = T))
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
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                                labels = paste0("c", k, "beta", p - 1, 1:nTICs))
          }
          BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[k]][[3]][nGF, ],
                                labels = paste0("c", k, "beta", "g", 1:nTICs))
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:3],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                              labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                      "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                     mxPath(from = "eta0", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                     mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L3", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                         name = paste0("c", k, "Y_alpha0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02, c", k, "Y_psi0g), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12, c", k, "Y_psi1g), ",
                                                                "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22, c", k, "Y_psi2g), ",
                                                                "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psi2g, c", k, "Y_psigg))"),
                                                         name = paste0("c", k, "Y_psi_r")),
                                     mxMatrix("Full", 4, length(growth_TIC), free = TRUE, values = starts[[k]][[3]],
                                              labels = c(paste0("c", k, "beta0", 1:length(growth_TIC)),
                                                         paste0("c", k, "beta1", 1:length(growth_TIC)),
                                                         paste0("c", k, "beta2", 1:length(growth_TIC)),
                                                         paste0("c", k, "betag", 1:length(growth_TIC))),
                                              byrow = T, name = paste0("c", k, "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]][[2]][[1]],
                                              labels = paste0("c", k, "mux", 1:length(growth_TIC)),
                                              byrow = F, name = paste0("c", k, "mux")),
                                     mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                     BETA, GF_loadings[[k]], AddPara[[k]], mxFitFunctionML(vector = T))
        }
        else if (is.null(growth_TIC)){
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:3],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                              labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                      "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                     mxPath(from = "eta0", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                     mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L3", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02, c", k, "Y_psi0g), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12, c", k, "Y_psi1g), ",
                                                                "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22, c", k, "Y_psi2g), ",
                                                                "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psi2g, c", k, "Y_psigg))"),
                                                         name = paste0("c", k, "Y_psi0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]],
                                     mxFitFunctionML(vector = T))
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
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                                labels = paste0("c", k, "beta", p - 1, 1:nTICs))
          }
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:3],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                              labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                     mxPath(from = "eta0", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                         name = paste0("c", k, "Y_alpha0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                         name = paste0("c", k, "Y_psi_r")),
                                     mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[k]][[3]][1:3, ],
                                              labels = c(paste0("c", k, "beta0", 1:length(growth_TIC)),
                                                         paste0("c", k, "beta1", 1:length(growth_TIC)),
                                                         paste0("c", k, "beta2", 1:length(growth_TIC))),
                                              byrow = T, name = paste0("c", k, "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]][[2]][[1]],
                                              labels = paste0("c", k, "mux", 1:length(growth_TIC)),
                                              byrow = F, name = paste0("c", k, "mux")),
                                     mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:3, ] + c", k, "beta %*% c", k, "mux"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], TIC_MEAN[[k]], TIC_VAR[[k]],
                                     BETA, GF_loadings[[k]], AddPara[[k]], mxFitFunctionML(vector = T))
        }
        else if (is.null(growth_TIC)){
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:3],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                              labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                     mxPath(from = "eta0", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta2", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                         name = paste0("c", k, "Y_psi0")),
                                     PATH_L[[k]], PATH_SLP[[k]], PATH_AUTO[[k]], GF_loadings[[k]], AddPara[[k]],
                                     mxFitFunctionML(vector = T))
        }
      }
    }
  }
  ## Define the output of the function
  return(class.list)
}

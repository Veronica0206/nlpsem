#' @title Construct An Object of mxModel for Latent Growth Curve Model with Time-invariant Covariates (If Any) To Be
#' Evaluated
#'
#' @description his function builds up an object of mxModel for a latent growth curve model with user-specified
#' functional form (including whether intrinsically nonlinear) with time-invariant covariates (if any).
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions, and time-invariant covariates (TICs) if any. It takes the value passed from \code{getLGCM()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' It takes the value passed from \code{getLGCM()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getLGCM()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for latent growth curve
#' models are: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}).
#' It takes the value passed from \code{getLGCM()}.
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. It takes the
#' value passed from \code{getLGCM()}.
#' @param records A numeric vector specifying indices of the study waves. It takes the value passed from \code{getLGCM()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing to the
#' variability of growth factors if any. It takes the value passed from \code{getLGCM()}.
#' @param starts A list containing initial values for the free parameters, either takes the value passed from \code{getLGCM()} or
#' derived by the helper function \code{getUNI.initial()}.
#'
#' @return A pre-optimized mxModel for a Latent Growth Curve Model.
#'
#' @keywords internal
#'
#' @importFrom OpenMx mxPath mxModel mxAlgebra mxMatrix mxData
#'
getLGCM.mxModel <- function(dat, t_var, y_var, curveFun, intrinsic, records, growth_TIC, starts){
  # Define manifest variables
  manifests <- paste0(y_var, records)
  # Define paths for adding growth TICs if any
  TIC_MEAN <- TIC_VAR <- BETA <- list()
  if (!is.null(growth_TIC)){
    nTICs <- length(growth_TIC)
    # Mean values of TIC(s)
    TIC_MEAN <- mxPath(from = "one", to = growth_TIC, arrows = 1, free = TRUE, values = starts[[2]][[1]],
                       labels = paste0("mux", 1:nTICs))
    # Var-cov of TIC(s)
    TIC_VAR <- mxPath(from = growth_TIC, to = growth_TIC, connect = "unique.pairs", arrows = 2, free = TRUE,
                      values = starts[[2]][[2]], labels = paste0("phi", 1:(nTICs * (nTICs + 1)/2)))
  }
  # Obtain factor loadigs for the specified functional form
  GF_loadings <- getUNI.loadings(y_model = "LGCM", t_var = t_var, y_var = y_var, curveFun = curveFun,
                                 intrinsic = intrinsic, records = records)
  # Define mxModel for a LGCM with linear functional form and TICs (if any)
  if (curveFun %in% c("linear", "LIN")){
    latents <- c("eta0", "eta1")
    if (!is.null(growth_TIC)){
      manifests <- c(manifests, growth_TIC)
      nGF <- length(latents)
      nTICs <- length(growth_TIC)
      for (p in 1:nGF){
        BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[3]][p, ],
                            labels = paste0("beta", p - 1, 1:nTICs))
      }
      model_mx <- mxModel(name = "LGCM_TICs, Linear Functional Form", type = "RAM",
                          mxData(observed = dat, type = "raw"),
                          manifestVars = manifests, latentVars = latents,
                          mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                 labels = c("Y_mueta0", "Y_mueta1")),
                          mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                 values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                          mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", records, "[1,1]")),
                          mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                 values = starts[[1]][[3]], labels = "Y_residuals"),
                          mxAlgebra(rbind(Y_mueta0, Y_mueta1), name = "Y_alpha0"),
                          mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                          cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                          mxMatrix("Full", 2, length(growth_TIC), free = TRUE, values = starts[[3]],
                                   labels = c(paste0("beta0", 1:length(growth_TIC)),
                                              paste0("beta1", 1:length(growth_TIC))), byrow = T, name = "beta"),
                          mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[2]][[1]],
                                   labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                          mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                          GF_loadings, TIC_MEAN, TIC_VAR, BETA)
    }
    else if (is.null(growth_TIC)){
      model_mx <- mxModel(name = "LGCM, Linear Functional Form", type = "RAM",
                          mxData(observed = dat, type = "raw"),
                          manifestVars = manifests, latentVars = latents,
                          mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                 labels = c("Y_mueta0", "Y_mueta1")),
                          mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                 values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                          mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", records, "[1,1]")),
                          mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                 values = starts[[1]][[3]], labels = "Y_residuals"),
                          mxAlgebra(rbind(Y_mueta0, Y_mueta1), name = "Y_mean0"),
                          mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                          cbind(Y_psi01, Y_psi11)), name = "Y_psi0"),
                          GF_loadings)
    }
  }
  # Define mxModel for a LGCM with quadratic functional form and TICs (if any)
  else if (curveFun %in% c("quadratic", "QUAD")){
    latents <- c("eta0", "eta1", "eta2")
    if (!is.null(growth_TIC)){
      manifests <- c(manifests, growth_TIC)
      nGF <- length(latents)
      nTICs <- length(growth_TIC)
      for (p in 1:nGF){
        BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[3]][p, ],
                            labels = paste0("beta", p - 1, 1:nTICs))
      }
      model_mx <- mxModel(name = "LGCM_TICs, Quadratic Functional Form", type = "RAM",
                          mxData(observed = dat, type = "raw"),
                          manifestVars = manifests, latentVars = latents,
                          mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                 labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                          mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                 values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                          mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", records, "[1,1]")),
                          mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L2", records, "[1,1]")),
                          mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                 values = starts[[1]][[3]], labels = "Y_residuals"),
                          mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2), name = "Y_alpha0"),
                          mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                          cbind(Y_psi01, Y_psi11, Y_psi12),
                                          cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                          mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[3]],
                                   labels = c(paste0("beta0", 1:length(growth_TIC)),
                                              paste0("beta1", 1:length(growth_TIC)),
                                              paste0("beta2", 1:length(growth_TIC))), byrow = T, name = "beta"),
                          mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[2]][[1]],
                                   labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                          mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                          GF_loadings, TIC_MEAN, TIC_VAR, BETA)
    }
    else if (is.null(growth_TIC)){
      model_mx <- mxModel(name = "LGCM, Quadratic Functional Form", type = "RAM",
                          mxData(observed = dat, type = "raw"),
                          manifestVars = manifests, latentVars = latents,
                          mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                 labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                          mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                 values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                          mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                          mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L1", records, "[1,1]")),
                          mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                 labels = paste0("L2", records, "[1,1]")),
                          mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                 values = starts[[1]][[3]], labels = "Y_residuals"),
                          mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2), name = "Y_mean0"),
                          mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                          cbind(Y_psi01, Y_psi11, Y_psi12),
                                          cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi0"),
                          GF_loadings)
    }
  }
  # Define mxModel for a LGCM with negative exponential functional form and TICs (if any)
  else if (curveFun %in% c("negative exponential", "EXP")){
    # Random ratio of growth rate
    if (intrinsic){
      latents <- c("eta0", "eta1", "deltag")
      if (!is.null(growth_TIC)){
        manifests <- c(manifests, growth_TIC)
        nGF <- length(latents)
        nTICs <- length(growth_TIC)
        for (p in 1:(nGF - 1)){
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[3]][p, ],
                              labels = paste0("beta", p - 1, 1:nTICs))
        }
        BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[3]][nGF, ],
                              labels = paste0("beta", "g", 1:nTICs))
        model_mx <- mxModel(name = "LGCM_TICs, Negative Exponential Functional Form, Random Ratio", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                   labels = c("Y_mueta0", "Y_mueta1")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio", name = "Y_mug"),
                            mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg")),
                            mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L2", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                            mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi0g),
                                            cbind(Y_psi01, Y_psi11, Y_psi1g),
                                            cbind(Y_psi0g, Y_psi1g, Y_psigg)), name = "Y_psi_r"),
                            mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[3]],
                                     labels = c(paste0("beta0", 1:length(growth_TIC)),
                                                paste0("beta1", 1:length(growth_TIC)),
                                                paste0("betag", 1:length(growth_TIC))), byrow = T, name = "beta"),
                            mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[2]][[1]],
                                     labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                            mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                            GF_loadings, TIC_MEAN, TIC_VAR, BETA)
      }
      else if (is.null(growth_TIC)){
        model_mx <- mxModel(name = "LGCM, Negative Exponential Functional Form, Random Ratio", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                   labels = c("Y_mueta0", "Y_mueta1")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio", name = "Y_mug"),
                            mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg")),
                            mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L2", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_mean0"),
                            mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi0g),
                                            cbind(Y_psi01, Y_psi11, Y_psi1g),
                                            cbind(Y_psi0g, Y_psi1g, Y_psigg)), name = "Y_psi0"),
                            GF_loadings)
      }
    }
    # fixed ratio of growth rate
    else if (!intrinsic){
      latents <- c("eta0", "eta1")
      if (!is.null(growth_TIC)){
        manifests <- c(manifests, growth_TIC)
        nGF <- length(latents)
        nTICs <- length(growth_TIC)
        for (p in 1:nGF){
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[3]][p, ],
                              labels = paste0("beta", p - 1, 1:nTICs))
        }
        model_mx <- mxModel(name = "LGCM_TICs, Negative Exponential Functional Form, Fixed Ratio", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                   labels = c("Y_mueta0", "Y_mueta1")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio", name = "Y_mug"),
                            mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]][c(1:2, 4)], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                            mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                            mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                            cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                            mxMatrix("Full", 2, length(growth_TIC), free = TRUE, values = starts[[3]][1:2, ],
                                     labels = c(paste0("beta0", 1:length(growth_TIC)),
                                                paste0("beta1", 1:length(growth_TIC))), byrow = T, name = "beta"),
                            mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[2]][[1]],
                                     labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                            mxAlgebra(Y_alpha0[1:2, ] + beta %*% mux, name = "Y_mean0"),
                            GF_loadings, TIC_MEAN, TIC_VAR, BETA)
      }
      else if (is.null(growth_TIC)){
        model_mx <- mxModel(name = "LGCM, Negative Exponential Functional Form, Fixed Ratio", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                   labels = c("Y_mueta0", "Y_mueta1")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio", name = "Y_mug"),
                            mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]][c(1:2, 4)], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                            mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_mean0"),
                            mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                            cbind(Y_psi01, Y_psi11)), name = "Y_psi0"),
                            GF_loadings)
      }
    }
  }
  # Define mxModel for a LGCM with Jenss-Bayley functional form and TICs (if any)
  else if (curveFun %in% c("Jenss-Bayley", "JB")){
    # Random ratio of growth acceleration
    if (intrinsic){
      latents <- c("eta0", "eta1", "eta2", "deltag")
      if (!is.null(growth_TIC)){
        manifests <- c(manifests, growth_TIC)
        nGF <- length(latents)
        nTICs <- length(growth_TIC)
        for (p in 1:(nGF - 1)){
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[3]][p, ],
                              labels = paste0("beta", p - 1, 1:nTICs))
        }
        BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[3]][nGF, ],
                              labels = paste0("beta", "g", 1:nTICs))
        model_mx <- mxModel(name = "LGCM_TICs, Jenss Bayley Functional Form, Random Ratio", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                   labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio", name = "Y_mug"),
                            mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                         "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg")),
                            mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L2", records, "[1,1]")),
                            mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L3", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                            mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02, Y_psi0g),
                                            cbind(Y_psi01, Y_psi11, Y_psi12, Y_psi1g),
                                            cbind(Y_psi02, Y_psi12, Y_psi22, Y_psi2g),
                                            cbind(Y_psi0g, Y_psi1g, Y_psi2g, Y_psigg)), name = "Y_psi_r"),
                            mxMatrix("Full", 4, length(growth_TIC), free = TRUE, values = starts[[3]],
                                     labels = c(paste0("beta0", 1:length(growth_TIC)),
                                                paste0("beta1", 1:length(growth_TIC)),
                                                paste0("beta2", 1:length(growth_TIC)),
                                                paste0("betag", 1:length(growth_TIC))), byrow = T, name = "beta"),
                            mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[2]][[1]],
                                     labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                            mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                            GF_loadings, TIC_MEAN, TIC_VAR, BETA)
      }
      else if (is.null(growth_TIC)){
        model_mx <- mxModel(name = "LGCM, Jenss Bayley Functional Form, Random Ratio", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                   labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio", name = "Y_mug"),
                            mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                         "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg")),
                            mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L2", records, "[1,1]")),
                            mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L3", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_mean0"),
                            mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02, Y_psi0g),
                                            cbind(Y_psi01, Y_psi11, Y_psi12, Y_psi1g),
                                            cbind(Y_psi02, Y_psi12, Y_psi22, Y_psi2g),
                                            cbind(Y_psi0g, Y_psi1g, Y_psi2g, Y_psigg)), name = "Y_psi0"),
                            GF_loadings)
      }
    }
    # Fixed ratio of growth acceleration
    else if (!intrinsic){
      latents <- c("eta0", "eta1", "eta2")
      if (!is.null(growth_TIC)){
        manifests <- c(manifests, growth_TIC)
        nGF <- length(latents)
        nTICs <- length(growth_TIC)
        for (p in 1:nGF){
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[3]][p, ],
                              labels = paste0("beta", p - 1, 1:nTICs))
        }
        model_mx <- mxModel(name = "LGCM_TICs, Jenss Bayley Functional Form, Fixed Ratio", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                   labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio", name = "Y_mug"),
                            mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                   labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                            mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L2", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                            mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                            cbind(Y_psi01, Y_psi11, Y_psi12),
                                            cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                            mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[3]][1:3, ],
                                     labels = c(paste0("beta0", 1:length(growth_TIC)),
                                                paste0("beta1", 1:length(growth_TIC)),
                                                paste0("beta2", 1:length(growth_TIC))), byrow = T, name = "beta"),
                            mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[2]][[1]],
                                     labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                            mxAlgebra(Y_alpha0[1:3, ] + beta %*% mux, name = "Y_mean0"),
                            GF_loadings, TIC_MEAN, TIC_VAR, BETA)
      }
      else if (is.null(growth_TIC)){
        model_mx <- mxModel(name = "LGCM, Jenss Bayley Functional Form, Fixed Ratio", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                   labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio", name = "Y_mug"),
                            mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                   labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                            mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L2", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_mean0"),
                            mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                            cbind(Y_psi01, Y_psi11, Y_psi12),
                                            cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi0"),
                            GF_loadings)
      }
    }
  }
  # Define mxModel for a LGCM with bilinear spline functional form and TICs (if any)
  else if (curveFun %in% c("bilinear spline", "BLS")){
    # Random knot
    if (intrinsic){
      latents <- c("eta0s", "eta1s", "eta2s", "deltag")
      if (!is.null(growth_TIC)){
        manifests <- c(manifests, growth_TIC)
        nGF <- length(latents)
        nTICs <- length(growth_TIC)
        for (p in 1:(nGF - 1)){
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[3]][p, ],
                              labels = paste0("beta", p - 1, 1:nTICs, "s"))
        }
        BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[3]][nGF, ],
                              labels = paste0("beta", "g", 1:nTICs, "s"))
        model_mx <- mxModel(name = "LGCM_TICs, Bilinear Spline Functional Form, Random knot", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                   labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "knot", name = "Y_mug"),
                            mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]], labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                         "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs")),
                            mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L2", records, "[1,1]")),
                            mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L3", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(cbind(1, -knot, knot, 0),
                                            cbind(0, 1, -1, 0),
                                            cbind(0, 1, 1, 0),
                                            cbind(0, 0, 0, 1)), name = "func"),
                            mxAlgebra(rbind(cbind(1, -knot, knot, 0),
                                            cbind(0, 1, -1, 0),
                                            cbind(0, 1, 1, 0),
                                            cbind(0, 0, 0, 1)), name = "grad"),
                            mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_alpha_s"),
                            mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s, Y_psi0gs),
                                            cbind(Y_psi01s, Y_psi11s, Y_psi12s, Y_psi1gs),
                                            cbind(Y_psi02s, Y_psi12s, Y_psi22s, Y_psi2gs),
                                            cbind(Y_psi0gs, Y_psi1gs, Y_psi2gs, Y_psiggs)), name = "Y_psi_s"),
                            mxMatrix("Full", 4, length(growth_TIC), free = TRUE, values = starts[[3]],
                                     labels = c(paste0("beta0", 1:length(growth_TIC), "s"),
                                                paste0("beta1", 1:length(growth_TIC), "s"),
                                                paste0("beta2", 1:length(growth_TIC), "s"),
                                                paste0("betag", 1:length(growth_TIC), "s")),
                                     byrow = T, name = "beta_s"),
                            mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[2]][[1]],
                                     labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                            mxAlgebra(rbind(func[1:3, 1:3] %*% Y_alpha_s, knot), name = "Y_alpha0"),
                            mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi_r"),
                            mxAlgebra(grad %*% beta_s, name = "beta"),
                            mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                            GF_loadings, TIC_MEAN, TIC_VAR, BETA)
      }
      else if (is.null(growth_TIC)){
        model_mx <- mxModel(name = "LGCM, Bilinear Spline Functional Form, Random knot", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                   labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "knot", name = "Y_mug"),
                            mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]], labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                         "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs")),
                            mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L2", records, "[1,1]")),
                            mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L3", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(cbind(1, -knot, knot, 0),
                                            cbind(0, 1, -1, 0),
                                            cbind(0, 1, 1, 0),
                                            cbind(0, 0, 0, 1)), name = "func"),
                            mxAlgebra(rbind(cbind(1, -knot, knot, 0),
                                            cbind(0, 1, -1, 0),
                                            cbind(0, 1, 1, 0),
                                            cbind(0, 0, 0, 1)), name = "grad"),
                            mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_mean_s"),
                            mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s, Y_psi0gs),
                                            cbind(Y_psi01s, Y_psi11s, Y_psi12s, Y_psi1gs),
                                            cbind(Y_psi02s, Y_psi12s, Y_psi22s, Y_psi2gs),
                                            cbind(Y_psi0gs, Y_psi1gs, Y_psi2gs, Y_psiggs)), name = "Y_psi_s"),
                            mxAlgebra(rbind(func[1:3, 1:3] %*% Y_mean_s, knot), name = "Y_mean0"),
                            mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi0"),
                            GF_loadings)
      }
    }
    # Fixed knot
    else if (!intrinsic){
      latents <- c("eta0s", "eta1s", "eta2s")
      if (!is.null(growth_TIC)){
        manifests <- c(manifests, growth_TIC)
        nGF <- length(latents)
        nTICs <- length(growth_TIC)
        for (p in 1:nGF){
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[3]][p, ],
                              labels = paste0("beta", p - 1, 1:nTICs, "s"))
        }
        model_mx <- mxModel(name = "LGCM_TICs, Bilinear Spline Functional Form, Fixed knot", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                   labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "knot", name = "Y_mug"),
                            mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                   labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi11s", "Y_psi12s", "Y_psi22s")),
                            mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L2", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(cbind(1, -knot, knot),
                                            cbind(0, 1, -1),
                                            cbind(0, 1, 1)), name = "func"),
                            mxAlgebra(rbind(cbind(1, -knot, knot),
                                            cbind(0, 1, -1),
                                            cbind(0, 1, 1)), name = "grad"),
                            mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_alpha_s"),
                            mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s),
                                            cbind(Y_psi01s, Y_psi11s, Y_psi12s),
                                            cbind(Y_psi02s, Y_psi12s, Y_psi22s)), name = "Y_psi_s"),
                            mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[3]][1:3, ],
                                     labels = c(paste0("beta0", 1:length(growth_TIC), "s"),
                                                paste0("beta1", 1:length(growth_TIC), "s"),
                                                paste0("beta2", 1:length(growth_TIC), "s")), byrow = T, name = "beta_s"),
                            mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[2]][[1]],
                                     labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                            mxAlgebra(rbind(func %*% Y_alpha_s, knot), name = "Y_alpha0"),
                            mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi_r"),
                            mxAlgebra(grad %*% beta_s, name = "beta"),
                            mxAlgebra(Y_alpha0[1:3, 1] + beta %*% mux, name = "Y_mean0"),
                            GF_loadings, TIC_MEAN, TIC_VAR, BETA)
      }
      else if (is.null(growth_TIC)){
        model_mx <- mxModel(name = "LGCM, Bilinear Spline Functional Form, Fixed knot", type = "RAM",
                            mxData(observed = dat, type = "raw"),
                            manifestVars = manifests, latentVars = latents,
                            mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                   labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                            mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "knot", name = "Y_mug"),
                            mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                   labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi11s", "Y_psi12s", "Y_psi22s")),
                            mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                            mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L1", records, "[1,1]")),
                            mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                   labels = paste0("L2", records, "[1,1]")),
                            mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2, free = TRUE,
                                   values = starts[[1]][[3]], labels = "Y_residuals"),
                            mxAlgebra(rbind(cbind(1, -knot, knot),
                                            cbind(0, 1, -1),
                                            cbind(0, 1, 1)), name = "func"),
                            mxAlgebra(rbind(cbind(1, -knot, knot),
                                            cbind(0, 1, -1),
                                            cbind(0, 1, 1)), name = "grad"),
                            mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_mean_s"),
                            mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s),
                                            cbind(Y_psi01s, Y_psi11s, Y_psi12s),
                                            cbind(Y_psi02s, Y_psi12s, Y_psi22s)), name = "Y_psi_s"),
                            mxAlgebra(rbind(func %*% Y_mean_s, knot), name = "Y_mean0"),
                            mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi0"),
                            GF_loadings)
      }
    }
  }
  return(model_mx)
}

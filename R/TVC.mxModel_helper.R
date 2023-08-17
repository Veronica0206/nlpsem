#' @title Construct An Object of mxModel for Latent Growth Curve Models or Latent Change Score Models with a Time
#' Varying Covariate and Time-invariant Covariates (If Any) To Be Evaluated
#'
#' @description This function builds up an object of mxModel for a latent growth curve model or latent change score model
#' with user-specified functional form (including whether intrinsically nonlinear), time-varying covariate,
#' and with time-invariant covariates (if any).
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements (for the longitudinal outcome and time-varying covariates), occasions, and time-invariant covariates
#' (TICs) if any. It takes the value passed from \code{getTVCmodel()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' It takes the value passed from \code{getTVCmodel()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getTVCmodel()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the
#' value passed from \code{getTVCmodel()}.
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. It takes the
#' value passed from \code{getTVCmodel()}.
#' @param records A numeric vector specifying the indices of the observed study waves. It takes the value passed from
#' \code{getTVCmodel()}.
#' @param y_model A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' It takes the value passed from \code{getTVCmodel()}.
#' @param TVC A string specifying the prefix of the column names corresponding to the time-varying covariate at each study wave.
#' It takes the value passed from \code{getTVCmodel()}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). It takes the
#' value passed from \code{getTVCmodel()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) that account for the
#' variability of growth factors, if any. It takes the value passed from \code{getTVCmodel()}.
#' @param starts A list of initial values for the parameters, either takes the value passed from \code{getTVCmodel()}
#' or derived by the helper function \code{getTVC.initial()}.
#'
#' @return A pre-optimized mxModel for a latent growth curve model or a latent change score model with a
#' time-varying covariate and time-invariant covariates (if any).
#'
#' @keywords internal
#'
#' @importFrom OpenMx mxPath mxModel mxAlgebra mxMatrix
#'
getTVC.mxModel <- function(dat, t_var, y_var, curveFun, intrinsic, records, y_model, TVC, decompose,
                           growth_TIC, starts){
  # Define manifest variables
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
  TVC_info <- getTVC.info(y_var = y_var, records = records, growth_TIC = growth_TIC, TVC = TVC, decompose = decompose,
                          starts = starts)
  BETA <- list()
  # Obtain factor loadigs for the specified functional form
  GF_loadings <- getUNI.loadings(y_model = y_model, t_var = t_var, y_var = y_var, curveFun = curveFun,
                                 intrinsic = intrinsic, records = records)
  if (y_model == "LGCM"){
    if (curveFun %in% c("linear", "LIN")){
      latents <- c("eta0Y", "eta1Y")
      Y_nGF <- length(latents)
      if (decompose == 0){
        if (!is.null(growth_TIC)){
          nTICs <- length(growth_TIC)
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
          }
          model_mx <- mxModel(name = "LGCM with TVC and TICs, Linear Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1")),
                              mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]],
                                     labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                              mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L1", records, "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                              cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                              mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                     values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                              mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                     values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                              mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                       values = starts[[4]][1:2, 1:length(growth_TIC)],
                                       labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                  paste0("beta1TIC", 1:length(growth_TIC))), byrow = T, name = "beta"),
                              mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                       labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                              mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                              TVC_info, GF_loadings, BETA)
        }
        else if (is.null(growth_TIC)){
          model_mx <- mxModel(name = "LGCM with TVC, Linear Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1")),
                              mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                              mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L1", records, "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1), name = "Y_mean0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                              cbind(Y_psi01, Y_psi11)), name = "Y_psi0"),
                              mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                     values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                              mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                     values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                              TVC_info, GF_loadings)
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        if (!is.null(growth_TIC)){
          nTICs <- length(growth_TIC)
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
          }
          model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Linear Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1")),
                              mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                              mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L1", records, "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                              cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                              mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                              mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                              cbind(X_psi01, X_psi11)), name = "X_psi0"),
                              mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                       values = starts[[4]][1:2, 1:length(growth_TIC)],
                                       labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                  paste0("beta1TIC", 1:length(growth_TIC))), byrow = T, name = "beta_TIC"),
                              mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC)), name = "beta"),
                              mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                       labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                              mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                              mxAlgebra(Y_alpha0 + beta %*% BL_mean, name = "Y_mean0"),
                              TVC_info, GF_loadings, BETA)
        }
        else if (is.null(growth_TIC)){
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                labels = paste0("beta", p - 1, "TVC"))
          }
          model_mx <- mxModel(name = "LGCM with TVC, Linear Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1")),
                              mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                              mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L1", records, "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                              cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                              mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                              mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                              cbind(X_psi01, X_psi11)), name = "X_psi0"),
                              mxAlgebra(rbind(beta0TVC, beta1TVC), name = "beta"),
                              mxAlgebra(Y_alpha0 + beta %*% X_mueta0, name = "Y_mean0"),
                              TVC_info, GF_loadings, BETA)
        }
      }
    }
    else if (curveFun %in% c("quadratic", "QUAD")){
      latents <- c("eta0Y", "eta1Y", "eta2Y")
      Y_nGF <- length(latents)
      if (decompose == 0){
        if (!is.null(growth_TIC)){
          nTICs <- length(growth_TIC)
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
          }
          model_mx <- mxModel(name = "LGCM with TVC and TICs, Quadratic Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                              mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]],
                                     labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                              mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L1", records, "[1,1]")),
                              mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L2", records, "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                              cbind(Y_psi01, Y_psi11, Y_psi12),
                                              cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                              mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                     values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                              mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                     values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                              mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                       values = starts[[4]][1:3, 1:length(growth_TIC)],
                                       labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                  paste0("beta1TIC", 1:length(growth_TIC)),
                                                  paste0("beta2TIC", 1:length(growth_TIC))), byrow = T, name = "beta"),
                              mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                       labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                              mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                              TVC_info, GF_loadings, BETA)
        }
        else if (is.null(growth_TIC)){
          model_mx <- mxModel(name = "LGCM with TVC, Quadratic Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                              mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]],
                                     labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                              mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L1", records, "[1,1]")),
                              mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L2", records, "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2), name = "Y_mean0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                              cbind(Y_psi01, Y_psi11, Y_psi12),
                                              cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi0"),
                              mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                     values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                              mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                     values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                              TVC_info, GF_loadings)
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        if (!is.null(growth_TIC)){
          nTICs <- length(growth_TIC)
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
          }
          model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Quadratic Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                              mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]],
                                     labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                              mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L1", records, "[1,1]")),
                              mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L2", records, "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                              cbind(Y_psi01, Y_psi11, Y_psi12),
                                              cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                              mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                              mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                              cbind(X_psi01, X_psi11)), name = "X_psi0"),
                              mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                       values = starts[[4]][1:3, 1:length(growth_TIC)],
                                       labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                  paste0("beta1TIC", 1:length(growth_TIC)),
                                                  paste0("beta2TIC", 1:length(growth_TIC))), byrow = T, name = "beta_TIC"),
                              mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC, beta2TVC)), name = "beta"),
                              mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                       labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                              mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                              mxAlgebra(Y_alpha0 + beta %*% BL_mean, name = "Y_mean0"),
                              TVC_info, GF_loadings, BETA)
        }
        else if (is.null(growth_TIC)){
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                labels = paste0("beta", p - 1, "TVC"))
          }
          model_mx <- mxModel(name = "LGCM with TVC, Quadratic Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                              mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]],
                                     labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                              mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L1", records, "[1,1]")),
                              mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L2", records, "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                              cbind(Y_psi01, Y_psi11, Y_psi12),
                                              cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                              mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                              mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                              cbind(X_psi01, X_psi11)), name = "X_psi0"),
                              mxAlgebra(rbind(beta0TVC, beta1TVC, beta2TVC), name = "beta"),
                              mxAlgebra(Y_alpha0 + beta %*% X_mueta0, name = "Y_mean0"),
                              TVC_info, GF_loadings, BETA)
        }
      }
    }
    else if (curveFun %in% c("negative exponential", "EXP")){
      if (intrinsic){
        latents <- c("eta0Y", "eta1Y", "deltag")
        Y_nGF <- length(latents)
        if (decompose == 0){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", c(paste0("TIC", 1:length(growth_TIC)))))
            model_mx <- mxModel(name = "LGCM with TVC and TICs, Negative Exponential Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                       free = TRUE, values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi1g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psigg)), name = "Y_psi_r"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:3, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("betagTIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            model_mx <- mxModel(name = "LGCM with TVC, Negative Exponential Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                       free = TRUE, values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_mean0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi1g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psigg)), name = "Y_psi0"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                TVC_info, GF_loadings)
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Negative Exponential Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                       free = TRUE, values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi1g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psigg)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:3, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("betagTIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_TIC"),
                                mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC, betagTVC)), name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                                mxAlgebra(Y_alpha0 + beta %*% BL_mean, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                  labels = paste0("beta", p - 1, "TVC"))
            }
            BETA[[Y_nGF]] <- mxPath(from = "lx1", to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF],
                                    labels = paste0("beta", "g", "TVC"))
            model_mx <- mxModel(name = "LGCM with TVC, Negative Exponential Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                       free = TRUE, values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi1g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psigg)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxAlgebra(rbind(beta0TVC, beta1TVC, betagTVC), name = "beta"),
                                mxAlgebra(Y_alpha0 + beta %*% X_mueta0, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
        }
      }
      else if (!intrinsic){
        latents <- c("eta0Y", "eta1Y")
        Y_nGF <- length(latents)
        if (decompose == 0){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            model_mx <- mxModel(name = "LGCM with TVC and TICs, Negative Exponential Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:2, 4)], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                                cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:2, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(Y_alpha0[1:2, ] + beta %*% mux, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            model_mx <- mxModel(name = "LGCM with TVC, Negative Exponential Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:2, 4)], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_mean0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                                cbind(Y_psi01, Y_psi11)), name = "Y_psi0"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                TVC_info, GF_loadings)
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Negative Exponential Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:2, 4)], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                                cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:2, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_TIC"),
                                mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC)), name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                                mxAlgebra(Y_alpha0[1:2, ] + beta %*% BL_mean, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                  labels = paste0("beta", p - 1, "TVC"))
            }
            model_mx <- mxModel(name = "LGCM with TVC, Negative Exponential Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:2, 4)], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                                cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxAlgebra(rbind(beta0TVC, beta1TVC), name = "beta"),
                                mxAlgebra(Y_alpha0[1:2, ] + beta %*% X_mueta0, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
        }
      }
    }
    else if (curveFun %in% c("Jenss-Bayley", "JB")){
      if (intrinsic){
        latents <- c("eta0Y", "eta1Y", "eta2Y", "deltag")
        Y_nGF <- length(latents)
        if (decompose == 0){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", c(paste0("TIC", 1:length(growth_TIC)))))
            model_mx <- mxModel(name = "LGCM with TVC and TICs, Jenss Bayley Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                       free = TRUE,values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                  "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi12, Y_psi1g),
                                                cbind(Y_psi02, Y_psi12, Y_psi22, Y_psi2g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psi2g, Y_psigg)), name = "Y_psi_r"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:4, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC)),
                                                    paste0("betagTIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            model_mx <- mxModel(name = "LGCM with TICs, Jenss Bayley Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                       free = TRUE,values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                  "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_mean0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi12, Y_psi1g),
                                                cbind(Y_psi02, Y_psi12, Y_psi22, Y_psi2g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psi2g, Y_psigg)), name = "Y_psi0"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                TVC_info, GF_loadings)
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Jenss Bayley Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                       free = TRUE,values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                  "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi12, Y_psi1g),
                                                cbind(Y_psi02, Y_psi12, Y_psi22, Y_psi2g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psi2g, Y_psigg)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:4, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC)),
                                                    paste0("betagTIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_TIC"),
                                mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC, beta2TVC, betagTVC)), name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                                mxAlgebra(Y_alpha0 + beta %*% BL_mean, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                  labels = paste0("beta", p - 1, "TVC"))
            }
            BETA[[Y_nGF]] <- mxPath(from = "lx1", to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF],
                                    labels = paste0("beta", "g", "TVC"))
            model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Jenss Bayley Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                       free = TRUE,values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                  "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi12, Y_psi1g),
                                                cbind(Y_psi02, Y_psi12, Y_psi22, Y_psi2g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psi2g, Y_psigg)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxAlgebra(rbind(beta0TVC, beta1TVC, beta2TVC, betagTVC), name = "beta"),
                                mxAlgebra(Y_alpha0 + beta %*% X_mueta0, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
        }
      }
      else if (!intrinsic){
        latents <- c("eta0Y", "eta1Y", "eta2Y")
        Y_nGF <- length(latents)
        if (decompose == 0){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Jenss Bayley Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                                cbind(Y_psi01, Y_psi11, Y_psi12),
                                                cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:3, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(Y_alpha0[1:3, ] + beta %*% mux, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            model_mx <- mxModel(name = "LGCM with TVC, Jenss Bayley Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_mean0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                                cbind(Y_psi01, Y_psi11, Y_psi12),
                                                cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi0"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                TVC_info, GF_loadings)
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Jenss Bayley Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                                cbind(Y_psi01, Y_psi11, Y_psi12),
                                                cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:3, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_TIC"),
                                mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC, beta2TVC)), name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                                mxAlgebra(Y_alpha0[1:3, ] + beta %*% BL_mean, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                  labels = paste0("beta", p - 1, "TVC"))
            }
            model_mx <- mxModel(name = "LGCM with TVC, Jenss Bayley Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                                mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                                cbind(Y_psi01, Y_psi11, Y_psi12),
                                                cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxAlgebra(rbind(beta0TVC, beta1TVC, beta2TVC), name = "beta"),
                                mxAlgebra(Y_alpha0[1:3, ] + beta %*% X_mueta0, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
        }
      }
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      if (intrinsic){
        latents <- c("eta0sY", "eta1sY", "eta2sY", "deltag")
        Y_nGF <- length(latents)
        if (decompose == 0){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", c(paste0("TIC", 1:length(growth_TIC)))))
            model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Bilinear Spline Functional Form (random knot)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_knot",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]],
                                       labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                  "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs")),
                                mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot, 0),
                                                cbind(0, 1, -1, 0),
                                                cbind(0, 1, 1, 0),
                                                cbind(0, 0, 0, 1)), name = "func"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot, 0),
                                                cbind(0, 1, -1, 0),
                                                cbind(0, 1, 1, 0),
                                                cbind(0, 0, 0, 1)), name = "grad"),
                                mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_alpha_s"),
                                mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s, Y_psi0gs),
                                                cbind(Y_psi01s, Y_psi11s, Y_psi12s, Y_psi1gs),
                                                cbind(Y_psi02s, Y_psi12s, Y_psi22s, Y_psi2gs),
                                                cbind(Y_psi0gs, Y_psi1gs, Y_psi2gs, Y_psiggs)), name = "Y_psi_s"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:4, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC)),
                                                    paste0("betagTIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_s"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(func[1:3, 1:3] %*% Y_alpha_s, Y_knot), name = "Y_alpha0"),
                                mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi_r"),
                                mxAlgebra(grad %*% beta_s, name = "beta"),
                                mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            model_mx <- mxModel(name = "LGCM with TVC, Bilinear Spline Functional Form (random knot)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_knot",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]],
                                       labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                  "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs")),
                                mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot, 0),
                                                cbind(0, 1, -1, 0),
                                                cbind(0, 1, 1, 0),
                                                cbind(0, 0, 0, 1)), name = "func"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot, 0),
                                                cbind(0, 1, -1, 0),
                                                cbind(0, 1, 1, 0),
                                                cbind(0, 0, 0, 1)), name = "grad"),
                                mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_mean_s"),
                                mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s, Y_psi0gs),
                                                cbind(Y_psi01s, Y_psi11s, Y_psi12s, Y_psi1gs),
                                                cbind(Y_psi02s, Y_psi12s, Y_psi22s, Y_psi2gs),
                                                cbind(Y_psi0gs, Y_psi1gs, Y_psi2gs, Y_psiggs)), name = "Y_psi_s"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxAlgebra(rbind(func[1:3, 1:3] %*% Y_mean_s, Y_knot), name = "Y_mean0"),
                                mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi0"),
                                TVC_info, GF_loadings)
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Bilinear Spline Functional Form (random knot)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_knot",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]],
                                       labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                  "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs")),
                                mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot, 0),
                                                cbind(0, 1, -1, 0),
                                                cbind(0, 1, 1, 0),
                                                cbind(0, 0, 0, 1)), name = "func"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot, 0),
                                                cbind(0, 1, -1, 0),
                                                cbind(0, 1, 1, 0),
                                                cbind(0, 0, 0, 1)), name = "grad"),
                                mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_alpha_s"),
                                mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s, Y_psi0gs),
                                                cbind(Y_psi01s, Y_psi11s, Y_psi12s, Y_psi1gs),
                                                cbind(Y_psi02s, Y_psi12s, Y_psi22s, Y_psi2gs),
                                                cbind(Y_psi0gs, Y_psi1gs, Y_psi2gs, Y_psiggs)), name = "Y_psi_s"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:4, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC)),
                                                    paste0("betagTIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_TIC"),
                                mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC, beta2TVC, betagTVC)), name = "beta_s"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                                mxAlgebra(rbind(func[1:3, 1:3] %*% Y_alpha_s, Y_knot), name = "Y_alpha0"),
                                mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi_r"),
                                mxAlgebra(grad %*% beta_s, name = "beta"),
                                mxAlgebra(Y_alpha0 + beta %*% BL_mean, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, "TVC"))
            }
            BETA[[Y_nGF]] <- mxPath(from = "lx1", to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", "TVC"))
            model_mx <- mxModel(name = "LGCM with TVC, Bilinear Spline Functional Form (random knot)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_knot",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]],
                                       labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                  "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs")),
                                mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot, 0),
                                                cbind(0, 1, -1, 0),
                                                cbind(0, 1, 1, 0),
                                                cbind(0, 0, 0, 1)), name = "func"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot, 0),
                                                cbind(0, 1, -1, 0),
                                                cbind(0, 1, 1, 0),
                                                cbind(0, 0, 0, 1)), name = "grad"),
                                mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_alpha_s"),
                                mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s, Y_psi0gs),
                                                cbind(Y_psi01s, Y_psi11s, Y_psi12s, Y_psi1gs),
                                                cbind(Y_psi02s, Y_psi12s, Y_psi22s, Y_psi2gs),
                                                cbind(Y_psi0gs, Y_psi1gs, Y_psi2gs, Y_psiggs)), name = "Y_psi_s"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxAlgebra(rbind(beta0TVC, beta1TVC, beta2TVC, betagTVC), name = "beta_s"),
                                mxAlgebra(rbind(func[1:3, 1:3] %*% Y_alpha_s, Y_knot), name = "Y_alpha0"),
                                mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi_r"),
                                mxAlgebra(grad %*% beta_s, name = "beta"),
                                mxAlgebra(Y_alpha0 + beta %*% X_mueta0, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
        }
      }
      else if (!intrinsic){
        latents <- c("eta0sY", "eta1sY", "eta2sY")
        Y_nGF <- length(latents)
        if (decompose == 0){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            model_mx <- mxModel(name = "LGCM with TVC and TICs, Bilinear Spline Functional Form (fixed knot)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_knot",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi11s", "Y_psi12s", "Y_psi22s")),
                                mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot),
                                                cbind(0, 1, -1),
                                                cbind(0, 1, 1)), name = "func"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot),
                                                cbind(0, 1, -1),
                                                cbind(0, 1, 1)), name = "grad"),
                                mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_alpha_s"),
                                mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s),
                                                cbind(Y_psi01s, Y_psi11s, Y_psi12s),
                                                cbind(Y_psi02s, Y_psi12s, Y_psi22s)), name = "Y_psi_s"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:3, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_s"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(func %*% Y_alpha_s, Y_knot), name = "Y_alpha0"),
                                mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi_r"),
                                mxAlgebra(grad %*% beta_s, name = "beta"),
                                mxAlgebra(Y_alpha0[1:3, 1] + beta %*% mux, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            model_mx <- mxModel(name = "LGCM with TVC, Bilinear Spline Functional Form (fixed knot)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_knot",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi11s", "Y_psi12s", "Y_psi22s")),
                                mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot),
                                                cbind(0, 1, -1),
                                                cbind(0, 1, 1)), name = "func"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot),
                                                cbind(0, 1, -1),
                                                cbind(0, 1, 1)), name = "grad"),
                                mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_mean_s"),
                                mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s),
                                                cbind(Y_psi01s, Y_psi11s, Y_psi12s),
                                                cbind(Y_psi02s, Y_psi12s, Y_psi22s)), name = "Y_psi_s"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxAlgebra(rbind(func %*% Y_mean_s, Y_knot), name = "Y_mean0"),
                                mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi0"),
                                TVC_info, GF_loadings)
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            model_mx <- mxModel(name = "LGCM with decomposed TVC and TICs, Bilinear Spline Functional Form (fixed knot)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_knot",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi11s", "Y_psi12s", "Y_psi22s")),
                                mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot),
                                                cbind(0, 1, -1),
                                                cbind(0, 1, 1)), name = "func"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot),
                                                cbind(0, 1, -1),
                                                cbind(0, 1, 1)), name = "grad"),
                                mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_alpha_s"),
                                mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s),
                                                cbind(Y_psi01s, Y_psi11s, Y_psi12s),
                                                cbind(Y_psi02s, Y_psi12s, Y_psi22s)), name = "Y_psi_s"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:3, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_TIC"),
                                mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC, beta2TVC)), name = "beta_s"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                                mxAlgebra(rbind(func %*% Y_alpha_s, Y_knot), name = "Y_alpha0"),
                                mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi_r"),
                                mxAlgebra(grad %*% beta_s, name = "beta"),
                                mxAlgebra(Y_alpha0[1:3, 1] + beta %*% BL_mean, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                  labels = paste0("beta", p - 1, "TVC"))
            }
            model_mx <- mxModel(name = "LGCM with TVC, Bilinear Spline Functional Form (fixed knot)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_knot",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi11s", "Y_psi12s", "Y_psi22s")),
                                mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records, "[1,1]")),
                                mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records, "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot),
                                                cbind(0, 1, -1),
                                                cbind(0, 1, 1)), name = "func"),
                                mxAlgebra(rbind(cbind(1, -Y_knot, Y_knot),
                                                cbind(0, 1, -1),
                                                cbind(0, 1, 1)), name = "grad"),
                                mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_alpha_s"),
                                mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s),
                                                cbind(Y_psi01s, Y_psi11s, Y_psi12s),
                                                cbind(Y_psi02s, Y_psi12s, Y_psi22s)), name = "Y_psi_s"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxAlgebra(rbind(beta0TVC, beta1TVC, beta2TVC), name = "beta_s"),
                                mxAlgebra(rbind(func %*% Y_alpha_s, Y_knot), name = "Y_alpha0"),
                                mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_psi_r"),
                                mxAlgebra(grad %*% beta_s, name = "beta"),
                                mxAlgebra(Y_alpha0[1:3, 1] + beta %*% X_mueta0, name = "Y_mean0"),
                                TVC_info, GF_loadings, BETA)
          }
        }
      }
    }
  }
  if (y_model == "LCSM"){
    ## Define paths from latent true scores to observed scores
    Y_PATH_L <- mxPath(from = paste0("ly", records), to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1)
    ## Define paths from latent instantaneous rate of change at each measurement to true scores
    Y_PATH_SLP <- mxPath(from = paste0("dy", records[-1]), to = paste0("ly", records[-1]), arrows = 1, free = FALSE,
                         values = 0, labels = paste0("lag", records[-1], "[1,1]"))
    #### Define autoregressive paths
    Y_PATH_AUTO <- mxPath(from = paste0("ly", records[-length(records)]), to = paste0("ly", records[-1]), arrows = 1,
                          free = FALSE, values = 1)
    # Obtain additional parameters (derived) for the specified functional form
    AddPara <- getUNI.addpara(dat = dat, curveFun = curveFun, intrinsic = intrinsic, t_var = t_var,
                              records = records, growth_TIC = growth_TIC, decompose = decompose,
                              starts = starts)
    if (curveFun %in% c("nonparametric", "NonP")){
      latents <- c("eta0Y", "eta1Y", paste0("dy", records[-1]), paste0("ly", records))
      Y_nGF <- length(latents) - (length(records) * 2 - 1)
      if (decompose == 0){
        if (!is.null(growth_TIC)){
          nTICs <- length(growth_TIC)
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
          }
          model_mx <- mxModel(name = "LCSM with TVC and TICs, Nonparametric Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1")),
                              mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                              mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                     free = c(F, rep(T, length(records) - 2)), values = c(1, starts[[1]][[4]][-1]),
                                     labels = paste0("Y_rel_rate", 1:(length(records) - 1))),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                              cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                              mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                     values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                              mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                     values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                              mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                       values = starts[[4]][1:2, 1:length(growth_TIC)],
                                       labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                  paste0("beta1TIC", 1:length(growth_TIC))), byrow = T, name = "beta"),
                              mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                       labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                              mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                              Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
        }
        else if (is.null(growth_TIC)){
          model_mx <- mxModel(name = "LCSM with TVC, Nonparametric Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1")),
                              mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                              mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                     free = c(F, rep(T, length(records) - 2)), values = c(1, starts[[1]][[4]][-1]),
                                     labels = paste0("Y_rel_rate", 1:(length(records) - 1))),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1), name = "Y_mean0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                              cbind(Y_psi01, Y_psi11)), name = "Y_psi0"),
                              mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                     values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                              mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                     values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                              Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara)
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        if (!is.null(growth_TIC)){
          nTICs <- length(growth_TIC)
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
          }
          model_mx <- mxModel(name = "LCSM with decomposed TVC and TICs, Nonparametric Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1")),
                              mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                              mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                     free = c(F, rep(T, length(records) - 2)), values = c(1, starts[[1]][[4]][-1]),
                                     labels = paste0("Y_rel_rate", 1:(length(records) - 1))),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                              cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                              mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                              mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                              cbind(X_psi01, X_psi11)), name = "X_psi0"),
                              mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                       values = starts[[4]][1:2, 1:length(growth_TIC)],
                                       labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                  paste0("beta1TIC", 1:length(growth_TIC))), byrow = T, name = "beta_TIC"),
                              mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC)), name = "beta"),
                              mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                       labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                              mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                              mxAlgebra(Y_alpha0 + beta %*% BL_mean, name = "Y_mean0"),
                              Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
        }
        else if (is.null(growth_TIC)){
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                labels = paste0("beta", p - 1, "TVC"))
          }
          model_mx <- mxModel(name = "LCSM with TVC, Nonparametric Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1")),
                              mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                              mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                     free = c(F, rep(T, length(records) - 2)), values = c(1, starts[[1]][[4]][-1]),
                                     labels = paste0("Y_rel_rate", 1:(length(records) - 1))),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                              cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                              mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                              mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                              cbind(X_psi01, X_psi11)), name = "X_psi0"),
                              mxAlgebra(rbind(beta0TVC, beta1TVC), name = "beta"),
                              mxAlgebra(Y_alpha0 + beta %*% X_mueta0, name = "Y_mean0"),
                              Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
        }
      }
    }
    else if (curveFun %in% c("quadratic", "QUAD")){
      latents <- c("eta0Y", "eta1Y", "eta2Y", paste0("dy", records[-1]), paste0("ly", records))
      Y_nGF <- length(latents) - (length(records) * 2 - 1)
      if (decompose == 0){
        if (!is.null(growth_TIC)){
          nTICs <- length(growth_TIC)
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
          }
          model_mx <- mxModel(name = "LCSM with TVC and TICs, Quadratic Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                              mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]],
                                     labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                              mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L2", records[-1], "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                              cbind(Y_psi01, Y_psi11, Y_psi12),
                                              cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                              mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                     values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                              mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                     values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                              mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                       values = starts[[4]][1:3, 1:length(growth_TIC)],
                                       labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                  paste0("beta1TIC", 1:length(growth_TIC)),
                                                  paste0("beta2TIC", 1:length(growth_TIC))), byrow = T, name = "beta"),
                              mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                       labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                              mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                              Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
        }
        else if (is.null(growth_TIC)){
          model_mx <- mxModel(name = "LCSM with TVC, Quadratic Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                              mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]],
                                     labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                              mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L2", records[-1], "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2), name = "Y_mean0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                              cbind(Y_psi01, Y_psi11, Y_psi12),
                                              cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi0"),
                              mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                     values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                              mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                     values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                              Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara)
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        if (!is.null(growth_TIC)){
          nTICs <- length(growth_TIC)
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
          }
          model_mx <- mxModel(name = "LCSM with decomposed TVC and TICs, Quadratic Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                              mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]],
                                     labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                              mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L2", records[-1], "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                              cbind(Y_psi01, Y_psi11, Y_psi12),
                                              cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                              mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                              mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                              cbind(X_psi01, X_psi11)), name = "X_psi0"),
                              mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                       values = starts[[4]][1:3, 1:length(growth_TIC)],
                                       labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                  paste0("beta1TIC", 1:length(growth_TIC)),
                                                  paste0("beta2TIC", 1:length(growth_TIC))), byrow = T, name = "beta_TIC"),
                              mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC, beta2TVC)), name = "beta"),
                              mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                       labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                              mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                              mxAlgebra(Y_alpha0 + beta %*% BL_mean, name = "Y_mean0"),
                              Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
        }
        else if (is.null(growth_TIC)){
          for (p in 1:Y_nGF){
            BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                labels = paste0("beta", p - 1, "TVC"))
          }
          model_mx <- mxModel(name = "LCSM with TVC, Quadratic Functional Form",
                              type = "RAM", mxData(observed = dat, type = "raw"),
                              manifestVars = manifests, latentVars = latents,
                              mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]],
                                     labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                              mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                     free = TRUE, values = starts[[1]][[2]],
                                     labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                              mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                     labels = paste0("L2", records[-1], "[1,1]")),
                              mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                     free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                              mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2), name = "Y_alpha0"),
                              mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                              cbind(Y_psi01, Y_psi11, Y_psi12),
                                              cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                              mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                              mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                              cbind(X_psi01, X_psi11)), name = "X_psi0"),
                              mxAlgebra(rbind(beta0TVC, beta1TVC, beta2TVC), name = "beta"),
                              mxAlgebra(Y_alpha0 + beta %*% X_mueta0, name = "Y_mean0"),
                              Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
        }
      }
    }
    else if (curveFun %in% c("negative exponential", "EXP")){
      if (intrinsic){
        latents <- c("eta0Y", "eta1Y", "deltag", paste0("dy", records[-1]), paste0("ly", records))
        Y_nGF <- length(latents) - (length(records) * 2 - 1)
        if (decompose == 0){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", c(paste0("TIC", 1:length(growth_TIC)))))
            model_mx <- mxModel(name = "LCSM with TVC and TICs, Negative Exponential Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                       free = TRUE, values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records[-1], "[1,1]")),
                                mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi1g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psigg)), name = "Y_psi_r"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:3, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("betagTIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
          else if (is.null(growth_TIC)){
            model_mx <- mxModel(name = "LCSM with TVC, Negative Exponential Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                       free = TRUE, values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records[-1], "[1,1]")),
                                mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_mean0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi1g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psigg)), name = "Y_psi0"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara)
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            model_mx <- mxModel(name = "LCSM with decomposed TVC and TICs, Negative Exponential Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                       free = TRUE, values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records[-1], "[1,1]")),
                                mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi1g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psigg)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:3, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("betagTIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_TIC"),
                                mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC, betagTVC)), name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                                mxAlgebra(Y_alpha0 + beta %*% BL_mean, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                  labels = paste0("beta", p - 1, "TVC"))
            }
            BETA[[Y_nGF]] <- mxPath(from = "lx1", to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF],
                                    labels = paste0("beta", "g", "TVC"))
            model_mx <- mxModel(name = "LCSM with TVC, Negative Exponential Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                       free = TRUE, values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records[-1], "[1,1]")),
                                mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi1g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psigg)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxAlgebra(rbind(beta0TVC, beta1TVC, betagTVC), name = "beta"),
                                mxAlgebra(Y_alpha0 + beta %*% X_mueta0, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
        }
      }
      else if (!intrinsic){
        latents <- c("eta0Y", "eta1Y", paste0("dy", records[-1]), paste0("ly", records))
        Y_nGF <- length(latents) - (length(records) * 2 - 1)
        if (decompose == 0){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            model_mx <- mxModel(name = "LCSM with TVC and TICs, Negative Exponential Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:2, 4)], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                                cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:2, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(Y_alpha0[1:2, ] + beta %*% mux, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
          else if (is.null(growth_TIC)){
            model_mx <- mxModel(name = "LCSM with TVC, Negative Exponential Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:2, 4)], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_mean0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                                cbind(Y_psi01, Y_psi11)), name = "Y_psi0"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara)
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            model_mx <- mxModel(name = "LCSM with decomposed TVC and TICs, Negative Exponential Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:2, 4)], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                                cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:2, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_TIC"),
                                mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC)), name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                                mxAlgebra(Y_alpha0[1:2, ] + beta %*% BL_mean, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                  labels = paste0("beta", p - 1, "TVC"))
            }
            model_mx <- mxModel(name = "LCSM with TVC, Negative Exponential Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:2],
                                       labels = c("Y_mueta0", "Y_mueta1")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][3], labels = "Y_slp_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:2, 4)], labels = c("Y_psi00", "Y_psi01", "Y_psi11")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L1", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_slp_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01),
                                                cbind(Y_psi01, Y_psi11)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxAlgebra(rbind(beta0TVC, beta1TVC), name = "beta"),
                                mxAlgebra(Y_alpha0[1:2, ] + beta %*% X_mueta0, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
        }
      }
    }
    else if (curveFun %in% c("Jenss-Bayley", "JB")){
      if (intrinsic){
        latents <- c("eta0Y", "eta1Y", "eta2Y", "deltag", paste0("dy", records[-1]), paste0("ly", records))
        Y_nGF <- length(latents) - (length(records) * 2 - 1)
        if (decompose == 0){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", c(paste0("TIC", 1:length(growth_TIC)))))
            model_mx <- mxModel(name = "LCSM with TVC and TICs, Jenss Bayley Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                       free = TRUE,values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                  "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi12, Y_psi1g),
                                                cbind(Y_psi02, Y_psi12, Y_psi22, Y_psi2g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psi2g, Y_psigg)), name = "Y_psi_r"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:4, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC)),
                                                    paste0("betagTIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(Y_alpha0 + beta %*% mux, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
          else if (is.null(growth_TIC)){
            model_mx <- mxModel(name = "LCSM with decomposed TVC and TICs, Jenss Bayley Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                       free = TRUE,values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                  "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_mean0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi12, Y_psi1g),
                                                cbind(Y_psi02, Y_psi12, Y_psi22, Y_psi2g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psi2g, Y_psigg)), name = "Y_psi0"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara)
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF, ],
                                    labels = paste0("beta", "g", c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            model_mx <- mxModel(name = "LCSM with decomposed TVC and TICs, Jenss Bayley Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                       free = TRUE,values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                  "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi12, Y_psi1g),
                                                cbind(Y_psi02, Y_psi12, Y_psi22, Y_psi2g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psi2g, Y_psigg)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:4, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC)),
                                                    paste0("betagTIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_TIC"),
                                mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC, beta2TVC, betagTVC)), name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                                mxAlgebra(Y_alpha0 + beta %*% BL_mean, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:(Y_nGF - 1)){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                  labels = paste0("beta", p - 1, "TVC"))
            }
            BETA[[Y_nGF]] <- mxPath(from = "lx1", to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[4]][Y_nGF],
                                    labels = paste0("beta", "g", "TVC"))
            model_mx <- mxModel(name = "LCSM with decomposed TVC and TICs, Jenss Bayley Functional Form (random ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                       free = TRUE,values = starts[[1]][[2]],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                  "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L3", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02, Y_psi0g),
                                                cbind(Y_psi01, Y_psi11, Y_psi12, Y_psi1g),
                                                cbind(Y_psi02, Y_psi12, Y_psi22, Y_psi2g),
                                                cbind(Y_psi0g, Y_psi1g, Y_psi2g, Y_psigg)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxAlgebra(rbind(beta0TVC, beta1TVC, beta2TVC, betagTVC), name = "beta"),
                                mxAlgebra(Y_alpha0 + beta %*% X_mueta0, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
        }
      }
      else if (!intrinsic){
        latents <- c("eta0Y", "eta1Y", "eta2Y", paste0("dy", records[-1]), paste0("ly", records))
        Y_nGF <- length(latents) - (length(records) * 2 - 1)
        if (decompose == 0){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            model_mx <- mxModel(name = "LCSM with decomposed TVC and TICs, Jenss Bayley Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                                cbind(Y_psi01, Y_psi11, Y_psi12),
                                                cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:3, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(Y_alpha0[1:3, ] + beta %*% mux, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
          else if (is.null(growth_TIC)){
            model_mx <- mxModel(name = "LCSM with TVC, Jenss Bayley Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_mean0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                                cbind(Y_psi01, Y_psi11, Y_psi12),
                                                cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi0"),
                                mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                       values = starts[[2]][[1]], labels = paste0("TVC_m", records)),
                                mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                       values = starts[[2]][[2]], labels = paste0("TVC_v", records)),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara)
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p, ],
                                  labels = paste0("beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            model_mx <- mxModel(name = "LCSM with decomposed TVC and TICs, Jenss Bayley Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                                cbind(Y_psi01, Y_psi11, Y_psi12),
                                                cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                         values = starts[[4]][1:3, 1:length(growth_TIC)],
                                         labels = c(paste0("beta0TIC", 1:length(growth_TIC)),
                                                    paste0("beta1TIC", 1:length(growth_TIC)),
                                                    paste0("beta2TIC", 1:length(growth_TIC))),
                                         byrow = T, name = "beta_TIC"),
                                mxAlgebra(cbind(beta_TIC, rbind(beta0TVC, beta1TVC, beta2TVC)), name = "beta"),
                                mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[3]][[1]][1:length(growth_TIC)],
                                         labels = c(paste0("mux", 1:length(growth_TIC))), byrow = F, name = "mux"),
                                mxAlgebra(rbind(mux, X_mueta0), name = "BL_mean"),
                                mxAlgebra(Y_alpha0[1:3, ] + beta %*% BL_mean, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[4]][p],
                                  labels = paste0("beta", p - 1, "TVC"))
            }
            model_mx <- mxModel(name = "LCSM with TVC, Jenss Bayley Functional Form (fixed ratio)",
                                type = "RAM", mxData(observed = dat, type = "raw"),
                                manifestVars = manifests, latentVars = latents,
                                mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[1]][[1]][1:3],
                                       labels = c("Y_mueta0", "Y_mueta1", "Y_mueta2")),
                                mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[1]][4], labels = "Y_acc_ratio",
                                         name = "Y_mug"),
                                mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                       values = starts[[1]][[2]][c(1:3, 5:6, 8)],
                                       labels = c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22")),
                                mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                       labels = paste0("L2", records[-1], "[1,1]")),
                                mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                       free = TRUE, values = starts[[1]][[3]], labels = "Y_residuals"),
                                mxAlgebra(rbind(Y_mueta0, Y_mueta1, Y_mueta2, Y_acc_ratio), name = "Y_alpha0"),
                                mxAlgebra(rbind(cbind(Y_psi00, Y_psi01, Y_psi02),
                                                cbind(Y_psi01, Y_psi11, Y_psi12),
                                                cbind(Y_psi02, Y_psi12, Y_psi22)), name = "Y_psi_r"),
                                mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean0"),
                                mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                                cbind(X_psi01, X_psi11)), name = "X_psi0"),
                                mxAlgebra(rbind(beta0TVC, beta1TVC, beta2TVC), name = "beta"),
                                mxAlgebra(Y_alpha0[1:3, ] + beta %*% X_mueta0, name = "Y_mean0"),
                                Y_PATH_L, Y_PATH_SLP, Y_PATH_AUTO, TVC_info, GF_loadings, AddPara, BETA)
          }
        }
      }
    }
  }
  return(model_mx)
}



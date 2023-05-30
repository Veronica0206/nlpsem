#' @title Define a Latent Growth Curve Model or Latent Change Score Model with a Time-varying Covariate as Class-specific Models
#' (Submodels) for a Longitudinal Mixture Model.
#'
#' @description This function defines a latent growth curve model or latent change score model with time-varying covariate as class-
#' specific models (submodels) for a longitudinal mixture model.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions for each longitudinal process, and time-invariant covariates (TICs) if any.
#' It takes the value passed from \code{getMIX()}.
#' @param nClass An integer specifying the number of latent classes for the mixture model. It takes the value passed from \code{getMIX()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' It takes the value passed from \code{getMIX()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getMIX()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the
#' value passed from \code{getMIX()}.
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. It takes the
#' value passed from \code{getMIX()}.
#' @param records A numeric vector specifying the indices of the observed study waves. It takes the value passed from
#' \code{getMIX()}.
#' @param y_model A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' It takes the value passed from \code{getMIX()}.
#' @param TVC A string specifying the prefix of the column names corresponding to the time-varying covariate at each study wave.
#' It takes the value passed from \code{getMIX()}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). It takes the
#' value passed from \code{getMIX()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) that account for the
#' variability of growth factors, if any. It takes the value passed from \code{getMIX()}.
#' @param starts A list of initial values for the parameters, either takes the value passed from \code{getMIX()}
#' or derived by the helper function \code{getMIX.initial()}.
#'
#' @return A list of manifest and latent variables and paths for an mxModel object.
#'
#' @importFrom OpenMx mxPath mxModel mxAlgebraFromString mxMatrix mxFitFunctionML
#'
getsub.TVC_l <- function(dat, nClass, t_var, records, y_var, curveFun, intrinsic, y_model, TVC, decompose,
                         growth_TIC, starts){
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
  TVC_info <- getMIX_TVC.info(nClass = nClass, y_var = y_var, records = records, growth_TIC = growth_TIC, TVC = TVC,
                              decompose = decompose, starts = starts)
  # Obtain factor loadigs for the specified functional form
  GF_loadings <- getMIX_UNI.loadings(nClass = nClass, y_model = y_model, t_var = t_var, y_var = y_var,
                                     curveFun = curveFun, intrinsic = intrinsic, records = records)
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
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                  labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 0,
                                              labels = paste0("c", k, "L1", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                           name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                       mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                values = starts[[k]][[4]][1:2, 1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta1TIC", 1:length(growth_TIC))),
                                                byrow = T, name = paste0("c", k, "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                byrow = F, name = paste0("c", k, "mux")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
          }
          else if (is.null(growth_TIC)){
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 0,
                                              labels = paste0("c", k, "L1", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                           name = paste0("c", k, "Y_psi0")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                       TVC_info[[k]], GF_loadings[[k]], mxFitFunctionML(vector = T))
          }
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                  labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 0,
                                              labels = paste0("c", k, "L1", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                           name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                           name = paste0("c", k, "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                  "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                           name = paste0("c", k, "X_psi0")),
                                       mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                values = starts[[k]][[4]][1:2, 1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta1TIC", 1:length(growth_TIC))),
                                                byrow = T, name = paste0("c", k, "beta_TIC")),
                                       mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                  "rbind(c", k, "beta0TVC, c", k, "beta1TVC))"),
                                                           name = paste0("c", k, "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                byrow = F, name = paste0("c", k, "mux")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                           name = paste0("c", k, "BL_mean")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "BL_mean"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                  labels = paste0("c", k, "beta", p - 1, "TVC"))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records),
                                              arrows = 1, free = FALSE, values = 0,
                                              labels = paste0("c", k, "L1", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                           name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                           name = paste0("c", k, "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                  "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                           name = paste0("c", k, "X_psi0")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC)"),
                                                           name = paste0("c", k, "beta")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "X_mueta0"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
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
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                  labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                       mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                              arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                              labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                  "Y_mueta2)"), name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                  "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                       mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta2TIC", 1:length(growth_TIC))),
                                                byrow = T, name = paste0("c", k, "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                byrow = F, name = paste0("c", k, "mux")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
          }
          else if (is.null(growth_TIC)){
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                       mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                              arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                              labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                  "Y_mueta2)"), name = paste0("c", k, "Y_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                  "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                           name = paste0("c", k, "Y_psi0")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                       TVC_info[[k]], GF_loadings[[k]], mxFitFunctionML(vector = T))
          }
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                  labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                       mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                              arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                              labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                  "Y_mueta2)"), name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                  "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                           name = paste0("c", k, "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                  "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                           name = paste0("c", k, "X_psi0")),
                                       mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta2TIC", 1:length(growth_TIC))),
                                                byrow = T, name = paste0("c", k, "beta_TIC")),
                                       mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                  "rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                  "beta2TVC))"),
                                                           name = paste0("c", k, "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                byrow = F, name = paste0("c", k, "mux")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                           name = paste0("c", k, "BL_mean")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "BL_mean"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                  labels = paste0("c", k, "beta", p - 1, "TVC"))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 1),
                                       mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                       mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                              values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                              arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                              labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                  "Y_mueta2)"), name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                  "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                           name = paste0("c", k, "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                  "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                           name = paste0("c", k, "X_psi0")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                  "beta2TVC)"), name = paste0("c", k, "beta")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "X_mueta0"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
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
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
              }
              BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", c(paste0("TIC", 1:length(growth_TIC)))))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"),
                                                  name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "betagTIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"),
                                                  name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], mxFitFunctionML(vector = T))
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (decompose == 1){
              Y_nGF <- length(latents) - (length(records) * 2 - 1) - 2
            }
            else if (I(decompose == 2 | decompose == 3)) {
              Y_nGF <- length(latents) - (length(records) * 2 - 1) - (length(records) - 1) - 2
            }
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              }
              BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"),
                                                  name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "betagTIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                    "rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "betagTVC))"),
                                                             name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                             name = paste0("c", k, "BL_mean")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "BL_mean"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                    labels = paste0("c", k, "beta", p - 1, "TVC"))
              }
              BETA[[Y_nGF]] <- mxPath(from = "lx1", to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF],
                                      labels = paste0("c", k, "beta", "g", "TVC"))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"),
                                                  name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "betagTVC)"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "X_mueta0"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
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
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:2, 4)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:2, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:2, ] + c", k, "beta %*% c", k, "mux"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:2, 4)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], mxFitFunctionML(vector = T))
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:2, 4)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:2, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                    "rbind(c", k, "beta0TVC, c", k, "beta1TVC))"),
                                                             name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                             name = paste0("c", k, "BL_mean")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:2, ] + c", k, "beta %*% c", k, "BL_mean"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                    labels = paste0("c", k, "beta", p - 1, "TVC"))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:2, 4)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC)"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:2, ] + c", k, "beta %*% c", k, "X_mueta0"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
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
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
              }
              BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", c(paste0("TIC", 1:length(growth_TIC)))))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22, c", k, "Y_psi2g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psi2g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:4, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "betagTIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22, c", k, "Y_psi2g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psi2g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], mxFitFunctionML(vector = T))
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              }
              BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22, c", k, "Y_psi2g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psi2g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:4, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "betagTIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                    "rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC, c", k, "betagTVC))"),
                                                             name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                             name = paste0("c", k, "BL_mean")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "BL_mean"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                    labels = paste0("c", k, "beta", p - 1, "TVC"))
              }
              BETA[[Y_nGF]] <- mxPath(from = "lx1", to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF],
                                      labels = paste0("c", k, "beta", "g", "TVC"))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22, c", k, "Y_psi2g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psi2g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC, c", k, "betagTVC)"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "X_mueta0"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
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
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:3, ] + c", k, "beta %*% c", k, "mux"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], mxFitFunctionML(vector = T))
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                    "rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC))"),
                                                             name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                             name = paste0("c", k, "BL_mean")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:3, ] + c", k, "beta %*% c", k, "BL_mean"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                    labels = paste0("c", k, "beta", p - 1, "TVC"))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 1),
                                         mxPath(from = "eta1Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2Y", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC)"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:3, ] + c", k, "beta %*% c", k, "X_mueta0"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
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
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
              }
              BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", c(paste0("TIC", 1:length(growth_TIC)))))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                          "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0("c", k, "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0("c", k, "grad")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                             name = paste0("c", k, "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s, c", k, "Y_psi0gs), ",
                                                                    "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s, c", k, "Y_psi1gs), ",
                                                                    "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s, c", k, "Y_psi2gs), ",
                                                                    "cbind(c", k, "Y_psi0gs, c", k, "Y_psi1gs, c", k, "Y_psi2gs, c", k, "Y_psiggs))"),
                                                             name = paste0("c", k, "Y_psi_s")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:4, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "betagTIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_s")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "func[1:3, 1:3] %*% c", k,
                                                                    "Y_alpha_s, c", k, "Y_knot)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "beta_s"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                          "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0("c", k, "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0("c", k, "grad")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                             name = paste0("c", k, "Y_mean_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s, c", k, "Y_psi0gs), ",
                                                                    "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s, c", k, "Y_psi1gs), ",
                                                                    "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s, c", k, "Y_psi2gs), ",
                                                                    "cbind(c", k, "Y_psi0gs, c", k, "Y_psi1gs, c", k, "Y_psi2gs, c", k, "Y_psiggs))"),
                                                             name = paste0("c", k, "Y_psi_s")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxAlgebraFromString(paste0("rbind(c", k, "func[1:3, 1:3] %*% c", k,
                                                                    "Y_mean_s, c", k, "Y_knot)"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         TVC_info[[k]], GF_loadings[[k]], mxFitFunctionML(vector = T))
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              }
              BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                          "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0("c", k, "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0("c", k, "grad")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                             name = paste0("c", k, "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s, c", k, "Y_psi0gs), ",
                                                                    "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s, c", k, "Y_psi1gs), ",
                                                                    "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s, c", k, "Y_psi2gs), ",
                                                                    "cbind(c", k, "Y_psi0gs, c", k, "Y_psi1gs, c", k, "Y_psi2gs, c", k, "Y_psiggs))"),
                                                             name = paste0("c", k, "Y_psi_s")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:4, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "betagTIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                    "rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC, c", k, "betagTVC))"),
                                                             name = paste0("c", k, "beta_s")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                             name = paste0("c", k, "BL_mean")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "func[1:3, 1:3] %*% c", k,
                                                                    "Y_alpha_s, c", k, "Y_knot)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "beta_s"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "BL_mean"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, "TVC"))
              }
              BETA[[Y_nGF]] <- mxPath(from = "lx1", to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", "TVC"))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                          "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L3", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0("c", k, "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                    "cbind(0, 1, -1, 0),",
                                                                    "cbind(0, 1, 1, 0),",
                                                                    "cbind(0, 0, 0, 1))"), name = paste0("c", k, "grad")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                             name = paste0("c", k, "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s, c", k, "Y_psi0gs), ",
                                                                    "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s, c", k, "Y_psi1gs), ",
                                                                    "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s, c", k, "Y_psi2gs), ",
                                                                    "cbind(c", k, "Y_psi0gs, c", k, "Y_psi1gs, c", k, "Y_psi2gs, c", k, "Y_psiggs))"),
                                                             name = paste0("c", k, "Y_psi_s")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC, c", k, "betagTVC)"),
                                                             name = paste0("c", k, "beta_s")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "func[1:3, 1:3] %*% c", k,
                                                                    "Y_alpha_s, c", k, "Y_knot)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "beta_s"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "X_mueta0"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
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
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s",
                                                                          "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0("c", k, "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0("c", k, "grad")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                             name = paste0("c", k, "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s), ",
                                                                    "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s), ",
                                                                    "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s))"),
                                                             name = paste0("c", k, "Y_psi_s")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_s")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "func %*% c", k,
                                                                    "Y_alpha_s, c", k, "Y_knot)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "beta_s"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:3, 1] + c", k, "beta %*% c", k,
                                                                    "mux"), name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s",
                                                                          "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0("c", k, "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0("c", k, "grad")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                             name = paste0("c", k, "Y_mean_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s), ",
                                                                    "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s), ",
                                                                    "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s))"),
                                                             name = paste0("c", k, "Y_psi_s")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxAlgebraFromString(paste0("rbind(c", k, "func %*% c", k,
                                                                    "Y_mean_s, c", k, "Y_knot)"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         TVC_info[[k]], GF_loadings[[k]], mxFitFunctionML(vector = T))
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s",
                                                                          "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0("c", k, "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0("c", k, "grad")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                             name = paste0("c", k, "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s), ",
                                                                    "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s), ",
                                                                    "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s))"),
                                                             name = paste0("c", k, "Y_psi_s")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                    "rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC))"),
                                                             name = paste0("c", k, "beta_s")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                             name = paste0("c", k, "BL_mean")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "func %*% c", k,
                                                                    "Y_alpha_s, c", k, "Y_knot)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "beta_s"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:3, 1] + c", k, "beta %*% c", k,
                                                                    "BL_mean"), name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                    labels = paste0("c", k, "beta", p - 1, "TVC"))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s",
                                                                          "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                         mxPath(from = "eta0sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records, "[1,1]")),
                                         mxPath(from = "eta2sY", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records, "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0("c", k, "func")),
                                         mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                    "cbind(0, 1, -1),",
                                                                    "cbind(0, 1, 1))"), name = paste0("c", k, "grad")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                             name = paste0("c", k, "Y_alpha_s")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s), ",
                                                                    "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s), ",
                                                                    "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s))"),
                                                             name = paste0("c", k, "Y_psi_s")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k, "beta2TVC)"),
                                                             name = paste0("c", k, "beta_s")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "func %*% c", k,
                                                                    "Y_alpha_s, c", k, "Y_knot)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "beta_s"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 [1:3, 1]+ c", k, "beta %*% c", k,
                                                                    "X_mueta0"), name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], BETA, mxFitFunctionML(vector = T))
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
    AddPara <- getMIX_UNI.addpara(dat = dat, nClass = nClass, curveFun = curveFun, intrinsic = intrinsic,
                                  t_var = t_var, records = records, growth_TIC = growth_TIC, decompose = decompose,
                                  starts = starts)
    if (curveFun %in% c("nonparametric", "NonP")){
      latents <- c("eta0Y", "eta1Y", paste0("dy", records[-1]), paste0("ly", records))
      Y_nGF <- length(latents) - (length(records) * 2 - 1)
      if (decompose == 0){
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                  labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                              free = c(F, rep(T, length(records) - 2)), values = c(1, starts[[k]][[1]][[4]][-1]),
                                              labels = paste0("c", k, "Y_rel_rate", 1:(length(records) - 1))),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]],
                                              labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                           name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                       mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                values = starts[[k]][[4]][1:2, 1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta1TIC", 1:length(growth_TIC))),
                                                byrow = T, name = paste0("c", k, "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                byrow = F, name = paste0("c", k, "mux")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                       mxFitFunctionML(vector = T))
          }
          else if (is.null(growth_TIC)){
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                              free = c(F, rep(T, length(records) - 2)), values = c(1, starts[[k]][[1]][[4]][-1]),
                                              labels = paste0("c", k, "Y_rel_rate", 1:(length(records) - 1))),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]],
                                              labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                           name = paste0("c", k, "Y_psi0")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                       mxFitFunctionML(vector = T))
          }
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                  labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                              free = c(F, rep(T, length(records) - 2)), values = c(1, starts[[k]][[1]][[4]][-1]),
                                              labels = paste0("c", k, "Y_rel_rate", 1:(length(records) - 1))),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]],
                                              labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                           name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                           name = paste0("c", k, "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                  "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                           name = paste0("c", k, "X_psi0")),
                                       mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                values = starts[[k]][[4]][1:2, 1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta1TIC", 1:length(growth_TIC))),
                                                byrow = T, name = paste0("c", k, "beta_TIC")),
                                       mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                  "rbind(c", k, "beta0TVC, c", k, "beta1TVC))"),
                                                           name = paste0("c", k, "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                byrow = F, name = paste0("c", k, "mux")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                           name = paste0("c", k, "BL_mean")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "BL_mean"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                       mxFitFunctionML(vector = T))
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                  labels = paste0("c", k, "beta", p - 1, "TVC"))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                       mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                       mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1,
                                              free = c(F, rep(T, length(records) - 2)), values = c(1, starts[[k]][[1]][[4]][-1]),
                                              labels = paste0("c", k, "Y_rel_rate", 1:(length(records) - 1))),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]],
                                              labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                           name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                           name = paste0("c", k, "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                  "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                           name = paste0("c", k, "X_psi0")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC)"),
                                                           name = paste0("c", k, "beta")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "X_mueta0"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                       mxFitFunctionML(vector = T))
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
              BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                  labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                              labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                  "Y_mueta2)"), name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                  "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                       mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta2TIC", 1:length(growth_TIC))),
                                                byrow = T, name = paste0("c", k, "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                byrow = F, name = paste0("c", k, "mux")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                       mxFitFunctionML(vector = T))
          }
          else if (is.null(growth_TIC)){
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                              labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                  "Y_mueta2)"), name = paste0("c", k, "Y_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                  "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                           name = paste0("c", k, "Y_psi0")),
                                       mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                              values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                       mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                              values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                       mxFitFunctionML(vector = T))
          }
        }
      }
      else if (decompose != 0){
        latents <- c(latents, X_latents)
        for (k in 1:nClass){
          if (!is.null(growth_TIC)){
            nTICs <- length(growth_TIC)
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                  labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                              labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                  "Y_mueta2)"), name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                  "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                           name = paste0("c", k, "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                  "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                           name = paste0("c", k, "X_psi0")),
                                       mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                           paste0("c", k, "beta2TIC", 1:length(growth_TIC))),
                                                byrow = T, name = paste0("c", k, "beta_TIC")),
                                       mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                  "rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                  "beta2TVC))"),
                                                           name = paste0("c", k, "beta")),
                                       mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                byrow = F, name = paste0("c", k, "mux")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                           name = paste0("c", k, "BL_mean")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "BL_mean"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                       mxFitFunctionML(vector = T))
          }
          else if (is.null(growth_TIC)){
            for (p in 1:Y_nGF){
              BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                  labels = paste0("c", k, "beta", p - 1, "TVC"))
            }
            class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                       manifestVars = manifests, latentVars = latents,
                                       mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]],
                                              labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                       mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                              free = TRUE, values = starts[[k]][[1]][[2]],
                                              labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                        "Y_psi11", "Y_psi12", "Y_psi22"))),
                                       mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                       mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                              labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                       mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                              free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                  "Y_mueta2)"), name = paste0("c", k, "Y_alpha0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                  "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                  "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                           name = paste0("c", k, "Y_psi_r")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                           name = paste0("c", k, "X_mean0")),
                                       mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                  "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                           name = paste0("c", k, "X_psi0")),
                                       mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                  "beta2TVC)"),
                                                           name = paste0("c", k, "beta")),
                                       mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "X_mueta0"),
                                                           name = paste0("c", k, "Y_mean0")),
                                       TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                       Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                       mxFitFunctionML(vector = T))
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
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
              }
              BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", c(paste0("TIC", 1:length(growth_TIC)))))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"),
                                                  name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "betagTIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"),
                                                  name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              }
              BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"),
                                                  name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "betagTIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                    "rbind(c", k, "beta0TVC, c", k, "beta1TVC, c",
                                                                    k, "betagTVC))"),
                                                             name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                             name = paste0("c", k, "BL_mean")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "BL_mean"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                    labels = paste0("c", k, "beta", p - 1, "TVC"))
              }
              BETA[[Y_nGF]] <- mxPath(from = "lx1", to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF],
                                      labels = paste0("c", k, "beta", "g", "TVC"))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"),
                                                  name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                                free = TRUE, values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g",
                                                                          "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE,
                                                values = 0, labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC, c",
                                                                    k, "betagTVC)"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "X_mueta0"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
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
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:2, 4)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:2, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:2, ] + c", k, "beta %*% c", k, "mux"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:2, 4)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:2, 4)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxMatrix("Full", 2, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:2, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                    "rbind(c", k, "beta0TVC, c", k, "beta1TVC))"),
                                                             name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                             name = paste0("c", k, "BL_mean")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:2, ] + c", k, "beta %*% c", k, "BL_mean"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                    labels = paste0("c", k, "beta", p - 1, "TVC"))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:2],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                                  labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:2, 4)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L1", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]],
                                                labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k,
                                                                    "Y_slp_ratio)"), name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC)"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:2, ] + c", k, "beta %*% c", k, "X_mueta0"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
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
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
              }
              BETA[[Y_nGF]] <- mxPath(from = growth_TIC, to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", c(paste0("TIC", 1:length(growth_TIC)))))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L3", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22, c", k, "Y_psi2g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psi2g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:4, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "betagTIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L3", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22, c", k, "Y_psi2g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psi2g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              }
              BETA[[Y_nGF]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF, ],
                                      labels = paste0("c", k, "beta", "g", c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L3", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22, c", k, "Y_psi2g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psi2g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxMatrix("Full", 4, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:4, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "betagTIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                    "rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC, c", k, "betagTVC))"),
                                                             name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                             name = paste0("c", k, "BL_mean")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "BL_mean"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              for (p in 1:(Y_nGF - 1)){
                BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                    labels = paste0("c", k, "beta", p - 1, "TVC"))
              }
              BETA[[Y_nGF]] <- mxPath(from = "lx1", to = latents[Y_nGF], arrows = 1, free = TRUE, values = starts[[k]][[4]][Y_nGF],
                                      labels = paste0("c", k, "beta", "g", "TVC"))
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                                free = TRUE,values = starts[[k]][[1]][[2]],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi0g", "Y_psi11",
                                                                          "Y_psi12", "Y_psi1g", "Y_psi22", "Y_psi2g", "Y_psigg"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = "deltag", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L3", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02, c", k, "Y_psi0g), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12, c", k, "Y_psi1g), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22, c", k, "Y_psi2g), ",
                                                                    "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psi2g, c", k, "Y_psigg))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC, c", k, "betagTVC)"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "X_mueta0"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
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
                BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)))))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:3, ] + c", k, "beta %*% c", k, "mux"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                             name = paste0("c", k, "Y_psi0")),
                                         mxPath(from = "one", to = paste0(TVC, records), arrows = 1, free = TRUE,
                                                values = starts[[k]][[2]][[1]], labels = paste0("c", k, "TVC_m", records)),
                                         mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2, free = TRUE,
                                                values = starts[[k]][[2]][[2]], labels = paste0("c", k, "TVC_v", records)),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]],
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
          }
        }
        else if (decompose != 0){
          latents <- c(latents, X_latents)
          for (k in 1:nClass){
            if (!is.null(growth_TIC)){
              nTICs <- length(growth_TIC)
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = c(growth_TIC, "lx1"), to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p, ],
                                    labels = paste0("c", k, "beta", p - 1, c(paste0("TIC", 1:length(growth_TIC)), "TVC")))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxMatrix("Full", 3, length(growth_TIC), free = TRUE,
                                                  values = starts[[k]][[4]][1:3, 1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "beta0TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta1TIC", 1:length(growth_TIC)),
                                                             paste0("c", k, "beta2TIC", 1:length(growth_TIC))),
                                                  byrow = T, name = paste0("c", k, "beta_TIC")),
                                         mxAlgebraFromString(paste0("cbind(c", k, "beta_TIC,",
                                                                    "rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC))"),
                                                             name = paste0("c", k, "beta")),
                                         mxMatrix("Full", length(growth_TIC), 1, free = TRUE, starts[[k]][[3]][[1]][1:length(growth_TIC)],
                                                  labels = c(paste0("c", k, "mux", 1:length(growth_TIC))),
                                                  byrow = F, name = paste0("c", k, "mux")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "mux, c", k, "X_mueta0)"),
                                                             name = paste0("c", k, "BL_mean")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:3, ] + c", k, "beta %*% c", k, "BL_mean"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
            else if (is.null(growth_TIC)){
              for (p in 1:Y_nGF){
                BETA[[p]] <- mxPath(from = "lx1", to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[4]][p],
                                    labels = paste0("c", k, "beta", p - 1, "TVC"))
              }
              class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                                         manifestVars = manifests, latentVars = latents,
                                         mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE, values = starts[[k]][[1]][[1]][1:3],
                                                labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                         mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                                  labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                         mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs", free = TRUE,
                                                values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                                labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02",
                                                                          "Y_psi11", "Y_psi12", "Y_psi22"))),
                                         mxPath(from = "eta0Y", to = "ly1", arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta1Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 1),
                                         mxPath(from = "eta2Y", to = paste0("dy", records[-1]), arrows = 1, free = FALSE, values = 0,
                                                labels = paste0("c", k, "L2", records[-1], "[1,1]")),
                                         mxPath(from = paste0(y_var, records), to = paste0(y_var, records), arrows = 2,
                                                free = TRUE, values = starts[[k]][[1]][[3]], labels = paste0("c", k, "Y_residuals")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                             name = paste0("c", k, "Y_alpha0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                    "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                    "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                             name = paste0("c", k, "Y_psi_r")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"),
                                                             name = paste0("c", k, "X_mean0")),
                                         mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                                    "cbind(c", k, "X_psi01, c", k, "X_psi11))"),
                                                             name = paste0("c", k, "X_psi0")),
                                         mxAlgebraFromString(paste0("rbind(c", k, "beta0TVC, c", k, "beta1TVC, c", k,
                                                                    "beta2TVC)"),
                                                             name = paste0("c", k, "beta")),
                                         mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:3, ] + c", k, "beta %*% c", k, "X_mueta0"),
                                                             name = paste0("c", k, "Y_mean0")),
                                         TVC_info[[k]], GF_loadings[[k]], AddPara[[k]], BETA,
                                         Y_PATH_L_L[[k]], Y_PATH_SLP_L[[k]], Y_PATH_AUTO_L[[k]],
                                         mxFitFunctionML(vector = T))
            }
          }
        }
      }
    }
  }
  return(class.list)
}



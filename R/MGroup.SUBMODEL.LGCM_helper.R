#' @title Define Latent Growth Curve Models as Class-specific Models (Submodels) for a Longitudinal Multiple Group Model
#'
#' @description This function defines latent growth curve models as class-specific models (submodels) for a longitudinal multiple group model.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions for each longitudinal process, and time-invariant covariates (TICs) if any.
#' It takes the value passed from \code{getMGroup()}.
#' @param nClass An integer specifying the number of manifested classes for the multiple group model. It takes the value passed from \code{getMGroup()}.
#' @param grp_var A string specifying the column that indicates manifested classes. It takes the value passed from \code{getMGroup()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' It takes the value passed from \code{getMGroup()}.
#' @param records A numeric vector specifying indices of the study waves. It takes the value passed from \code{getMGroup()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getMGroup()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for latent growth curve
#' models are: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). It takes the
#' value passed from \code{getMGroup()}.
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. It takes the value
#' passed from \code{getMGroup()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing to the
#' variability of growth factors if any. It takes the value passed from \code{getMGroup()}.
#' @param starts A list of initial values for the parameters, either takes the value passed from \code{getMGroup()} or
#' derived by the helper function \code{getMGroup.initial()}.
#'
#' @return A list of manifest and latent variables and paths for an mxModel object.
#'
#' @keywords internal
#'
#' @importFrom OpenMx mxPath mxModel mxAlgebraFromString mxMatrix
#'
getsub.LGCM_m <- function(dat, nClass, grp_var, t_var, y_var, curveFun, intrinsic, records, growth_TIC, starts){
  ## Define manifest variables
  manifests <- paste0(y_var, records)
  ## Define paths for adding growth TICs if any
  TIC_mean <- TIC_VAR <- BETA <- list()
  if (!is.null(growth_TIC)){
    for (k in 1:nClass){
      nTICs <- length(growth_TIC)
      ### Y_mean values of TIC(s)
      TIC_mean[[k]] <- mxPath(from = "one", to = growth_TIC, arrows = 1, free = TRUE, values = starts[[k]][[2]][[1]],
                              labels = paste0("c", k, "mux", 1:nTICs))
      ### Var-cov of TIC(s)
      TIC_VAR[[k]] <- mxPath(from = growth_TIC, to = growth_TIC, connect = "unique.pairs", arrows = 2, free = TRUE,
                             values = starts[[k]][[2]][[2]],
                             labels = paste0("c", k, "phi", 1:(nTICs * (nTICs + 1)/2)))
    }
  }
  GF_loadings <- getMIX_UNI.loadings(nClass = nClass, y_model = "LGCM", t_var = t_var, y_var = y_var,
                                     curveFun = curveFun, intrinsic = intrinsic, records = records)
  class.list <- list()
  ## Define latent variables, growth factor loadings, paths of the longitudinal outcome
  if (curveFun %in% c("linear", "LIN")){
    latents <- c("eta0", "eta1")
    for (k in 1:nClass){
      if (!is.null(growth_TIC)){
        nGF <- length(latents)
        nTICs <- length(growth_TIC)
        for (p in 1:nGF){
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                              labels = paste0("c", k, "beta", p - 1, 1:nTICs))
        }
        subdat <- dat[dat[, grp_var] == k, ]
        class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                   manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                   mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                          values = starts[[k]][[1]][[1]],
                                          labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                   mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]][[1]][[2]],
                                          labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                   mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                          values = 1),
                                   mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                          values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
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
                                   GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
      }
      else if (is.null(growth_TIC)){
        subdat <- dat[dat[, grp_var] == k, ]
        class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                   manifestVars = manifests, latentVars = latents,
                                   mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                          values = starts[[k]][[1]][[1]],
                                          labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                   mxPath(from = latents[1:2], to = latents[1:2], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]][[1]][[2]],
                                          labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi11"))),
                                   mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                          values = 1),
                                   mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE,
                                          values = 0, labels = paste0("c", k, "L1", records, "[1,1]")),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                          labels = paste0("c", k, "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1)"),
                                                       name = paste0("c", k, "Y_mean0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                              "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                       name = paste0("c", k, "Y_psi0")),
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
          BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                              labels = paste0("c", k, "beta", p - 1, 1:nTICs))
        }
        subdat <- dat[dat[, grp_var] == k, ]
        class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                   manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                   mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                          values = starts[[k]][[1]][[1]],
                                          labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                   mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]][[1]][[2]],
                                          labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                   mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0("c", k, "L1", records, "[1,1]")),
                                   mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0("c", k, "L2", records, "[1,1]")),
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
                                   GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
      }
      else if (is.null(growth_TIC)){
        subdat <- dat[dat[, grp_var] == k, ]
        class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                   manifestVars = manifests, latentVars = latents,
                                   mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                          values = starts[[k]][[1]][[1]],
                                          labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                   mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                          free = TRUE, values = starts[[k]][[1]][[2]],
                                          labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                   mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                   mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0("c", k, "L1", records, "[1,1]")),
                                   mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                          labels = paste0("c", k, "L2", records, "[1,1]")),
                                   mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                          arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                          labels = paste0("c", k, "Y_residuals")),
                                   mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2)"),
                                                       name = paste0("c", k, "Y_mean0")),
                                   mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                              "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                              "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                       name = paste0("c", k, "Y_psi0")),
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
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                                labels = paste0("c", k, "beta", p - 1, 1:nTICs))
          }
          BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[k]][[3]][nGF, ],
                                labels = paste0("c", k, "beta", "g", 1:nTICs))
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:2],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                              labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records, "[1,1]")),
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
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:2], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:2],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][3],
                                              labels = paste0("c", k, "Y_slp_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi0g", "Y_psi11", "Y_psi1g", "Y_psigg"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_slp_ratio)"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi0g), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi1g), ",
                                                                "cbind(c", k, "Y_psi0g, c", k, "Y_psi1g, c", k, "Y_psigg))"),
                                                         name = paste0("c", k, "Y_psi0")),
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
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                                labels = paste0("c", k, "beta", p - 1, 1:nTICs))
          }
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
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
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
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
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
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
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_slp_ratio)"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11))"),
                                                         name = paste0("c", k, "Y_psi0")),
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
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                                labels = paste0("c", k, "beta", p - 1, 1:nTICs))
          }
          BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[k]][[3]][nGF, ],
                                labels = paste0("c", k, "beta", "g", 1:nTICs))
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
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
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L3", records, "[1,1]")),
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
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
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
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L3", records, "[1,1]")),
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
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                                labels = paste0("c", k, "beta", p - 1, 1:nTICs))
          }
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:3],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                              labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records, "[1,1]")),
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
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:3],
                                            labels = paste0("c", k, c("Y_mueta0", "Y_mueta1", "Y_mueta2"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                              labels = paste0("c", k, "Y_acc_ratio"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                            labels = paste0("c", k, c("Y_psi00", "Y_psi01", "Y_psi02", "Y_psi11", "Y_psi12", "Y_psi22"))),
                                     mxPath(from = "eta0", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = "eta2", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0, c", k, "Y_mueta1, c", k, "Y_mueta2, c", k, "Y_acc_ratio)"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00, c", k, "Y_psi01, c", k, "Y_psi02), ",
                                                                "cbind(c", k, "Y_psi01, c", k, "Y_psi11, c", k, "Y_psi12), ",
                                                                "cbind(c", k, "Y_psi02, c", k, "Y_psi12, c", k, "Y_psi22))"),
                                                         name = paste0("c", k, "Y_psi0")),
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
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                                labels = paste0("c", k, "beta", p - 1, 1:nTICs))
          }
          BETA[[nGF]] <- mxPath(from = growth_TIC, to = latents[nGF], arrows = 1, free = TRUE, values = starts[[k]][[3]][nGF, ],
                                labels = paste0("c", k, "beta", "g", 1:nTICs))
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:3],
                                            labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                              labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]],
                                            labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                      "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                     mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L3", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                         name = paste0("c", k, "Y_alpha_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s, c", k, "Y_psi0gs), ",
                                                                "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s, c", k, "Y_psi1gs), ",
                                                                "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s, c", k, "Y_psi2gs), ",
                                                                "cbind(c", k, "Y_psi0gs, c", k, "Y_psi1gs, c", k, "Y_psi2gs, c", k, "Y_psiggs))"),
                                                         name = paste0("c", k, "Y_psi_s")),
                                     mxMatrix("Full", 4, length(growth_TIC), free = TRUE, values = starts[[k]][[3]],
                                              labels = c(paste0("c", k, "beta0", 1:length(growth_TIC)),
                                                         paste0("c", k, "beta1", 1:length(growth_TIC)),
                                                         paste0("c", k, "beta2", 1:length(growth_TIC)),
                                                         paste0("c", k, "betag", 1:length(growth_TIC))),
                                              byrow = T, name = paste0("c", k, "beta_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                "cbind(0, 1, -1, 0),",
                                                                "cbind(0, 1, 1, 0),",
                                                                "cbind(0, 0, 0, 1))"), name = paste0("c", k, "func")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                "cbind(0, 1, -1, 0),",
                                                                "cbind(0, 1, 1, 0),",
                                                                "cbind(0, 0, 0, 1))"), name = paste0("c", k, "grad")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "func[1:3, 1:3] %*% c", k, "Y_alpha_s, c", k, "Y_knot)"),
                                                         name = paste0("c", k, "Y_alpha0")),
                                     mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                         name = paste0("c", k, "Y_psi_r")),
                                     mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "beta_s"), name = paste0("c", k, "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]][[2]][[1]],
                                              labels = paste0("c", k, "mux", 1:length(growth_TIC)),
                                              byrow = F, name = paste0("c", k, "mux")),
                                     mxAlgebraFromString(paste0("c", k, "Y_alpha0 + c", k, "beta %*% c", k, "mux"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:3],
                                            labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                              labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:4], to = latents[1:4], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]],
                                            labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi0gs", "Y_psi11s",
                                                                      "Y_psi12s", "Y_psi1gs", "Y_psi22s", "Y_psi2gs", "Y_psiggs"))),
                                     mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records, "[1,1]")),
                                     mxPath(from = "deltag", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L3", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                         name = paste0("c", k, "Y_mean_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s, c", k, "Y_psi0gs), ",
                                                                "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s, c", k, "Y_psi1gs), ",
                                                                "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s, c", k, "Y_psi2gs), ",
                                                                "cbind(c", k, "Y_psi0gs, c", k, "Y_psi1gs, c", k, "Y_psi2gs, c", k, "Y_psiggs))"),
                                                         name = paste0("c", k, "Y_psi_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                "cbind(0, 1, -1, 0),",
                                                                "cbind(0, 1, 1, 0),",
                                                                "cbind(0, 0, 0, 1))"), name = paste0("c", k, "func")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot, 0),",
                                                                "cbind(0, 1, -1, 0),",
                                                                "cbind(0, 1, 1, 0),",
                                                                "cbind(0, 0, 0, 1))"), name = paste0("c", k, "grad")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "func[1:3, 1:3] %*% c", k, "Y_mean_s, c", k, "Y_knot)"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                         name = paste0("c", k, "Y_psi0")),
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
            BETA[[p]] <- mxPath(from = growth_TIC, to = latents[p], arrows = 1, free = TRUE, values = starts[[k]][[3]][p, ],
                                labels = paste0("c", k, "beta", p - 1, 1:nTICs))
          }
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                     manifestVars = c(manifests, growth_TIC), latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:3],
                                            labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                              labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                            labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                     mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                         name = paste0("c", k, "Y_alpha_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s), ",
                                                                "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s), ",
                                                                "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s))"),
                                                         name = paste0("c", k, "Y_psi_s")),
                                     mxMatrix("Full", 3, length(growth_TIC), free = TRUE, values = starts[[k]][[3]][1:3, ],
                                              labels = c(paste0("c", k, "beta0", 1:length(growth_TIC)),
                                                         paste0("c", k, "beta1", 1:length(growth_TIC)),
                                                         paste0("c", k, "beta2", 1:length(growth_TIC))),
                                              byrow = T, name = paste0("c", k, "beta_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                "cbind(0, 1, -1),",
                                                                "cbind(0, 1, 1))"), name = paste0("c", k, "func")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                "cbind(0, 1, -1),",
                                                                "cbind(0, 1, 1))"), name = paste0("c", k, "grad")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "func %*% c", k, "Y_alpha_s, c", k, "Y_knot)"),
                                                         name = paste0("c", k, "Y_alpha0")),
                                     mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                         name = paste0("c", k, "Y_psi_r")),
                                     mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "beta_s"), name = paste0("c", k, "beta")),
                                     mxMatrix("Full", length(growth_TIC), 1, free = TRUE, values = starts[[k]][[2]][[1]],
                                              labels = paste0("c", k, "mux", 1:length(growth_TIC)),
                                              byrow = F, name = paste0("c", k, "mux")),
                                     mxAlgebraFromString(paste0("c", k, "Y_alpha0[1:3, 1] + c", k, "beta %*% c", k, "mux"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     GF_loadings[[k]], TIC_mean[[k]], TIC_VAR[[k]], BETA)
        }
        else if (is.null(growth_TIC)){
          subdat <- dat[dat[, grp_var] == k, ]
          class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", mxData(observed = subdat, type = "raw"),
                                     manifestVars = manifests, latentVars = latents,
                                     mxPath(from = "one", to = latents[1:3], arrows = 1, free = TRUE,
                                            values = starts[[k]][[1]][[1]][1:3],
                                            labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                                     mxMatrix("Full", 1, 1, free = TRUE, values = starts[[k]][[1]][[1]][4],
                                              labels = paste0("c", k, "Y_knot"), name = paste0("c", k, "Y_mug")),
                                     mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                                            free = TRUE, values = starts[[k]][[1]][[2]][c(1:3, 5:6, 8)],
                                            labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s", "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                                     mxPath(from = "eta0s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 1),
                                     mxPath(from = "eta1s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L1", records, "[1,1]")),
                                     mxPath(from = "eta2s", to = paste0(y_var, records), arrows = 1, free = FALSE, values = 0,
                                            labels = paste0("c", k, "L2", records, "[1,1]")),
                                     mxPath(from = paste0(y_var, records), to = paste0(y_var, records),
                                            arrows = 2, free = TRUE, values = starts[[k]][[1]][[3]],
                                            labels = paste0("c", k, "Y_residuals")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"),
                                                         name = paste0("c", k, "Y_mean_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s), ",
                                                                "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s), ",
                                                                "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s))"),
                                                         name = paste0("c", k, "Y_psi_s")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                "cbind(0, 1, -1),",
                                                                "cbind(0, 1, 1))"), name = paste0("c", k, "func")),
                                     mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "Y_knot, c", k, "Y_knot),",
                                                                "cbind(0, 1, -1),",
                                                                "cbind(0, 1, 1))"), name = paste0("c", k, "grad")),
                                     mxAlgebraFromString(paste0("rbind(c", k, "func %*% c", k, "Y_mean_s, c", k, "Y_knot)"),
                                                         name = paste0("c", k, "Y_mean0")),
                                     mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"),
                                                         name = paste0("c", k, "Y_psi0")),
                                     GF_loadings[[k]])
        }
      }
    }
  }
  ## Define the output of the function
  return(class.list)
}

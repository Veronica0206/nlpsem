#' @title Construct An Object of mxModel for Multivariate Latent Growth Curve Models or Multivariate Latent Change Score
#' Models To Be Evaluated
#'
#' @description This function builds up an object of mxModel for a multivariate latent growth curve model or a multivariate latent
#' change score model with user-specified functional form (including whether intrinsically nonlinear).
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables
#' with repeated measurements and occasions for multiple longitudinal outcomes. It takes the value passed from \code{getMGM()}.
#' @param t_var A vector of strings, with each element representing the prefix for column names related to the time
#' variable for the corresponding outcome variable at each study wave. It takes the value passed from \code{getMGM()}.
#' @param y_var A vector of strings, with each element representing the prefix for column names corresponding to a
#' particular outcome variable at each study wave. It takes the value passed from \code{getMGM()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the
#' value passed from \code{getMGM()}.
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. It takes the
#' value passed from \code{getMGM()}.
#' @param records A list of numeric vectors, with each vector specifying the indices of the observed study waves for
#' the corresponding outcome variable. It takes the value passed from \code{getMGM()}.
#' @param y_model A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' It takes the value passed from \code{getMGM()}.
#' @param starts A list of initial values for the parameters, either takes the value passed from \code{getMGM()} or derived
#' by the helper function \code{getMULTI.initial()}.
#'
#' @return A pre-optimized mxModel for a multivariate latent growth curve model or a multivariate latent change
#' score model.
#'
getMGM.mxModel <- function(dat, t_var, y_var, curveFun, intrinsic, records, y_model, starts){
  y_var_L <- tolower(y_var)
  # Define manifest variables
  traj_list <- list()
  for (traj in 1:length(y_var)){
    traj_list[[length(traj_list) + 1]] <- paste0(y_var[traj], records[[traj]])
  }
  manifests <- unlist(traj_list)
  # Define variances and covariance(s) of residuals
  RES <- COV <- list()
  for (traj in 1:length(y_var)){
    RES[[traj]] <- mxPath(from = traj_list[[traj]], to = traj_list[[traj]], arrows = 2, free = TRUE,
                          values = starts[[3]][traj, traj],
                          labels = paste0(y_var[traj], "_residuals"))
  }
  for (traj_i in 1:(length(y_var) - 1)){
    for (traj_j in traj_i:(length(y_var) - 1)){
      if (setequal(substr(traj_list[[traj_i]], 2, 2), substr(traj_list[[traj_j + 1]], 2, 2))){
        COV[[traj_i + traj_j - 1]] <- mxPath(from = traj_list[[traj_i]], to = traj_list[[traj_j + 1]],
                                             arrows = 2, free = TRUE, values = starts[[3]][traj_i, traj_j + 1],
                                             labels = paste0(y_var[traj_i], y_var[traj_j + 1], "_RES"))
      }
      else{
        T_common <- Reduce(intersect, list(substr(traj_list[[traj_i]], 2, 2), substr(traj_list[[traj_j + 1]], 2, 2)))
        COV[[traj_i + traj_j - 1]] <- mxPath(from = paste0(y_var[traj_i], T_common),
                                             to = paste0(y_var[traj_j + 1], T_common),
                                             arrows = 2, free = TRUE, values = starts[[3]][traj_i, traj_j + 1],
                                             labels = paste0(y_var[traj_i], y_var[traj_j + 1], "_RES"))
      }
    }
  }
  # Obtain factor loadigs for the specified trajectory model (either "LGCM" or "LCSM") and functional form
  GF_loadings <- getMULTI.loadings(y_model = y_model, t_var = t_var, y_var = y_var, curveFun = curveFun,
                                   intrinsic = intrinsic, records = records)
  if (y_model == "LGCM"){
    if (curveFun %in% c("linear", "LIN")){
      latents <- paste0(rep(c("eta0", "eta1"), length(y_var)), rep(y_var, each = 2))
      multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
      GF_MEAN <- GF_LOADINGS <- list()
      mean_L <- psi_L <- list()
      gf_var_label <- gf_cov_label <- list()
      psi_btw_L <- list()
      for (traj in 1:length(y_var)){
        GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                  arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:2],
                                  labels = paste0(y_var[traj], c("_mueta0", "_mueta1")))
        GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                           arrows = 1, free = FALSE, values = 1),
                                    mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                           arrows = 1, free = FALSE,
                                           labels = paste0("L1", records[[traj]], y_var[traj], "[1,1]")))
        mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj], "_mueta1)"),
                                              name = paste0(y_var[traj], "_mean0"))
        psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj], "_psi01), ",
                                                    "cbind(", y_var[traj], "_psi01, ", y_var[traj], "_psi11))"),
                                             name = paste0(y_var[traj], "_psi0"))
        gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi", c("00", "01", "01", "11")),
                                       byrow = T, nrow = 2, ncol = 2)
      }
      for (traj_i in 1:(length(y_var) - 1)){
        for (traj_j in traj_i:(length(y_var) - 1)){
          psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ", y_var[traj_i], y_var[traj_j + 1], "_psi01), ",
                                                                         "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ", y_var[traj_i], y_var[traj_j + 1], "_psi11))"),
                                                                  name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
          gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi", c("00", "01", "10", "11")),
                                                        byrow = T, nrow = 2, ncol = 2)
        }
      }
      for (traj in 1:length(y_var)){
        multi_label[((traj - 1) * 2 + 1):((traj - 1) * 2 + 2), ((traj - 1) * 2 + 1):((traj - 1) * 2 + 2)] <- gf_var_label[[traj]]
      }
      for (traj_i in 1:(length(y_var) - 1)){
        for (traj_j in traj_i:(length(y_var) - 1)){
          multi_label[((traj_i - 1) * 2 + 1):((traj_i - 1) * 2 + 2), (traj_j * 2 + 1):(traj_j * 2 + 2)] <- gf_cov_label[[traj_i + traj_j - 1]]
        }
      }
      model_mx <- mxModel("Multivariate Latent Growth Curve Model, Linear", type = "RAM",
                          manifestVars = manifests, latentVars = latents,
                          mxData(observed = dat, type = "raw"),
                          mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                 free = TRUE, values = t(starts[[2]])[row(t(starts[[2]])) >= col(t(starts[[2]]))],
                                 labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                          GF_loadings, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L, psi_btw_L)
    }
    else if (curveFun %in% c("quadratic", "QUAD")){
      latents <- paste0(rep(c("eta0", "eta1", "eta2"), length(y_var)), rep(y_var, each = 3))
      multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
      GF_MEAN <- GF_LOADINGS <- list()
      mean_L <- psi_L <- list()
      gf_var_label <- gf_cov_label <- list()
      psi_btw_L <- list()
      for (traj in 1:length(y_var)){
        GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                  arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:3],
                                  labels = paste0(y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
        GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                           arrows = 1, free = FALSE, values = 1),
                                    mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                           arrows = 1, free = FALSE,
                                           labels = paste0("L1", records[[traj]], y_var[traj], "[1,1]")),
                                    mxPath(from = paste0("eta2", y_var[traj]), to = traj_list[[traj]],
                                           arrows = 1, free = FALSE,
                                           labels = paste0("L2", records[[traj]], y_var[traj], "[1,1]")))
        mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj],
                                                     "_mueta1, ", y_var[traj], "_mueta2)"),
                                              name = paste0(y_var[traj], "_mean0"))
        psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj],
                                                    "_psi01, ", y_var[traj], "_psi02), ",
                                                    "cbind(", y_var[traj], "_psi01, ", y_var[traj],
                                                    "_psi11, ", y_var[traj], "_psi12), ",
                                                    "cbind(", y_var[traj], "_psi02, ", y_var[traj],
                                                    "_psi12, ", y_var[traj], "_psi22))"),
                                             name = paste0(y_var[traj], "_psi0"))
        gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi",
                                              c("00", "01", "02", "01", "11", "12", "02", "12", "22")),
                                       byrow = T, nrow = 3, ncol = 3)
      }
      for (traj_i in 1:(length(y_var) - 1)){
        for (traj_j in traj_i:(length(y_var) - 1)){
          psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi01, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi02),",
                                                                         "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi12), ",
                                                                         "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi20, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi21, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi22))"),
                                                                  name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
          gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                               c("00", "01", "02", "10", "11", "12", "20", "21", "22")),
                                                        byrow = T, nrow = 3, ncol = 3)
        }
      }
      for (traj in 1:length(y_var)){
        multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <- gf_var_label[[traj]]
      }
      for (traj_i in 1:(length(y_var) - 1)){
        for (traj_j in traj_i:(length(y_var) - 1)){
          multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <- gf_cov_label[[traj_i + traj_j - 1]]
        }
      }
      model_mx <- mxModel("Multivariate Latent Growth Curve Model, Quadratic", type = "RAM",
                          manifestVars = manifests, latentVars = latents,
                          mxData(observed = dat, type = "raw"),
                          mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                 free = TRUE, values = t(starts[[2]])[row(t(starts[[2]])) >= col(t(starts[[2]]))],
                                 labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                          GF_loadings, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L, psi_btw_L)
    }
    else if (curveFun %in% c("negative exponential", "EXP")){
      if (intrinsic){
        latents <- paste0(rep(c("eta0", "eta1", "deltag"), length(y_var)), rep(y_var, each = 3))
        multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
        GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in 1:length(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:2],
                                    labels = paste0(y_var[traj], c("_mueta0", "_mueta1")))
          GAMMA[[traj]]  <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[traj]][3],
                                     labels = paste0(y_var[traj], "_slp_ratio"),
                                     name = paste0(y_var[traj], "_mug"))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj], "_mueta1, ",
                                                       y_var[traj], "_slp_ratio)"),
                                                name = paste0(y_var[traj], "_mean0"))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE,
                                             labels = paste0("L1", records[[traj]], y_var[traj], "[1,1]")),
                                      mxPath(from = paste0("deltag", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L2", records[[traj]], y_var[traj], "[1,1]")))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj],
                                                      "_psi01, ", y_var[traj], "_psi0g), ",
                                                      "cbind(", y_var[traj], "_psi01, ", y_var[traj],
                                                      "_psi11, ", y_var[traj], "_psi1g), ",
                                                      "cbind(", y_var[traj], "_psi0g, ", y_var[traj],
                                                      "_psi1g, ", y_var[traj], "_psigg))"),
                                               name = paste0(y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi",
                                                c("00", "01", "0g", "01", "11", "1g", "0g", "1g", "gg")),
                                         byrow = T, nrow = 3, ncol = 3)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                           y_var[traj_i], y_var[traj_j + 1],
                                                                           "_psi01, ", y_var[traj_i], y_var[traj_j + 1], "_psi0g),",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi1g), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psig0, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psig1, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psigg))"),
                                                                    name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                 c("00", "01", "0g", "10", "11", "1g", "g0", "g1", "gg")),
                                                          byrow = T, nrow = 3, ncol = 3)
          }
        }
        for (traj in 1:length(y_var)){
          multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
            gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
              gf_cov_label[[traj_i + traj_j - 1]]
          }
        }
        model_mx <- mxModel("Multivariate Latent Growth Curve Model, Negative Exponential (random ratio)", type = "RAM",
                            manifestVars = manifests, latentVars = latents,
                            mxData(observed = dat, type = "raw"),
                            mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                   free = TRUE, values = t(starts[[2]])[row(t(starts[[2]])) >= col(t(starts[[2]]))],
                                   labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                            GF_loadings, GAMMA, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L, psi_btw_L)
      }
      else if (!intrinsic){
        latents <- paste0(rep(c("eta0", "eta1"), length(y_var)), rep(y_var, each = 2))
        multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
        GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in 1:length(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:2],
                                    labels = paste0(y_var[traj], c("_mueta0", "_mueta1")))
          GAMMA[[traj]]  <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[traj]][3],
                                     labels = paste0(y_var[traj], "_slp_ratio"),
                                     name = paste0(y_var[traj], "_mug"))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj], "_mueta1, ",
                                                       y_var[traj], "_slp_ratio)"),
                                                name = paste0(y_var[traj], "_mean0"))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE,
                                             labels = paste0("L1", records[[traj]], y_var[traj], "[1,1]")))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj], "_psi01), ",
                                                      "cbind(", y_var[traj], "_psi01, ", y_var[traj], "_psi11))"),
                                               name = paste0(y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi", c("00", "01", "01", "11")),
                                         byrow = T, nrow = 2, ncol = 2)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi01), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi11))"),
                                                                    name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi", c("00", "01", "10", "11")),
                                                          byrow = T, nrow = 2, ncol = 2)
          }
        }
        for (traj in 1:length(y_var)){
          multi_label[((traj - 1) * 2 + 1):((traj - 1) * 2 + 2), ((traj - 1) * 2 + 1):((traj - 1) * 2 + 2)] <-
            gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 2 + 1):((traj_i - 1) * 2 + 2), (traj_j * 2 + 1):(traj_j * 2 + 2)] <-
              gf_cov_label[[traj_i + traj_j - 1]]
          }
        }
        model_mx <- mxModel("Multivariate Latent Growth Curve Model, Negative Exponential (fixed ratio)", type = "RAM",
                            manifestVars = manifests, latentVars = latents,
                            mxData(observed = dat, type = "raw"),
                            mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = t(starts[[2]][c(1:2, 4:5), c(1:2, 4:5)])[row(t(starts[[2]][c(1:2, 4:5), c(1:2, 4:5)])) >=
                                                                                       col(t(starts[[2]][c(1:2, 4:5), c(1:2, 4:5)]))],
                                   labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                            GF_loadings, GAMMA, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L, psi_btw_L)
      }
    }
    else if (curveFun %in% c("Jenss-Bayley", "JB")){
      if (intrinsic){
        latents <- paste0(rep(c("eta0", "eta1", "eta2", "deltag"), length(y_var)), rep(y_var, each = 4))
        multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
        GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in 1:length(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:3],
                                    labels = paste0(y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
          GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[traj]][4],
                                    labels = paste0(y_var[traj], "_acc_ratio"),
                                    name = paste0(y_var[traj], "_mug"))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj], "_mueta1, ",
                                                       y_var[traj], "_mueta2, ", y_var[traj], "_acc_ratio)"),
                                                name = paste0(y_var[traj], "_mean0"))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L1", records[[traj]], y_var[traj], "[1,1]")),
                                      mxPath(from = paste0("eta2", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L2", records[[traj]], y_var[traj], "[1,1]")),
                                      mxPath(from = paste0("deltag", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L3", records[[traj]], y_var[traj], "[1,1]")))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj], "_psi01, ",
                                                      y_var[traj], "_psi02, ", y_var[traj], "_psi0g), ",
                                                      "cbind(", y_var[traj], "_psi01, ", y_var[traj], "_psi11, ",
                                                      y_var[traj], "_psi12, ", y_var[traj], "_psi1g), ",
                                                      "cbind(", y_var[traj], "_psi02, ", y_var[traj], "_psi12, ",
                                                      y_var[traj], "_psi22, ", y_var[traj], "_psi2g), ",
                                                      "cbind(", y_var[traj], "_psi0g, ", y_var[traj], "_psi1g, ",
                                                      y_var[traj], "_psi2g, ", y_var[traj], "_psigg))"),
                                               name = paste0(y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi", c("00", "01", "02", "0g", "01", "11", "12", "1g",
                                                                       "02", "12", "22", "2g", "0g", "1g", "2g", "gg")),
                                         byrow = T, nrow = 4, ncol = 4)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi01, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi02, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi0g), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi12, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi1g), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi20, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi21, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi22, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi2g), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psig0, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psig1, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psig2, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psigg))"),
                                                                    name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                 c("00", "01", "02", "0g", "10", "11", "12", "1g",
                                                                   "20", "21", "22", "2g", "g0", "g1", "g2", "gg")),
                                                          byrow = T, nrow = 4, ncol = 4)
          }
        }
        for (traj in 1:length(y_var)){
          multi_label[((traj - 1) * 4 + 1):((traj - 1) * 4 + 4), ((traj - 1) * 4 + 1):((traj - 1) * 4 + 4)] <-
            gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 4 + 1):((traj_i - 1) * 4 + 4), (traj_j * 4 + 1):(traj_j * 4 + 4)] <-
              gf_cov_label[[traj_i + traj_j - 1]]
          }
        }
        model_mx <- mxModel("Multivariate Latent Growth Curve Model, Jenss Bayley (random ratio)", type = "RAM",
                            manifestVars = manifests, latentVars = latents,
                            mxData(observed = dat, type = "raw"),
                            mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                   free = TRUE, values = t(starts[[2]])[row(t(starts[[2]])) >= col(t(starts[[2]]))],
                                   labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                            GF_loadings, GAMMA, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L, psi_btw_L)
      }
      else if (!intrinsic){
        latents <- paste0(rep(c("eta0", "eta1", "eta2"), length(y_var)), rep(y_var, each = 3))
        multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
        GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in 1:length(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:3],
                                    labels = paste0(y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
          GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[traj]][4],
                                    labels = paste0(y_var[traj], "_acc_ratio"),
                                    name = paste0(y_var[traj], "_mug"))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj], "_mueta1, ",
                                                       y_var[traj], "_mueta2, ", y_var[traj], "_acc_ratio)"),
                                                name = paste0(y_var[traj], "_mean0"))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L1", records[[traj]], y_var[traj], "[1,1]")),
                                      mxPath(from = paste0("eta2", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L2", records[[traj]], y_var[traj], "[1,1]")))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj], "_psi01, ",
                                                      y_var[traj], "_psi02), ",
                                                      "cbind(", y_var[traj], "_psi01, ", y_var[traj], "_psi11, ",
                                                      y_var[traj], "_psi12), ",
                                                      "cbind(", y_var[traj], "_psi02, ", y_var[traj], "_psi12, ",
                                                      y_var[traj], "_psi22))"),
                                               name = paste0(y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi", c("00", "01", "02", "01", "11", "12",
                                                                       "02", "12", "22")),
                                         byrow = T, nrow = 3, ncol = 3)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi01, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi02), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi12), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi20, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi21, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi22))"),
                                                                    name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                 c("00", "01", "02", "10", "11", "12",
                                                                   "20", "21", "22")),
                                                          byrow = T, nrow = 3, ncol = 3)
          }
        }
        for (traj in 1:length(y_var)){
          multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
            gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
              gf_cov_label[[traj_i + traj_j - 1]]
          }
        }
        model_mx <- mxModel("Multivariate Latent Growth Curve Model, Jenss Bayley (fixed ratio)", type = "RAM",
                            manifestVars = manifests, latentVars = latents,
                            mxData(observed = dat, type = "raw"),
                            mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = t(starts[[2]][c(1:3, 5:7), c(1:3, 5:7)])[row(t(starts[[2]][c(1:3, 5:7), c(1:3, 5:7)])) >=
                                                                                       col(t(starts[[2]][c(1:3, 5:7), c(1:3, 5:7)]))],
                                   labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                            GF_loadings, GAMMA, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L, psi_btw_L)
      }
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      if (intrinsic){
        latents <- paste0(rep(c("eta0s", "eta1s", "eta2s", "deltag"), length(y_var)), rep(y_var, each = 4))
        multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
        GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
        func_L <- grad_L <- mean_s_L <- psi_s_L <- mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_s_L <- psi_btw_L <- list()
        for (traj in 1:length(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0s", "eta1s", "eta2s"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:3],
                                    labels = paste0(y_var[traj], c("_mueta0s", "_mueta1s", "_mueta2s")))
          GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[traj]][4],
                                    labels = paste0(y_var[traj], "_knot"),
                                    name = paste0(y_var[traj], "_mug"))
          mean_s_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0s, ", y_var[traj], "_mueta1s, ",
                                                         y_var[traj], "_mueta2s, ", y_var[traj], "_knot)"),
                                                  name = paste0(y_var[traj], "_mean_s"))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0s", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1s", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L1", records[[traj]], y_var[traj], "[1,1]")),
                                      mxPath(from = paste0("eta2s", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L2", records[[traj]], y_var[traj], "[1,1]")),
                                      mxPath(from = paste0("deltag", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L3", records[[traj]], y_var[traj], "[1,1]")))
          func_L[[length(func_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -", y_var[traj], "_mug, ", y_var[traj], "_mug, 0), ",
                                                                     "cbind(0, 1, -1, 0), ", "cbind(0, 1, 1, 0)", ", cbind(0, 0, 0, 1))"),
                                                              name = paste0("func", y_var[traj]))
          grad_L[[length(grad_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -", y_var[traj], "_mug, ", y_var[traj], "_mug, 0), ",
                                                                     "cbind(0, 1, -1, 0), ", "cbind(0, 1, 1, 0)", ", cbind(0, 0, 0, 1))"),
                                                              name = paste0("grad", y_var[traj]))
          psi_s_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00s, ", y_var[traj], "_psi01s, ",
                                                        y_var[traj], "_psi02s, ", y_var[traj], "_psi0gs), ",
                                                        "cbind(", y_var[traj], "_psi01s, ", y_var[traj], "_psi11s, ",
                                                        y_var[traj], "_psi12s, ", y_var[traj], "_psi1gs), ",
                                                        "cbind(", y_var[traj], "_psi02s, ", y_var[traj], "_psi12s, ",
                                                        y_var[traj], "_psi22s, ", y_var[traj], "_psi2gs), ",
                                                        "cbind(", y_var[traj], "_psi0gs, ", y_var[traj], "_psi1gs, ",
                                                        y_var[traj], "_psi2gs, ", y_var[traj], "_psiggs))"),
                                                 name = paste0(y_var[traj], "_psi_s"))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(func", y_var[traj], " %*% ", y_var[traj],
                                                       "_mean_s)"), name = paste0(y_var[traj], "_mean0"))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("grad", y_var[traj], " %*% ", y_var[traj], "_psi_s",
                                                      " %*% t(grad", y_var[traj], ")"),
                                               name = paste0(y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi", c("00s", "01s", "02s", "0gs", "01s", "11s", "12s", "1gs",
                                                                       "02s", "12s", "22s", "2gs", "0gs", "1gs", "2gs", "ggs")),
                                         byrow = T, nrow = 4, ncol = 4)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_s_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi01s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi02s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi0gs), ",
                                                                             "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi11s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi12s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi1gs), ",
                                                                             "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi20s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi21s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi22s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi2gs), ",
                                                                             "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psig0s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psig1s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psig2s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psiggs))"),
                                                                      name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi_s"))
            psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("grad", y_var[traj_i], " %*% ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi_s",
                                                                           " %*% t(grad", y_var[traj_j + 1], ")"),
                                                                    name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                 c("00s", "01s", "02s", "0gs", "10s", "11s", "12s", "1gs",
                                                                   "20s", "21s", "22s", "2gs", "g0s", "g1s", "g2s", "ggs")),
                                                          byrow = T, nrow = 4, ncol = 4)
          }
        }
        for (traj in 1:length(y_var)){
          multi_label[((traj - 1) * 4 + 1):((traj - 1) * 4 + 4), ((traj - 1) * 4 + 1):((traj - 1) * 4 + 4)] <-
            gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 4 + 1):((traj_i - 1) * 4 + 4), (traj_j * 4 + 1):(traj_j * 4 + 4)] <-
              gf_cov_label[[traj_i + traj_j - 1]]
          }
        }
        model_mx <- mxModel("Multivariate Latent Growth Curve Model, Bilinear Spline (random knot)", type = "RAM",
                            manifestVars = manifests, latentVars = latents,
                            mxData(observed = dat, type = "raw"),
                            mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                   free = TRUE, values = t(starts[[2]])[row(t(starts[[2]])) >= col(t(starts[[2]]))],
                                   labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                            GF_loadings, GAMMA, func_L, grad_L, mean_s_L, psi_s_L, psi_btw_s_L,
                            GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L, psi_btw_L)
      }
      else if (!intrinsic){
        latents <- paste0(rep(c("eta0s", "eta1s", "eta2s"), length(y_var)), rep(y_var, each = 3))
        multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
        GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
        func_L <- grad_L <- mean_s_L <- psi_s_L <- mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_s_L <- psi_btw_L <- list()
        for (traj in 1:length(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0s", "eta1s", "eta2s"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:3],
                                    labels = paste0(y_var[traj], c("_mueta0s", "_mueta1s", "_mueta2s")))
          GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[traj]][4],
                                    labels = paste0(y_var[traj], "_knot"),
                                    name = paste0(y_var[traj], "_mug"))
          mean_s_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0s, ", y_var[traj], "_mueta1s, ",
                                                         y_var[traj], "_mueta2s)"),
                                                  name = paste0(y_var[traj], "_mean_s"))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0s", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1s", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L1", records[[traj]], y_var[traj], "[1,1]")),
                                      mxPath(from = paste0("eta2s", y_var[traj]), to = traj_list[[traj]],
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0("L2", records[[traj]], y_var[traj], "[1,1]")))
          func_L[[length(func_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -", y_var[traj], "_mug, ", y_var[traj], "_mug), ",
                                                                     "cbind(0, 1, -1), ", "cbind(0, 1, 1))"),
                                                              name = paste0("func", y_var[traj]))
          grad_L[[length(grad_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -", y_var[traj], "_mug, ", y_var[traj], "_mug), ",
                                                                     "cbind(0, 1, -1), ", "cbind(0, 1, 1))"),
                                                              name = paste0("grad", y_var[traj]))
          psi_s_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00s, ", y_var[traj], "_psi01s, ",
                                                        y_var[traj], "_psi02s), ",
                                                        "cbind(", y_var[traj], "_psi01s, ", y_var[traj], "_psi11s, ",
                                                        y_var[traj], "_psi12s), ",
                                                        "cbind(", y_var[traj], "_psi02s, ", y_var[traj], "_psi12s, ",
                                                        y_var[traj], "_psi22s))"),
                                                 name = paste0(y_var[traj], "_psi_s"))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(func", y_var[traj], " %*% ", y_var[traj],
                                                       "_mean_s, ", y_var[traj], "_knot)"),
                                                name = paste0(y_var[traj], "_mean0"))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("grad", y_var[traj], " %*% ", y_var[traj], "_psi_s",
                                                      " %*% t(grad", y_var[traj], ")"),
                                               name = paste0(y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi", c("00s", "01s", "02s", "01s", "11s", "12s",
                                                                       "02s", "12s", "22s")),
                                         byrow = T, nrow = 3, ncol = 3)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_s_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi01s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi02s), ",
                                                                             "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi11s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi12s), ",
                                                                             "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi20s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi21s, ",
                                                                             y_var[traj_i], y_var[traj_j + 1], "_psi22s))"),
                                                                      name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi_s"))
            psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("grad", y_var[traj_i], " %*% ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi_s",
                                                                           " %*% t(grad", y_var[traj_j + 1], ")"),
                                                                    name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                 c("00s", "01s", "02s", "10s", "11s", "12s", "20s", "21s", "22s")),
                                                          byrow = T, nrow = 3, ncol = 3)
          }
        }
        for (traj in 1:length(y_var)){
          multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
            gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
              gf_cov_label[[traj_i + traj_j - 1]]
          }
        }
        model_mx <- mxModel("Multivariate Latent Growth Curve Model, Bilinear Spline (fixed knot)", type = "RAM",
                            manifestVars = manifests, latentVars = latents,
                            mxData(observed = dat, type = "raw"),
                            mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = t(starts[[2]][c(1:3, 5:7), c(1:3, 5:7)])[row(t(starts[[2]][c(1:3, 5:7), c(1:3, 5:7)])) >=
                                                                                       col(t(starts[[2]][c(1:3, 5:7), c(1:3, 5:7)]))],
                                   labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                            GF_loadings, GAMMA, func_L, grad_L, mean_s_L, psi_s_L, psi_btw_s_L,
                            GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L, psi_btw_L)
      }
    }
  }
  else if (y_model == "LCSM"){
    # Obtain additional parameters (derived) for the specified functional form
    AddPara <- getMULTI.addpara(dat = dat, y_var = y_var, curveFun = curveFun, intrinsic = intrinsic,
                                t_var = t_var, records = records, starts = starts)
    # Define latent variables in addition to growth factors
    L_list <- D_list <- list()
    for (traj in 1:length(y_var)){
      L_list[[traj]] <- paste0("l", y_var_L[traj], records[[traj]])
      D_list[[traj]] <- paste0("d", y_var_L[traj], records[[traj]][-1])
    }
    PATH_L_L <- PATH_SLP_L <- PATH_AUTO_L <- list()
    for (traj in 1:length(y_var)){
      # Define paths from latent true scores to observed scores
      PATH_L_L[[traj]] <- mxPath(from = paste0("l", y_var_L[traj], records[[traj]]),
                                 to = paste0(y_var[[traj]], records[[traj]]), arrows = 1, free = FALSE, values = 1)
      # Define paths from latent instantaneous rate of change at each measurement to true scores
      PATH_SLP_L[[traj]] <- mxPath(from = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                   to = paste0("l", y_var_L[traj], records[[traj]][-1]),
                                   arrows = 1, free = FALSE, values = 0,
                                   labels = paste0(y_var[traj], "lag", records[[traj]][-1], "[1,1]"))
      # Define autoregressive paths
      PATH_AUTO_L[[traj]] <- mxPath(from = paste0("l", y_var_L[traj], records[[traj]][-length(records[[traj]])]),
                                    to = paste0("l", y_var_L[traj], records[[traj]][-1]), arrows = 1,
                                    free = FALSE, values = 1)
    }
    if (curveFun %in% c("nonparametric", "NonP")){
      latents <- c(paste0(rep(c("eta0", "eta1"), length(y_var)), rep(y_var, each = 2)),
                   unlist(L_list), unlist(D_list))
      multi_label <- matrix(NA, nrow = 2 * length(y_var), ncol = 2 * length(y_var))
      GF_MEAN <- GF_LOADINGS <- list()
      mean_L <- psi_L <- list()
      gf_var_label <- gf_cov_label <- list()
      psi_btw_L <- list()
      for (traj in 1:length(y_var)){
        GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                  arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:2],
                                  labels = paste0(y_var[traj], c("_mueta0", "_mueta1")))
        GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                           to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                           arrows = 1, free = FALSE, values = 1),
                                    mxPath(from = paste0("eta1", y_var[traj]),
                                           to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                           arrows = 1, free = c(F, rep(T, length(records[[traj]]) - 2)),
                                           values = c(1, starts[[4]][[traj]][-1]),
                                           labels = paste0(y_var[traj], "_rel_rate",
                                                           1:(length(records[[traj]]) - 1))))
        mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj], "_mueta1)"),
                                              name = paste0(y_var[traj], "_mean0"))
        psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj], "_psi01), ",
                                                    "cbind(", y_var[traj], "_psi01, ", y_var[traj], "_psi11))"),
                                             name = paste0(y_var[traj], "_psi0"))
        gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi", c("00", "01", "01", "11")),
                                       byrow = T, nrow = 2, ncol = 2)
      }
      for (traj_i in 1:(length(y_var) - 1)){
        for (traj_j in traj_i:(length(y_var) - 1)){
          psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ", y_var[traj_i], y_var[traj_j + 1], "_psi01), ",
                                                                         "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ", y_var[traj_i], y_var[traj_j + 1], "_psi11))"),
                                                                  name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
          gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi", c("00", "01", "10", "11")),
                                                        byrow = T, nrow = 2, ncol = 2)
        }
      }
      for (traj in 1:length(y_var)){
        multi_label[((traj - 1) * 2 + 1):((traj - 1) * 2 + 2), ((traj - 1) * 2 + 1):((traj - 1) * 2 + 2)] <- gf_var_label[[traj]]
      }
      for (traj_i in 1:(length(y_var) - 1)){
        for (traj_j in traj_i:(length(y_var) - 1)){
          multi_label[((traj_i - 1) * 2 + 1):((traj_i - 1) * 2 + 2), (traj_j * 2 + 1):(traj_j * 2 + 2)] <- gf_cov_label[[traj_i + traj_j - 1]]
        }
      }
      model_mx <- mxModel("Multivariate Latent Change Score Model, Nonparametric", type = "RAM",
                          manifestVars = manifests, latentVars = latents,
                          mxData(observed = dat, type = "raw"),
                          mxPath(from = latents[1:(length(y_var) * 2)], to = latents[1:(length(y_var) * 2)],
                                 arrows = 2, connect = "unique.pairs", free = TRUE, values = t(starts[[2]])[row(t(starts[[2]])) >= col(t(starts[[2]]))],
                                 labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                          GF_loadings, AddPara, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L, psi_btw_L,
                          PATH_L_L, PATH_SLP_L, PATH_AUTO_L)
    }
    else if (curveFun %in% c("quadratic", "QUAD")){
      latents <- c(paste0(rep(c("eta0", "eta1", "eta2"), length(y_var)), rep(y_var, each = 3)),
                   unlist(L_list), unlist(D_list))
      multi_label <- matrix(NA, nrow = 3 * length(y_var), ncol = 3 * length(y_var))
      GF_MEAN <- GF_LOADINGS <- list()
      mean_L <- psi_L <- list()
      gf_var_label <- gf_cov_label <- list()
      psi_btw_L <- list()
      for (traj in 1:length(y_var)){
        GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                  arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:3],
                                  labels = paste0(y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
        GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                           to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                           arrows = 1, free = FALSE, values = 1),
                                    mxPath(from = paste0("eta1", y_var[traj]),
                                           to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                           arrows = 1, free = FALSE, values = 1),
                                    mxPath(from = paste0("eta2", y_var[traj]),
                                           to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                           arrows = 1, free = FALSE, values = 0,
                                           labels = paste0(y_var[traj], "L2", records[[traj]][-1], "[1,1]")))
        mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj],
                                                     "_mueta1, ", y_var[traj], "_mueta2)"),
                                              name = paste0(y_var[traj], "_mean0"))
        psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj],
                                                    "_psi01, ", y_var[traj], "_psi02), ",
                                                    "cbind(", y_var[traj], "_psi01, ", y_var[traj],
                                                    "_psi11, ", y_var[traj], "_psi12), ",
                                                    "cbind(", y_var[traj], "_psi02, ", y_var[traj],
                                                    "_psi12, ", y_var[traj], "_psi22))"),
                                             name = paste0(y_var[traj], "_psi0"))
        gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi",
                                              c("00", "01", "02", "01", "11", "12", "02", "12", "22")),
                                       byrow = T, nrow = 3, ncol = 3)
      }
      for (traj_i in 1:(length(y_var) - 1)){
        for (traj_j in traj_i:(length(y_var) - 1)){
          psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi01, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi02),",
                                                                         "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi12), ",
                                                                         "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi20, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi21, ",
                                                                         y_var[traj_i], y_var[traj_j + 1], "_psi22))"),
                                                                  name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
          gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                               c("00", "01", "02", "10", "11", "12", "20", "21", "22")),
                                                        byrow = T, nrow = 3, ncol = 3)
        }
      }
      for (traj in 1:length(y_var)){
        multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <- gf_var_label[[traj]]
      }
      for (traj_i in 1:(length(y_var) - 1)){
        for (traj_j in traj_i:(length(y_var) - 1)){
          multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <- gf_cov_label[[traj_i + traj_j - 1]]
        }
      }
      model_mx <- mxModel("Multivariate Latent Change Score Model, Quadratic", type = "RAM",
                          manifestVars = manifests, latentVars = latents,
                          mxData(observed = dat, type = "raw"),
                          mxPath(from = latents[1:(length(y_var) * 3)], to = latents[1:(length(y_var) * 3)],
                                 arrows = 2, connect = "unique.pairs", free = TRUE, values = t(starts[[2]])[row(t(starts[[2]])) >= col(t(starts[[2]]))],
                                 labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                          GF_loadings, AddPara, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L, psi_btw_L,
                          PATH_L_L, PATH_SLP_L, PATH_AUTO_L)
    }
    else if (curveFun %in% c("negative exponential", "EXP")){
      if (intrinsic){
        latents <- c(paste0(rep(c("eta0", "eta1", "deltag"), length(y_var)), rep(y_var, each = 3)),
                     unlist(L_list), unlist(D_list))
        multi_label <- matrix(NA, nrow = 3 * length(y_var), ncol = 3 * length(y_var))
        GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in 1:length(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:2],
                                    labels = paste0(y_var[traj], c("_mueta0", "_mueta1")))
          GAMMA[[traj]]  <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[traj]][3],
                                     labels = paste0(y_var[traj], "_slp_ratio"),
                                     name = paste0(y_var[traj], "_mug"))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj], "_mueta1, ",
                                                       y_var[traj], "_slp_ratio)"),
                                                name = paste0(y_var[traj], "_mean0"))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                             to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]),
                                             to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0(y_var[traj], "L1", records[[traj]][-1], "[1,1]")),
                                      mxPath(from = paste0("deltag", y_var[traj]),
                                             to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0(y_var[traj], "L2", records[[traj]][-1], "[1,1]")))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj],
                                                      "_psi01, ", y_var[traj], "_psi0g), ",
                                                      "cbind(", y_var[traj], "_psi01, ", y_var[traj],
                                                      "_psi11, ", y_var[traj], "_psi1g), ",
                                                      "cbind(", y_var[traj], "_psi0g, ", y_var[traj],
                                                      "_psi1g, ", y_var[traj], "_psigg))"),
                                               name = paste0(y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi",
                                                c("00", "01", "0g", "01", "11", "1g", "0g", "1g", "gg")),
                                         byrow = T, nrow = 3, ncol = 3)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                           y_var[traj_i], y_var[traj_j + 1],
                                                                           "_psi01, ", y_var[traj_i], y_var[traj_j + 1], "_psi0g),",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi1g), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psig0, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psig1, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psigg))"),
                                                                    name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                 c("00", "01", "0g", "10", "11", "1g", "g0", "g1", "gg")),
                                                          byrow = T, nrow = 3, ncol = 3)
          }
        }
        for (traj in 1:length(y_var)){
          multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
            gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
              gf_cov_label[[traj_i + traj_j - 1]]
          }
        }
        model_mx <- mxModel("Multivariate Latent Change Score Model, Negative Exponential", type = "RAM",
                            manifestVars = manifests, latentVars = latents,
                            mxData(observed = dat, type = "raw"),
                            mxPath(from = latents[1:(length(y_var) * 3)], to = latents[1:(length(y_var) * 3)],
                                   arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = t(starts[[2]])[row(t(starts[[2]])) >= col(t(starts[[2]]))],
                                   labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                            GF_loadings, AddPara, GAMMA, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L,
                            psi_btw_L, PATH_L_L, PATH_SLP_L, PATH_AUTO_L)
      }
      else if (!intrinsic){
        latents <- c(paste0(rep(c("eta0", "eta1"), length(y_var)), rep(y_var, each = 2)),
                     unlist(L_list), unlist(D_list))
        multi_label <- matrix(NA, nrow = 2 * length(y_var), ncol = 2 * length(y_var))
        GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in 1:length(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:2],
                                    labels = paste0(y_var[traj], c("_mueta0", "_mueta1")))
          GAMMA[[traj]]  <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[traj]][3],
                                     labels = paste0(y_var[traj], "_slp_ratio"),
                                     name = paste0(y_var[traj], "_mug"))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj], "_mueta1, ",
                                                       y_var[traj], "_slp_ratio)"),
                                                name = paste0(y_var[traj], "_mean0"))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                             to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]),
                                             to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0(y_var[traj], "L1", records[[traj]][-1], "[1,1]")))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj], "_psi01), ",
                                                      "cbind(", y_var[traj], "_psi01, ", y_var[traj], "_psi11))"),
                                               name = paste0(y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi", c("00", "01", "01", "11")),
                                         byrow = T, nrow = 2, ncol = 2)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi01), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi11))"),
                                                                    name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi", c("00", "01", "10", "11")),
                                                          byrow = T, nrow = 2, ncol = 2)
          }
        }
        for (traj in 1:length(y_var)){
          multi_label[((traj - 1) * 2 + 1):((traj - 1) * 2 + 2), ((traj - 1) * 2 + 1):((traj - 1) * 2 + 2)] <-
            gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 2 + 1):((traj_i - 1) * 2 + 2), (traj_j * 2 + 1):(traj_j * 2 + 2)] <-
              gf_cov_label[[traj_i + traj_j - 1]]
          }
        }
        model_mx <- mxModel("Multivariate Latent Change Score Model, Negative Exponential", type = "RAM",
                            manifestVars = manifests, latentVars = latents,
                            mxData(observed = dat, type = "raw"),
                            mxPath(from = latents[1:(length(y_var) * 2)], to = latents[1:(length(y_var) * 2)],
                                   arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = t(starts[[2]][c(1:2, 4:5), c(1:2, 4:5)])[row(t(starts[[2]][c(1:2, 4:5), c(1:2, 4:5)])) >=
                                                                                       col(t(starts[[2]][c(1:2, 4:5), c(1:2, 4:5)]))],
                                   labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                            GF_loadings, AddPara, GAMMA, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L,
                            psi_btw_L, PATH_L_L, PATH_SLP_L, PATH_AUTO_L)
      }
    }
    else if (curveFun %in% c("Jenss-Bayley", "JB")){
      if (intrinsic){
        latents <- c(paste0(rep(c("eta0", "eta1", "eta2", "deltag"), length(y_var)), rep(y_var, each = 4)),
                     unlist(L_list), unlist(D_list))
        multi_label <- matrix(NA, nrow = 4 * length(y_var), ncol = 4 * length(y_var))
        GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in 1:length(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:3],
                                    labels = paste0(y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
          GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[traj]][4],
                                    labels = paste0(y_var[traj], "_acc_ratio"),
                                    name = paste0(y_var[traj], "_mug"))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj], "_mueta1, ",
                                                       y_var[traj], "_mueta2, ", y_var[traj], "_acc_ratio)"),
                                                name = paste0(y_var[traj], "_mean0"))
          GF_LOADINGS[[traj]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                             to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta1", y_var[traj]),
                                             to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                             arrows = 1, free = FALSE, values = 1),
                                      mxPath(from = paste0("eta2", y_var[traj]),
                                             to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                             arrows = 1, free = FALSE, values = 0,
                                             labels = paste0(y_var[traj], "L2", records[[traj]][-1], "[1,1]")),
                                      mxPath(from = paste0("deltag", y_var[traj]),
                                             to = paste0("d", y_var_L[traj], records[[traj]][-1]), arrows = 1, free = FALSE, values = 0,
                                             labels = paste0(y_var[traj], "L3", records[[traj]][-1], "[1,1]")))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj], "_psi01, ",
                                                      y_var[traj], "_psi02, ", y_var[traj], "_psi0g), ",
                                                      "cbind(", y_var[traj], "_psi01, ", y_var[traj], "_psi11, ",
                                                      y_var[traj], "_psi12, ", y_var[traj], "_psi1g), ",
                                                      "cbind(", y_var[traj], "_psi02, ", y_var[traj], "_psi12, ",
                                                      y_var[traj], "_psi22, ", y_var[traj], "_psi2g), ",
                                                      "cbind(", y_var[traj], "_psi0g, ", y_var[traj], "_psi1g, ",
                                                      y_var[traj], "_psi2g, ", y_var[traj], "_psigg))"),
                                               name = paste0(y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi", c("00", "01", "02", "0g", "01", "11", "12", "1g",
                                                                       "02", "12", "22", "2g", "0g", "1g", "2g", "gg")),
                                         byrow = T, nrow = 4, ncol = 4)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi01, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi02, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi0g), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi12, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi1g), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi20, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi21, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi22, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi2g), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psig0, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psig1, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psig2, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psigg))"),
                                                                    name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                 c("00", "01", "02", "0g", "10", "11", "12", "1g",
                                                                   "20", "21", "22", "2g", "g0", "g1", "g2", "gg")),
                                                          byrow = T, nrow = 4, ncol = 4)
          }
        }
        for (traj in 1:length(y_var)){
          multi_label[((traj - 1) * 4 + 1):((traj - 1) * 4 + 4), ((traj - 1) * 4 + 1):((traj - 1) * 4 + 4)] <-
            gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 4 + 1):((traj_i - 1) * 4 + 4), (traj_j * 4 + 1):(traj_j * 4 + 4)] <-
              gf_cov_label[[traj_i + traj_j - 1]]
          }
        }
        model_mx <- mxModel("Multivariate Latent Change Score Model, Jenss Bayley (random ratio)", type = "RAM",
                            manifestVars = manifests, latentVars = latents,
                            mxData(observed = dat, type = "raw"),
                            mxPath(from = latents[1:(length(y_var) * 4)], to = latents[1:(length(y_var) * 4)],
                                   arrows = 2, connect = "unique.pairs",
                                   free = TRUE, values = t(starts[[2]])[row(t(starts[[2]])) >= col(t(starts[[2]]))],
                                   labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                            GF_loadings, AddPara, GAMMA, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L,
                            psi_btw_L, PATH_L_L, PATH_SLP_L, PATH_AUTO_L)
      }
      else if (!intrinsic){
        latents <- c(paste0(rep(c("eta0", "eta1", "eta2"), length(y_var)), rep(y_var, each = 3)),
                     unlist(L_list), unlist(D_list))
        multi_label <- matrix(NA, nrow = 3 * length(y_var), ncol = 3 * length(y_var))
        GF_MEAN <- GF_LOADINGS <- GAMMA <- list()
        mean_L <- psi_L <- list()
        gf_var_label <- gf_cov_label <- list()
        psi_btw_L <- list()
        for (traj in 1:length(y_var)){
          GF_MEAN[[traj]] <- mxPath(from = "one", to = paste0(c("eta0", "eta1", "eta2"), y_var[traj]),
                                    arrows = 1, free = TRUE, values = starts[[1]][[traj]][1:3],
                                    labels = paste0(y_var[traj], c("_mueta0", "_mueta1", "_mueta2")))
          GAMMA[[traj]] <- mxMatrix("Full", 1, 1, free = TRUE, values = starts[[1]][[traj]][4],
                                    labels = paste0(y_var[traj], "_acc_ratio"),
                                    name = paste0(y_var[traj], "_mug"))
          mean_L[[traj]] <- mxAlgebraFromString(paste0("rbind(", y_var[traj], "_mueta0, ", y_var[traj], "_mueta1, ",
                                                       y_var[traj], "_mueta2, ", y_var[traj], "_acc_ratio)"),
                                                name = paste0(y_var[traj], "_mean0"))
          GF_LOADINGS[[length(GF_LOADINGS) + 1]] <- list(mxPath(from = paste0("eta0", y_var[traj]),
                                                                to = paste0("l", y_var_L[traj], records[[traj]][1]),
                                                                arrows = 1, free = FALSE, values = 1),
                                                         mxPath(from = paste0("eta1", y_var[traj]),
                                                                to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                                                arrows = 1, free = FALSE, values = 1),
                                                         mxPath(from = paste0("eta2", y_var[traj]),
                                                                to = paste0("d", y_var_L[traj], records[[traj]][-1]),
                                                                arrows = 1, free = FALSE, values = 0,
                                                                labels = paste0(y_var[traj], "L2", records[[traj]][-1], "[1,1]")))
          psi_L[[traj]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj], "_psi00, ", y_var[traj], "_psi01, ",
                                                      y_var[traj], "_psi02), ",
                                                      "cbind(", y_var[traj], "_psi01, ", y_var[traj], "_psi11, ",
                                                      y_var[traj], "_psi12), ",
                                                      "cbind(", y_var[traj], "_psi02, ", y_var[traj], "_psi12, ",
                                                      y_var[traj], "_psi22))"),
                                               name = paste0(y_var[traj], "_psi0"))
          gf_var_label[[traj]] <- matrix(paste0(y_var[traj], "_psi", c("00", "01", "02", "01", "11", "12",
                                                                       "02", "12", "22")),
                                         byrow = T, nrow = 3, ncol = 3)
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            psi_btw_L[[traj_i + traj_j - 1]] <- mxAlgebraFromString(paste0("rbind(cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi00, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi01, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi02), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi10, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi11, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi12), ",
                                                                           "cbind(", y_var[traj_i], y_var[traj_j + 1], "_psi20, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi21, ",
                                                                           y_var[traj_i], y_var[traj_j + 1], "_psi22))"),
                                                                    name = paste0(y_var[traj_i], y_var[traj_j + 1], "_psi"))
            gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0(y_var[traj_i], y_var[traj_j + 1], "_psi",
                                                                 c("00", "01", "02", "10", "11", "12",
                                                                   "20", "21", "22")),
                                                          byrow = T, nrow = 3, ncol = 3)
          }
        }
        for (traj in 1:length(y_var)){
          multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
            gf_var_label[[traj]]
        }
        for (traj_i in 1:(length(y_var) - 1)){
          for (traj_j in traj_i:(length(y_var) - 1)){
            multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
              gf_cov_label[[traj_i + traj_j - 1]]
          }
        }
        model_mx <- mxModel("Multivariate Latent Change Score Model, Jenss Bayley (fixed ratio)", type = "RAM",
                            manifestVars = manifests, latentVars = latents,
                            mxData(observed = dat, type = "raw"),
                            mxPath(from = latents[1:(length(y_var) * 3)], to = latents[1:(length(y_var) * 3)],
                                   arrows = 2, connect = "unique.pairs", free = TRUE,
                                   values = t(starts[[2]][c(1:3, 5:7), c(1:3, 5:7)])[row(t(starts[[2]][c(1:3, 5:7), c(1:3, 5:7)])) >=
                                                                                       col(t(starts[[2]][c(1:3, 5:7), c(1:3, 5:7)]))],
                                   labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                            GF_loadings, AddPara, GAMMA, GF_MEAN, GF_LOADINGS, RES, COV, mean_L, psi_L,
                            psi_btw_L, PATH_L_L, PATH_SLP_L, PATH_AUTO_L)
      }
    }
  }
  return(model_mx)
}

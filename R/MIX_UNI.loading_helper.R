#' @title Get Factor Loadings for a Mixture Model or Multiple Group Model with Univariate Longitudinal Outcome
#'
#' @description This function specifies the factor loadings for a mixture model with a univariate longitudinal outcome. The longitudinal
#' outcome is fit by a Latent Growth Curve Model or a Latent Change Score Model.
#'
#' @param nClass An integer specifying the number of classes for the mixture model or multiple group model. It takes the value
#' passed from \code{getMIX()} or \code{getMGroup()}.
#' @param y_model A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and
#' \code{"LCSM"}. It takes the value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable at each study wave.
#' It takes the value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the
#' value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param intrinsic A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. It takes the
#' value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param records A numeric vector specifying indices of the study waves. It takes the value passed from \code{getMIX()} or
#' \code{getMGroup()}.
#'
#' @return A list containing the specification of definition variables (i.e., individual time points for the latent
#' growth curve models, and individual time points and individual time lags (intervals) between adjacent time points
#' for latent change score models) and factor loadings of a univariate longitudinal outcome.
#'
#' @keywords internal
#'
getMIX_UNI.loadings <- function(nClass, y_model, t_var, records, y_var, curveFun, intrinsic){
  # Define `outPoint` to store individual time points from raw data
  output_CL <- list()
  for (k in 1:nClass){
    outPoint <- list()
    for (j in records){
      outPoint[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var, j),
                                name = paste0("t", j))
    }
    if (y_model == "LGCM"){
      if (curveFun %in% c("linear", "LIN")){
        outLoads1 <- list()
        for(j in records){
          outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j), name = paste0("c", k, "L1", j))
        }
        outLoads <- list(outLoads1)
      }
      else if (curveFun %in% c("quadratic", "QUAD")){
        outLoads1 <- outLoads2 <- list()
        for(j in records){
          outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j), name = paste0("c", k, "L1", j))
          outLoads2[[j]] <- mxAlgebraFromString(paste0("t", j, "^2"), name = paste0("c", k, "L2", j))
        }
        outLoads <- list(outLoads1, outLoads2)
      }
      else if (curveFun %in% c("negative exponential", "EXP")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("1 - exp(-c", k, "Y_mug * t", j, ")"),
                                                  name = paste0("c", k, "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("c", k, "Y_mueta1 * exp(-c", k, "Y_mug * t", j, ") * t", j),
                                                  name = paste0("c", k, "L2", j))
          }
          outLoads <- list(outLoads1, outLoads2)
        }
        else if (!intrinsic){
          outLoads1 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("1 - exp(-c", k, "Y_mug * t", j, ")"),
                                                  name = paste0("c", k, "L1", j))
          }
          outLoads <- list(outLoads1)
        }
      }
      else if (curveFun %in% c("Jenss-Bayley", "JB")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- outLoads3 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j), name = paste0("c", k, "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("exp(c", k, "Y_mug * t", j, ") - 1"),
                                                  name = paste0("c", k, "L2", j))
            outLoads3[[j]] <- mxAlgebraFromString(paste0("c", k, "Y_mueta2 * exp(c", k, "Y_mug * t", j, ") * t", j),
                                                  name = paste0("c", k, "L3", j))
          }
          outLoads <- list(outLoads1, outLoads2, outLoads3)
        }
        else if (!intrinsic){
          outLoads1 <- outLoads2 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j), name = paste0("c", k, "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("exp(c", k, "Y_mug * t", j, ") - 1"),
                                                  name = paste0("c", k, "L2", j))
          }
          outLoads <- list(outLoads1, outLoads2)
        }
      }
      else if (curveFun %in% c("bilinear spline", "BLS")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- outLoads3 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " - c", k, "Y_mug"),
                                                  name = paste0("c", k, "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " - c", k, "Y_mug)"),
                                                  name = paste0("c", k, "L2", j))
            outLoads3[[j]] <- mxAlgebraFromString(paste0("-c", k, "Y_mueta2s * (t", j, " - c", k, "Y_mug)/abs(t", j,
                                                         " - c", k, "Y_mug) - c", k, "Y_mueta2s"),
                                                  name = paste0("c", k, "L3", j))
          }
          outLoads <- list(outLoads1, outLoads2, outLoads3)
        }
        else if (!intrinsic){
          outLoads1 <- outLoads2 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " - c", k, "Y_mug"),
                                                  name = paste0("c", k, "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " - c", k, "Y_mug)"),
                                                  name = paste0("c", k, "L2", j))
          }
          outLoads <- list(outLoads1, outLoads2)
        }
      }
      output_CL[[k]] <- list(outPoint, outLoads)
    }
    else if (y_model == "LCSM"){
      # Define `outLag` to store individual time lags from raw data
      outLag <- list()
      outLag[[1]] <- mxAlgebra(0, name = paste0("lag", 1))
      for(j in records[-1]){
        outLag[[j]] <- mxAlgebraFromString(paste0("t", j , " -  t", j - 1), name = paste0("lag", j))
      }
      # Define `midTime` to store individual midway time point through each interval from raw data, which help
      # define factor loadings of latent change score models.
      midTime <- list()
      midTime[[1]] <- mxAlgebra(0, name = paste0("mid_t", 1))
      for(j in records[-1]){
        midTime[[j]] <- mxAlgebraFromString(paste0("(t", j , " +  t", j - 1, ")/2"),
                                            name = paste0("mid_t", j))
      }
      if (curveFun %in% c("nonparametric", "NonP")){
        outLoads <- list()
      }
      else if (curveFun %in% c("quadratic", "QUAD")){
        outLoads2 <- list()
        outLoads2[[1]] <- mxAlgebra(0, name = paste0("c", k, "L2", 1))
        for(j in records[-1]){
          outLoads2[[j]] <- mxAlgebraFromString(paste0("2 * mid_t", j), name = paste0("c", k, "L2", j))
        }
        outLoads <- outLoads2
      }
      else if (curveFun %in% c("negative exponential", "EXP")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- list()
          outLoads1[[1]] <- mxAlgebra(0, name = paste0("c", k, "L1", 1))
          outLoads2[[1]] <- mxAlgebra(0, name = paste0("c", k, "L2", 1))
          for (j in records[-1]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("c", k, "Y_mug * exp(-c", k, "Y_mug * mid_t", j, ")"),
                                                  name = paste0("c", k, "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("c", k, "Y_mueta1 * exp(-c", k, "Y_mug * mid_t", j,
                                                         ") * (1 - c", k, "Y_mug * mid_t", j, ")"),
                                                  name = paste0("c", k, "L2", j))
          }
          outLoads <- list(outLoads1, outLoads2)
        }
        else if (!intrinsic){
          outLoads1 <- list()
          outLoads1[[1]] <- mxAlgebra(0, name = paste0("c", k, "L1", 1))
          for (j in records[-1]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("c", k, "Y_mug * exp(-c", k, "Y_mug * mid_t", j, ")"),
                                                  name = paste0("c", k, "L1", j))
          }
          outLoads <- list(outLoads1)
        }
      }
      else if (curveFun %in% c("Jenss-Bayley", "JB")){
        if (intrinsic){
          outLoads2 <- outLoads3 <- list()
          outLoads2[[1]] <- mxAlgebra(0, name = paste0("c", k, "L2", 1))
          outLoads3[[1]] <- mxAlgebra(0, name = paste0("c", k, "L3", 1))
          for (j in records[-1]){
            outLoads2[[j]] <- mxAlgebraFromString(paste0("c", k, "Y_mug * exp(c", k, "Y_mug * mid_t", j, ")"),
                                                  name = paste0("c", k, "L2", j))
            outLoads3[[j]] <- mxAlgebraFromString(paste0("c", k, "Y_mueta2 * (1 + c", k, "Y_mug *  mid_t", j, ")",
                                                         " * exp(c", k, "Y_mug * mid_t", j, ")"),
                                                  name = paste0("c", k, "L3", j))
          }
          outLoads <- list(outLoads2, outLoads3)
        }
        else if (!intrinsic){
          outLoads2 <- list()
          outLoads2[[1]] <- mxAlgebra(0, name = paste0("c", k, "L2", 1))
          for (j in records[-1]){
            outLoads2[[j]] <- mxAlgebraFromString(paste0("c", k, "Y_mug * exp(c", k, "Y_mug * mid_t", j, ")"),
                                                  name = paste0("c", k, "L2", j))
          }
          outLoads <- list(outLoads2)
        }
      }
      output_CL[[k]] <- list(outPoint, outLag, midTime, outLoads)
    }
  }
  return(output_CL)
}


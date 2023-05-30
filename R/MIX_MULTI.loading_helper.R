#' @title Get Factor Loadings for a Mixture Model with MGM as Submodels
#'
#' @description This function specifies the factor loadings for a mixture model with MGM as submodels. The longitudinal outcomes are fit
#' by Latent Growth Curve Models or a Latent Change Score Models.
#'
#' @param nClass A numeric value to indicate the number of latent classes.
#' @param y_model A character string specifying how to fit the longitudinal outcome. Supported values are "LGCM" and
#' "LCSM". It takes the value passed from `getMIX()`.
#' @param t_var t_var A character vector where each element represents a prefix for column names corresponding to the
#' time variables for the respective longitudinal process. It takes the value passed from `getMIX()`.
#' @param y_var A character vector where each element represents a prefix for column names corresponding to the outcome
#' variables for the respective longitudinal process. It takes the value passed from `getMIX()`.
#' @param curveFun The specified functional form of the growth curve. Supported options include: "linear" (or "LIN"),
#' "quadratic" (or "QUAD"), "negative exponential" (or "EXP"), "Jenss-Bayley" (or "JB"), "bilinear spline" (or "BLS"),
#' and "nonparametric" (or "NonP"). It takes the value passed from `getMIX()`.
#' @param intrinsic A boolean flag for whether an intrinsically nonlinear longitudinal model is built up. It takes the
#' value passed from `getMIX()`.
#' @param records A list of numeric vectors where each vector contains indices of the observed time points for the
#' respective longitudinal processes. It takes the value passed from `getMIX()`.
#'
#' @return A list containing the specification of definition variables (i.e., individual measurement occasions) and factor
#' loadings of multivariate longitudinal outcomes.
#'
getMIX_MULTI.loadings <- function(nClass, y_model, t_var, y_var, curveFun, intrinsic = NULL, records){
  output_CL <- list()
  for (k in 1:nClass){
    outPoint_L <- outLag_L <- midTime_L <- outLoads_L <- list()
    for (traj in 1:length(y_var)){
      outPoint <- list()
      for (j in records[[traj]]){
        outPoint[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[traj], j),
                                  name = paste0(y_var[traj], "t", j))
      }
      outPoint_L[[traj]] <- outPoint
      if (y_model == "LGCM"){
        if (curveFun %in% c("linear", "LIN")){
          outLoads1 <- list()
          for (j in records[[traj]]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                  name = paste0("c", k, "L1", j, y_var[traj]))
          }
          outLoads_L[[traj]] <- list(outLoads1)
        }
        else if (curveFun %in% c("quadratic", "QUAD")){
          outLoads1 <- outLoads2 <- list()
          for(j in records[[traj]]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                  name = paste0("c", k, "L1", j, y_var[traj]))
            outLoads2[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j, "^2"),
                                                  name = paste0("c", k, "L2", j, y_var[traj]))
          }
          outLoads_L[[traj]] <- list(outLoads1, outLoads2)
        }
        else if (curveFun %in% c("negative exponential", "EXP")){
          if (intrinsic){
            outLoads1 <- outLoads2 <- list()
            for (j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0("1 - exp(-c", k, y_var[traj], "_mug",
                                                           " * ", y_var[traj], "t", j, ")"),
                                                    name = paste0("c", k, "L1", j, y_var[traj]))
              outLoads2[[j]] <- mxAlgebraFromString(paste0("c", k, y_var[traj], "_mueta1 * exp(-c", k, y_var[traj], "_mug",
                                                           " * ", y_var[traj], "t", j, ") * ", y_var[traj], "t", j),
                                                    name = paste0("c", k, "L2", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2)
          }
          else if (!intrinsic){
            outLoads1 <- list()
            for (j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0("1 - exp(-c", k, y_var[traj], "_mug",
                                                           " * ", y_var[traj], "t", j, ")"),
                                                    name = paste0("c", k, "L1", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1)
          }
        }
        else if (curveFun %in% c("Jenss-Bayley", "JB")){
          if (intrinsic){
            outLoads1 <- outLoads2 <- outLoads3 <- list()
            for(j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                    name = paste0("c", k, "L1", j, y_var[traj]))
              outLoads2[[j]] <- mxAlgebraFromString(paste0("exp(c", k, y_var[traj], "_mug * ", y_var[traj], "t", j, ") - 1"),
                                                    name = paste0("c", k, "L2", j, y_var[traj]))
              outLoads3[[j]] <- mxAlgebraFromString(paste0("c", k, y_var[traj], "_mueta2 * exp(c", k, y_var[traj], "_mug * ",
                                                           y_var[traj], "t", j, ") * ", y_var[traj], "t", j),
                                                    name = paste0("c", k, "L3", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2, outLoads3)
          }
          else if (!intrinsic){
            outLoads1 <- outLoads2 <- list()
            for(j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                    name = paste0("c", k, "L1", j, y_var[traj]))
              outLoads2[[j]] <- mxAlgebraFromString(paste0("exp(c", k, y_var[traj], "_mug * ", y_var[traj], "t", j, ") - 1"),
                                                    name = paste0("c", k, "L2", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2)
          }
        }
        else if (curveFun %in% c("bilinear spline", "BLS")){
          if (intrinsic){
            outLoads1 <- outLoads2 <- outLoads3 <- list()
            for(j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j, " - c", k, y_var[traj], "_mug"),
                                                    name = paste0("c", k, "L1", j, y_var[traj]))
              outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(", y_var[traj], "t", j, " - c", k, y_var[traj], "_mug)"),
                                                    name = paste0("c", k, "L2", j, y_var[traj]))
              outLoads3[[j]] <- mxAlgebraFromString(paste0("-c", k, y_var[traj], "_mueta2s * (", y_var[traj], "t", j, " - c", k, y_var[traj],
                                                           "_mug)/abs(", y_var[traj], "t", j,
                                                           " - c", k, y_var[traj], "_mug) - c", k, y_var[traj], "_mueta2s"),
                                                    name = paste0("c", k, "L3", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2, outLoads3)
          }
          else if (!intrinsic){
            outLoads1 <- outLoads2 <- list()
            for(j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j, " - c", k, y_var[traj], "_mug"),
                                                    name = paste0("c", k, "L1", j, y_var[traj]))
              outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(", y_var[traj], "t", j, " - c", k, y_var[traj], "_mug)"),
                                                    name = paste0("c", k, "L2", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2)
          }
        }
      }
      else if (y_model == "LCSM"){
        outLag <- list()
        for(j in records[[traj]][-1]){
          outLag[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j , " - ", y_var[traj], "t", j - 1),
                                             name = paste0(y_var[traj], "lag", j))
        }
        outLag_L[[traj]] <- outLag
        # Define `midTime` to store individual midway time point through each interval from raw data
        midTime <- list()
        midTime[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "mid_t", 1))
        for(j in records[[traj]][-1]){
          midTime[[j]] <- mxAlgebraFromString(paste0("(", y_var[traj], "t", j , " + ", y_var[traj], "t", j - 1, ")/2"),
                                              name = paste0(y_var[traj], "mid_t", j))
        }
        midTime_L[[traj]] <- midTime
        if (curveFun %in% c("nonparametric", "NonP")){
          outLoads_L[[traj]] <- list()
        }
        else if (curveFun %in% c("quadratic", "QUAD")){
          outLoads2 <- list()
          outLoads2[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L2", 1))
          for(j in records[[traj]][-1]){
            outLoads2[[j]] <- mxAlgebraFromString(paste0("2 * ", y_var[traj], "mid_t", j),
                                                  name = paste0("c", k, y_var[traj], "L2", j))
          }
          outLoads_L[[traj]] <- outLoads2
        }
        else if (curveFun %in% c("negative exponential", "EXP")){
          if (intrinsic){
            outLoads1 <- outLoads2 <- list()
            outLoads1[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L1", 1))
            outLoads2[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L2", 1))
            for (j in records[[traj]][-1]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0("c", k, y_var[traj], "_mug * exp(-c", k, y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0("c", k, y_var[traj], "L1", j))
              outLoads2[[j]] <- mxAlgebraFromString(paste0("c", k, y_var[traj], "_mueta1 * exp(-c", k, y_var[traj], "_mug * ", y_var[traj], "mid_t", j,
                                                           ") * (1 - c", k, y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0("c", k, y_var[traj], "L2", j))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2)
          }
          else if (!intrinsic){
            outLoads1 <- list()
            outLoads1[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L1", 1))
            for (j in records[[traj]][-1]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0("c", k, y_var[traj], "_mug * exp(-c", k, y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0("c", k, y_var[traj], "L1", j))
            }
            outLoads_L[[traj]] <- list(outLoads1)
          }
        }
        else if (curveFun %in% c("Jenss-Bayley", "JB")){
          if (intrinsic){
            outLoads2 <- outLoads3 <- list()
            outLoads2[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L2", 1))
            outLoads3[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L3", 1))
            for (j in records[[traj]][-1]){
              outLoads2[[j]] <- mxAlgebraFromString(paste0("c", k, y_var[traj], "_mug * exp(c", k, y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0("c", k, y_var[traj], "L2", j))
              outLoads3[[j]] <- mxAlgebraFromString(paste0("c", k, y_var[traj], "_mueta2 * (1 + c", k, y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")",
                                                           " * exp(c", k, y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0("c", k, y_var[traj], "L3", j))
            }
            outLoads_L[[traj]] <- list(outLoads2, outLoads3)
          }
          else if (!intrinsic){
            outLoads2 <- list()
            outLoads2[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L2", 1))
            for (j in records[[traj]][-1]){
              outLoads2[[j]] <- mxAlgebraFromString(paste0("c", k, y_var[traj], "_mug * exp(c", k, y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0("c", k, y_var[traj], "L2", j))
            }
            outLoads_L[[traj]] <- list(outLoads2)
          }
        }
      }
    }
    if (y_model == "LGCM"){
      output_CL[[k]] <- list(outPoint_L, outLoads_L)
    }
    else if (y_model == "LCSM"){
      output_CL[[k]] <- list(outPoint_L, outLag_L, midTime_L, outLoads_L)
    }
  }
  return(output_CL)
}

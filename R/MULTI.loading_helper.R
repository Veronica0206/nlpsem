#' @title Get Factor Loadings for a Multivariate Longitudinal Outcomes with Specified Functional Curves
#'
#' @description This function specifies the factor loadings for a multivariate longitudinal outcomes with a given functional form.
#' The longitudinal outcomes are fit by Latent Growth Curve Models or a Latent Change Score Models.
#'
#' @param y_model A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' It takes the value passed from \code{getMGM()}.
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
#'
#' @return A list containing the specification of definition variables (i.e., individual time points for the latent
#' growth curve models, and individual time points and individual time lags (intervals) between adjacent time points
#' for latent change score models) and factor loadings of a multivariate longitudinal outcomes.
#'
#' @keywords internal
#'
#' @importFrom OpenMx mxMatrix mxAlgebraFromString mxAlgebra
#'
getMULTI.loadings <- function(y_model, t_var, y_var, curveFun, intrinsic, records){
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
                                                name = paste0("L1", j, y_var[traj]))
        }
        outLoads_L[[traj]] <- list(outLoads1)
      }
      else if (curveFun %in% c("quadratic", "QUAD")){
        outLoads1 <- outLoads2 <- list()
        for(j in records[[traj]]){
          outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                name = paste0("L1", j, y_var[traj]))
          outLoads2[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j, "^2"),
                                                name = paste0("L2", j, y_var[traj]))
        }
        outLoads_L[[traj]] <- list(outLoads1, outLoads2)
      }
      else if (curveFun %in% c("negative exponential", "EXP")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- list()
          for (j in records[[traj]]){
            outPoint[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[traj], j),
                                    name = paste0(y_var[traj], "t", j))
            outLoads1[[j]] <- mxAlgebraFromString(paste0("1 - exp(-", y_var[traj], "_mug",
                                                         " * ", y_var[traj], "t", j, ")"),
                                                  name = paste0("L1", j, y_var[traj]))
            outLoads2[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "_mueta1 * exp(-", y_var[traj], "_mug",
                                                         " * ", y_var[traj], "t", j, ") * ", y_var[traj], "t", j),
                                                  name = paste0("L2", j, y_var[traj]))
          }
          outLoads_L[[traj]] <- list(outLoads1, outLoads2)
        }
        else if (!intrinsic){
          outLoads1 <- list()
          for (j in records[[traj]]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("1 - exp(-", y_var[traj], "_mug",
                                                         " * ", y_var[traj], "t", j, ")"),
                                                  name = paste0("L1", j, y_var[traj]))
          }
          outLoads_L[[traj]] <- list(outLoads1)
        }
      }
      else if (curveFun %in% c("Jenss-Bayley", "JB")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- outLoads3 <- list()
          for(j in records[[traj]]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                  name = paste0("L1", j, y_var[traj]))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("exp(", y_var[traj], "_mug * ", y_var[traj], "t", j, ") - 1"),
                                                  name = paste0("L2", j, y_var[traj]))
            outLoads3[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "_mueta2 * exp(", y_var[traj], "_mug * ",
                                                         y_var[traj], "t", j, ") * ", y_var[traj], "t", j),
                                                  name = paste0("L3", j, y_var[traj]))
          }
          outLoads_L[[traj]] <- list(outLoads1, outLoads2, outLoads3)
        }
        else if (!intrinsic){
          outLoads1 <- outLoads2 <- list()
          for(j in records[[traj]]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                  name = paste0("L1", j, y_var[traj]))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("exp(", y_var[traj], "_mug * ", y_var[traj], "t", j, ") - 1"),
                                                  name = paste0("L2", j, y_var[traj]))
          }
          outLoads_L[[traj]] <- list(outLoads1, outLoads2)
        }
      }
      else if (curveFun %in% c("bilinear spline", "BLS")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- outLoads3 <- list()
          for(j in records[[traj]]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j, " - ", y_var[traj], "_mug"),
                                                  name = paste0("L1", j, y_var[traj]))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(", y_var[traj], "t", j, " - ", y_var[traj], "_mug)"),
                                                  name = paste0("L2", j, y_var[traj]))
            outLoads3[[j]] <- mxAlgebraFromString(paste0("-", y_var[traj], "_mueta2s * (", y_var[traj], "t", j, " - ", y_var[traj],
                                                         "_mug)/abs(", y_var[traj], "t", j,
                                                         " - ", y_var[traj], "_mug) - ", y_var[traj], "_mueta2s"),
                                                  name = paste0("L3", j, y_var[traj]))
          }
          outLoads_L[[traj]] <- list(outLoads1, outLoads2, outLoads3)
        }
        else if (!intrinsic){
          outLoads1 <- outLoads2 <- list()
          for(j in records[[traj]]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j, " - ", y_var[traj], "_mug"),
                                                  name = paste0("L1", j, y_var[traj]))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(", y_var[traj], "t", j, " - ", y_var[traj], "_mug)"),
                                                  name = paste0("L2", j, y_var[traj]))
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
                                                name = paste0(y_var[traj], "L2", j))
        }
        outLoads_L[[traj]] <- outLoads2
      }
      else if (curveFun %in% c("negative exponential", "EXP")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- list()
          outLoads1[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L1", 1))
          outLoads2[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L2", 1))
          for (j in records[[traj]][-1]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "_mug * exp(-", y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                  name = paste0(y_var[traj], "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "_mueta1 * exp(-", y_var[traj], "_mug * ", y_var[traj], "mid_t", j,
                                                         ") * (1 - ", y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                  name = paste0(y_var[traj], "L2", j))
          }
          outLoads_L[[traj]] <- list(outLoads1, outLoads2)
        }
        else if (!intrinsic){
          outLoads1 <- list()
          outLoads1[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L1", 1))
          for (j in records[[traj]][-1]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "_mug * exp(-", y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                  name = paste0(y_var[traj], "L1", j))
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
            outLoads2[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "_mug * exp(", y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                  name = paste0(y_var[traj], "L2", j))
            outLoads3[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "_mueta2 * (1 + ", y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")",
                                                         " * exp(", y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                  name = paste0(y_var[traj], "L3", j))
          }
          outLoads_L[[traj]] <- list(outLoads2, outLoads3)
        }
        else if (!intrinsic){
          outLoads2 <- list()
          outLoads2[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L2", 1))
          for (j in records[[traj]][-1]){
            outLoads2[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "_mug * exp(", y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                  name = paste0(y_var[traj], "L2", j))
          }
          outLoads_L[[traj]] <- list(outLoads2)
        }
      }
    }
  }
  if (y_model == "LGCM"){
    return(list(outPoint_L, outLoads_L))
  }
  else if (y_model == "LCSM"){
    output <- list(outPoint_L, outLag_L, midTime_L, outLoads_L)
  }
}

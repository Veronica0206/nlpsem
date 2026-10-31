#' @title Growth-factor loadings of multivariate outcomes, one set per prefix
#' @description \code{pfx = ""} gives the single-group objects; \code{pfx = c("c1", "c2", ...)} the class-specific
#' ones used by the multiple-group and mixture models. The names of matrices and algebras carry the prefix.
#' @keywords internal
#' @noRd
#'
#' @importFrom OpenMx mxMatrix mxAlgebraFromString mxAlgebra
#'
.loadings_MULTI_classes <- function(pfx, y_model, t_var, y_var, curveFun, intrinsic = NULL, records){
  output_CL <- list()
  for (k in seq_along(pfx)){
    outPoint_L <- outLag_L <- midTime_L <- outLoads_L <- list()
    for (traj in seq_along(y_var)){
      outPoint <- list()
      for (j in records[[traj]]){
        outPoint[[j]] <- mxMatrix("Full", 1, 1, free = FALSE, labels = paste0("data.", t_var[traj], j),
                                  name = paste0(y_var[traj], "t", j))
      }
      outPoint_L[[traj]] <- outPoint
      if (y_model == "LGCM"){
        if (curveFun %in% c("linear", "LIN")){
          outLoads1 <- list()
          for (j in records[[traj]]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                  name = paste0(pfx[k], "L1", j, y_var[traj]))
          }
          outLoads_L[[traj]] <- list(outLoads1)
        }
        else if (curveFun %in% c("quadratic", "QUAD")){
          outLoads1 <- outLoads2 <- list()
          for(j in records[[traj]]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                  name = paste0(pfx[k], "L1", j, y_var[traj]))
            outLoads2[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j, "^2"),
                                                  name = paste0(pfx[k], "L2", j, y_var[traj]))
          }
          outLoads_L[[traj]] <- list(outLoads1, outLoads2)
        }
        else if (curveFun %in% c("negative exponential", "EXP")){
          if (intrinsic){
            outLoads1 <- outLoads2 <- list()
            for (j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0("1 - exp(-", pfx[k], y_var[traj], "_mug",
                                                           " * ", y_var[traj], "t", j, ")"),
                                                    name = paste0(pfx[k], "L1", j, y_var[traj]))
              outLoads2[[j]] <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "_mueta1 * exp(-", pfx[k], y_var[traj], "_mug",
                                                           " * ", y_var[traj], "t", j, ") * ", y_var[traj], "t", j),
                                                    name = paste0(pfx[k], "L2", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2)
          }
          else if (!intrinsic){
            outLoads1 <- list()
            for (j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0("1 - exp(-", pfx[k], y_var[traj], "_mug",
                                                           " * ", y_var[traj], "t", j, ")"),
                                                    name = paste0(pfx[k], "L1", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1)
          }
        }
        else if (curveFun %in% c("Jenss-Bayley", "JB")){
          if (intrinsic){
            outLoads1 <- outLoads2 <- outLoads3 <- list()
            for(j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                    name = paste0(pfx[k], "L1", j, y_var[traj]))
              outLoads2[[j]] <- mxAlgebraFromString(paste0("exp(", pfx[k], y_var[traj], "_mug * ", y_var[traj], "t", j, ") - 1"),
                                                    name = paste0(pfx[k], "L2", j, y_var[traj]))
              outLoads3[[j]] <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "_mueta2 * exp(", pfx[k], y_var[traj], "_mug * ",
                                                           y_var[traj], "t", j, ") * ", y_var[traj], "t", j),
                                                    name = paste0(pfx[k], "L3", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2, outLoads3)
          }
          else if (!intrinsic){
            outLoads1 <- outLoads2 <- list()
            for(j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j),
                                                    name = paste0(pfx[k], "L1", j, y_var[traj]))
              outLoads2[[j]] <- mxAlgebraFromString(paste0("exp(", pfx[k], y_var[traj], "_mug * ", y_var[traj], "t", j, ") - 1"),
                                                    name = paste0(pfx[k], "L2", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2)
          }
        }
        else if (curveFun %in% c("bilinear spline", "BLS")){
          if (intrinsic){
            outLoads1 <- outLoads2 <- outLoads3 <- list()
            for(j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j, " - ", pfx[k], y_var[traj], "_mug"),
                                                    name = paste0(pfx[k], "L1", j, y_var[traj]))
              outLoads2[[j]] <- mxAlgebraFromString(paste0("sqrt((", y_var[traj], "t", j, " - ", pfx[k], y_var[traj], "_mug)^2 + 1e-6)"),
                                                    name = paste0(pfx[k], "L2", j, y_var[traj]))
              outLoads3[[j]] <- mxAlgebraFromString(paste0("-", pfx[k], y_var[traj], "_mueta2s * (", y_var[traj], "t", j, " - ", pfx[k], y_var[traj],
                                                           "_mug)/sqrt((", y_var[traj], "t", j, " - ", pfx[k], y_var[traj], "_mug)^2 + 1e-6) - ", pfx[k], y_var[traj], "_mueta1s"),
                                                    name = paste0(pfx[k], "L3", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2, outLoads3)
          }
          else if (!intrinsic){
            outLoads1 <- outLoads2 <- list()
            for(j in records[[traj]]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j, " - ", pfx[k], y_var[traj], "_mug"),
                                                    name = paste0(pfx[k], "L1", j, y_var[traj]))
              outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(", y_var[traj], "t", j, " - ", pfx[k], y_var[traj], "_mug)"),
                                                    name = paste0(pfx[k], "L2", j, y_var[traj]))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2)
          }
        }
      }
      else if (y_model == "LCSM"){
        outLag <- list()
        for (pos in seq_along(records[[traj]])[-1]){
          j <- records[[traj]][pos]
          outLag[[j]] <- mxAlgebraFromString(paste0(y_var[traj], "t", j , " - ", y_var[traj], "t", records[[traj]][pos - 1]),
                                             name = paste0(y_var[traj], "lag", j))
        }
        outLag_L[[traj]] <- outLag
        # Define `midTime` to store individual midway time point through each interval from raw data
        midTime <- list()
        midTime[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "mid_t", records[[traj]][1]))
        for (pos in seq_along(records[[traj]])[-1]){
          j <- records[[traj]][pos]
          midTime[[j]] <- mxAlgebraFromString(paste0("(", y_var[traj], "t", j , " + ", y_var[traj], "t", records[[traj]][pos - 1], ")/2"),
                                              name = paste0(y_var[traj], "mid_t", j))
        }
        midTime_L[[traj]] <- midTime
        if (curveFun %in% c("nonparametric", "NonP")){
          outLoads_L[[traj]] <- list()
        }
        else if (curveFun %in% c("quadratic", "QUAD")){
          outLoads2 <- list()
          outLoads2[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L2", records[[traj]][1]))
          for(j in records[[traj]][-1]){
            outLoads2[[j]] <- mxAlgebraFromString(paste0("2 * ", y_var[traj], "mid_t", j),
                                                  name = paste0(pfx[k], y_var[traj], "L2", j))
          }
          outLoads_L[[traj]] <- outLoads2
        }
        else if (curveFun %in% c("negative exponential", "EXP")){
          if (intrinsic){
            outLoads1 <- outLoads2 <- list()
            outLoads1[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L1", records[[traj]][1]))
            outLoads2[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L2", records[[traj]][1]))
            for (j in records[[traj]][-1]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "_mug * exp(-", pfx[k], y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0(pfx[k], y_var[traj], "L1", j))
              outLoads2[[j]] <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "_mueta1 * exp(-", pfx[k], y_var[traj], "_mug * ", y_var[traj], "mid_t", j,
                                                           ") * (1 - ", pfx[k], y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0(pfx[k], y_var[traj], "L2", j))
            }
            outLoads_L[[traj]] <- list(outLoads1, outLoads2)
          }
          else if (!intrinsic){
            outLoads1 <- list()
            outLoads1[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L1", records[[traj]][1]))
            for (j in records[[traj]][-1]){
              outLoads1[[j]] <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "_mug * exp(-", pfx[k], y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0(pfx[k], y_var[traj], "L1", j))
            }
            outLoads_L[[traj]] <- list(outLoads1)
          }
        }
        else if (curveFun %in% c("Jenss-Bayley", "JB")){
          if (intrinsic){
            outLoads2 <- outLoads3 <- list()
            outLoads2[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L2", records[[traj]][1]))
            outLoads3[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L3", records[[traj]][1]))
            for (j in records[[traj]][-1]){
              outLoads2[[j]] <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "_mug * exp(", pfx[k], y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0(pfx[k], y_var[traj], "L2", j))
              outLoads3[[j]] <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "_mueta2 * (1 + ", pfx[k], y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")",
                                                           " * exp(", pfx[k], y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0(pfx[k], y_var[traj], "L3", j))
            }
            outLoads_L[[traj]] <- list(outLoads2, outLoads3)
          }
          else if (!intrinsic){
            outLoads2 <- list()
            outLoads2[[1]] <- mxAlgebra(0, name = paste0(y_var[traj], "L2", records[[traj]][1]))
            for (j in records[[traj]][-1]){
              outLoads2[[j]] <- mxAlgebraFromString(paste0(pfx[k], y_var[traj], "_mug * exp(", pfx[k], y_var[traj], "_mug * ", y_var[traj], "mid_t", j, ")"),
                                                    name = paste0(pfx[k], y_var[traj], "L2", j))
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

#' @keywords internal
#' @noRd
getMULTI.loadings <- function(y_model, t_var, y_var, curveFun, intrinsic, records){
  .loadings_MULTI_classes(pfx = "", y_model = y_model, t_var = t_var, y_var = y_var, curveFun = curveFun,
                          intrinsic = intrinsic, records = records)[[1]]
}
getMIX_MULTI.loadings <- function(nClass, y_model, t_var, y_var, curveFun, intrinsic = NULL, records){
  .loadings_MULTI_classes(pfx = paste0("c", seq_len(nClass)), y_model = y_model, t_var = t_var, y_var = y_var,
                          curveFun = curveFun, intrinsic = intrinsic, records = records)
}

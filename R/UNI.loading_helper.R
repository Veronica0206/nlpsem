#' @title Growth-factor loadings of a univariate outcome, one set per prefix
#' @description \code{pfx = ""} gives the single-group objects; \code{pfx = c("c1", "c2", ...)} the class-specific
#' ones used by the multiple-group and mixture models. The names of matrices and algebras carry the prefix.
#' @keywords internal
#' @noRd
#'
#' @importFrom OpenMx mxMatrix mxAlgebraFromString mxAlgebra
#'
.loadings_UNI_classes <- function(pfx, y_model, t_var, records, y_var, curveFun, intrinsic){
  # Define `outPoint` to store individual time points from raw data
  output_CL <- list()
  for (k in seq_along(pfx)){
    outPoint <- list()
    for (j in records){
      outPoint[[j]] <- mxMatrix("Full", 1, 1, free = FALSE, labels = paste0("data.", t_var, j),
                                name = paste0("t", j))
    }
    if (y_model == "LGCM"){
      if (curveFun %in% c("linear", "LIN")){
        outLoads1 <- list()
        for(j in records){
          outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j), name = paste0(pfx[k], "L1", j))
        }
        outLoads <- list(outLoads1)
      }
      else if (curveFun %in% c("quadratic", "QUAD")){
        outLoads1 <- outLoads2 <- list()
        for(j in records){
          outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j), name = paste0(pfx[k], "L1", j))
          outLoads2[[j]] <- mxAlgebraFromString(paste0("t", j, "^2"), name = paste0(pfx[k], "L2", j))
        }
        outLoads <- list(outLoads1, outLoads2)
      }
      else if (curveFun %in% c("negative exponential", "EXP")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("1 - exp(-", pfx[k], "Y_mug * t", j, ")"),
                                                  name = paste0(pfx[k], "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0(pfx[k], "Y_mueta1 * exp(-", pfx[k], "Y_mug * t", j, ") * t", j),
                                                  name = paste0(pfx[k], "L2", j))
          }
          outLoads <- list(outLoads1, outLoads2)
        }
        else if (!intrinsic){
          outLoads1 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("1 - exp(-", pfx[k], "Y_mug * t", j, ")"),
                                                  name = paste0(pfx[k], "L1", j))
          }
          outLoads <- list(outLoads1)
        }
      }
      else if (curveFun %in% c("Jenss-Bayley", "JB")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- outLoads3 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j), name = paste0(pfx[k], "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("exp(", pfx[k], "Y_mug * t", j, ") - 1"),
                                                  name = paste0(pfx[k], "L2", j))
            outLoads3[[j]] <- mxAlgebraFromString(paste0(pfx[k], "Y_mueta2 * exp(", pfx[k], "Y_mug * t", j, ") * t", j),
                                                  name = paste0(pfx[k], "L3", j))
          }
          outLoads <- list(outLoads1, outLoads2, outLoads3)
        }
        else if (!intrinsic){
          outLoads1 <- outLoads2 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j), name = paste0(pfx[k], "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("exp(", pfx[k], "Y_mug * t", j, ") - 1"),
                                                  name = paste0(pfx[k], "L2", j))
          }
          outLoads <- list(outLoads1, outLoads2)
        }
      }
      else if (curveFun %in% c("bilinear spline", "BLS")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- outLoads3 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " - ", pfx[k], "Y_mug"),
                                                  name = paste0(pfx[k], "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("sqrt((t", j, " - ", pfx[k], "Y_mug)^2 + 1e-6)"),
                                                  name = paste0(pfx[k], "L2", j))
            outLoads3[[j]] <- mxAlgebraFromString(paste0("-", pfx[k], "Y_mueta2s * (t", j, " - ", pfx[k], "Y_mug)/sqrt((t", j,
                                                         " - ", pfx[k], "Y_mug)^2 + 1e-6) - ", pfx[k], "Y_mueta1s"),
                                                  name = paste0(pfx[k], "L3", j))
          }
          outLoads <- list(outLoads1, outLoads2, outLoads3)
        }
        else if (!intrinsic){
          outLoads1 <- outLoads2 <- list()
          for(j in records){
            outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " - ", pfx[k], "Y_mug"),
                                                  name = paste0(pfx[k], "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " - ", pfx[k], "Y_mug)"),
                                                  name = paste0(pfx[k], "L2", j))
          }
          outLoads <- list(outLoads1, outLoads2)
        }
      }
      output_CL[[k]] <- list(outPoint, outLoads)
    }
    else if (y_model == "LCSM"){
      # Define `outLag` to store individual time lags from raw data
      outLag <- list()
      outLag[[1]] <- mxAlgebra(0, name = paste0("lag", records[1]))
      for (pos in seq_along(records)[-1]){
        j <- records[pos]
        outLag[[j]] <- mxAlgebraFromString(paste0("t", j , " -  t", records[pos - 1]), name = paste0("lag", j))
      }
      # Define `midTime` to store individual midway time point through each interval from raw data, which help
      # define factor loadings of latent change score models.
      midTime <- list()
      midTime[[1]] <- mxAlgebra(0, name = paste0("mid_t", records[1]))
      for (pos in seq_along(records)[-1]){
        j <- records[pos]
        midTime[[j]] <- mxAlgebraFromString(paste0("(t", j , " +  t", records[pos - 1], ")/2"),
                                            name = paste0("mid_t", j))
      }
      if (curveFun %in% c("nonparametric", "NonP")){
        outLoads <- list()
      }
      else if (curveFun %in% c("quadratic", "QUAD")){
        outLoads2 <- list()
        outLoads2[[1]] <- mxAlgebra(0, name = paste0(pfx[k], "L2", records[1]))
        for(j in records[-1]){
          outLoads2[[j]] <- mxAlgebraFromString(paste0("2 * mid_t", j), name = paste0(pfx[k], "L2", j))
        }
        outLoads <- outLoads2
      }
      else if (curveFun %in% c("negative exponential", "EXP")){
        if (intrinsic){
          outLoads1 <- outLoads2 <- list()
          outLoads1[[1]] <- mxAlgebra(0, name = paste0(pfx[k], "L1", records[1]))
          outLoads2[[1]] <- mxAlgebra(0, name = paste0(pfx[k], "L2", records[1]))
          for (j in records[-1]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(pfx[k], "Y_mug * exp(-", pfx[k], "Y_mug * mid_t", j, ")"),
                                                  name = paste0(pfx[k], "L1", j))
            outLoads2[[j]] <- mxAlgebraFromString(paste0(pfx[k], "Y_mueta1 * exp(-", pfx[k], "Y_mug * mid_t", j,
                                                         ") * (1 - ", pfx[k], "Y_mug * mid_t", j, ")"),
                                                  name = paste0(pfx[k], "L2", j))
          }
          outLoads <- list(outLoads1, outLoads2)
        }
        else if (!intrinsic){
          outLoads1 <- list()
          outLoads1[[1]] <- mxAlgebra(0, name = paste0(pfx[k], "L1", records[1]))
          for (j in records[-1]){
            outLoads1[[j]] <- mxAlgebraFromString(paste0(pfx[k], "Y_mug * exp(-", pfx[k], "Y_mug * mid_t", j, ")"),
                                                  name = paste0(pfx[k], "L1", j))
          }
          outLoads <- list(outLoads1)
        }
      }
      else if (curveFun %in% c("Jenss-Bayley", "JB")){
        if (intrinsic){
          outLoads2 <- outLoads3 <- list()
          outLoads2[[1]] <- mxAlgebra(0, name = paste0(pfx[k], "L2", records[1]))
          outLoads3[[1]] <- mxAlgebra(0, name = paste0(pfx[k], "L3", records[1]))
          for (j in records[-1]){
            outLoads2[[j]] <- mxAlgebraFromString(paste0(pfx[k], "Y_mug * exp(", pfx[k], "Y_mug * mid_t", j, ")"),
                                                  name = paste0(pfx[k], "L2", j))
            outLoads3[[j]] <- mxAlgebraFromString(paste0(pfx[k], "Y_mueta2 * (1 + ", pfx[k], "Y_mug *  mid_t", j, ")",
                                                         " * exp(", pfx[k], "Y_mug * mid_t", j, ")"),
                                                  name = paste0(pfx[k], "L3", j))
          }
          outLoads <- list(outLoads2, outLoads3)
        }
        else if (!intrinsic){
          outLoads2 <- list()
          outLoads2[[1]] <- mxAlgebra(0, name = paste0(pfx[k], "L2", records[1]))
          for (j in records[-1]){
            outLoads2[[j]] <- mxAlgebraFromString(paste0(pfx[k], "Y_mug * exp(", pfx[k], "Y_mug * mid_t", j, ")"),
                                                  name = paste0(pfx[k], "L2", j))
          }
          outLoads <- list(outLoads2)
        }
      }
      output_CL[[k]] <- list(outPoint, outLag, midTime, outLoads)
    }
  }
  return(output_CL)
}

#' @keywords internal
#' @noRd
getUNI.loadings <- function(y_model, t_var, y_var, curveFun, intrinsic, records){
  .loadings_UNI_classes(pfx = "", y_model = y_model, t_var = t_var, records = records, y_var = y_var,
                        curveFun = curveFun, intrinsic = intrinsic)[[1]]
}
getMIX_UNI.loadings <- function(nClass, y_model, t_var, records, y_var, curveFun, intrinsic){
  .loadings_UNI_classes(pfx = paste0("c", seq_len(nClass)), y_model = y_model, t_var = t_var, records = records,
                        y_var = y_var, curveFun = curveFun, intrinsic = intrinsic)
}

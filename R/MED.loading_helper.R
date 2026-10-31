#' @title Growth-factor loadings of a longitudinal mediation model, one set per prefix
#' @description \code{pfx = ""} gives the single-group objects; \code{pfx = c("c1", "c2", ...)} the class-specific
#' ones used by the multiple-group and mixture models. The names of matrices and algebras carry the prefix.
#' @keywords internal
#' @noRd
#'

#'
.loadings_MED_classes <- function(pfx, t_var, y_var, m_var, x_var, x_type, curveFun, y_records, m_records, x_records = NULL){
  output_CL <- list()
  for (k in seq_along(pfx)){
    outPointM <- outPointY <- list()
    for (j in m_records){
      outPointM[[j]] <- mxMatrix("Full", 1, 1, free = FALSE, labels = paste0("data.", t_var[2], j),
                                 name = paste0("t", m_var, j))
    }
    for (j in y_records){
      outPointY[[j]] <- mxMatrix("Full", 1, 1, free = FALSE, labels = paste0("data.", t_var[1], j),
                                 name = paste0("t", y_var, j))
    }
    if (x_type == "baseline"){
      if (curveFun %in% c("linear", "LIN")){
        outLoadsY1 <- outLoadsM1 <- list()
        for (j in m_records){
          outLoadsM1[[j]] <- mxAlgebraFromString(paste0("t", m_var, j), name = paste0(pfx[k], "L1", j, "M"))
        }
        for (j in y_records){
          outLoadsY1[[j]] <- mxAlgebraFromString(paste0("t", y_var, j), name = paste0(pfx[k], "L1", j, "Y"))
        }
        outLoads <- list(outLoadsY1, outLoadsM1)
      }
      else if (curveFun %in% c("bilinear spline", "BLS")){
        outLoadsY1 <- outLoadsY2 <- outLoadsM1 <- outLoadsM2 <- list()
        for (j in m_records){
          outLoadsM1[[j]] <- mxAlgebraFromString(paste0("min(0, t", m_var, j, " - ", pfx[k], "M_mug)"),
                                                 name = paste0(pfx[k], "L1", j, "M"))
          outLoadsM2[[j]] <- mxAlgebraFromString(paste0("max(0, t", m_var, j, " - ", pfx[k], "M_mug)"),
                                                 name = paste0(pfx[k], "L2", j, "M"))
        }
        for (j in y_records){
          outLoadsY1[[j]] <- mxAlgebraFromString(paste0("min(0, t", y_var, j, " - ", pfx[k], "Y_mug)"),
                                                 name = paste0(pfx[k], "L1", j, "Y"))
          outLoadsY2[[j]] <- mxAlgebraFromString(paste0("max(0, t", y_var, j, " - ", pfx[k], "Y_mug)"),
                                                 name = paste0(pfx[k], "L2", j, "Y"))
        }
        outLoads <- list(outLoadsY1, outLoadsY2, outLoadsM1, outLoadsM2)
      }
      outPoint <- list(outPointY, outPointM)
    }
    else if (x_type == "longitudinal"){
      outPointX <- list()
      for (j in x_records){
        outPointX[[j]] <- mxMatrix("Full", 1, 1, free = FALSE, labels = paste0("data.", t_var[3], j),
                                   name = paste0("t", x_var, j))
      }
      if (curveFun %in% c("linear", "LIN")){
        outLoadsY1 <- outLoadsM1 <- outLoadsX1 <- list()
        for (j in x_records){
          outLoadsX1[[j]] <- mxAlgebraFromString(paste0("t", x_var, j), name = paste0(pfx[k], "L1", j, "X"))
        }
        for (j in m_records){
          outLoadsM1[[j]] <- mxAlgebraFromString(paste0("t", m_var, j), name = paste0(pfx[k], "L1", j, "M"))
        }
        for (j in y_records){
          outLoadsY1[[j]] <- mxAlgebraFromString(paste0("t", y_var, j), name = paste0(pfx[k], "L1", j, "Y"))
        }
        outLoads <- list(outLoadsY1, outLoadsM1, outLoadsX1)
      }
      else if (curveFun %in% c("bilinear spline", "BLS")){
        outLoadsY1 <- outLoadsY2 <- outLoadsM1 <- outLoadsM2 <- outLoadsX1 <- outLoadsX2 <- list()
        for (j in x_records){
          outLoadsX1[[j]] <- mxAlgebraFromString(paste0("min(0, t", x_var, j, " - ", pfx[k], "X_mug)"),
                                                 name = paste0(pfx[k], "L1", j, "X"))
          outLoadsX2[[j]] <- mxAlgebraFromString(paste0("max(0, t", x_var, j, " - ", pfx[k], "X_mug)"),
                                                 name = paste0(pfx[k], "L2", j, "X"))
        }
        for (j in m_records){
          outLoadsM1[[j]] <- mxAlgebraFromString(paste0("min(0, t", m_var, j, " - ", pfx[k], "M_mug)"),
                                                 name = paste0(pfx[k], "L1", j, "M"))
          outLoadsM2[[j]] <- mxAlgebraFromString(paste0("max(0, t", m_var, j, " - ", pfx[k], "M_mug)"),
                                                 name = paste0(pfx[k], "L2", j, "M"))
        }
        for (j in y_records){
          outLoadsY1[[j]] <- mxAlgebraFromString(paste0("min(0, t", y_var, j, " - ", pfx[k], "Y_mug)"),
                                                 name = paste0(pfx[k], "L1", j, "Y"))
          outLoadsY2[[j]] <- mxAlgebraFromString(paste0("max(0, t", y_var, j, " - ", pfx[k], "Y_mug)"),
                                                 name = paste0(pfx[k], "L2", j, "Y"))
        }
        outLoads <- list(outLoadsY1, outLoadsY2, outLoadsM1, outLoadsM2, outLoadsX1, outLoadsX2)
      }
      outPoint <- list(outPointY, outPointM, outPointX)
    }
    output_CL[[k]] <- list(outPoint, outLoads)
  }
  return(output_CL)
}


#' @keywords internal
#' @noRd
getMED.loadings <- function(t_var, y_var, m_var, x_type, x_var, curveFun, y_records, m_records, x_records = NULL){
  .loadings_MED_classes(pfx = "", t_var = t_var, y_var = y_var, m_var = m_var, x_var = x_var, x_type = x_type,
                        curveFun = curveFun, y_records = y_records, m_records = m_records, x_records = x_records)[[1]]
}
getMIX_MED.loadings <- function(nClass, t_var, y_var, m_var, x_var, x_type, curveFun, y_records, m_records,
                                x_records = NULL){
  .loadings_MED_classes(pfx = paste0("c", seq_len(nClass)), t_var = t_var, y_var = y_var, m_var = m_var, x_var = x_var,
                        x_type = x_type, curveFun = curveFun, y_records = y_records, m_records = m_records,
                        x_records = x_records)
}

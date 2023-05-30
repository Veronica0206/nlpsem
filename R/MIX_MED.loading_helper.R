#' @title Get Factor Loadings for a Mixture Model or Multiple Group Model with Longitudinal Mediation Model with Specified Functional Curves as Submodels
#'
#' @description This function specifies the factor loadings for a mixture model with longitudinal mediation model with given functional
#' form as submodels
#'
#' @param nClass An integer specifying the number of classes for the mixture model or multiple group model. It takes the value passed
#' from \code{getMIX()} or \code{getMGroup()}.
#' @param t_var A vector of strings, with each element representing the prefix for column names related to the time variable for the
#' corresponding longitudinal variable at each study wave. It takes the value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave. It takes the value
#' passed from \code{getMIX()} or \code{getMGroup()}.
#' @param m_var A string specifying the prefix of the column names corresponding to the mediator variable at each study wave.
#' It takes the value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. It takes the value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. It takes the value passed from \code{getMIX()}
#' or \code{getMGroup()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options include: "linear" (or "LIN"),
#' and "bilinear spline" (or "BLS"). It takes the value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param y_records A numeric vector specifying indices of the study waves for the outcome variable. It is the first vector
#' in \code{records} specified in \code{getMIX()} or \code{getMGroup()}.
#' @param m_records A numeric vector specifying indices of the study waves for the mediator variable. It is the second vector
#' in \code{records} specified in \code{getMIX()} or \code{getMGroup()}.
#' @param x_records A numeric vector specifying indices of the study waves for the predictor variable. Default value is \code{NULL},
#' which is for \code{x_type = "baseline"}. If \code{x_type = "longitudinal"}, it is the third vector in \code{records} specified in
#' \code{getMIX()} or \code{getMGroup()}.
#'
#' @return A list containing the specification of definition variables (i.e., individual measurement occasions) and factor
#' loadings of a longitudinal mediation model.
#'
getMIX_MED.loadings <- function(nClass, t_var, y_var, m_var, x_var, x_type, curveFun, y_records, m_records,
                                x_records = NULL){
  output_CL <- list()
  for (k in 1:nClass){
    outPointM <- outPointY <- list()
    for (j in m_records){
      outPointM[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[1], j),
                                 name = paste0("t", m_var, j))
    }
    for (j in y_records){
      outPointY[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[2], j),
                                 name = paste0("t", y_var, j))
    }
    if (x_type == "baseline"){
      if (curveFun %in% c("linear", "LIN")){
        outLoadsY1 <- outLoadsM1 <- list()
        for (j in m_records){
          outLoadsM1[[j]] <- mxAlgebraFromString(paste0("t", m_var, j), name = paste0("c", k, "L1", j, "M"))
        }
        for (j in y_records){
          outLoadsY1[[j]] <- mxAlgebraFromString(paste0("t", y_var, j), name = paste0("c", k, "L1", j, "Y"))
        }
        outLoads <- list(outLoadsY1, outLoadsM1)
      }
      else if (curveFun %in% c("bilinear spline", "BLS")){
        outLoadsY1 <- outLoadsY2 <- outLoadsM1 <- outLoadsM2 <- list()
        for (j in m_records){
          outLoadsM1[[j]] <- mxAlgebraFromString(paste0("min(0, t", m_var, j, " - c", k, "M_mug)"),
                                                 name = paste0("c", k, "L1", j, "M"))
          outLoadsM2[[j]] <- mxAlgebraFromString(paste0("max(0, t", m_var, j, " - c", k, "M_mug)"),
                                                 name = paste0("c", k, "L2", j, "M"))
        }
        for (j in y_records){
          outLoadsY1[[j]] <- mxAlgebraFromString(paste0("min(0, t", y_var, j, " - c", k, "Y_mug)"),
                                                 name = paste0("c", k, "L1", j, "Y"))
          outLoadsY2[[j]] <- mxAlgebraFromString(paste0("max(0, t", y_var, j, " - c", k, "Y_mug)"),
                                                 name = paste0("c", k, "L2", j, "Y"))
        }
        outLoads <- list(outLoadsY1, outLoadsY2, outLoadsM1, outLoadsM2)
      }
      outPoint <- list(outPointY, outPointM)
    }
    else if (x_type == "longitudinal"){
      outPointX <- list()
      for (j in x_records){
        outPointX[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[1], j),
                                   name = paste0("t", x_var, j))
      }
      if (curveFun %in% c("linear", "LIN")){
        outLoadsY1 <- outLoadsM1 <- outLoadsX1 <- list()
        for (j in x_records){
          outLoadsX1[[j]] <- mxAlgebraFromString(paste0("t", x_var, j), name = paste0("c", k, "L1", j, "X"))
        }
        for (j in m_records){
          outLoadsM1[[j]] <- mxAlgebraFromString(paste0("t", m_var, j), name = paste0("c", k, "L1", j, "M"))
        }
        for (j in y_records){
          outLoadsY1[[j]] <- mxAlgebraFromString(paste0("t", y_var, j), name = paste0("c", k, "L1", j, "Y"))
        }
        outLoads <- list(outLoadsY1, outLoadsM1, outLoadsX1)
      }
      else if (curveFun %in% c("bilinear spline", "BLS")){
        outLoadsY1 <- outLoadsY2 <- outLoadsM1 <- outLoadsM2 <- outLoadsX1 <- outLoadsX2 <- list()
        for (j in x_records){
          outLoadsX1[[j]] <- mxAlgebraFromString(paste0("min(0, t", x_var, j, " - c", k, "X_mug)"),
                                                 name = paste0("c", k, "L1", j, "X"))
          outLoadsX2[[j]] <- mxAlgebraFromString(paste0("max(0, t", x_var, j, " - c", k, "X_mug)"),
                                                 name = paste0("c", k, "L2", j, "X"))
        }
        for (j in m_records){
          outLoadsM1[[j]] <- mxAlgebraFromString(paste0("min(0, t", m_var, j, " - c", k, "M_mug)"),
                                                 name = paste0("c", k, "L1", j, "M"))
          outLoadsM2[[j]] <- mxAlgebraFromString(paste0("max(0, t", m_var, j, " - c", k, "M_mug)"),
                                                 name = paste0("c", k, "L2", j, "M"))
        }
        for (j in y_records){
          outLoadsY1[[j]] <- mxAlgebraFromString(paste0("min(0, t", y_var, j, " - c", k, "Y_mug)"),
                                                 name = paste0("c", k, "L1", j, "Y"))
          outLoadsY2[[j]] <- mxAlgebraFromString(paste0("max(0, t", y_var, j, " - c", k, "Y_mug)"),
                                                 name = paste0("c", k, "L2", j, "Y"))
        }
        outLoads <- list(outLoadsY1, outLoadsY2, outLoadsM1, outLoadsM2, outLoadsX1, outLoadsX2)
      }
      outPoint <- list(outPointY, outPointM, outPointX)
    }
    output_CL[[k]] <- list(outPoint, outLoads)
  }
  return(output_CL)
}


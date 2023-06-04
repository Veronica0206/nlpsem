#' @title Get Factor Loadings for a Longitudinal Mediation Model with Specified Functional Curves
#'
#' @description This function specifies the factor loadings for a longitudinal mediation model with given functional form.
#'
#' @param t_var A vector of strings, with each element representing the prefix for column names related to the time
#' variable for the corresponding longitudinal variable at each study wave. It takes the value passed from \code{getMediation()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getMediation()}.
#' @param m_var A string specifying the prefix of the column names corresponding to the mediator variable at each study wave.
#' It takes the value passed from \code{getMediation()}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. It takes the value passed from \code{getMediation()}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. It takes the value passed from
#' \code{getMediation()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options include: \code{"linear"}
#' (or \code{"LIN"}), and \code{"bilinear spline"} (or \code{"BLS"}). It takes the value passed from \code{getMediation()}.
#' @param y_records A numeric vector specifying indices of the study waves for the outcome variable. It is the first vector
#' in \code{records} specified in \code{getMediation()}.
#' @param m_records A numeric vector specifying indices of the study waves for the mediator variable. It is the second vector
#' in \code{records} specified in \code{getMediation()}.
#' @param x_records A numeric vector specifying indices of the study waves for the predictor variable. Default value is \code{NULL},
#' which is for \code{x_type = "baseline"}. If \code{x_type = "longitudinal"}, it is the third vector in \code{records} specified in
#' \code{getMediation()}.
#'
#' @return A list containing the specification of definition variables (i.e., individual measurement occasions) and factor
#' loadings of a longitudinal mediation model.
#'
#' @keywords internal
#'
getMED.loadings <- function(t_var, y_var, m_var, x_type, x_var, curveFun, y_records, m_records,
                            x_records = NULL){
  outDefM <- outDefY <- list()
  for (j in m_records){
    outDefM[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[1], j),
                             name = paste0("t", m_var, j))
  }
  for (j in y_records){
    outDefY[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[2], j),
                             name = paste0("t", y_var, j))
  }
  if (x_type == "baseline"){
    if (curveFun %in% c("linear", "LIN")){
      outLoadsY1 <- outLoadsM1 <- list()
      for (j in m_records){
        outLoadsM1[[j]] <- mxAlgebraFromString(paste0("t", m_var, j), name = paste0("L1", j, "M"))
      }
      for (j in y_records){
        outLoadsY1[[j]] <- mxAlgebraFromString(paste0("t", y_var, j), name = paste0("L1", j, "Y"))
      }
      outLoads <- list(outLoadsY1, outLoadsM1)
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      outLoadsY1 <- outLoadsY2 <- outLoadsM1 <- outLoadsM2 <- list()
      for (j in m_records){
        outLoadsM1[[j]] <- mxAlgebraFromString(paste0("min(0, t", m_var, j, " - M_mug)"),
                                               name = paste0("L1", j, "M"))
        outLoadsM2[[j]] <- mxAlgebraFromString(paste0("max(0, t", m_var, j, " - M_mug)"),
                                               name = paste0("L2", j, "M"))
      }
      for (j in y_records){
        outLoadsY1[[j]] <- mxAlgebraFromString(paste0("min(0, t", y_var, j, " - Y_mug)"),
                                               name = paste0("L1", j, "Y"))
        outLoadsY2[[j]] <- mxAlgebraFromString(paste0("max(0, t", y_var, j, " - Y_mug)"),
                                               name = paste0("L2", j, "Y"))
      }
      outLoads <- list(outLoadsY1, outLoadsY2, outLoadsM1, outLoadsM2)
    }
    outDef <- list(outDefY, outDefM)
  }
  else if (x_type == "longitudinal"){
    outDefX <- list()
    for (j in x_records){
      outDefX[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[1], j),
                               name = paste0("t", x_var, j))
    }
    if (curveFun %in% c("linear", "LIN")){
      outLoadsY1 <- outLoadsM1 <- outLoadsX1 <- list()
      for (j in x_records){
        outLoadsX1[[j]] <- mxAlgebraFromString(paste0("t", x_var, j), name = paste0("L1", j, "X"))
      }
      for (j in m_records){
        outLoadsM1[[j]] <- mxAlgebraFromString(paste0("t", m_var, j), name = paste0("L1", j, "M"))
      }
      for (j in y_records){
        outLoadsY1[[j]] <- mxAlgebraFromString(paste0("t", y_var, j), name = paste0("L1", j, "Y"))
      }
      outLoads <- list(outLoadsY1, outLoadsM1, outLoadsX1)
    }
    else if (curveFun %in% c("bilinear spline", "BLS")){
      outLoadsY1 <- outLoadsY2 <- outLoadsM1 <- outLoadsM2 <- outLoadsX1 <- outLoadsX2 <- list()
      for (j in x_records){
        outLoadsX1[[j]] <- mxAlgebraFromString(paste0("min(0, t", x_var, j, " - X_mug)"),
                                               name = paste0("L1", j, "X"))
        outLoadsX2[[j]] <- mxAlgebraFromString(paste0("max(0, t", x_var, j, " - X_mug)"),
                                               name = paste0("L2", j, "X"))
      }
      for (j in m_records){
        outLoadsM1[[j]] <- mxAlgebraFromString(paste0("min(0, t", m_var, j, " - M_mug)"),
                                               name = paste0("L1", j, "M"))
        outLoadsM2[[j]] <- mxAlgebraFromString(paste0("max(0, t", m_var, j, " - M_mug)"),
                                               name = paste0("L2", j, "M"))
      }
      for (j in y_records){
        outLoadsY1[[j]] <- mxAlgebraFromString(paste0("min(0, t", y_var, j, " - Y_mug)"),
                                               name = paste0("L1", j, "Y"))
        outLoadsY2[[j]] <- mxAlgebraFromString(paste0("max(0, t", y_var, j, " - Y_mug)"),
                                               name = paste0("L2", j, "Y"))
      }
      outLoads <- list(outLoadsY1, outLoadsY2, outLoadsM1, outLoadsM2, outLoadsX1, outLoadsX2)
    }
    outDef <- list(outDefY, outDefM, outDefX)
  }
  return(list(outDef, outLoads))
}


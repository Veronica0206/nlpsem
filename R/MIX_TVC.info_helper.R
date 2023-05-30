#' @title Get the Time-Varying Covariate (TVC) Information for a Mixture Model or Multiple Group Model with a Time-varying Covariate
#'
#' @description This function constructs the OpenMx model paths and parameters for a TVC and its relationship with the parameters
#' related to growth factors of a longitudinal outcome.
#'
#' @param nClass An integer specifying the number of classes for the mixture model or multiple group model. It takes the value
#' passed from \code{getMIX()} or \code{getMGroup()}.
#' @param y_var A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
#' It takes the value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param records A numeric vector specifying the indices of the observed study waves. It takes the value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param growth_TIC A string or character vector specifying the column name(s) of time-invariant covariate(s) that account for the
#' variability of growth factors, if any. Default is \code{NULL}, indicating no growth TICs present in the model. It takes the value
#' passed from \code{getMIX()} or \code{getMGroup()}.
#' @param TVC A string specifying the prefix of the column names corresponding to the time-varying covariate at each study wave. It
#' takes the value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). It takes the
#' value passed from \code{getMIX()} or \code{getMGroup()}.
#' @param starts A list of initial values for the parameters, either takes the value passed from \code{getMIX()} or \code{getMGroup()}, or
#' derived by the helper function \code{getMIX.initial()} or \code{getMGroup.initial()}.
#'
#' @return A list with each sub-list containing a list of OpenMx paths and parameters for the TVC and an OpenMx path for the
#' state effect of the TVC on the corresponding longitudinal outcome.
#'
#' @importFrom OpenMx mxPath mxAlgebra mxAlgebraFromString
#'
getMIX_TVC.info <- function(nClass, y_var, records, growth_TIC, TVC, decompose, starts){
  output_CL <- list()
  # Case when no decomposition of TVC is required
  if (decompose == 0){
    # Loop over all classes
    for (k in 1:nClass){
      # If growth_TIC is provided
      if (!is.null(growth_TIC)){
        nTICs <- length(growth_TIC)
        # Means of TICs
        X_MEAN <- mxPath(from = "one", to = growth_TIC, arrows = 1, free = TRUE, values = starts[[k]][[3]][[1]],
                         labels = c(paste0("c", k, "mux", 1:length(growth_TIC))))
        # Var-cov matrix of TICs
        X_VAR <- mxPath(from = growth_TIC, to = growth_TIC, arrows = 2, connect = "unique.pairs", free = TRUE,
                        values = starts[[k]][[3]][[2]], labels = paste0("c", k, "phi", 1:(nTICs * (nTICs + 1)/2)))
        # TVC effect on the corresponding longitudinal outcome
        KAPPA <- mxPath(from = paste0(TVC, records), to = paste0(y_var, records), arrows = 1, free = TRUE,
                        values = starts[[k]][[5]], labels = paste0("c", k, "kappa"))
        output_CL[[k]] <- list(KAPPA, X_MEAN, X_VAR)
      }
      # If growth_TIC is not provided
      else if (is.null(growth_TIC)){
        # TVC effect on the corresponding longitudinal outcome
        KAPPA <- mxPath(from = paste0(TVC, records), to = paste0(y_var, records), arrows = 1, free = TRUE,
                        values = starts[[k]][[3]], labels = paste0("c", k, "kappa"))
        output_CL[[k]] <- list(KAPPA)
      }
    }
  }
  # Case when decomposition of TVC is required
  else if (decompose != 0){
    # Loop over all classes
    for (k in 1:nClass){
      # Define additional parameters for the TVC (i.e., interval-specific slopes)
      abs_rate <- outLag <- list()
      outLag[[1]] <- mxAlgebra(0, name = paste0("lag", 1))
      abs_rate[[1]] <- mxAlgebra(0, name = paste0("c", k, "abs_rate", 1))
      for (j in records[-1]){
        abs_rate[[j]] <- mxAlgebraFromString(paste0("c", k, "X_mueta1 * c", k, "X_rel_rate", j - 1),
                                             name = paste0("c", k, "abs_rate", j))
        outLag[[j]] <- mxAlgebraFromString(paste0("t", j , " -  t", j - 1),
                                           name = paste0("lag", j))
      }
      # Define growth factor loadings of the TVC, paths of the TVC and TICs (if any)
      if (!is.null(growth_TIC)){
        nTICs <- length(growth_TIC)
        start_X <- matrix(0, nrow = nTICs + 2, ncol = nTICs + 2)
        start_X[1:(nTICs + 1), 1:(nTICs + 1)] <- starts[[k]][[3]][[2]]
        start_X[nTICs + 2, nTICs + 1] <- start_X[nTICs + 1, nTICs + 2] <- starts[[k]][[2]][[2]][2]
        start_X[nTICs + 2, nTICs + 2] <- starts[[k]][[2]][[2]][3]
        name_X <- matrix(NA, nrow = nTICs + 2, ncol = nTICs + 2)
        name_X[(nrow(name_X) - 1):nrow(name_X), (ncol(name_X) - 1):ncol(name_X)] <-
          matrix(paste0("c", k, c("X_psi00", "X_psi01", "X_psi01", "X_psi11")), byrow = T, nrow = 2)
        for (i in 1:nTICs){
          for (j in 1:nTICs){
            name_X[i, j] <- paste0("c", k, "phi", j, i)
          }
        }
        name_X[1:nTICs, nTICs + 1] <- name_X[nTICs + 1, 1:nTICs] <- paste0("c", k, "covBL", 1:nTICs)
        X_MEAN <- mxPath(from = "one", to = c(growth_TIC, "eta0x", "eta1x"), arrows = 1, free = TRUE,
                         values = c(starts[[k]][[3]][[1]][1:length(growth_TIC)], starts[[k]][[2]][[1]]),
                         labels = c(paste0("c", k, "mux", 1:length(growth_TIC)),
                                    paste0("c", k, c("X_mueta0", "X_mueta1"))))
        X_VAR <- mxPath(from = c(growth_TIC, "eta0x", "eta1x"), to = c(growth_TIC, "eta0x", "eta1x"), arrows = 2,
                        connect = "unique.pairs", free = !(start_X[row(start_X) >= col(start_X)] == 0),
                        values = start_X[row(start_X) >= col(start_X)],
                        labels = name_X[row(name_X) >= col(name_X)])
        X_GF_LOADINGS <- list(mxPath(from = "eta0x", to = "lx1", arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1x", to = paste0("dx", records[-1]), arrows = 1,
                                     free = c(F, rep(T, length(records) - 2)), values = c(1, starts[[k]][[2]][[4]][-1]),
                                     labels = paste0("c", k, "X_rel_rate", 1:(length(records) - 1))))
      }
      else if (is.null(growth_TIC)){
        X_MEAN <- mxPath(from = "one", to = c("eta0x", "eta1x"), arrows = 1, free = TRUE,
                         values = starts[[k]][[2]][[1]], labels = paste0("c", k, c("X_mueta0", "X_mueta1")))
        X_VAR <- mxPath(from = c("eta0x", "eta1x"), to = c("eta0x", "eta1x"), arrows = 2, connect = "unique.pairs",
                        free = TRUE, values = starts[[k]][[2]][[2]],
                        labels = paste0("c", k, c("X_psi00", "X_psi01", "X_psi11")))
        X_GF_LOADINGS <- list(mxPath(from = "eta0x", to = "lx1", arrows = 1, free = FALSE, values = 1),
                              mxPath(from = "eta1x", to = paste0("dx", records[-1]), arrows = 1,
                                     free = c(F, rep(T, length(records) - 2)), values = c(1, starts[[k]][[2]][[4]][-1]),
                                     labels = paste0("c", k, "X_rel_rate", 1:(length(records) - 1))))
      }
      X_RES <- mxPath(from = paste0(TVC, records), to = paste0(TVC, records), arrows = 2,
                      free = TRUE, values = starts[[k]][[2]][[3]], labels = paste0("c", k, "X_residuals"))
      XY_COV <- mxPath(from = paste0(TVC, records), to = paste0(y_var, records), arrows = 2,
                       free = TRUE, values = starts[[k]][[6]], labels = paste0("c", k, "Cov_XYres"))
      X_PATH_L <- mxPath(from = paste0("lx", records), to = paste0(TVC, records), arrows = 1,
                         free = FALSE, values = 1)

      # Decomposition method 1: true score of baseline and interval-specific slopes
      if (decompose == 1){
        X_PATH_SLP <- mxPath(from = paste0("dx", records[-1]), to = paste0("lx", records[-1]), arrows = 1, free = FALSE,
                             values = 0, labels = paste0("lag", records[-1], "[1,1]"))
        X_PATH_AUTO <- mxPath(from = paste0("lx", records[-length(records)]), to = paste0("lx", records[-1]),
                              arrows = 1, free = FALSE, values = 1)
        KAPPA <- mxPath(from = paste0("dx", records[-1]), to = paste0(y_var, records[-1]), arrows = 1, free = TRUE,
                        values = starts[[k]][[5]], labels = paste0("c", k, "kappa"))
        X_PATH_SLP_L <- list(X_PATH_SLP)
        X_PATH_AUTO_L <- list(X_PATH_AUTO)
      }
      # Decomposition method 1: true score of baseline and interval-specific changes
      else if (decompose == 2){
        X_PATH_SLP1 <- mxPath(from = paste0("dx", records[-1]), to = paste0("deltax", records[-1]), arrows = 1, free = FALSE,
                              values = 0, labels = paste0("lag", records[-1], "[1,1]"))
        X_PATH_SLP2 <- mxPath(from = paste0("deltax", records[-1]), to = paste0("lx",records[-1]), arrows = 1,
                              free = FALSE, values = 1)
        X_PATH_AUTO <- mxPath(from = paste0("lx", records[-length(records)]), to = paste0("lx", records[-1]),
                              arrows = 1, free = FALSE, values = 1)
        KAPPA <- mxPath(from = paste0("deltax", records[-1]), to = paste0(y_var, records[-1]), arrows = 1, free = TRUE,
                        values = starts[[k]][[5]], labels = paste0("c", k, "kappa"))
        X_PATH_SLP_L <- list(X_PATH_SLP1, X_PATH_SLP2)
        X_PATH_AUTO_L <- list(X_PATH_AUTO)
      }
      # Decomposition method 1: true score of baseline and change from baseline
      else if (decompose == 3){
        #### Define path from latent instantaneous rate of change at each measurement to true scores
        X_PATH_SLP1 <- mxPath(from = paste0("dx", records[-1]), to = paste0("Deltax", records[-1]), arrows = 1, free = FALSE,
                              values = 0, labels = paste0("lag", records[-1], "[1,1]"))
        X_PATH_SLP2 <- mxPath(from = paste0("Deltax", records[-1]), to = paste0("lx", records[-1]), arrows = 1,
                              free = FALSE, values = 1)
        X_PATH_AUTO1 <- mxPath(from = paste0("Deltax", records[-c(1, length(records))]),
                               to = paste0("Deltax", records[-c(1, 2)]), arrows = 1, free = FALSE, values = 1)
        X_PATH_AUTO2 <- mxPath(from = "lx1", to = paste0("lx", records[-1]), arrows = 1, free = FALSE, values = 1)
        KAPPA <- mxPath(from = paste0("Deltax", records[-1]), to = paste0(y_var, records[-1]), arrows = 1, free = TRUE,
                        values = starts[[k]][[5]], labels = paste0("c", k, "kappa"))
        X_PATH_SLP_L <- list(X_PATH_SLP1, X_PATH_SLP2)
        X_PATH_AUTO_L <- list(X_PATH_AUTO1, X_PATH_AUTO2)
      }
      # Add the results to the output list
      output_CL[[k]] <- list(KAPPA, X_MEAN, X_VAR, X_GF_LOADINGS, X_RES, XY_COV, X_PATH_L,
                             X_PATH_SLP_L, X_PATH_AUTO_L, abs_rate, outLag)
    }
  }
  return(output_CL)
}

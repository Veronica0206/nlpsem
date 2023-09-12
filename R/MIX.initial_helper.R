#' @title Compute Initial Values for Parameters of Mixture Models
#'
#' @description This function computes the initial values for the parameters for a longitudinal mixture model. The supported submodels
#' (i.e., class-specific models) include (1) latent growth curve models, (2) latent change score models, (3) latent growth
#' curve models or latent change score models with a time varying covariate, (4) multivariate latent growth curve models or
#' multivariate latent change score models, (5) longitudinal mediation models. For the first three submodels, time-invariant
#' covariates are allowed.
#'
#' @param dat A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
#' repeated measurements and occasions for each longitudinal process, and time-invariant covariates (TICs) if any. It takes the value
#' passed from \code{getMIX()}.
#' @param nClass An integer specifying the number of latent classes for the mixture model. It takes the value passed from \code{getMIX()}.
#' @param prop_starts A numeric vector of user-specified initial component proportions of latent classes. It takes the value passed from
#' \code{getMIX()}.
#' @param sub_Model A string that specifies the sub-model for latent classes. Supported sub-models include \code{"LGCM"} (for latent
#' growth curve models), \code{"LCSM"} (for latent change score models), \code{"TVC"} (for latent growth curve models or latent change
#' score models with a time-varying covariate), \code{"MGM"} (for multivariate latent growth curve models or latent change score models),
#' and \code{"MED"} (for longitudinal mediation models). It takes the value passed from \code{getMIX()}.
#' @param cluster_TIC A string or character vector representing the column name(s) for time-invariant covariate(s) indicating cluster
#' formations. It takes the value passed from \code{getMIX()}.
#' @param t_var A string specifying the prefix of the column names corresponding to the time variable for each study wave. This applies when
#' \code{sub_Model} is \code{"LGCM"}, \code{"LCSM"} or \code{"TVC"}. For \code{sub_Model} being \code{"MGM"} or \code{"MED"}, \code{t_var} should
#' be a string vector where each element corresponds to the time variable prefix for each respective longitudinal process. It takes the value
#' passed from \code{getMIX()}.
#' @param records A numeric vector denoting the indices of the observed study waves. This applies when \code{sub_Model} is \code{"LGCM"},
#' \code{"LCSM"} or \code{"TVC"}. For \code{sub_Model} being \code{"MGM"} or \code{"MED"}, \code{records} should be a list of numeric vectors,
#' where each vector provides the indices of the observed study waves for each longitudinal process. It takes the value passed from \code{getMIX()}.
#' @param y_var A string defining the prefix of the column names corresponding to the outcome variable for each study wave. This is applicable
#' when \code{sub_Model} is not \code{"MGM"}. For \code{sub_Model} being \code{"MGM"}, \code{y_var} should be a string vector where each element
#' corresponds to the prefix of the column names for each outcome variable across the study waves. It takes the value passed from \code{getMIX()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model = "LGCM"} include:
#' \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported options for
#' \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the value passed from \code{getMIX()}.
#' @param m_var A string that specifies the prefix of the column names corresponding to the mediator variable at each study wave.
#' It takes the value passed from \code{getMIX()}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. It takes the value passed from \code{getMIX()}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. It takes the value passed from \code{getMIX()}.
#' @param TVC A string that specifies the prefix of the column names corresponding to the time-varying covariate at each time
#' point. It takes the value passed from \code{getMIX()}.
#' @param decompose An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
#' decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
#' specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). It takes the value passed
#' from \code{getMIX()}.
#' @param growth_TIC A string or character vector of column names of time-invariant covariate(s) accounting for the variability
#' of growth factors if any. It takes the value passed from \code{getMIX()}.
#' @param res_scale A list where each element is a (vector of) numeric scaling factor(s) for residual variance to calculate the
#' corresponding initial value for a latent class, between \code{0} and \code{1} exclusive. It takes the value passed from \code{getMIX()}.
#' @param res_cor A list where each element is a (vector of) numeric initial value(s) for residual correlation in each class. It
#' needs to be specified if the sub_Model is \code{"TVC"} (when \code{decompose != 0}), \code{"MGM"}, or \code{"MED"}. It takes the value
#' passed from \code{getMIX()}.
#'
#' @return
#' A list containing initial values for each class in the specified model.
#'
#' @keywords internal
#'
#' @importFrom stats var cov lm na.exclude kmeans
#'
getMIX.initial <- function(dat, nClass, prop_starts, sub_Model, cluster_TIC, t_var, records, y_var, curveFun,
                           m_var, x_var, x_type, TVC, decompose, growth_TIC, res_scale, res_cor){
  # Initialize an empty list for starts
  starts <- list()
  # Create an ID variable based on the number of rows in dat
  ID <- 1:nrow(dat)
  # Case 1: LGCMs or LCSMs
  if (sub_Model %in% c("LGCM", "LCSM")){
    # Extract the measurement values and occasions
    dat_traj <- dat[, paste0(y_var, records)]
    dat_time <- dat[, paste0(t_var, records)]
    # Calculate the number of time points
    nT <- length(records)
    uni_GF <- getUNI.GF(dat_traj = dat_traj, dat_time = dat_time, nT = nT, curveFun = curveFun)
    if (curveFun %in% c("Jenss-Bayley", "JB")){
      uni_GF[, 3] <- abs(uni_GF[, 3])
    }
    out_Kmeans <- kmeans(x = uni_GF, centers = nClass, iter.max = 100)
    order_indices <- order(apply(out_Kmeans$centers, 1, mean))
    dat$label <- match(out_Kmeans$cluster, order_indices)
    for (k in 1:nClass){
      starts[[k]] <- getUNI.initial(dat = dat[dat$label == k, ], t_var = t_var, y_var = y_var, curveFun = curveFun,
                                    records = records, growth_TIC = growth_TIC, res_scale = res_scale[[k]])
    }
  }
  # Case 2: TVC model
  else if (sub_Model == "TVC"){
    # Extract the measurement values and occasions
    dat_traj <- dat[, paste0(y_var, records)]
    dat_time <- dat[, paste0(t_var, records)]
    # Calculate the number of time points
    nT <- length(records)
    growth_factor <- getUNI.GF(dat_traj = dat_traj, dat_time = dat_time, nT = nT, curveFun = curveFun)
    out_Kmeans <- kmeans(x = growth_factor, centers = nClass, iter.max = 100)
    order_indices <- order(out_Kmeans$centers[, 1])
    order_indices <- order(apply(out_Kmeans$centers, 1, mean))
    dat$label <- match(out_Kmeans$cluster, order_indices)
    for (k in 1:nClass){
      if (decompose == 0){
        starts[[k]] <- getTVC.initial(dat = dat[dat$label == k, ], t_var = t_var, y_var = y_var, curveFun = curveFun,
                                      records = records, growth_TIC = growth_TIC, TVC = TVC, decompose = decompose,
                                      res_scale = res_scale[[k]])
      }
      else if (decompose != 0){
        starts[[k]] <- getTVC.initial(dat = dat[dat$label == k, ], t_var = t_var, y_var = y_var, curveFun = curveFun,
                                      records = records, growth_TIC = growth_TIC, TVC = TVC, decompose = decompose,
                                      res_scale = res_scale[[k]], res_cor = res_cor[[k]])
      }
    }
  }
  # Case 3: MGM model
  else if (sub_Model == "MGM"){
    uni_GF <- list()
    for (traj in 1:length(y_var)){
      nT <- length(records[[traj]])
      dat_traj <- dat[, paste0(y_var[traj], records[[traj]])]
      dat_time <- dat[, paste0(t_var[traj], records[[traj]])]
      # Calculate individual growth factors for each longitudinal outcome using `getUNI.GF()` function
      uni_GF[[traj]] <- getUNI.GF(dat_traj = dat_traj, dat_time = dat_time, nT = nT, curveFun = curveFun)
    }
    uni_GF_dat <- do.call(cbind, uni_GF)
    out_Kmeans <- kmeans(x = uni_GF_dat, centers = nClass, iter.max = 100)
    order_indices <- order(apply(out_Kmeans$centers, 1, mean))
    dat$label <- match(out_Kmeans$cluster, order_indices)
    for (k in 1:nClass){
      starts[[k]] <- getMULTI.initial(dat = dat[dat$label == k, ], t_var = t_var, y_var = y_var,
                                      curveFun = curveFun, records = records, res_scale = res_scale[[k]],
                                      res_cor = res_cor[[k]])
    }
  }
  # Case 4: MED model
  else if (sub_Model == "MED"){
    if (x_type == "baseline"){
      traj_var <- c(y_var, m_var)
    }
    else if (x_type == "longitudinal"){
      traj_var <- c(y_var, m_var, x_var)
    }
    uni_GF <- list()
    for (traj in 1:length(traj_var)){
      nT <- length(records[[traj]])
      dat_traj <- dat[, paste0(traj_var[traj], records[[traj]])]
      dat_time <- dat[, paste0(t_var[traj], records[[traj]])]
      # Calculate individual growth factors for each longitudinal outcome using `getUNI.GF()` function
      uni_GF[[traj]] <- getUNI.GF(dat_traj = dat_traj, dat_time = dat_time, nT = nT, curveFun = curveFun)
    }
    uni_GF_dat <- do.call(cbind, uni_GF)
    out_Kmeans <- kmeans(x = uni_GF_dat, centers = nClass, iter.max = 100)
    order_indices <- order(apply(out_Kmeans$centers, 1, mean))
    dat$label <- match(out_Kmeans$cluster, order_indices)
    for (k in 1:nClass){
      starts[[k]] <- getMED.initial(dat = dat[dat$label == k, ], t_var = t_var, y_var = y_var, m_var = m_var,
                                    x_var = x_var, x_type = x_type, curveFun = curveFun, records = records,
                                    res_scale = res_scale[[k]], res_cor = res_cor[[k]])
    }
  }
  # Check if cluster_TIC is not null
  if (!is.null(cluster_TIC)){
    dat_nnet <- dat[dat$label != 0, c("label", cluster_TIC)]
    mod_nnet <- nnet::multinom(label ~ ., data = dat_nnet)
    logit_starts <- (as.matrix(rbind(rep(0, length(cluster_TIC) + 1), summary(mod_nnet)$coefficient)))
    rownames(logit_starts) <- colnames(logit_starts) <- NULL
    starts[[length(starts) + 1]] <- logit_starts
  }
  # If cluster_TIC is null
  else if (is.null(cluster_TIC)){
    # Return the starts list
    starts[[length(starts) + 1]] <- prop_starts/prop_starts[1]
  }
  return(starts)
}

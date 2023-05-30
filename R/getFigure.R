#' @title Generate Visualization for Fitted Model
#'
#' @description This function generates visualizations for the output of a fitted model. When a Latent Growth Curve Model
#' (LGCM) is fitted for the longitudinal process, it provides (class-specific) estimated growth status with
#' 95% confidence intervals. When a Latent Change Score Model (LCSM) is fitted for the longitudinal process,
#' it provides (class-specific) estimated growth rate with 95% confidence intervals and change from baseline
#' with 95% confidence intervals. These visualizations are particularly useful for understanding the results
#' and trajectories of different classes or groups within the model.
#'
#' @param model A fitted mxModel object. This is the output from one of the estimation functions in this
#' package.
#' @param nClass An integer specifying the number of classes for the mixture model or multiple group model. Default
#' is \code{NULL}, indicating a single-group model.
#' @param cluster_TIC A string or character vector representing the column name(s) for time-invariant covariate(s)
#' indicating cluster formations. Default is \code{NULL}, indicating no such time-invariant covariates are present
#' in the model.
#' @param grp_var A string specifying the column that indicates manifested classes when applicable.
#' @param sub_Model A string that specifies the sub-model for latent classes. Supported sub-models include \code{"LGCM"} (for latent
#' growth curve models), \code{"LCSM"} (for latent change score models), \code{"TVC"} (for latent growth curve models or latent change
#' score models with a time-varying covariate), \code{"MGM"} (for multivariate latent growth curve models or latent change score models),
#' and \code{"MED"} (for longitudinal mediation models).
#' @param y_var A string or character vector representing the prefix of the column names for the outcome variable(s)
#' at each study wave.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model = "LGCM"} include:
#' \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported options for
#' \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}).
#' @param y_model A string that specifies how to fit longitudinal outcomes. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' By default, this is \code{NULL} as this argument only requires when \code{sub_Model} is \code{"TVC"} or \code{"MGM"}.
#' @param t_var A string representing the prefix of the column names corresponding to the time variable at each study
#' wave.
#' @param records A numeric vector representing the indices of the study waves.
#' @param m_var A string that specifies the prefix of the column names corresponding to the mediator variable at each time point.
#' Default is \code{NULL} as this argument only requires when \code{sub_Model} is \code{"MED"}.
#' @param x_var A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
#' corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. Default is \code{NULL} as this
#' argument only requires when \code{sub_Model} is \code{"MED"}.
#' @param x_type A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
#' and \code{"longitudinal"}. Default is \code{NULL} as this argument only requires when \code{sub_Model} is \code{"MED"}.
#' @param xstarts A numeric value to indicate the starting time of the longitudinal process.
#' @param xlab A string representing the time unit (e.g., "Week", "Month", or "Year") for the x-axis. Default is
#' "Time".
#' @param outcome A string or character vector representing the name(s) of the longitudinal process(es) under examination.
#'
#' @return A ggplot object or a list of ggplot objects, each representing a figure for the fitted model. If a list of
#' ggplot objects is returned, it can be visualized using the \code{print} function.
#'
#' @export
#'
#' @examples
#' OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' # Load ECLS-K (2011) data
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' # Re-baseline the data so that the estimated initial status is for the starting point of the study
#' baseT <- RMS_dat0$T1
#' RMS_dat0$T1 <- RMS_dat0$T1 - baseT
#' RMS_dat0$T2 <- RMS_dat0$T2 - baseT
#' RMS_dat0$T3 <- RMS_dat0$T3 - baseT
#' RMS_dat0$T4 <- RMS_dat0$T4 - baseT
#' RMS_dat0$T5 <- RMS_dat0$T5 - baseT
#' RMS_dat0$T6 <- RMS_dat0$T6 - baseT
#' RMS_dat0$T7 <- RMS_dat0$T7 - baseT
#' RMS_dat0$T8 <- RMS_dat0$T8 - baseT
#' RMS_dat0$T9 <- RMS_dat0$T9 - baseT
#' xstarts <- mean(baseT)
#'
#' \donttest{
#' # Plot single group LGCM model
#' set.seed(20191029)
#' BLS_LGCM1 <- getLGCM(dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS",
#'                      intrinsic = FALSE, records = 1:9, res_scale = 0.1)
#' Figure1 <- getFigure(
#'   model = BLS_LGCM1, nClass = NULL, cluster_TIC = NULL, sub_Model = "LGCM",
#'   y_var = "M", curveFun = "BLS", y_model = "LGCM", t_var = "T", records = 1:9,
#'   m_var = NULL, x_var = NULL, x_type = NULL, xstarts = xstarts, xlab = "Month",
#'   outcome = "Mathematics"
#' )
#' print(Figure1)
#' # Plot mixture LGCM model
#' BLS_LGCM2 <-  getMIX(
#'   dat = RMS_dat0, prop_starts = c(0.45, 0.55), sub_Model = "LGCM",
#'   cluster_TIC = NULL, y_var = "M", t_var = "T", records = 1:9,
#'   curveFun = "BLS", intrinsic = FALSE, res_scale = list(0.3, 0.3)
#' )
#' Figure2 <- getFigure(
#'   model = BLS_LGCM2, nClass = 2, cluster_TIC = NULL, sub_Model = "LGCM",
#'   y_var = "M", curveFun = "BLS", y_model = "LGCM", t_var = "T", records = 1:9,
#'   m_var = NULL, x_var = NULL, x_type = NULL, xstarts = xstarts, xlab = "Month",
#'   outcome = "Mathematics"
#' )
#' print(Figure2)
#' }
#'
getFigure <- function(model, nClass = NULL, cluster_TIC = NULL, grp_var = NULL, sub_Model, y_var,
                      curveFun, y_model = NULL, t_var, records, m_var = NULL, x_var = NULL, x_type = NULL,
                      xstarts, xlab = "Time", outcome = "Process"){
  if (I(sub_Model == "MGM" & length(y_var) != length(outcome))){
    stop("Please ensure provide a name for each longitudinal process!")
  }
  if (sub_Model == "MED"){
    if (I(x_type == "baseline" & length(outcome) != 2)){
      stop("Please ensure provide a name for both longitudinal mediator and longitudinal outcome!")
    }
    else if (I(x_type == "longitudinal" & length(outcome) != 3)){
      stop("Please ensure provide a name for longitudinal predictor, longitudinal mediator and longitudinal outcome!")
    }
  }
  if (is.null(nClass)){
    if (sub_Model %in% c("LGCM", "LCSM")){
      figures <- getFitFig(model = model, nClass = nClass, cluster_TIC = cluster_TIC, grp_var = grp_var,
                           sub_Model = sub_Model, t_var = t_var, records = records, y_var = y_var,
                           curveFun = curveFun, y_model = sub_Model, xstarts = xstarts, xlab = xlab,
                           outcome = outcome)
      return(figures)
    }
    else if (sub_Model == "TVC"){
      figures <- getFitFig(model = model, nClass = nClass, cluster_TIC = cluster_TIC, grp_var = grp_var,
                            sub_Model = sub_Model, t_var = t_var, records = records,
                           y_var = y_var, curveFun = curveFun, y_model = y_model, xstarts = xstarts, xlab = xlab,
                           outcome = outcome)
      return(figures)
    }
    else if (sub_Model == "MGM"){
      figures_L <- list()
      for (traj in 1:length(y_var)){
        figures_L[[traj]] <- getFitFig(model = model, nClass = nClass, cluster_TIC = cluster_TIC,
                                       grp_var = grp_var, sub_Model = sub_Model, t_var = t_var[traj],
                                       records = records[[traj]], y_var = y_var[traj], curveFun = curveFun,
                                       y_model = y_model, xstarts = xstarts, xlab = xlab, outcome = outcome)
      }
      return(figures_L)
    }
    else if (sub_Model == "MED"){
      figures_L <- list()
      if (x_type == "baseline"){
        traj_var <- c(y_var, m_var)
        for (traj in 1:length(y_var)){
          figures_L[[traj]] <- getFitFig(model = model, nClass = nClass, cluster_TIC = cluster_TIC,
                                         grp_var = grp_var, sub_Model = sub_Model, t_var = t_var[traj],
                                         records = records[[traj]], y_var = traj_var[traj], curveFun = curveFun,
                                         y_model = "LGCM", xstarts = xstarts, xlab = xlab, outcome = outcome)
        }
      }
      else if (x_type == "longitudinal"){
        traj_var <- c(y_var, m_var, x_var)
        for (traj in 1:length(y_var)){
          figures_L[[traj]] <- getFitFig(model = model, nClass = nClass, cluster_TIC = cluster_TIC,
                                         grp_var = grp_var, sub_Model = sub_Model, t_var = t_var[traj],
                                         records = records[[traj]], y_var = traj_var[traj], curveFun = curveFun,
                                         y_model = "LGCM", xstarts = xstarts, xlab = xlab, outcome = outcome)
        }
      }
    }
  }
  else if (!is.null(nClass)){
    if (sub_Model %in% c("LGCM", "LCSM")){
      figures <- getFitFig(model = model, nClass = nClass, cluster_TIC = cluster_TIC, grp_var = grp_var,
                            sub_Model = sub_Model, t_var = t_var, records = records,
                           y_var = y_var, curveFun = curveFun, y_model = sub_Model, xstarts = xstarts, xlab = xlab,
                           outcome = outcome)
      return(figures)
    }
    else if (sub_Model == "TVC"){
      figures <- getFitFig(model = model, nClass = nClass, cluster_TIC = cluster_TIC, grp_var = grp_var,
                           sub_Model = sub_Model, t_var = t_var, records = records, y_var = y_var, curveFun = curveFun,
                           y_model = y_model, xstarts = xstarts, xlab = xlab, outcome = outcome)
      return(figures)
    }
    else if (sub_Model == "MGM"){
      figures_L <- list()
      for (traj in 1:length(y_var)){
        figures_L[[traj]] <- getFitFig(model = model, nClass = nClass, cluster_TIC = cluster_TIC, grp_var = grp_var,
                                       sub_Model = sub_Model, t_var = t_var[traj],
                                       records = records[[traj]], y_var = y_var[traj], curveFun = curveFun,
                                       y_model = y_model, xstarts = xstarts, xlab = xlab, outcome = outcome)
      }
      return(figures_L)
    }
    else if (sub_Model == "MED"){
      figures_L <- list()
      if (x_type == "baseline"){
        traj_var <- c(y_var, m_var)
        for (traj in 1:length(y_var)){
          figures_L[[traj]] <- getFitFig(model = model, nClass = nClass, cluster_TIC = cluster_TIC, grp_var = grp_var,
                                         sub_Model = sub_Model, t_var = t_var[traj], records = records[[traj]],
                                         y_var = traj_var[traj], curveFun = curveFun, y_model = "LGCM", xstarts = xstarts,
                                         xlab = xlab, outcome = outcome)
        }
      }
      else if (x_type == "longitudinal"){
        traj_var <- c(y_var, m_var, x_var)
        for (traj in 1:length(y_var)){
          figures_L[[traj]] <- getFitFig(model = model, nClass = nClass, cluster_TIC = cluster_TIC, grp_var = grp_var,
                                         sub_Model = sub_Model, t_var = t_var[traj], records = records[[traj]],
                                         y_var = traj_var[traj], curveFun = curveFun, y_model = "LGCM", xstarts = xstarts,
                                         xlab = xlab, outcome = outcome)
        }
      }
    }
  }
  return(figures_L)
}


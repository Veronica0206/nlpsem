#' @title Helper Function to Generate Visualization for a Fitted Model
#'
#' @description
#' This is an internal function that generates a ggplot object for a fitted model. It is called by the \code{getFigure}
#' function.
#'
#' @param model A fitted mxModel object. This is the output from one of the estimation functions in this package.
#' It takes value passed from \code{getFigure()}.
#' @param nClass An integer specifying the number of classes for the mixture model or multiple group model. It
#' takes value passed from \code{getFigure()}.
#' @param cluster_TIC A string or character vector representing the column name(s) for time-invariant covariate(s)
#' indicating cluster formations. It takes value passed from \code{getFigure()}.
#' @param grp_var A string specifying the column that indicates manifested classes when applicable. It takes the value
#' passed from \code{getFigure()}.
#' @param sub_Model A string that specifies the sub-model for latent classes. Supported sub-models include \code{"LGCM"} (for latent
#' growth curve models), \code{"LCSM"} (for latent change score models), \code{"TVC"} (for latent growth curve models or latent change
#' score models with a time-varying covariate), \code{"MGM"} (for multivariate latent growth curve models or latent change score models),
#' and \code{"MED"} (for longitudinal mediation models). It takes value passed from \code{getFigure()}.
#' @param t_var A string representing the prefix of the column names corresponding to the time variable at each study
#' wave. It takes value passed from \code{getFigure()}.
#' @param records A numeric vector representing the indices of the study waves. It takes value passed from \code{getFigure()}.
#' @param y_var A string or character vector representing the prefix of the column names for the outcome variable(s)
#' at each study wave. It takes value passed from \code{getFigure()}.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model = "LGCM"} include:
#' \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported options for
#' \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}),
#' \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes value passed from \code{getFigure()}.
#' @param y_model A string that specifies how to fit longitudinal outcomes. Supported values are \code{"LGCM"} and \code{"LCSM"}.
#' It takes value passed from \code{getFigure()}.
#' @param xstarts A numeric value to indicate the starting time of the longitudinal process. It takes value passed from \code{getFigure()}.
#' @param xlab A string representing the time unit (e.g., "Week", "Month", or "Year") for the x-axis. Default is
#' "Time". It takes value passed from \code{getFigure()}.
#' @param outcome A string or character vector representing the name(s) of the longitudinal process(es) under examination.
#' It takes value passed from \code{getFigure()}.
#'
#' @return A ggplot object or a list of ggplot objects.
#'
#' @keywords internal
#'
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_smooth labs scale_linetype_manual scale_x_continuous guides guide_legend theme_bw theme element_text element_rect geom_step scale_color_manual margin
#' @importFrom grDevices colorRampPalette
#'
getFitFig <- function(model, nClass, cluster_TIC, grp_var, sub_Model, t_var, records, y_var, curveFun, y_model,
                      xstarts, xlab, outcome){
  if (!is.null(grp_var)){
    model_names <- names(model@submodels)
    sub_dat <- lapply(model_names, function(name) {
      as.data.frame(model[[name]]$data$observed)
    })
    dat <- do.call(rbind, sub_dat)
  }
  else if (is.null(grp_var)){
    dat <- as.data.frame(model@data$observed)
  }
  dat$ID <- 1:nrow(dat)
  if (!is.null(nClass)){
    if (is.null(grp_var)){
      dat$Class <- getPosterior(model = model, nClass = nClass, cluster_TIC = cluster_TIC)[[2]]
      class_props <- table(dat$Class)/length(dat$Class)
      class_labels <- paste0("Class ", seq_len(nClass), " (", round(class_props * 100, 1), "%)")
    }
    else if (!is.null(grp_var)){
      class_props <- table(dat[, grp_var])/length(dat[, grp_var])
      class_labels <- paste0("Group ", seq_len(nClass), " (", round(class_props * 100, 1), "%)")
    }
  }
  if (is.null(nClass)){
    if (y_model == "LGCM"){
      dat_traj <- dat[, c("ID", paste0(y_var, records))]
      if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
        y_var <- "Y"
      }
      dat_time <- dat[, c("ID", paste0(t_var, records))]
      dat_traj_l <- pivot_longer(dat_traj, cols = -ID, names_to = "Wave", values_to = "value")
      dat_traj_l$Wave <- substring(dat_traj_l$Wave, 2)
      dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
      dat_time_l$time <- dat_time_l$time + xstarts
      dat_time_l$Wave <- substring(dat_time_l$Wave, 2)
      dat_raw <- merge(dat_traj_l, dat_time_l, by = c("ID", "Wave"))
      dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
      t_seq <- apply(dat[, paste0(t_var, records)], 2, mean)
      Y_mean0 <- c(mxEvalByName(paste0(y_var, "_mean0"), model = model))
      Y_mean0_se <- c(mxSE(paste0(y_var, "_mean0"), model, forceName = T))
      Y_mean0_u <- Y_mean0 + 1.96 * Y_mean0_se
      Y_mean0_l <- Y_mean0 - 1.96 * Y_mean0_se
      if (curveFun %in% c("linear", "LIN")){
        Y_hat <- Y_mean0[1] + Y_mean0[2] * t_seq
        Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * t_seq
        Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * t_seq
      }
      else if (curveFun %in% c("quadratic", "QUAD")){
        Y_hat <- Y_mean0[1] + Y_mean0[2] * t_seq + Y_mean0[3] * t_seq^2
        Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * t_seq + Y_mean0_u[3] * t_seq^2
        Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * t_seq + Y_mean0_l[3] * t_seq^2
      }
      else if (curveFun %in% c("exponential", "EXP")){
        Y_mean0 <- if (length(Y_mean0) == 3){
          Y_mean0
        } else {
          c(Y_mean0, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[3])
        }
        Y_mean0_u <- if (length(Y_mean0_u) == 3){
          Y_mean0_u
        } else {
          c(Y_mean0_u, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[3] + 1.96 * mxSE(paste0(y_var, "_alpha0"), model, forceName = T)[3])
        }
        Y_mean0_l <- if (length(Y_mean0_l) == 3){
          Y_mean0_l
        } else {
          c(Y_mean0_l, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[3] - 1.96 * mxSE(paste0(y_var, "_alpha0"), model, forceName = T)[3])
        }
        Y_hat <- Y_mean0[1] + Y_mean0[2] * (1 - exp(-Y_mean0[3] * t_seq))
        Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * (1 - exp(-Y_mean0_u[3] * t_seq))
        Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * (1 - exp(-Y_mean0_l[3] * t_seq))
      }
      else if (curveFun %in% c("Jenss-Bayley", "JB")){
        Y_mean0 <- if (length(Y_mean0) == 4){
          Y_mean0
        } else {
          c(Y_mean0, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[4])
        }
        Y_mean0_u <- if (length(Y_mean0_u) == 4){
          Y_mean0_u
        } else {
          c(Y_mean0_u, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[4] + 1.96 * mxSE(paste0(y_var, "_alpha0"), model, forceName = T)[4])
        }
        Y_mean0_l <- if (length(Y_mean0_l) == 4){
          Y_mean0_l
        } else {
          c(Y_mean0_l, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[4] - 1.96 * mxSE(paste0(y_var, "_alpha0"), model, forceName = T)[4])
        }
        Y_hat <- Y_mean0[1] + Y_mean0[2] * t_seq + Y_mean0[3] * (exp(Y_mean0[4] * t_seq) - 1)
        Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * t_seq + Y_mean0_u[3] * (exp(Y_mean0_u[4] * t_seq) - 1)
        Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * t_seq + Y_mean0_l[3] * (exp(Y_mean0_l[4] * t_seq) - 1)
      }
      else if (curveFun %in% c("bilinear spline", "BLS")){
        Y_mean0 <- if (length(Y_mean0) == 4){
          Y_mean0
        } else {
          c(Y_mean0, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[4])
        }
        Y_mean0_u <- if (length(Y_mean0_u) == 4){
          Y_mean0_u
        } else {
          c(Y_mean0_u, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[4] + 1.96 * mxSE(paste0(y_var, "_alpha0"), model, forceName = T)[4])
        }
        Y_mean0_l <- if (length(Y_mean0_l) == 4){
          Y_mean0_l
        } else {
          c(Y_mean0_l, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[4] - 1.96 * mxSE(paste0(y_var, "_alpha0"), model, forceName = T)[4])
        }
        Y_hat <- ifelse(t_seq <= Y_mean0[4], Y_mean0[1] + Y_mean0[2] * t_seq,
                        Y_mean0[1] + Y_mean0[2] * Y_mean0[4] + Y_mean0[3] * (t_seq - Y_mean0[4]))
        Y_hat_u <- ifelse(t_seq <= Y_mean0_u[4], Y_mean0_u[1] + Y_mean0_u[2] * t_seq,
                          Y_mean0_u[1] + Y_mean0_u[2] * Y_mean0_u[4] + Y_mean0_u[3] * (t_seq - Y_mean0_u[4]))
        Y_hat_l <- ifelse(t_seq <= Y_mean0_l[4], Y_mean0_l[1] + Y_mean0_l[2] * t_seq,
                          Y_mean0_l[1] + Y_mean0_l[2] * Y_mean0_l[4] + Y_mean0_l[3] * (t_seq - Y_mean0_l[4]))
      }
      dat_est <- data.frame(time = t_seq, Y_hat = Y_hat, Y_hat_u = Y_hat_u, Y_hat_l = Y_hat_l)
      dat_est$time <- dat_est$time + xstarts
      lines <- c("Model Implied Growth Status" = "solid",
                 "Smooth Line of Observed Growth Status" = "twodash",
                 "95% Confidence Interval on Model Implied Growth Status" = "dotdash")
      fig.status <- ggplot(data = dat_est, aes(x = time, y = Y_hat, group = 1)) +
        geom_line(aes(linetype = "Model Implied Growth Status"),  color = "black", size = 1) +
        geom_line(aes(x = time, y = Y_hat_l, linetype = "95% Confidence Interval on Model Implied Growth Status"), size = 1) +
        geom_line(aes(x = time, y = Y_hat_u, linetype = "95% Confidence Interval on Model Implied Growth Status"), size = 1) +
        geom_smooth(data = dat_raw, aes(x = time, y = value, group = 1,
                                                          linetype = "Smooth Line of Observed Growth Status"),
                             color = "black", linewidth = 1, se = F) +
        labs(x = xlab, y = paste0("Measurement of ", outcome)) +
        scale_linetype_manual(values = c("dotdash", "solid", "twodash")) +
        scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
        guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3)) +
        theme_bw() +
        theme(strip.text.y = element_text(size = 4),
                       strip.background = element_rect(color = "white", fill = "white"),
                       strip.placement = "outside",
                       legend.position = "bottom",
                       legend.box = "vertical")
      figure <- list(fig.status)
    }
    else if (y_model == "LCSM"){
      dat_chg <- cbind(dat[, "ID"],
                       dat[, c(paste0(y_var, records[-1]))] - dat[, c(paste0(y_var, records[1]))])
      if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
        y_var <- "Y"
      }
      names(dat_chg) <- c("ID", paste0("Delta", names(dat_chg)[-1]))
      dat_chg_l <- pivot_longer(dat_chg, cols = -ID, names_to = "Wave", values_to = "value")
      dat_chg_l$Wave <- substring(dat_chg_l$Wave, 7)
      dat_time <- dat[, c("ID", paste0(t_var, records[-1]))]
      dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
      dat_time_l$time <- dat_time_l$time + xstarts
      dat_time_l$Wave <- substring(dat_time_l$Wave, 2)
      dat_raw <- merge(dat_chg_l, dat_time_l, by = c("ID", "Wave"))
      dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
      t_seq <- apply(dat[, paste0(t_var, records[-1])], 2, mean)
      Y_mean0 <- c(mxEvalByName(paste0(y_var, "_mean0"), model = model))
      Y_mean0_se <- c(mxSE(paste0(y_var, "_mean0"), model, forceName = T))
      Y_mean0_u <- Y_mean0 + 1.96 * Y_mean0_se
      Y_mean0_l <- Y_mean0 - 1.96 * Y_mean0_se
      if (curveFun %in% c("nonparametric", "NonP")){
        rel_rate <- c(1, model@output$estimate[grep(paste0(y_var, "_rel_rate"), names(model@output$estimate))])
        dY_hat <- Y_mean0[2] * rel_rate
        dY_hat_u <- Y_mean0_u[2] * rel_rate
        dY_hat_l <- Y_mean0_l[2] * rel_rate
        cY_hat <- mxEvalByName(paste0(y_var, "chg_bl_m"), model = model)
        cY_hat_u <- mxEvalByName(paste0(y_var, "chg_bl_m"), model = model) +
          1.96 * mxSE(paste0(y_var, "chg_bl_m"), model, forceName = T)
        cY_hat_l <- mxEvalByName(paste0(y_var, "chg_bl_m"), model = model) -
          1.96 * mxSE(paste0(y_var, "chg_bl_m"), model, forceName = T)
      }
      else if (curveFun %in% c("quadratic", "QUAD")){
        dY_hat <- Y_mean0[2] + 2 * Y_mean0[3] * t_seq
        dY_hat_u <- Y_mean0_u[2] + 2 * Y_mean0_u[3] * t_seq
        dY_hat_l <- Y_mean0_l[2] + 2 * Y_mean0_l[3] * t_seq
        cY_hat <- Y_mean0[2] * t_seq + Y_mean0[3] * t_seq^2
        cY_hat_u <- Y_mean0_u[2] * t_seq + Y_mean0_u[3] * t_seq^2
        cY_hat_l <- Y_mean0_l[2] * t_seq + Y_mean0_l[3] * t_seq^2
      }
      else if (curveFun %in% c("exponential", "EXP")){
        Y_mean0 <- if (length(Y_mean0) == 3){
          Y_mean0
        } else {
          c(Y_mean0, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[3])
        }
        Y_mean0_u <- if (length(Y_mean0_u) == 3){
          Y_mean0_u
        } else {
          c(Y_mean0_u, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[3] + 1.96 * mxSE(paste0(y_var, "_alpha0"), model, forceName = T)[3])
        }
        Y_mean0_l <- if (length(Y_mean0_l) == 3){
          Y_mean0_l
        } else {
          c(Y_mean0_l, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[3] - 1.96 * mxSE(paste0(y_var, "_alpha0"), model, forceName = T)[3])
        }
        dY_hat <- Y_mean0[2] * Y_mean0[3] * exp(-Y_mean0[3] * t_seq)
        dY_hat_u <- Y_mean0_u[2] * Y_mean0_u[3] * exp(-Y_mean0_u[3] * t_seq)
        dY_hat_l <- Y_mean0_l[2] * Y_mean0_l[3] * exp(-Y_mean0_l[3] * t_seq)
        cY_hat <- Y_mean0[2] * (1 - exp(-Y_mean0[3] * t_seq))
        cY_hat_u <- Y_mean0_u[2] * (1 - exp(-Y_mean0_u[3] * t_seq))
        cY_hat_l <- Y_mean0_l[2] * (1 - exp(-Y_mean0_l[3] * t_seq))
      }
      else if (curveFun %in% c("Jenss-Bayley", "JB")){
        Y_mean0 <- if (length(Y_mean0) == 4){
          Y_mean0
        } else {
          c(Y_mean0, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[4])
        }
        Y_mean0_u <- if (length(Y_mean0_u) == 4){
          Y_mean0_u
        } else {
          c(Y_mean0_u, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[4] + 1.96 * mxSE(paste0(y_var, "_alpha0"), model, forceName = T)[4])
        }
        Y_mean0_l <- if (length(Y_mean0_l) == 4){
          Y_mean0_l
        } else {
          c(Y_mean0_l, mxEvalByName(paste0(y_var, "_alpha0"), model = model)[4] - 1.96 * mxSE(paste0(y_var, "_alpha0"), model, forceName = T)[4])
        }
        dY_hat <- Y_mean0[2] + Y_mean0[3] * Y_mean0[4] * exp(Y_mean0[4] * t_seq)
        dY_hat_u <- Y_mean0_u[2] + Y_mean0_u[3] * Y_mean0_u[4] * exp(Y_mean0_u[4] * t_seq)
        dY_hat_l <- Y_mean0_l[2] + Y_mean0_l[3] * Y_mean0_l[4] * exp(Y_mean0_l[4] * t_seq)
        cY_hat <- Y_mean0[2] * t_seq + Y_mean0[3] * (exp(Y_mean0[4] * t_seq) - 1)
        cY_hat_u <- Y_mean0_u[2] * t_seq + Y_mean0_u[3] * (exp(Y_mean0_u[4] * t_seq) - 1)
        cY_hat_l <- Y_mean0_l[2] * t_seq + Y_mean0_l[3] * (exp(Y_mean0_l[4] * t_seq) - 1)
      }
      dat_est <- data.frame(time = t_seq, dY_hat = dY_hat, dY_hat_l = dY_hat_l, dY_hat_u = dY_hat_u,
                            cY_hat = cY_hat, cY_hat_l = cY_hat_l, cY_hat_u = cY_hat_u)
      dat_est$time <- dat_est$time + xstarts
      lines <- c("Model Implied Change from Baseline" = "solid",
                 "Smooth Line of Observed Change from Baseline" = "twodash",
                 "95% Confidence Interval on Model Implied Change from Baseline" = "dotdash")
      fig.CHG_BL <- ggplot(data = dat_est, aes(x = time, y = cY_hat, group = 1)) +
        geom_line(aes(linetype = "Model Implied Change from Baseline"),  color = "black", size = 1) +
        geom_line(aes(x = time, y = cY_hat_l, linetype = "95% Confidence Interval on Model Implied Change from Baseline"), size = 1) +
        geom_line(aes(x = time, y = cY_hat_u, linetype = "95% Confidence Interval on Model Implied Change from Baseline"), size = 1) +
        geom_smooth(data = dat_raw, aes(x = time, y = value, group = 1,
                                                          linetype = "Smooth Line of Observed Change from Baseline"),
                             color = "black", size = 1, se = F) +
        labs(x = xlab, y = paste0("Change from baseline of ", outcome)) +
        scale_linetype_manual(values = c("dotdash", "solid", "twodash")) +
        scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
        guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3)) +
        theme_bw() +
        theme(strip.text.y = element_text(size = 4),
                       strip.background = element_rect(color = "white", fill = "white"),
                       strip.placement = "outside",
                       legend.position = "bottom",
                       legend.box = "vertical")
      lines <- c("Model Implied Change from Baseline" = "solid",
                 "95% Confidence Interval on Model Implied Change from Baseline" = "dotdash")
      if (curveFun %in% c("nonparametric", "NonP")){
        fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = 1)) +
          geom_step(aes(linetype = "Model Implied Growth Rate"),  color = "black", size = 1) +
          geom_step(aes(x = time, y = dY_hat_l, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
          geom_step(aes(x = time, y = dY_hat_u, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
          labs(x = xlab, y = paste0("Growth Rate of ", outcome)) +
          scale_linetype_manual(values = c("dotdash", "solid")) +
          scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
          guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2)) +
          theme_bw() +
          theme(strip.text.y = element_text(size = 4),
                strip.background = element_rect(color = "white", fill = "white"),
                strip.placement = "outside",
                legend.position = "bottom",
                legend.box = "vertical")
      }
      else{
        fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = 1)) +
          geom_line(aes(linetype = "Model Implied Growth Rate"),  color = "black", size = 1) +
          geom_line(aes(x = time, y = dY_hat_l, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
          geom_line(aes(x = time, y = dY_hat_u, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
          labs(x = xlab, y = paste0("Growth Rate of ", outcome)) +
          scale_linetype_manual(values = c("dotdash", "solid")) +
          scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
          guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2)) +
          theme_bw() +
          theme(strip.text.y = element_text(size = 4),
                strip.background = element_rect(color = "white", fill = "white"),
                strip.placement = "outside",
                legend.position = "bottom",
                legend.box = "vertical")
      }
      figure <- list(fig.CHG_BL, fig.SLP)
    }
  }
  else if (!is.null(nClass)){
    if (is.null(grp_var)){
        if (y_model == "LGCM"){
          dat$Class <- getPosterior(model = model, nClass = nClass, cluster_TIC = cluster_TIC)[[2]]
          dat_traj <- dat[, c("ID", "Class", paste0(y_var, records))]
          if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
            y_var <- "Y"
          }
          dat_time <- dat[, c("ID", paste0(t_var, records))]
          dat_traj_l <- pivot_longer(dat_traj, cols = -c(ID, Class), names_to = "Wave", values_to = "value")
          dat_traj_l$Wave <- substring(dat_traj_l$Wave, 2)
          dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
          dat_time_l$time <- dat_time_l$time + xstarts
          dat_time_l$Wave <- substring(dat_time_l$Wave, 2)
          dat_raw <- merge(dat_traj_l, dat_time_l, by = c("ID", "Wave"))
          dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
          t_seq <- apply(dat[, paste0(t_var, records)], 2, mean)
          dat_est_L <- list()
          for (k in 1:nClass){
            Y_mean0 <- mxEvalByName(paste0("c", k, y_var, "_mean0"), model = model@submodels[[k]])
            Y_mean0_se <- mxSE(paste0("Class", k, ".c", k, y_var, "_mean0"), model, forceName = T)
            Y_mean0_u <- Y_mean0 + 1.96 * Y_mean0_se
            Y_mean0_l <- Y_mean0 - 1.96 * Y_mean0_se
            if (curveFun %in% c("linear", "LIN")){
              Y_hat <- Y_mean0[1] + Y_mean0[2] * t_seq
              Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * t_seq
              Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * t_seq
            }
            else if (curveFun %in% c("quadratic", "QUAD")){
              Y_hat <- Y_mean0[1] + Y_mean0[2] * t_seq + Y_mean0[3] * t_seq^2
              Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * t_seq + Y_mean0_u[3] * t_seq^2
              Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * t_seq + Y_mean0_l[3] * t_seq^2
            }
            else if (curveFun %in% c("exponential", "EXP")){
              Y_mean0 <- if (length(Y_mean0) == 3){
                Y_mean0
              } else {
                c(Y_mean0, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3])
              }
              Y_mean0_u <- if (length(Y_mean0_u) == 3){
                Y_mean0_u
              } else {
                c(Y_mean0_u, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3] +
                    1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[3])
              }
              Y_mean0_l <- if (length(Y_mean0_l) == 3){
                Y_mean0_l
              } else {
                c(Y_mean0_l, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3] -
                    1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[3])
              }
              Y_hat <- Y_mean0[1] + Y_mean0[2] * (1 - exp(-Y_mean0[3] * t_seq))
              Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * (1 - exp(-Y_mean0_u[3] * t_seq))
              Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * (1 - exp(-Y_mean0_l[3] * t_seq))
            }
            else if (curveFun %in% c("Jenss-Bayley", "JB")){
              Y_mean0 <- if (length(Y_mean0) == 4){
                Y_mean0
              } else {
                c(Y_mean0, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4])
              }
              Y_mean0_u <- if (length(Y_mean0_u) == 4){
                Y_mean0_u
              } else {
                c(Y_mean0_u, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] +
                    1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
              }
              Y_mean0_l <- if (length(Y_mean0_l) == 4){
                Y_mean0_l
              } else {
                c(Y_mean0_l, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] -
                    1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
              }
              Y_hat <- Y_mean0[1] + Y_mean0[2] * t_seq + Y_mean0[3] * (exp(Y_mean0[4] * t_seq) - 1)
              Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * t_seq + Y_mean0_u[3] * (exp(Y_mean0_u[4] * t_seq) - 1)
              Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * t_seq + Y_mean0_l[3] * (exp(Y_mean0_l[4] * t_seq) - 1)
            }
            else if (curveFun %in% c("bilinear spline", "BLS")){
              Y_mean0 <- if (length(Y_mean0) == 4){
                Y_mean0
              } else {
                c(Y_mean0, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4])
              }
              Y_mean0_u <- if (length(Y_mean0_u) == 4){
                Y_mean0_u
              } else {
                c(Y_mean0_u, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] +
                    1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
              }
              Y_mean0_l <- if (length(Y_mean0_l) == 4){
                Y_mean0_l
              } else {
                c(Y_mean0_l, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] -
                    1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
              }
              Y_hat <- ifelse(t_seq <= Y_mean0[4], Y_mean0[1] + Y_mean0[2] * t_seq,
                              Y_mean0[1] + Y_mean0[2] * Y_mean0[4] + Y_mean0[3] * (t_seq - Y_mean0[4]))
              Y_hat_u <- ifelse(t_seq <= Y_mean0_u[4], Y_mean0_u[1] + Y_mean0_u[2] * t_seq,
                                Y_mean0_u[1] + Y_mean0_u[2] * Y_mean0_u[4] + Y_mean0_u[3] * (t_seq - Y_mean0_u[4]))
              Y_hat_l <- ifelse(t_seq <= Y_mean0_l[4], Y_mean0_l[1] + Y_mean0_l[2] * t_seq,
                                Y_mean0_l[1] + Y_mean0_l[2] * Y_mean0_l[4] + Y_mean0_l[3] * (t_seq - Y_mean0_l[4]))
            }
            dat_est_L[[k]] <- data.frame(time = t_seq, Y_hat = Y_hat, Y_hat_u = Y_hat_u, Y_hat_l = Y_hat_l, Class = k)
          }
          dat_est <- do.call(rbind, dat_est_L)
          dat_est$time <- dat_est$time + xstarts
          lines <- c("Model Implied Growth Status" = "solid",
                     "Smooth Line of Observed Growth Status" = "twodash",
                     "95% Confidence Interval on Model Implied Growth Status" = "dotdash")
          fig.status <- ggplot(data = dat_est, aes(x = time, y = Y_hat, group = Class, color = as.factor(Class))) +
            geom_line(aes(linetype = "Model Implied Growth Status"), size = 1) +
            geom_line(aes(x = time, y = Y_hat_l, linetype = "95% Confidence Interval on Model Implied Growth Status"), size = 1) +
            geom_line(aes(x = time, y = Y_hat_u, linetype = "95% Confidence Interval on Model Implied Growth Status"), size = 1) +
            geom_smooth(data = dat_raw, aes(x = time, y = value, group = Class, linetype = "Smooth Line of Observed Growth Status"),
                        linewidth = 1, se = F) +
            labs(x = xlab, y = paste0("Growth Status of ", outcome), color = "Class") +
            scale_linetype_manual(values = c("dotdash", "solid", "twodash")) +
            scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
            scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
            guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3, order = 1),
                   color = guide_legend(title = "Class", order = 2)) +
            theme_bw() +
            theme(strip.text.y = element_text(size = 4),
                  strip.background = element_rect(color = "white", fill = "white"),
                  strip.placement = "outside",
                  legend.position = "bottom",
                  legend.box = "vertical",
                  legend.margin = margin(-3, 0, -3, 0))
          figure <- list(fig.status)
        }
        else if (y_model == "LCSM"){
          dat$Class <- getPosterior(model = model, nClass = nClass, cluster_TIC = cluster_TIC)[[2]]
          dat_chg <- cbind(dat[, c("ID", "Class")],
                           dat[, c(paste0(y_var, records[-1]))] - dat[, c(paste0(y_var, records[1]))])
          if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
            y_var <- "Y"
          }
          names(dat_chg) <- c("ID", "Class", paste0("Delta", names(dat_chg)[-c(1, 2)]))
          dat_chg_l <- pivot_longer(dat_chg, cols = -c(ID, Class), names_to = "Wave", values_to = "value")
          dat_chg_l$Wave <- substring(dat_chg_l$Wave, 7)
          dat_time <- dat[, c("ID", paste0(t_var, records[-1]))]
          dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
          dat_time_l$time <- dat_time_l$time + xstarts
          dat_time_l$Wave <- substring(dat_time_l$Wave, 2)
          dat_raw <- merge(dat_chg_l, dat_time_l, by = c("ID", "Wave"))
          dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
          t_seq <- apply(dat[, paste0(t_var, records[-1])], 2, mean)
          dat_est_L <- list()
          for (k in 1:nClass){
            Y_mean0 <- mxEvalByName(paste0("c", k, y_var, "_mean0"), model = model@submodels[[k]])
            Y_mean0_se <- mxSE(paste0("Class", k, ".c", k, y_var, "_mean0"), model, forceName = T)
            Y_mean0_u <- Y_mean0 + 1.96 * Y_mean0_se
            Y_mean0_l <- Y_mean0 - 1.96 * Y_mean0_se
            if (curveFun %in% c("nonparametric", "NonP")){
              rel_rate <- c(1, model@output$estimate[grep(paste0("c", k, y_var, "_rel_rate"),
                                                          names(model@output$estimate))])
              dY_hat <- Y_mean0[2] * rel_rate
              dY_hat_u <- Y_mean0_u[2] * rel_rate
              dY_hat_l <- Y_mean0_l[2] * rel_rate
              cY_hat <- mxEvalByName(paste0("c", k, y_var, "chg_bl_m"), model = model@submodels[[k]])
              cY_hat_u <- mxEvalByName(paste0("c", k, y_var, "chg_bl_m"), model = model@submodels[[k]]) +
                1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "chg_bl_m"), model, forceName = T)
              cY_hat_l <- mxEvalByName(paste0("c", k, y_var, "chg_bl_m"), model = model@submodels[[k]]) -
                1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "chg_bl_m"), model, forceName = T)
            }
            else if (curveFun %in% c("quadratic", "QUAD")){
              dY_hat <- Y_mean0[2] + 2 * Y_mean0[3] * t_seq
              dY_hat_u <- Y_mean0_u[2] + 2 * Y_mean0_u[3] * t_seq
              dY_hat_l <- Y_mean0_l[2] + 2 * Y_mean0_l[3] * t_seq
              cY_hat <- Y_mean0[2] * t_seq + Y_mean0[3] * t_seq^2
              cY_hat_u <- Y_mean0_u[2] * t_seq + Y_mean0_u[3] * t_seq^2
              cY_hat_l <- Y_mean0_l[2] * t_seq + Y_mean0_l[3] * t_seq^2
            }
            else if (curveFun %in% c("exponential", "EXP")){
              Y_mean0 <- if (length(Y_mean0) == 3){
                Y_mean0
              } else {
                c(Y_mean0, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3])
              }
              Y_mean0_u <- if (length(Y_mean0_u) == 3){
                Y_mean0_u
              } else {
                c(Y_mean0_u, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3] +
                    1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[3])
              }
              Y_mean0_l <- if (length(Y_mean0_l) == 3){
                Y_mean0_l
              } else {
                c(Y_mean0_l, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3] -
                    1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[3])
              }
              dY_hat <- Y_mean0[2] * Y_mean0[3] * exp(-Y_mean0[3] * t_seq)
              dY_hat_u <- Y_mean0_u[2] * Y_mean0_u[3] * exp(-Y_mean0_u[3] * t_seq)
              dY_hat_l <- Y_mean0_l[2] * Y_mean0_l[3] * exp(-Y_mean0_l[3] * t_seq)
              cY_hat <- Y_mean0[2] * (1 - exp(-Y_mean0[3] * t_seq))
              cY_hat_u <- Y_mean0_u[2] * (1 - exp(-Y_mean0_u[3] * t_seq))
              cY_hat_l <- Y_mean0_l[2] * (1 - exp(-Y_mean0_l[3] * t_seq))
            }
            else if (curveFun %in% c("Jenss-Bayley", "JB")){
              Y_mean0 <- if (length(Y_mean0) == 4){
                Y_mean0
              } else {
                c(Y_mean0, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4])
              }
              Y_mean0_u <- if (length(Y_mean0_u) == 4){
                Y_mean0_u
              } else {
                c(Y_mean0_u, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] +
                    1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
              }
              Y_mean0_l <- if (length(Y_mean0_l) == 4){
                Y_mean0_l
              } else {
                c(Y_mean0_l, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] -
                    1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
              }
              dY_hat <- Y_mean0[2] + Y_mean0[3] * Y_mean0[4] * exp(Y_mean0[4] * t_seq)
              dY_hat_u <- Y_mean0_u[2] + Y_mean0_u[3] * Y_mean0_u[4] * exp(Y_mean0_u[4] * t_seq)
              dY_hat_l <- Y_mean0_l[2] + Y_mean0_l[3] * Y_mean0_l[4] * exp(Y_mean0_l[4] * t_seq)
              cY_hat <- Y_mean0[2] * t_seq + Y_mean0[3] * (exp(Y_mean0[4] * t_seq) - 1)
              cY_hat_u <- Y_mean0_u[2] * t_seq + Y_mean0_u[3] * (exp(Y_mean0_u[4] * t_seq) - 1)
              cY_hat_l <- Y_mean0_l[2] * t_seq + Y_mean0_l[3] * (exp(Y_mean0_l[4] * t_seq) - 1)
            }
            dat_est_L[[k]] <- data.frame(time = t_seq, dY_hat = dY_hat, dY_hat_l = dY_hat_l, dY_hat_u = dY_hat_u,
                                         cY_hat = cY_hat, cY_hat_l = cY_hat_l, cY_hat_u = cY_hat_u, Class = k)
          }
          dat_est <- do.call(rbind, dat_est_L)
          dat_est$time <- dat_est$time + xstarts
          lines <- c("Model Implied Change from Baseline" = "solid",
                     "Smooth Line of Observed Change from Baseline" = "twodash",
                     "95% Confidence Interval on Model Implied Change from Baseline" = "dotdash")
          fig.CHG_BL <- ggplot(data = dat_est, aes(x = time, y = cY_hat, group = Class, color = as.factor(Class))) +
            geom_line(aes(linetype = "Model Implied Change from Baseline"),  size = 1) +
            geom_line(aes(x = time, y = cY_hat_l, linetype = "95% Confidence Interval on Model Implied Change from Baseline"), size = 1) +
            geom_line(aes(x = time, y = cY_hat_u, linetype = "95% Confidence Interval on Model Implied Change from Baseline"), size = 1) +
            geom_smooth(data = dat_raw, aes(x = time, y = value, group = Class,
                                            linetype = "Smooth Line of Observed Change from Baseline"),
                        size = 1, se = F) +
            labs(x = xlab, y = paste0("Change from baseline of ", outcome), color = "Class") +
            scale_linetype_manual(values = c("dotdash", "solid", "twodash")) +
            scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
            scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
            guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3, order = 1),
                   color = guide_legend(title = "Class", order = 2)) +
            theme_bw() +
            theme(strip.text.y = element_text(size = 4),
                  strip.background = element_rect(color = "white", fill = "white"),
                  strip.placement = "outside",
                  legend.position = "bottom",
                  legend.box = "vertical",
                  legend.margin = margin(-3, 0, -3, 0))
          lines <- c("Model Implied Change from Baseline" = "solid",
                     "95% Confidence Interval on Model Implied Change from Baseline" = "dotdash")
          if (curveFun %in% c("nonparametric", "NonP")){
            fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = Class, color = as.factor(Class))) +
              geom_step(aes(linetype = "Model Implied Growth Rate"),  size = 1) +
              geom_step(aes(x = time, y = dY_hat_l, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
              geom_step(aes(x = time, y = dY_hat_u, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
              labs(x = xlab, y = paste0("Growth Rate of ", outcome), color = "Class") +
              scale_linetype_manual(values = c("dotdash", "solid")) +
              scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
              scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
              guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2, order = 1),
                     color = guide_legend(title = "Class", order = 2)) +
              theme_bw() +
              theme(strip.text.y = element_text(size = 4),
                    strip.background = element_rect(color = "white", fill = "white"),
                    strip.placement = "outside",
                    legend.position = "bottom",
                    legend.box = "vertical",
                    legend.margin = margin(-3, 0, -3, 0))
          }
          else{
            fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = Class, color = as.factor(Class))) +
              geom_line(aes(linetype = "Model Implied Growth Rate"),  size = 1) +
              geom_line(aes(x = time, y = dY_hat_l, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
              geom_line(aes(x = time, y = dY_hat_u, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
              labs(x = xlab, y = paste0("Growth Rate of ", outcome), color = "Class") +
              scale_linetype_manual(values = c("dotdash", "solid")) +
              scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
              scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
              guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2, order = 1),
                     color = guide_legend(title = "Class", order = 2)) +
              theme_bw() +
              theme(strip.text.y = element_text(size = 4),
                    strip.background = element_rect(color = "white", fill = "white"),
                    strip.placement = "outside",
                    legend.position = "bottom",
                    legend.box = "vertical",
                    legend.margin = margin(-3, 0, -3, 0))
          }
          figure <- list(fig.CHG_BL, fig.SLP)
        }
    }
    else if (!is.null(grp_var)){
      if (y_model == "LGCM"){
        dat_traj <- dat[, c("ID", grp_var, paste0(y_var, records))]
        if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
          y_var <- "Y"
        }
        names(dat_traj)[2] <- "Class"
        dat_time <- dat[, c("ID", grp_var, paste0(t_var, records))]
        names(dat_time)[2] <- "Class"
        dat_traj_l <- pivot_longer(dat_traj, cols = -c(ID, Class), names_to = "Wave", values_to = "value")
        dat_traj_l$Wave <- substring(dat_traj_l$Wave, 2)
        dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
        dat_time_l$time <- dat_time_l$time + xstarts
        dat_time_l$Wave <- substring(dat_time_l$Wave, 2)
        dat_raw <- merge(dat_traj_l, dat_time_l, by = c("ID", "Wave"))
        dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
        t_seq <- apply(dat[, paste0(t_var, records)], 2, mean)
        dat_est_L <- list()
        for (k in 1:nClass){
          Y_mean0 <- mxEvalByName(paste0("c", k, y_var, "_mean0"), model = model@submodels[[k]])
          Y_mean0_se <- mxSE(paste0("Class", k, ".c", k, y_var, "_mean0"), model, forceName = T)
          Y_mean0_u <- Y_mean0 + 1.96 * Y_mean0_se
          Y_mean0_l <- Y_mean0 - 1.96 * Y_mean0_se
          if (curveFun %in% c("linear", "LIN")){
            Y_hat <- Y_mean0[1] + Y_mean0[2] * t_seq
            Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * t_seq
            Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * t_seq
          }
          else if (curveFun %in% c("quadratic", "QUAD")){
            Y_hat <- Y_mean0[1] + Y_mean0[2] * t_seq + Y_mean0[3] * t_seq^2
            Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * t_seq + Y_mean0_u[3] * t_seq^2
            Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * t_seq + Y_mean0_l[3] * t_seq^2
          }
          else if (curveFun %in% c("exponential", "EXP")){
            Y_mean0 <- if (length(Y_mean0) == 3){
              Y_mean0
            } else {
              c(Y_mean0, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3])
            }
            Y_mean0_u <- if (length(Y_mean0_u) == 3){
              Y_mean0_u
            } else {
              c(Y_mean0_u, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3] +
                  1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[3])
            }
            Y_mean0_l <- if (length(Y_mean0_l) == 3){
              Y_mean0_l
            } else {
              c(Y_mean0_l, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3] -
                  1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[3])
            }
            Y_hat <- Y_mean0[1] + Y_mean0[2] * (1 - exp(-Y_mean0[3] * t_seq))
            Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * (1 - exp(-Y_mean0_u[3] * t_seq))
            Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * (1 - exp(-Y_mean0_l[3] * t_seq))
          }
          else if (curveFun %in% c("Jenss-Bayley", "JB")){
            Y_mean0 <- if (length(Y_mean0) == 4){
              Y_mean0
            } else {
              c(Y_mean0, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4])
            }
            Y_mean0_u <- if (length(Y_mean0_u) == 4){
              Y_mean0_u
            } else {
              c(Y_mean0_u, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] +
                  1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
            }
            Y_mean0_l <- if (length(Y_mean0_l) == 4){
              Y_mean0_l
            } else {
              c(Y_mean0_l, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] -
                  1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
            }
            Y_hat <- Y_mean0[1] + Y_mean0[2] * t_seq + Y_mean0[3] * (exp(Y_mean0[4] * t_seq) - 1)
            Y_hat_u <- Y_mean0_u[1] + Y_mean0_u[2] * t_seq + Y_mean0_u[3] * (exp(Y_mean0_u[4] * t_seq) - 1)
            Y_hat_l <- Y_mean0_l[1] + Y_mean0_l[2] * t_seq + Y_mean0_l[3] * (exp(Y_mean0_l[4] * t_seq) - 1)
          }
          else if (curveFun %in% c("bilinear spline", "BLS")){
            Y_mean0 <- if (length(Y_mean0) == 4){
              Y_mean0
            } else {
              c(Y_mean0, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4])
            }
            Y_mean0_u <- if (length(Y_mean0_u) == 4){
              Y_mean0_u
            } else {
              c(Y_mean0_u, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] +
                  1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
            }
            Y_mean0_l <- if (length(Y_mean0_l) == 4){
              Y_mean0_l
            } else {
              c(Y_mean0_l, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] -
                  1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
            }
            Y_hat <- ifelse(t_seq <= Y_mean0[4], Y_mean0[1] + Y_mean0[2] * t_seq,
                            Y_mean0[1] + Y_mean0[2] * Y_mean0[4] + Y_mean0[3] * (t_seq - Y_mean0[4]))
            Y_hat_u <- ifelse(t_seq <= Y_mean0_u[4], Y_mean0_u[1] + Y_mean0_u[2] * t_seq,
                              Y_mean0_u[1] + Y_mean0_u[2] * Y_mean0_u[4] + Y_mean0_u[3] * (t_seq - Y_mean0_u[4]))
            Y_hat_l <- ifelse(t_seq <= Y_mean0_l[4], Y_mean0_l[1] + Y_mean0_l[2] * t_seq,
                              Y_mean0_l[1] + Y_mean0_l[2] * Y_mean0_l[4] + Y_mean0_l[3] * (t_seq - Y_mean0_l[4]))
          }
          dat_est_L[[k]] <- data.frame(time = t_seq, Y_hat = Y_hat, Y_hat_u = Y_hat_u, Y_hat_l = Y_hat_l, Class = k)
        }
        dat_est <- do.call(rbind, dat_est_L)
        dat_est$time <- dat_est$time + xstarts
        lines <- c("Model Implied Growth Status" = "solid",
                   "Smooth Line of Observed Growth Status" = "twodash",
                   "95% Confidence Interval on Model Implied Growth Status" = "dotdash")
        fig.status <- ggplot(data = dat_est, aes(x = time, y = Y_hat, group = Class, color = as.factor(Class))) +
          geom_line(aes(linetype = "Model Implied Growth Status"), size = 1) +
          geom_line(aes(x = time, y = Y_hat_l, linetype = "95% Confidence Interval on Model Implied Growth Status"), size = 1) +
          geom_line(aes(x = time, y = Y_hat_u, linetype = "95% Confidence Interval on Model Implied Growth Status"), size = 1) +
          geom_smooth(data = dat_raw, aes(x = time, y = value, group = Class, linetype = "Smooth Line of Observed Growth Status"),
                      linewidth = 1, se = F) +
          labs(x = xlab, y = paste0("Growth Status of ", outcome), color = "Class") +
          scale_linetype_manual(values = c("dotdash", "solid", "twodash")) +
          scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
          scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
          guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3, order = 1),
                 color = guide_legend(title = grp_var, order = 2)) +
          theme_bw() +
          theme(strip.text.y = element_text(size = 4),
                strip.background = element_rect(color = "white", fill = "white"),
                strip.placement = "outside",
                legend.position = "bottom",
                legend.box = "vertical",
                legend.margin = margin(-3, 0, -3, 0))
        figure <- list(fig.status)
      }
      else if (y_model == "LCSM"){
        dat_chg <- cbind(dat[, c("ID", grp_var)],
                         dat[, c(paste0(y_var, records[-1]))] - dat[, c(paste0(y_var, records[1]))])
        if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
          y_var <- "Y"
        }
        names(dat_chg)[2] <- "Class"
        names(dat_chg) <- c("ID", "Class", paste0("Delta", names(dat_chg)[-c(1, 2)]))
        dat_chg_l <- pivot_longer(dat_chg, cols = -c(ID, Class), names_to = "Wave", values_to = "value")
        dat_chg_l$Wave <- substring(dat_chg_l$Wave, 7)
        dat_time <- dat[, c("ID", paste0(t_var, records[-1]))]
        dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
        dat_time_l$time <- dat_time_l$time + xstarts
        dat_time_l$Wave <- substring(dat_time_l$Wave, 2)
        dat_raw <- merge(dat_chg_l, dat_time_l, by = c("ID", "Wave"))
        dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
        t_seq <- apply(dat[, paste0(t_var, records[-1])], 2, mean)
        dat_est_L <- list()
        for (k in 1:nClass){
          Y_mean0 <- mxEvalByName(paste0("c", k, y_var, "_mean0"), model = model@submodels[[k]])
          Y_mean0_se <- mxSE(paste0("Class", k, ".c", k, y_var, "_mean0"), model, forceName = T)
          Y_mean0_u <- Y_mean0 + 1.96 * Y_mean0_se
          Y_mean0_l <- Y_mean0 - 1.96 * Y_mean0_se
          if (curveFun %in% c("nonparametric", "NonP")){
            rel_rate <- c(1, model@output$estimate[grep(paste0("c", k, y_var, "_rel_rate"),
                                                        names(model@output$estimate))])
            dY_hat <- Y_mean0[2] * rel_rate
            dY_hat_u <- Y_mean0_u[2] * rel_rate
            dY_hat_l <- Y_mean0_l[2] * rel_rate
            cY_hat <- mxEvalByName(paste0("c", k, y_var, "chg_bl_m"), model = model@submodels[[k]])
            cY_hat_u <- mxEvalByName(paste0("c", k, y_var, "chg_bl_m"), model = model@submodels[[k]]) +
              1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "chg_bl_m"), model, forceName = T)
            cY_hat_l <- mxEvalByName(paste0("c", k, y_var, "chg_bl_m"), model = model@submodels[[k]]) -
              1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "chg_bl_m"), model, forceName = T)
          }
          else if (curveFun %in% c("quadratic", "QUAD")){
            dY_hat <- Y_mean0[2] + 2 * Y_mean0[3] * t_seq
            dY_hat_u <- Y_mean0_u[2] + 2 * Y_mean0_u[3] * t_seq
            dY_hat_l <- Y_mean0_l[2] + 2 * Y_mean0_l[3] * t_seq
            cY_hat <- Y_mean0[2] * t_seq + Y_mean0[3] * t_seq^2
            cY_hat_u <- Y_mean0_u[2] * t_seq + Y_mean0_u[3] * t_seq^2
            cY_hat_l <- Y_mean0_l[2] * t_seq + Y_mean0_l[3] * t_seq^2
          }
          else if (curveFun %in% c("exponential", "EXP")){
            Y_mean0 <- if (length(Y_mean0) == 3){
              Y_mean0
            } else {
              c(Y_mean0, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3])
            }
            Y_mean0_u <- if (length(Y_mean0_u) == 3){
              Y_mean0_u
            } else {
              c(Y_mean0_u, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3] +
                  1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[3])
            }
            Y_mean0_l <- if (length(Y_mean0_l) == 3){
              Y_mean0_l
            } else {
              c(Y_mean0_l, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[3] -
                  1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[3])
            }
            dY_hat <- Y_mean0[2] * Y_mean0[3] * exp(-Y_mean0[3] * t_seq)
            dY_hat_u <- Y_mean0_u[2] * Y_mean0_u[3] * exp(-Y_mean0_u[3] * t_seq)
            dY_hat_l <- Y_mean0_l[2] * Y_mean0_l[3] * exp(-Y_mean0_l[3] * t_seq)
            cY_hat <- Y_mean0[2] * (1 - exp(-Y_mean0[3] * t_seq))
            cY_hat_u <- Y_mean0_u[2] * (1 - exp(-Y_mean0_u[3] * t_seq))
            cY_hat_l <- Y_mean0_l[2] * (1 - exp(-Y_mean0_l[3] * t_seq))
          }
          else if (curveFun %in% c("Jenss-Bayley", "JB")){
            Y_mean0 <- if (length(Y_mean0) == 4){
              Y_mean0
            } else {
              c(Y_mean0, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4])
            }
            Y_mean0_u <- if (length(Y_mean0_u) == 4){
              Y_mean0_u
            } else {
              c(Y_mean0_u, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] +
                  1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
            }
            Y_mean0_l <- if (length(Y_mean0_l) == 4){
              Y_mean0_l
            } else {
              c(Y_mean0_l, mxEvalByName(paste0("c", k, y_var, "_alpha0"), model = model@submodels[[k]])[4] -
                  1.96 * mxSE(paste0("Class", k, ".c", k, y_var, "_alpha0"), model, forceName = T)[4])
            }
            dY_hat <- Y_mean0[2] + Y_mean0[3] * Y_mean0[4] * exp(Y_mean0[4] * t_seq)
            dY_hat_u <- Y_mean0_u[2] + Y_mean0_u[3] * Y_mean0_u[4] * exp(Y_mean0_u[4] * t_seq)
            dY_hat_l <- Y_mean0_l[2] + Y_mean0_l[3] * Y_mean0_l[4] * exp(Y_mean0_l[4] * t_seq)
            cY_hat <- Y_mean0[2] * t_seq + Y_mean0[3] * (exp(Y_mean0[4] * t_seq) - 1)
            cY_hat_u <- Y_mean0_u[2] * t_seq + Y_mean0_u[3] * (exp(Y_mean0_u[4] * t_seq) - 1)
            cY_hat_l <- Y_mean0_l[2] * t_seq + Y_mean0_l[3] * (exp(Y_mean0_l[4] * t_seq) - 1)
          }
          dat_est_L[[k]] <- data.frame(time = t_seq, dY_hat = dY_hat, dY_hat_l = dY_hat_l, dY_hat_u = dY_hat_u,
                                       cY_hat = cY_hat, cY_hat_l = cY_hat_l, cY_hat_u = cY_hat_u, Class = k)
        }
        dat_est <- do.call(rbind, dat_est_L)
        dat_est$time <- dat_est$time + xstarts
        lines <- c("Model Implied Change from Baseline" = "solid",
                   "Smooth Line of Observed Change from Baseline" = "twodash",
                   "95% Confidence Interval on Model Implied Change from Baseline" = "dotdash")
        fig.CHG_BL <- ggplot(data = dat_est, aes(x = time, y = cY_hat, group = Class, color = as.factor(Class))) +
          geom_line(aes(linetype = "Model Implied Change from Baseline"),  size = 1) +
          geom_line(aes(x = time, y = cY_hat_l, linetype = "95% Confidence Interval on Model Implied Change from Baseline"), size = 1) +
          geom_line(aes(x = time, y = cY_hat_u, linetype = "95% Confidence Interval on Model Implied Change from Baseline"), size = 1) +
          geom_smooth(data = dat_raw, aes(x = time, y = value, group = Class,
                                          linetype = "Smooth Line of Observed Change from Baseline"),
                      size = 1, se = F) +
          labs(x = xlab, y = paste0("Change from baseline of ", outcome), color = "Class") +
          scale_linetype_manual(values = c("dotdash", "solid", "twodash")) +
          scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
          scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
          guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3, order = 1),
                 color = guide_legend(title = grp_var, order = 2)) +
          theme_bw() +
          theme(strip.text.y = element_text(size = 4),
                strip.background = element_rect(color = "white", fill = "white"),
                strip.placement = "outside",
                legend.position = "bottom",
                legend.box = "vertical",
                legend.margin = margin(-3, 0, -3, 0))
        lines <- c("Model Implied Change from Baseline" = "solid",
                   "95% Confidence Interval on Model Implied Change from Baseline" = "dotdash")
        if (curveFun %in% c("nonparametric", "NonP")){
          fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = Class, color = as.factor(Class))) +
            geom_step(aes(linetype = "Model Implied Growth Rate"),  size = 1) +
            geom_step(aes(x = time, y = dY_hat_l, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
            geom_step(aes(x = time, y = dY_hat_u, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
            labs(x = xlab, y = paste0("Growth Rate of ", outcome), color = "Class") +
            scale_linetype_manual(values = c("dotdash", "solid")) +
            scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
            scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
            guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2, order = 1),
                   color = guide_legend(title = grp_var, order = 2)) +
            theme_bw() +
            theme(strip.text.y = element_text(size = 4),
                  strip.background = element_rect(color = "white", fill = "white"),
                  strip.placement = "outside",
                  legend.position = "bottom",
                  legend.box = "vertical",
                  legend.margin = margin(-3, 0, -3, 0))
        }
        else{
          fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = Class, color = as.factor(Class))) +
            geom_line(aes(linetype = "Model Implied Growth Rate"),  size = 1) +
            geom_line(aes(x = time, y = dY_hat_l, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
            geom_line(aes(x = time, y = dY_hat_u, linetype = "95% Confidence Interval on Model Implied Growth Rate"), size = 1) +
            labs(x = xlab, y = paste0("Growth Rate of ", outcome), color = "Class") +
            scale_linetype_manual(values = c("dotdash", "solid")) +
            scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
            scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
            guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2, order = 1),
                   color = guide_legend(title = grp_var, order = 2)) +
            theme_bw() +
            theme(strip.text.y = element_text(size = 4),
                  strip.background = element_rect(color = "white", fill = "white"),
                  strip.placement = "outside",
                  legend.position = "bottom",
                  legend.box = "vertical",
                  legend.margin = margin(-3, 0, -3, 0))
        }
        figure <- list(fig.CHG_BL, fig.SLP)
      }

    }
  }
  return(figure)
}

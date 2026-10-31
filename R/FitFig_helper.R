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
#' @param med_role For \code{sub_Model = "MED"}, the role of the plotted process in the mediation model: \code{"Y"}
#' (outcome), \code{"M"} (mediator), or \code{"X"} (longitudinal predictor). Mediation builders name the growth-factor
#' algebras by role rather than by the user-supplied variable prefix. Default is \code{NULL}.
#'
#' @return A ggplot object or a list of ggplot objects.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_smooth stat_summary labs scale_linetype_manual scale_x_continuous guides guide_legend theme_bw theme element_text element_rect geom_step scale_color_manual margin
#' @importFrom grDevices colorRampPalette
#' @importFrom OpenMx mxGetExpected
#'
## Distribution of the growth-factor means in the parameterisation used by the plot.
## Asking mxSE for their joint covariance retains correlations between means, rates,
## knots and covariate-adjusted means, including mediation reparameterisations.
.plot_gf_distribution <- function(model, eval_model, eval_prefix, se_prefix, sub_Model, curveFun){
  suffix <- if (sub_Model == "MED") "_mean" else "_mean0"
  est <- c(mxEvalByName(paste0(eval_prefix, suffix), model = eval_model))
  base <- paste0(se_prefix, suffix)
  terms <- paste0(base, "[", seq_along(est), ",1]")
  if (sub_Model == "MED" && curveFun %in% c("bilinear spline", "BLS")){
    if (length(est) == 3L){
      est <- c(est, c(mxEvalByName(paste0(eval_prefix, "_mug"), model = eval_model)))
      terms <- c(terms, paste0(se_prefix, "_mug"))
    }
    terms <- c(paste0("(", terms[2], " - ", terms[1], " * ", terms[4], ")"), terms[c(1, 3, 4)])
    est <- c(est[2] - est[1] * est[4], est[c(1, 3, 4)])
  } else {
    needed <- if (curveFun %in% c("negative exponential", "exponential", "EXP")) 3L else
      if (curveFun %in% c("Jenss-Bayley", "JB", "bilinear spline", "BLS")) 4L else length(est)
    if (length(est) < needed){
      est <- c(est, c(mxEvalByName(paste0(eval_prefix, "_alpha0"), model = eval_model))[needed])
      terms <- c(terms, paste0(se_prefix, "_alpha0[", needed, ",1]"))
    }
  }
  expression <- paste0("rbind(", paste(terms, collapse = ","), ")")
  covariance <- mxSE(expression, model, details = TRUE, silent = TRUE)$Cov
  list(est = est, cov = covariance)
}

## Pointwise delta-method intervals. No display-range truncation is applied to bounds.
.plot_delta_ci <- function(estimate, gradient, covariance){
  variance <- rowSums((gradient %*% covariance) * gradient)
  # Remove only round-off at zero; materially negative variances remain undefined.
  tol <- 100 * .Machine$double.eps * rowSums(abs(gradient %*% covariance) * abs(gradient))
  tiny_negative <- !is.na(variance) & variance < 0 & variance >= -tol
  variance[tiny_negative] <- 0
  se <- suppressWarnings(sqrt(variance))
  list(est = as.numeric(estimate), lower = as.numeric(estimate - stats::qnorm(0.975) * se),
       upper = as.numeric(estimate + stats::qnorm(0.975) * se))
}

.plot_curve_delta <- function(distribution, times, curveFun, quantity = c("status", "rate", "change"), baseline_time = 0){
  quantity <- match.arg(quantity)
  p <- distribution$est
  plot_times <- as.numeric(times)
  tt <- if (quantity == "change") c(plot_times, baseline_time) else plot_times
  z <- rep(0, length(tt)); one <- rep(1, length(tt))
  if (curveFun %in% c("linear", "LIN")){
    estimate <- p[1] + p[2] * tt
    gradient <- cbind(one, tt)
    if (quantity == "rate") { estimate <- rep(p[2], length(tt)); gradient <- cbind(z, one) }
  } else if (curveFun %in% c("quadratic", "QUAD")){
    estimate <- p[1] + p[2] * tt + p[3] * tt^2
    gradient <- cbind(one, tt, tt^2)
    if (quantity == "rate") { estimate <- p[2] + 2 * p[3] * tt; gradient <- cbind(z, one, 2 * tt) }
  } else if (curveFun %in% c("negative exponential", "exponential", "EXP")){
    ee <- exp(-p[3] * tt)
    estimate <- p[1] + p[2] * (1 - ee)
    gradient <- cbind(one, 1 - ee, p[2] * tt * ee)
    if (quantity == "rate"){
      estimate <- p[2] * p[3] * ee
      gradient <- cbind(z, p[3] * ee, p[2] * ee * (1 - p[3] * tt))
    }
  } else if (curveFun %in% c("Jenss-Bayley", "JB")){
    ee <- exp(p[4] * tt)
    estimate <- p[1] + p[2] * tt + p[3] * (ee - 1)
    gradient <- cbind(one, tt, ee - 1, p[3] * tt * ee)
    if (quantity == "rate"){
      estimate <- p[2] + p[3] * p[4] * ee
      gradient <- cbind(z, one, p[4] * ee, p[3] * ee * (1 + p[4] * tt))
    }
  } else if (curveFun %in% c("bilinear spline", "BLS")){
    after <- tt > p[4]
    estimate <- p[1] + p[2] * pmin(tt, p[4]) + p[3] * pmax(0, tt - p[4])
    gradient <- cbind(one, pmin(tt, p[4]), pmax(0, tt - p[4]), (p[2] - p[3]) * after)
    if (quantity == "rate") stop("Bilinear-spline rate plots are not supported.")
  } else stop("Unsupported curve function for trajectory intervals: ", curveFun)
  if (quantity == "change"){
    baseline <- length(estimate)
    estimate <- estimate[-baseline] - estimate[baseline]
    gradient <- sweep(gradient[-baseline, , drop = FALSE], 2, gradient[baseline, ], "-")
  }
  out <- .plot_delta_ci(estimate, gradient, distribution$cov)
  if (curveFun %in% c("bilinear spline", "BLS") &&
      is.finite(distribution$cov[4, 4]) && distribution$cov[4, 4] > 0 && p[2] != p[3]){
    knot_tolerance <- sqrt(.Machine$double.eps) * max(1, abs(p[4]))
    at_knot <- abs(plot_times - p[4]) <= knot_tolerance
    if (quantity == "change" && abs(baseline_time - p[4]) <= knot_tolerance) at_knot[] <- TRUE
    # The derivative in the estimated knot is not unique exactly at the kink.
    out$lower[at_knot] <- out$upper[at_knot] <- NA_real_
  }
  out
}

## Model-implied mean of the outcome at each wave of a TVC model, including the effect of the TVC, at the
## mean occasions (and the means of the other definition variables) of observed, by default the model's
## own data: the occasions at which the growth-component curve is drawn.
.tvc_full_mean <- function(model, columns, observed = model@data@observed){
  num <- vapply(observed, is.numeric, logical(1))
  observed[1, num] <- lapply(observed[num], mean, na.rm = TRUE)
  model <- mxModel(model, mxData(observed, type = "raw"))
  expected <- mxGetExpected(model, "means", defvar.row = 1)
  as.numeric(expected[1, columns])
}

## The same for each class of a mixture, at the mean occasions of all persons, or for each group of a
## multiple-group TVC model, at the group's own data and mean occasions; as changes from the first wave when
## change = TRUE.
.tvc_class_full_means <- function(model, columns, time_columns, xstarts, change = FALSE){
  groups <- !length(model@data)
  do.call(rbind, lapply(seq_along(model@submodels), function(k){
    observed <- if (groups) model@submodels[[k]]@data$observed else model@data$observed
    times <- colMeans(observed[, time_columns, drop = FALSE], na.rm = TRUE) + xstarts
    value <- .tvc_full_mean(model@submodels[[k]], columns, observed)
    if (change) data.frame(time = times[-1], value = value[-1] - value[1], Class = k) else
      data.frame(time = times, value = value, Class = k)
  }))
}

## Observed curve: ggplot2's default smoother, unless that would be a GAM (1,000 or more rows in a group) on
## fewer than ten distinct occasions, which it cannot fit; the mean at each occasion is drawn then.
.observed_curve <- function(data, mapping, ...){
  # Match the rows available to the smoother: an occasion whose outcomes are all missing cannot
  # support a basis knot, even though its times are present in the fitted data.
  data <- data[is.finite(data$time) & is.finite(data$value), , drop = FALSE]
  group <- if ("Class" %in% names(data)) data$Class else rep(1L, nrow(data))
  rows <- if (nrow(data)) max(table(group)) else 0L
  occasions <- if (nrow(data)) min(vapply(split(data$time, group, drop = TRUE),
                                          function(x) length(unique(x)), integer(1))) else 0L
  if (occasions == 1L){
    stat_summary(data = data, mapping = mapping, fun = mean, fun.min = mean, fun.max = mean,
                 geom = "pointrange", shape = 4, na.rm = TRUE, ...)
  } else if (rows >= 1000 && occasions < 10){
    stat_summary(data = data, mapping = mapping, fun = mean, geom = "line", na.rm = TRUE, ...)
  } else{
    geom_smooth(data = data, mapping = mapping, se = FALSE, ...)
  }
}

.plot_algebra_ci <- function(model, eval_model, eval_name, se_name){
  estimate <- c(mxEvalByName(eval_name, model = eval_model))
  se <- c(mxSE(se_name, model, forceName = TRUE, silent = TRUE))
  list(est = estimate, lower = estimate - stats::qnorm(0.975) * se,
       upper = estimate + stats::qnorm(0.975) * se)
}

## Cumulative fitted recurrence at the model's original mean occasions. For random
## rates the last latent is a deviation from the expansion rate, including any
## covariate shift, rather than the rate itself. Keep this shift in the mean and
## in its joint uncertainty; a continuous growth curve is a different estimand.
.lcsm_plot_components <- function(eval_model, eval_prefix, se_prefix, curveFun, sub_Model){
  mean0 <- c(mxEvalByName(paste0(eval_prefix, "_mean0"), eval_model, compute = TRUE))
  loading_prefix <- if (sub_Model == "MGM") eval_prefix else substr(eval_prefix, 1, nchar(eval_prefix) - 1L)
  loading_se_prefix <- if (sub_Model == "MGM") se_prefix else substr(se_prefix, 1, nchar(se_prefix) - 1L)
  loads <- mxEvalByName(paste0(loading_prefix, "s_loads"), eval_model, compute = TRUE)
  term <- function(j) paste0(se_prefix, "_mean0[", j, ",1]")
  if (curveFun %in% c("nonparametric", "NonP")){
    means <- mean0[2]; terms <- term(2)
  } else if (curveFun %in% c("quadratic", "QUAD")){
    means <- mean0[2:3]; terms <- term(2:3)
  } else if (curveFun %in% c("negative exponential", "EXP", "Jenss-Bayley", "JB")){
    last <- if (curveFun %in% c("negative exponential", "EXP")) 3L else 4L
    means <- mean0[2:(last - 1L)]; terms <- term(2:(last - 1L))
    if (ncol(loads) == last - 1L){
      center <- c(mxEvalByName(paste0(eval_prefix, "_mug"), eval_model, compute = TRUE))
      means <- c(means, mean0[last] - center)
      terms <- c(terms, paste0("(", term(last), " - ", se_prefix, "_mug)"))
    }
  } else stop("Unsupported curve function for fitted change plots: ", curveFun)
  if (ncol(loads) != length(means)) stop("The fitted change loadings and growth means do not agree.")
  expression <- paste0(loading_se_prefix, "s_loads %*% rbind(", paste(terms, collapse = ","), ")")
  list(est = c(loads %*% means), expression = expression)
}

.lcsm_plot_distribution <- function(model, components){
  expression <- components$expression
  covariance <- mxSE(expression, model, details = TRUE, silent = TRUE)$Cov
  list(est = components$est, cov = covariance)
}

## Use fitted manifest order to map wave IDs, not sorted column names. Legacy
## models whose observed columns were explicitly aliased can be matched to the
## stored mean grid only when every requested occasion has one unique match.
.lcsm_plot_indices <- function(eval_model, sub_Model, raw_y_var, records, times, baseline_time, n){
  grid_name <- if (sub_Model == "MGM") paste0(raw_y_var, "M_TIME") else "M_TIME"
  grid <- c(mxEvalByName(grid_name, eval_model, compute = TRUE))
  selected_times <- c(baseline_time, times)
  if (length(grid) != n || any(!is.finite(grid)) || any(!is.finite(selected_times))){
    stop("Fitted change plots require a finite observation-time grid.", call. = FALSE)
  }
  indices <- NULL
  if (!is.null(records)){
    if (length(records) != length(selected_times) || anyNA(records) || anyDuplicated(records)){
      stop("Select at least two distinct fitted waves for a change plot.", call. = FALSE)
    }
    manifests <- eval_model@manifestVars
    suffix <- substring(manifests[startsWith(manifests, raw_y_var)], nchar(raw_y_var) + 1L)
    fitted_records <- as.numeric(suffix[grepl("^[0-9]+$", suffix)])
    if (length(fitted_records) >= n){
      indices <- match(records, fitted_records[seq_len(n)])
      if (anyNA(indices)) stop("Change plots can select only waves included in the fitted model.", call. = FALSE)
    }
  }
  if (is.null(indices)){
    indices <- vapply(selected_times, function(tt){
      candidates <- which(abs(grid - tt) <= 100 * .Machine$double.eps * pmax(1, abs(grid), abs(tt)))
      if (length(candidates) != 1L) stop("Requested occasions do not uniquely identify fitted waves.", call. = FALSE)
      candidates
    }, integer(1))
  }
  if (length(indices) < 2L || anyDuplicated(indices) ||
      !(all(diff(indices) > 0) || all(diff(indices) < 0))){
    stop("Select distinct waves in fitted order or its reverse for a change plot.", call. = FALSE)
  }
  tolerance <- 100 * .Machine$double.eps * pmax(1, abs(grid[indices]), abs(selected_times))
  if (any(abs(grid[indices] - selected_times) > tolerance) || any(diff(selected_times) == 0)){
    stop("Selected observation times do not match the fitted change-score grid.", call. = FALSE)
  }
  indices
}

.trajectory_plot_data <- function(model, eval_model, eval_prefix, se_prefix, sub_Model, curveFun, y_model, times,
                                  baseline_time = 0, records = NULL, raw_y_var = NULL){
  if (y_model == "LGCM"){
    distribution <- .plot_gf_distribution(model, eval_model, eval_prefix, se_prefix, sub_Model, curveFun)
    status <- .plot_curve_delta(distribution, times, curveFun, "status")
    return(data.frame(time = times, Y_hat = status$est, Y_hat_u = status$upper, Y_hat_l = status$lower))
  }
  # Mixture submodels inherit data from the parent. Supply those data only to
  # this evaluation copy; uncertainty still uses the parent fitted covariance.
  if (!length(eval_model@data) && length(model@data)){
    eval_model <- mxModel(eval_model, mxData(model@data$observed, type = "raw"))
  }
  components <- .lcsm_plot_components(eval_model, eval_prefix, se_prefix, curveFun, sub_Model)
  n <- length(components$est)
  indices <- .lcsm_plot_indices(eval_model, sub_Model, raw_y_var, records, times, baseline_time, n)
  distribution <- .lcsm_plot_distribution(model, components)
  selected <- diag(n)[indices, , drop = FALSE]
  change_contrast <- sweep(selected[-1, , drop = FALSE], 2, selected[1, ], "-")
  rate_contrast <- (selected[-1, , drop = FALSE] - selected[-nrow(selected), , drop = FALSE]) /
    diff(c(baseline_time, times))
  change <- .plot_delta_ci(change_contrast %*% distribution$est, change_contrast, distribution$cov)
  rate <- .plot_delta_ci(rate_contrast %*% distribution$est, rate_contrast, distribution$cov)
  data.frame(time = times, dY_hat = rate$est, dY_hat_l = rate$lower, dY_hat_u = rate$upper,
             cY_hat = change$est, cY_hat_l = change$lower, cY_hat_u = change$upper,
             interval_start = c(baseline_time, utils::head(times, -1L)))
}

getFitFig <- function(model, nClass, cluster_TIC, grp_var, sub_Model, t_var, records, y_var, curveFun, y_model,
                      xstarts, xlab, outcome, med_role = NULL){
  if (y_model == "LCSM" && (!is.numeric(records) || length(records) < 2L || anyNA(records) ||
                           any(!is.finite(records)) || any(records < 1 | records != trunc(records)) || anyDuplicated(records))){
    stop("Select at least two distinct fitted waves for a change plot.", call. = FALSE)
  }
  # In a TVC model the curves are the latent growth component: they omit the effect of the TVC on the
  # outcome, so they are labelled as such; the single-group plots add the full model-implied mean.
  tvc <- identical(sub_Model, "TVC")
  lt_labels <- function(x){
    if (tvc){
      x <- sub("Model Implied Growth Status", "Model Implied Growth Component", x, fixed = TRUE)
      x <- ifelse(grepl("Model Implied", x, fixed = TRUE) & !grepl("TVC Effect Included", x, fixed = TRUE),
                   paste0(x, " (TVC Effect Excluded)"), x)
    }
    if (y_model == "LCSM") x <- vapply(x, function(label) paste(strwrap(label, width = 65), collapse = "\n"), character(1))
    x
  }
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
  raw_y_var <- y_var
  # Preserve the manifested-group codes before reserving ID for plotting row identity.
  # The public grouping column itself may be named ID.
  group_codes <- if (!is.null(grp_var)) dat[[grp_var]] else NULL
  dat$ID <- 1:nrow(dat)
  if (!is.null(nClass)){
    if (is.null(grp_var)){
      dat$Class <- getPosterior(model = model, nClass = nClass, cluster_TIC = cluster_TIC)@membership
      class_props <- table(factor(dat$Class, levels = seq_len(nClass)))/length(dat$Class)
      class_labels <- paste0("Class ", seq_len(nClass), " (", round(class_props * 100, 1), "%)")
    }
    else if (!is.null(grp_var)){
      class_props <- table(group_codes)/length(group_codes)
      class_labels <- paste0("Group ", seq_len(nClass), " (", round(class_props * 100, 1), "%)")
    }
  }
  if (is.null(nClass)){
    if (y_model == "LGCM"){
      dat_traj <- dat[, c("ID", paste0(y_var, records))]
      if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
        y_var <- "Y"
      }
      else if (sub_Model == "MED"){
        y_var <- med_role
      }
      dat_time <- dat[, c("ID", paste0(t_var, records))]
      dat_traj_l <- pivot_longer(dat_traj, cols = -ID, names_to = "Wave", values_to = "value")
      dat_traj_l$Wave <- substring(dat_traj_l$Wave, nchar(raw_y_var) + 1L)
      dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
      dat_time_l$time <- dat_time_l$time + xstarts
      dat_time_l$Wave <- substring(dat_time_l$Wave, nchar(t_var) + 1L)
      dat_raw <- merge(dat_traj_l, dat_time_l, by = c("ID", "Wave"))
      dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
      t_seq0 <- apply(dat[, paste0(t_var, records)], 2, mean)
      t_seq <- seq(min(t_seq0), max(t_seq0), 0.1)
      dat_est <- .trajectory_plot_data(model, eval_model = model, eval_prefix = y_var,
                                        se_prefix = y_var, sub_Model = sub_Model, curveFun = curveFun,
                                        y_model = y_model, times = t_seq,
                                        baseline_time = mean(dat[[paste0(t_var, records[1])]]),
                                        records = records, raw_y_var = raw_y_var)

      dat_est$time <- dat_est$time + xstarts
      if ("interval_start" %in% names(dat_est)) dat_est$interval_start <- dat_est$interval_start + xstarts
      full <- if (tvc) data.frame(time = t_seq0 + xstarts, value = .tvc_full_mean(model, paste0(raw_y_var, records)))
      fig.status <- ggplot(data = dat_est, aes(x = time, y = Y_hat, group = 1)) +
        geom_line(aes(linetype = "Model Implied Growth Status"), color = "black", linewidth = 1) +
        geom_line(aes(x = time, y = Y_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Growth Status"), linewidth = 1) +
        geom_line(aes(x = time, y = Y_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Growth Status"), linewidth = 1) +
        .observed_curve(dat_raw, aes(x = time, y = value, group = 1,
                                        linetype = "Smooth Line of Observed Growth Status"),
                    color = "black", linewidth = 1) +
        (if (tvc) geom_line(data = full, aes(x = time, y = value, group = 1,
                                             linetype = "Model Implied Mean Outcome (TVC Effect Included)"),
                            color = "black", linewidth = 0.8)) +
        (if (tvc) geom_point(data = full, aes(x = time, y = value), inherit.aes = FALSE, color = "black", size = 1.5)) +
        labs(x = xlab, y = paste0("Measurement of ", outcome)) +
        scale_linetype_manual(values = c("Model Implied Growth Status" = "solid",
                                         "Smooth Line of Observed Growth Status" = "twodash",
                                         "95% Pointwise Confidence Interval on Model Implied Growth Status" = "dotdash",
                                         if (tvc) c("Model Implied Mean Outcome (TVC Effect Included)" = "dashed")),
                              labels = lt_labels) +
        scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
        guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3 + tvc)) +
        theme_bw() +
        theme(strip.text.y = element_text(size = 4),
              strip.background = element_rect(color = "white", fill = "white"),
              strip.placement = "outside",
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 12),
              legend.position = "bottom",
              legend.box = "vertical")
      figure <- list(fig.status)
    }
    else if (y_model == "LCSM"){
      dat_chg <- cbind(data.frame(ID = dat$ID),
                       dat[, paste0(y_var, records[-1]), drop = FALSE] - dat[[paste0(y_var, records[1])]])
      if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
        y_var <- "Y"
      }
      else if (sub_Model == "MED"){
        y_var <- med_role
      }
      names(dat_chg) <- c("ID", paste0("Delta", names(dat_chg)[-1]))
      dat_chg_l <- pivot_longer(dat_chg, cols = -ID, names_to = "Wave", values_to = "value")
      dat_chg_l$Wave <- substring(dat_chg_l$Wave, nchar(raw_y_var) + 6L)
      dat_time <- dat[, c("ID", paste0(t_var, records[-1]))]
      dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
      dat_time_l$time <- dat_time_l$time + xstarts
      dat_time_l$Wave <- substring(dat_time_l$Wave, nchar(t_var) + 1L)
      dat_raw <- merge(dat_chg_l, dat_time_l, by = c("ID", "Wave"))
      dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
      t_seq <- colMeans(dat[, paste0(t_var, records[-1]), drop = FALSE])
      dat_est <- .trajectory_plot_data(model, eval_model = model, eval_prefix = y_var,
                                        se_prefix = y_var, sub_Model = sub_Model, curveFun = curveFun,
                                        y_model = y_model, times = t_seq,
                                        baseline_time = mean(dat[[paste0(t_var, records[1])]]),
                                        records = records, raw_y_var = raw_y_var)


      dat_est$time <- dat_est$time + xstarts
      if ("interval_start" %in% names(dat_est)) dat_est$interval_start <- dat_est$interval_start + xstarts
      full <- if (tvc){
        full_mean <- .tvc_full_mean(model, paste0(raw_y_var, records))
        data.frame(time = dat_est$time, value = full_mean[-1] - full_mean[1])
      }
      fig.CHG_BL <- ggplot(data = dat_est, aes(x = time, y = cY_hat, group = 1)) +
        geom_line(aes(linetype = "Model Implied Change from Baseline"), color = "black", linewidth = 1) +
        geom_line(aes(x = time, y = cY_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Change from Baseline"), linewidth = 1) +
        geom_line(aes(x = time, y = cY_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Change from Baseline"), linewidth = 1) +
        .observed_curve(dat_raw, aes(x = time, y = value, group = 1,
                                        linetype = "Smooth Line of Observed Change from Baseline"),
                    color = "black", linewidth = 1) +
        (if (tvc) geom_line(data = full, aes(x = time, y = value, group = 1,
                                             linetype = "Model Implied Mean Change from Baseline (TVC Effect Included)"),
                            color = "black", linewidth = 0.8)) +
        (if (tvc) geom_point(data = full, aes(x = time, y = value), inherit.aes = FALSE, color = "black", size = 1.5)) +
        labs(x = xlab, y = paste0("Change from baseline of ", outcome)) +
        scale_linetype_manual(values = c("Model Implied Change from Baseline" = "solid",
                                         "Smooth Line of Observed Change from Baseline" = "twodash",
                                         "95% Pointwise Confidence Interval on Model Implied Change from Baseline" = "dotdash",
                                         if (tvc) c("Model Implied Mean Change from Baseline (TVC Effect Included)" = "dashed")),
                              labels = lt_labels) +
        scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
        guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3 + tvc)) +
        theme_bw() +
        theme(strip.text.y = element_text(size = 4),
              strip.background = element_rect(color = "white", fill = "white"),
              strip.placement = "outside",
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 12),
              legend.position = "bottom",
              legend.box = "vertical")
      if (curveFun %in% c("nonparametric", "NonP")){
        fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = 1)) +
          ggplot2::geom_segment(aes(x = interval_start, xend = time, yend = dY_hat, linetype = "Model Implied Interval-Average Growth Rate"), color = "black", linewidth = 1) +
          ggplot2::geom_segment(aes(x = interval_start, xend = time, y = dY_hat_l, yend = dY_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
          ggplot2::geom_segment(aes(x = interval_start, xend = time, y = dY_hat_u, yend = dY_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
          labs(x = xlab, y = paste0("Interval-Average Growth Rate of ", outcome)) +
          scale_linetype_manual(values = c("Model Implied Interval-Average Growth Rate" = "solid",
                                           "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate" = "dotdash"), labels = lt_labels) +
          scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
          guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2)) +
          theme_bw() +
          theme(strip.text.y = element_text(size = 4),
                strip.background = element_rect(color = "white", fill = "white"),
                strip.placement = "outside",
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 12),
                legend.position = "bottom",
                legend.box = "vertical")
      }
      else{
        fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = 1)) +
          geom_line(aes(linetype = "Model Implied Interval-Average Growth Rate"), color = "black", linewidth = 1) +
          geom_line(aes(x = time, y = dY_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
          geom_line(aes(x = time, y = dY_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
          labs(x = xlab, y = paste0("Interval-Average Growth Rate of ", outcome)) +
          scale_linetype_manual(values = c("Model Implied Interval-Average Growth Rate" = "solid",
                                           "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate" = "dotdash"), labels = lt_labels) +
          scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
          guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2)) +
          theme_bw() +
          theme(strip.text.y = element_text(size = 4),
                strip.background = element_rect(color = "white", fill = "white"),
                strip.placement = "outside",
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 12),
                legend.position = "bottom",
                legend.box = "vertical")
      }
      figure <- list(fig.CHG_BL, fig.SLP)
    }
  }
  else if (!is.null(nClass)){
    if (is.null(grp_var)){
        if (y_model == "LGCM"){
          dat$Class <- getPosterior(model = model, nClass = nClass, cluster_TIC = cluster_TIC)@membership
          dat_traj <- dat[, c("ID", "Class", paste0(y_var, records))]
          if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
            y_var <- "Y"
          }
          else if (sub_Model == "MED"){
            y_var <- med_role
          }
          dat_time <- dat[, c("ID", paste0(t_var, records))]
          dat_traj_l <- pivot_longer(dat_traj, cols = -c(ID, Class), names_to = "Wave", values_to = "value")
          dat_traj_l$Wave <- substring(dat_traj_l$Wave, nchar(raw_y_var) + 1L)
          dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
          dat_time_l$time <- dat_time_l$time + xstarts
          dat_time_l$Wave <- substring(dat_time_l$Wave, nchar(t_var) + 1L)
          dat_raw <- merge(dat_traj_l, dat_time_l, by = c("ID", "Wave"))
          dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
          t_seq0 <- apply(dat[, paste0(t_var, records)], 2, mean)
          t_seq <- seq(min(t_seq0), max(t_seq0), 0.1)
          dat_est_L <- list()
          for (k in 1:nClass){
            dat_est_L[[k]] <- .trajectory_plot_data(model, eval_model = model@submodels[[k]],
                                                   eval_prefix = paste0("c", k, y_var),
                                                   se_prefix = paste0("Class", k, ".c", k, y_var),
                                                   sub_Model = sub_Model, curveFun = curveFun,
                                                   y_model = y_model, times = t_seq,
                                        baseline_time = mean(dat[[paste0(t_var, records[1])]]),
                                        records = records, raw_y_var = raw_y_var)
            dat_est_L[[k]]$Class <- k
          }
          dat_est <- do.call(rbind, dat_est_L)

          dat_est$time <- dat_est$time + xstarts
      if ("interval_start" %in% names(dat_est)) dat_est$interval_start <- dat_est$interval_start + xstarts
          full <- if (tvc) .tvc_class_full_means(model, paste0(raw_y_var, records), paste0(t_var, records), xstarts, change = FALSE)
          fig.status <- ggplot(data = dat_est, aes(x = time, y = Y_hat, group = Class, color = as.factor(Class))) +
            geom_line(aes(linetype = "Model Implied Growth Status"), linewidth = 1) +
            geom_line(aes(x = time, y = Y_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Growth Status"), linewidth = 1) +
            geom_line(aes(x = time, y = Y_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Growth Status"), linewidth = 1) +
            .observed_curve(dat_raw, aes(x = time, y = value, group = Class, linetype = "Smooth Line of Observed Growth Status"),
                        linewidth = 1) +
            (if (tvc) geom_line(data = full, aes(x = time, y = value, group = Class, color = as.factor(Class),
                                                 linetype = "Model Implied Mean Outcome (TVC Effect Included)"), linewidth = 0.8)) +
            (if (tvc) geom_point(data = full, aes(x = time, y = value, color = as.factor(Class)), inherit.aes = FALSE, size = 1.5)) +
            labs(x = xlab, y = paste0("Growth Status of ", outcome), color = "Class") +
            scale_linetype_manual(values = c("Model Implied Growth Status" = "solid",
                                             "Smooth Line of Observed Growth Status" = "twodash",
                                             "95% Pointwise Confidence Interval on Model Implied Growth Status" = "dotdash",
                                             if (tvc) c("Model Implied Mean Outcome (TVC Effect Included)" = "dashed")), labels = lt_labels) +
            scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
            scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
            guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3 + tvc, order = 1),
                   color = guide_legend(title = "Class", order = 2)) +
            theme_bw() +
            theme(strip.text.y = element_text(size = 4),
                  strip.background = element_rect(color = "white", fill = "white"),
                  strip.placement = "outside",
                  legend.position = "bottom",
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12),
                  legend.box = "vertical",
                  legend.margin = margin(-3, 0, -3, 0))
          figure <- list(fig.status)
        }
        else if (y_model == "LCSM"){
          dat$Class <- getPosterior(model = model, nClass = nClass, cluster_TIC = cluster_TIC)@membership
          dat_chg <- cbind(dat[, c("ID", "Class")],
                           dat[, paste0(y_var, records[-1]), drop = FALSE] - dat[[paste0(y_var, records[1])]])
          if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
            y_var <- "Y"
          }
          else if (sub_Model == "MED"){
            y_var <- med_role
          }
          names(dat_chg) <- c("ID", "Class", paste0("Delta", names(dat_chg)[-c(1, 2)]))
          dat_chg_l <- pivot_longer(dat_chg, cols = -c(ID, Class), names_to = "Wave", values_to = "value")
          dat_chg_l$Wave <- substring(dat_chg_l$Wave, nchar(raw_y_var) + 6L)
          dat_time <- dat[, c("ID", paste0(t_var, records[-1]))]
          dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
          dat_time_l$time <- dat_time_l$time + xstarts
          dat_time_l$Wave <- substring(dat_time_l$Wave, nchar(t_var) + 1L)
          dat_raw <- merge(dat_chg_l, dat_time_l, by = c("ID", "Wave"))
          dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
          t_seq <- colMeans(dat[, paste0(t_var, records[-1]), drop = FALSE])
          dat_est_L <- list()
          for (k in 1:nClass){
            dat_est_L[[k]] <- .trajectory_plot_data(model, eval_model = model@submodels[[k]],
                                                   eval_prefix = paste0("c", k, y_var),
                                                   se_prefix = paste0("Class", k, ".c", k, y_var),
                                                   sub_Model = sub_Model, curveFun = curveFun,
                                                   y_model = y_model, times = t_seq,
                                        baseline_time = mean(dat[[paste0(t_var, records[1])]]),
                                        records = records, raw_y_var = raw_y_var)
            dat_est_L[[k]]$Class <- k
          }
          dat_est <- do.call(rbind, dat_est_L)


          dat_est$time <- dat_est$time + xstarts
      if ("interval_start" %in% names(dat_est)) dat_est$interval_start <- dat_est$interval_start + xstarts
          full <- if (tvc) .tvc_class_full_means(model, paste0(raw_y_var, records), paste0(t_var, records), xstarts, change = TRUE)
          fig.CHG_BL <- ggplot(data = dat_est, aes(x = time, y = cY_hat, group = Class, color = as.factor(Class))) +
            geom_line(aes(linetype = "Model Implied Change from Baseline"), linewidth = 1) +
            geom_line(aes(x = time, y = cY_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Change from Baseline"), linewidth = 1) +
            geom_line(aes(x = time, y = cY_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Change from Baseline"), linewidth = 1) +
            .observed_curve(dat_raw, aes(x = time, y = value, group = Class, linetype = "Smooth Line of Observed Change from Baseline"),
                        linewidth = 1) +
            (if (tvc) geom_line(data = full, aes(x = time, y = value, group = Class, color = as.factor(Class),
                                                 linetype = "Model Implied Mean Change from Baseline (TVC Effect Included)"), linewidth = 0.8)) +
            (if (tvc) geom_point(data = full, aes(x = time, y = value, color = as.factor(Class)), inherit.aes = FALSE, size = 1.5)) +
            labs(x = xlab, y = paste0("Change from baseline of ", outcome), color = "Class") +
            scale_linetype_manual(values = c("Model Implied Change from Baseline" = "solid",
                                             "Smooth Line of Observed Change from Baseline" = "twodash",
                                             "95% Pointwise Confidence Interval on Model Implied Change from Baseline" = "dotdash",
                                             if (tvc) c("Model Implied Mean Change from Baseline (TVC Effect Included)" = "dashed")), labels = lt_labels) +
            scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
            scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
            guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3 + tvc, order = 1),
                   color = guide_legend(title = "Class", order = 2)) +
            theme_bw() +
            theme(strip.text.y = element_text(size = 4),
                  strip.background = element_rect(color = "white", fill = "white"),
                  strip.placement = "outside",
                  legend.position = "bottom",
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12),
                  legend.box = "vertical",
                  legend.margin = margin(-3, 0, -3, 0))
          if (curveFun %in% c("nonparametric", "NonP")){
            fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = Class, color = as.factor(Class))) +
              ggplot2::geom_segment(aes(x = interval_start, xend = time, yend = dY_hat, linetype = "Model Implied Interval-Average Growth Rate"), linewidth = 1) +
              ggplot2::geom_segment(aes(x = interval_start, xend = time, y = dY_hat_l, yend = dY_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
              ggplot2::geom_segment(aes(x = interval_start, xend = time, y = dY_hat_u, yend = dY_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
              labs(x = xlab, y = paste0("Interval-Average Growth Rate of ", outcome), color = "Class") +
              scale_linetype_manual(values = c("Model Implied Interval-Average Growth Rate" = "solid",
                                               "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate" = "dotdash"), labels = lt_labels) +
              scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
              scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
              guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2, order = 1),
                     color = guide_legend(title = "Class", order = 2)) +
              theme_bw() +
              theme(strip.text.y = element_text(size = 4),
                    strip.background = element_rect(color = "white", fill = "white"),
                    strip.placement = "outside",
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 12),
                    legend.position = "bottom",
                    legend.box = "vertical",
                    legend.margin = margin(-3, 0, -3, 0))
          }
          else{
            fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = Class, color = as.factor(Class))) +
              geom_line(aes(linetype = "Model Implied Interval-Average Growth Rate"), linewidth = 1) +
              geom_line(aes(x = time, y = dY_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
              geom_line(aes(x = time, y = dY_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
              labs(x = xlab, y = paste0("Interval-Average Growth Rate of ", outcome), color = "Class") +
              scale_linetype_manual(values = c("Model Implied Interval-Average Growth Rate" = "solid",
                                               "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate" = "dotdash"), labels = lt_labels) +
              scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
              scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
              guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2, order = 1),
                     color = guide_legend(title = "Class", order = 2)) +
              theme_bw() +
              theme(strip.text.y = element_text(size = 4),
                    strip.background = element_rect(color = "white", fill = "white"),
                    strip.placement = "outside",
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 12),
                    legend.position = "bottom",
                    legend.box = "vertical",
                    legend.margin = margin(-3, 0, -3, 0))
          }
          figure <- list(fig.CHG_BL, fig.SLP)
        }
    }
    else if (!is.null(grp_var)){
      if (y_model == "LGCM"){
        dat_traj <- data.frame(ID = dat$ID, Class = group_codes,
                               dat[, paste0(y_var, records), drop = FALSE], check.names = FALSE)
        if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
          y_var <- "Y"
        }
        else if (sub_Model == "MED"){
          y_var <- med_role
        }
        dat_time <- dat[, c("ID", paste0(t_var, records))]
        dat_traj_l <- pivot_longer(dat_traj, cols = -c(ID, Class), names_to = "Wave", values_to = "value")
        dat_traj_l$Wave <- substring(dat_traj_l$Wave, nchar(raw_y_var) + 1L)
        dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
        dat_time_l$time <- dat_time_l$time + xstarts
        dat_time_l$Wave <- substring(dat_time_l$Wave, nchar(t_var) + 1L)
        dat_raw <- merge(dat_traj_l, dat_time_l, by = c("ID", "Wave"))
        dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
        t_seq0 <- apply(dat[, paste0(t_var, records)], 2, mean)
        t_seq <- seq(min(t_seq0), max(t_seq0), 0.1)
        dat_est_L <- list()
        for (k in 1:nClass){
          dat_est_L[[k]] <- .trajectory_plot_data(model, eval_model = model@submodels[[k]],
                                                 eval_prefix = paste0("c", k, y_var),
                                                 se_prefix = paste0("Class", k, ".c", k, y_var),
                                                 sub_Model = sub_Model, curveFun = curveFun,
                                                 y_model = y_model, times = t_seq,
                                        baseline_time = mean(dat[[paste0(t_var, records[1])]]),
                                        records = records, raw_y_var = raw_y_var)
          dat_est_L[[k]]$Class <- k
        }
        dat_est <- do.call(rbind, dat_est_L)

        dat_est$time <- dat_est$time + xstarts
      if ("interval_start" %in% names(dat_est)) dat_est$interval_start <- dat_est$interval_start + xstarts
        full <- if (tvc) .tvc_class_full_means(model, paste0(raw_y_var, records), paste0(t_var, records), xstarts, change = FALSE)
        fig.status <- ggplot(data = dat_est, aes(x = time, y = Y_hat, group = Class, color = as.factor(Class))) +
          geom_line(aes(linetype = "Model Implied Growth Status"), linewidth = 1) +
          geom_line(aes(x = time, y = Y_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Growth Status"), linewidth = 1) +
          geom_line(aes(x = time, y = Y_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Growth Status"), linewidth = 1) +
          .observed_curve(dat_raw, aes(x = time, y = value, group = Class, linetype = "Smooth Line of Observed Growth Status"),
                      linewidth = 1) +
          (if (tvc) geom_line(data = full, aes(x = time, y = value, group = Class, color = as.factor(Class),
                                               linetype = "Model Implied Mean Outcome (TVC Effect Included)"), linewidth = 0.8)) +
          (if (tvc) geom_point(data = full, aes(x = time, y = value, color = as.factor(Class)), inherit.aes = FALSE, size = 1.5)) +
          labs(x = xlab, y = paste0("Growth Status of ", outcome), color = "Class") +
          scale_linetype_manual(values = c("Model Implied Growth Status" = "solid",
                                           "Smooth Line of Observed Growth Status" = "twodash",
                                           "95% Pointwise Confidence Interval on Model Implied Growth Status" = "dotdash",
                                             if (tvc) c("Model Implied Mean Outcome (TVC Effect Included)" = "dashed")), labels = lt_labels) +
          scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
          scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
          guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3 + tvc, order = 1),
                 color = guide_legend(title = grp_var, order = 2)) +
          theme_bw() +
          theme(strip.text.y = element_text(size = 4),
                strip.background = element_rect(color = "white", fill = "white"),
                strip.placement = "outside",
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 12),
                legend.position = "bottom",
                legend.box = "vertical",
                legend.margin = margin(-3, 0, -3, 0))
        figure <- list(fig.status)
      }
      else if (y_model == "LCSM"){
        dat_chg <- cbind(data.frame(ID = dat$ID, Class = group_codes),
                         dat[, paste0(y_var, records[-1]), drop = FALSE] - dat[[paste0(y_var, records[1])]])
        if (sub_Model %in% c("LGCM", "LCSM", "TVC")){
          y_var <- "Y"
        }
        else if (sub_Model == "MED"){
          y_var <- med_role
        }
        names(dat_chg)[2] <- "Class"
        names(dat_chg) <- c("ID", "Class", paste0("Delta", names(dat_chg)[-c(1, 2)]))
        dat_chg_l <- pivot_longer(dat_chg, cols = -c(ID, Class), names_to = "Wave", values_to = "value")
        dat_chg_l$Wave <- substring(dat_chg_l$Wave, nchar(raw_y_var) + 6L)
        dat_time <- dat[, c("ID", paste0(t_var, records[-1]))]
        dat_time_l <- pivot_longer(dat_time, cols = -ID, names_to = "Wave", values_to = "time")
        dat_time_l$time <- dat_time_l$time + xstarts
        dat_time_l$Wave <- substring(dat_time_l$Wave, nchar(t_var) + 1L)
        dat_raw <- merge(dat_chg_l, dat_time_l, by = c("ID", "Wave"))
        dat_raw <- dat_raw[order(dat_raw$ID, as.numeric(dat_raw$Wave)), ]
        t_seq <- colMeans(dat[, paste0(t_var, records[-1]), drop = FALSE])
        dat_est_L <- list()
        for (k in 1:nClass){
          group_dat <- as.data.frame(model@submodels[[k]]$data$observed)
          t_seq <- colMeans(group_dat[, paste0(t_var, records[-1]), drop = FALSE])
          dat_est_L[[k]] <- .trajectory_plot_data(model, eval_model = model@submodels[[k]],
                                                 eval_prefix = paste0("c", k, y_var),
                                                 se_prefix = paste0("Class", k, ".c", k, y_var),
                                                 sub_Model = sub_Model, curveFun = curveFun,
                                                 y_model = y_model, times = t_seq,
                                        baseline_time = mean(group_dat[[paste0(t_var, records[1])]]),
                                        records = records, raw_y_var = raw_y_var)
          dat_est_L[[k]]$Class <- k
        }
        dat_est <- do.call(rbind, dat_est_L)


        dat_est$time <- dat_est$time + xstarts
      if ("interval_start" %in% names(dat_est)) dat_est$interval_start <- dat_est$interval_start + xstarts
        full <- if (tvc) .tvc_class_full_means(model, paste0(raw_y_var, records), paste0(t_var, records), xstarts, change = TRUE)
        fig.CHG_BL <- ggplot(data = dat_est, aes(x = time, y = cY_hat, group = Class, color = as.factor(Class))) +
          geom_line(aes(linetype = "Model Implied Change from Baseline"),  linewidth = 1) +
          geom_line(aes(x = time, y = cY_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Change from Baseline"), linewidth = 1) +
          geom_line(aes(x = time, y = cY_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Change from Baseline"), linewidth = 1) +
          .observed_curve(dat_raw, aes(x = time, y = value, group = Class,
                                          linetype = "Smooth Line of Observed Change from Baseline"),
                      linewidth = 1) +
          (if (tvc) geom_line(data = full, aes(x = time, y = value, group = Class, color = as.factor(Class),
                                               linetype = "Model Implied Mean Change from Baseline (TVC Effect Included)"), linewidth = 0.8)) +
          (if (tvc) geom_point(data = full, aes(x = time, y = value, color = as.factor(Class)), inherit.aes = FALSE, size = 1.5)) +
          labs(x = xlab, y = paste0("Change from baseline of ", outcome), color = "Class") +
          scale_linetype_manual(values = c("Model Implied Change from Baseline" = "solid",
                                           "Smooth Line of Observed Change from Baseline" = "twodash",
                                           "95% Pointwise Confidence Interval on Model Implied Change from Baseline" = "dotdash",
                                             if (tvc) c("Model Implied Mean Change from Baseline (TVC Effect Included)" = "dashed")), labels = lt_labels) +
          scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
          scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
          guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 3 + tvc, order = 1),
                 color = guide_legend(title = grp_var, order = 2)) +
          theme_bw() +
          theme(strip.text.y = element_text(size = 4),
                strip.background = element_rect(color = "white", fill = "white"),
                strip.placement = "outside",
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 12),
                legend.position = "bottom",
                legend.box = "vertical",
                legend.margin = margin(-3, 0, -3, 0))
        if (curveFun %in% c("nonparametric", "NonP")){
          fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = Class, color = as.factor(Class))) +
            ggplot2::geom_segment(aes(x = interval_start, xend = time, yend = dY_hat, linetype = "Model Implied Interval-Average Growth Rate"),  linewidth = 1) +
            ggplot2::geom_segment(aes(x = interval_start, xend = time, y = dY_hat_l, yend = dY_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
            ggplot2::geom_segment(aes(x = interval_start, xend = time, y = dY_hat_u, yend = dY_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
            labs(x = xlab, y = paste0("Interval-Average Growth Rate of ", outcome), color = "Class") +
            scale_linetype_manual(values = c("Model Implied Interval-Average Growth Rate" = "solid",
                                             "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate" = "dotdash"), labels = lt_labels) +
            scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
            scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
            guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2, order = 1),
                   color = guide_legend(title = grp_var, order = 2)) +
            theme_bw() +
            theme(strip.text.y = element_text(size = 4),
                  strip.background = element_rect(color = "white", fill = "white"),
                  strip.placement = "outside",
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12),
                  legend.position = "bottom",
                  legend.box = "vertical",
                  legend.margin = margin(-3, 0, -3, 0))
        }
        else{
          fig.SLP <- ggplot(data = dat_est, aes(x = time, y = dY_hat, group = Class, color = as.factor(Class))) +
            geom_line(aes(linetype = "Model Implied Interval-Average Growth Rate"),  linewidth = 1) +
            geom_line(aes(x = time, y = dY_hat_l, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
            geom_line(aes(x = time, y = dY_hat_u, linetype = "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate"), linewidth = 1) +
            labs(x = xlab, y = paste0("Interval-Average Growth Rate of ", outcome), color = "Class") +
            scale_linetype_manual(values = c("Model Implied Interval-Average Growth Rate" = "solid",
                                             "95% Pointwise Confidence Interval on Model Implied Interval-Average Growth Rate" = "dotdash"), labels = lt_labels) +
            scale_x_continuous(breaks = ceiling(seq(from = min(dat_est$time), to = max(dat_est$time), length.out = length(records)))) +
            scale_color_manual(values = colorRampPalette(c("grey20", "grey80"))(nClass), labels = class_labels) +
            guides(linetype = guide_legend(title = "Trajectory Type: ", nrow = 2, order = 1),
                   color = guide_legend(title = grp_var, order = 2)) +
            theme_bw() +
            theme(strip.text.y = element_text(size = 4),
                  strip.background = element_rect(color = "white", fill = "white"),
                  strip.placement = "outside",
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12),
                  legend.position = "bottom",
                  legend.box = "vertical",
                  legend.margin = margin(-3, 0, -3, 0))
        }
        figure <- list(fig.CHG_BL, fig.SLP)
      }
    }
  }
  if (y_model == "LCSM" && curveFun %in% c("nonparametric", "NonP")){
    rate_scale <- figure[[2]]$scales$get_scales("x")
    rate_scale$breaks <- pretty(range(c(figure[[2]]$data$interval_start, figure[[2]]$data$time)), n = length(records))
  }
  if (y_model == "LCSM" && length(records) == 2L){
    # One selected interval has one endpoint per class. Lines cannot display
    # that estimate or its uncertainty, so draw a point and its vertical CI.
    for (j in seq_along(figure)){
      mapping <- if (j == 1L) aes(ymin = cY_hat_l, ymax = cY_hat_u,
                                   linetype = "Model Implied Change from Baseline") else
        aes(ymin = dY_hat_l, ymax = dY_hat_u, linetype = "Model Implied Interval-Average Growth Rate")
      figure[[j]]$layers <- c(list(ggplot2::geom_pointrange(mapping = mapping, linewidth = 0.8)),
                               figure[[j]]$layers[-seq_len(3L)])
      x_scale <- figure[[j]]$scales$get_scales("x")
      x_scale$breaks <- sort(unique(figure[[j]]$data$time))
      figure[[j]]$labels$caption <- "Points and vertical bars: fitted means and 95% pointwise confidence intervals."
    }
  }
  return(figure)
}

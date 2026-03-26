#' @title Derive Individual Growth Factors for Latent Growth Curve Models or Latent Change Score Models with Time-Invariant Covariates (If Any)
#'
#' @description This function derives individual growth factors for the specified latent growth curve model or
#' latent change score model from raw data. These individual growth factors help further compute initial values
#' for parameters related to growth factors, time-invariant covariates (if any), and path coefficients (if any).
#'
#' @param dat_traj A data frame containing the records for the repeated measurements.
#' @param dat_time A data frame containing the records for measurement occasions associated with the repeated
#' measurements.
#' @param nT An integer representing the number of repeated measurements.
#' @param curveFun A string specifying the functional form of the growth curve. Supported options for \code{y_model =
#' "LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
#' options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
#' (or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the
#' value passed from \code{getUNI.initial()} or \code{getMULTI.initial()} or \code{getTVC.initial()} or \code{getMIX.initial()}.
#'
#' @return A data frame containing the derived individual growth factors from the raw data. The returned data frame
#' also carries an attribute \code{"mean_res_var"} — the median of per-person mean squared residuals from the
#' individual model fits — which can be used as a data-driven initial value for the residual variance parameter.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats coef nls nls.control residuals median
#'
.flatten_curve_data <- function(dat_traj, dat_time) {
  measures <- as.numeric(as.matrix(dat_traj))
  times <- as.numeric(as.matrix(dat_time))
  keep <- is.finite(measures) & is.finite(times)
  list(measures = measures[keep], times = times[keep])
}

.fit_basis_curve <- function(measures, times, basis_fun, param) {
  keep <- is.finite(measures) & is.finite(times)
  y <- as.numeric(measures[keep])
  t <- as.numeric(times[keep])
  if (length(y) < 2) {
    return(list(coef = numeric(0), rss = Inf, res_var = NA_real_, kappa = Inf))
  }

  basis <- basis_fun(t, param)
  if (is.null(dim(basis))) basis <- matrix(basis, ncol = 1)
  design <- cbind(Intercept = 1, basis)
  fit <- stats::lm.fit(x = design, y = y)
  coef <- fit$coefficients
  coef[!is.finite(coef)] <- 0
  fitted <- drop(design %*% coef)
  design_kappa <- suppressWarnings(kappa(design))
  if (!is.finite(design_kappa)) design_kappa <- Inf

  list(
    coef = as.numeric(coef),
    rss = sum((y - fitted)^2),
    res_var = mean((y - fitted)^2),
    kappa = design_kappa
  )
}

.grid_curve_fits <- function(measures, times, param_grid, basis_fun) {
  param_grid <- sort(unique(param_grid[is.finite(param_grid)]))
  if (length(param_grid) == 0) {
    return(list(list(param = NA_real_, coef = numeric(0), rss = Inf, res_var = NA_real_, kappa = Inf)))
  }

  lapply(param_grid, function(param) c(list(param = param), .fit_basis_curve(measures, times, basis_fun, param)))
}

.grid_search_curve <- function(measures, times, param_grid, basis_fun) {
  fits <- .grid_curve_fits(measures, times, param_grid, basis_fun)
  rss <- vapply(fits, `[[`, numeric(1), "rss")
  fits[[which.min(rss)]]
}

.select_stable_curve_fit <- function(fits, rss_tol = 0.001) {
  if (length(fits) == 0) {
    return(list(param = NA_real_, coef = numeric(0), rss = Inf, res_var = NA_real_, kappa = Inf))
  }

  rss <- vapply(fits, `[[`, numeric(1), "rss")
  params <- vapply(fits, `[[`, numeric(1), "param")
  kappas <- vapply(fits, `[[`, numeric(1), "kappa")
  min_rss <- min(rss)
  keep <- which(rss <= min_rss * (1 + rss_tol))
  if (length(keep) == 0) keep <- which.min(rss)
  keep <- keep[order(kappas[keep], -abs(params[keep]), rss[keep])]
  fits[[keep[1]]]
}

.shrink_curve_param <- function(anchor, param, anchor_weight = 0.75) {
  if (!is.finite(anchor)) return(param)
  if (!is.finite(param)) return(anchor)
  anchor_weight * anchor + (1 - anchor_weight) * param
}

.curve_span <- function(times) {
  rng <- range(times[is.finite(times)])
  span <- diff(rng)
  if (!is.finite(span) || span <= 0) 1 else span
}

.exp_param_grid <- function(times, anchor = NULL) {
  span <- .curve_span(times)
  lower <- max(1e-6, 1e-4 / span)
  upper <- max(lower * 25, 5 / span)

  if (is.null(anchor) || !is.finite(anchor) || anchor <= 0) {
    exp(seq(log(lower), log(upper), length.out = 25))
  } else {
    local_lower <- max(lower, anchor / 4)
    local_upper <- min(upper, anchor * 4)
    unique(c(anchor, exp(seq(log(local_lower), log(local_upper), length.out = 15))))
  }
}

.jb_param_grid <- function(times, anchor = NULL) {
  span <- .curve_span(times)
  max_abs <- max(0.05, 3 / span)
  min_abs <- max(1e-4 / span, max_abs / 100)

  if (is.null(anchor) || !is.finite(anchor) || abs(anchor) < min_abs) {
    pos <- exp(seq(log(min_abs), log(max_abs), length.out = 15))
    unique(c(-rev(pos), pos))
  } else {
    delta <- max(abs(anchor) * 1.5, max_abs / 4)
    local_grid <- seq(max(-max_abs, anchor - delta), min(max_abs, anchor + delta), length.out = 21)
    local_grid <- local_grid[abs(local_grid) >= min_abs]
    if (length(local_grid) == 0) local_grid <- c(sign(anchor) * min_abs, anchor)
    unique(local_grid)
  }
}

.bls_knot_grid <- function(times, anchor = NULL) {
  times <- sort(unique(times[is.finite(times)]))
  if (length(times) < 2) return(times)

  mids <- (times[-1] + times[-length(times)]) / 2
  if (length(mids) > 15) {
    mids <- unique(as.numeric(stats::quantile(mids, probs = seq(0.15, 0.85, length.out = 15), names = FALSE)))
  }
  if (!is.null(anchor) && is.finite(anchor) && anchor > min(times) && anchor < max(times)) {
    mids <- c(mids, anchor)
  }
  sort(unique(mids))
}

getUNI.GF <- function(dat_traj, dat_time, nT, curveFun){
  # Calculate individual growth factors eta0 (initial values) and eta1 (slopes) for linear trajectories.
  if (curveFun %in% c("linear", "LIN")){
    # Apply the linear model to each row
    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- ind[1:nT]
      times <- ind[-(1:nT)]
      ind_model <- lm(measures ~ times)
      c(eta0 = unname(coef(ind_model)[1]), eta1 = unname(coef(ind_model)[2]),
        res_var = mean(residuals(ind_model)^2))
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2])
    attr(growth_factor, "mean_res_var") <- median(coef_mat[, ncol(coef_mat)], na.rm = TRUE)
  }
  # Calculate individual growth factors eta0 (initial values), eta1 (linear slopes), and eta2 (quadratic slopes)
  # for quadratic trajectories.
  else if (curveFun %in% c("quadratic", "QUAD")){
    # Apply the quadratic model to each row
    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- ind[1:nT]
      times <- ind[-(1:nT)]
      times_sq <- times^2
      ind_model <- lm(measures ~ times + times_sq)
      c(eta0 = unname(coef(ind_model)[1]), eta1 = unname(coef(ind_model)[2]),
        eta2 = unname(coef(ind_model)[3]), res_var = mean(residuals(ind_model)^2))
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2], eta2 = coef_mat[, 3])
    attr(growth_factor, "mean_res_var") <- median(coef_mat[, ncol(coef_mat)], na.rm = TRUE)
  }
  # Calculate individual growth factors eta0 (initial values), eta1 (asymptotic change between final and initial
  # values), and gamma (growth rate parameter) for negative exponential trajectories.
  else if (curveFun %in% c("negative exponential", "EXP")){
    basis_fun <- function(t, param) 1 - exp(-param * t)
    pooled_data <- .flatten_curve_data(dat_traj, dat_time)
    pooled_fit <- .grid_search_curve(
      pooled_data$measures,
      pooled_data$times,
      .exp_param_grid(pooled_data$times),
      basis_fun
    )
    c_anchor <- pooled_fit$param

    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- as.numeric(ind[1:nT])
      times <- as.numeric(ind[-(1:nT)])
      person_fit <- .grid_search_curve(measures, times, .exp_param_grid(times, anchor = c_anchor), basis_fun)
      c(eta0 = person_fit$coef[1], eta1 = person_fit$coef[2], gamma = person_fit$param,
        res_var = person_fit$res_var)
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2], gamma = coef_mat[, 3])
    attr(growth_factor, "mean_res_var") <- median(coef_mat[, ncol(coef_mat)], na.rm = TRUE)
  }
  # Calculate individual growth factors eta0 (initial values), eta1 (slope of linear asymptote), eta2 (change
  # between intercept of linear asymptote and initial status), and gamma (growth acceleration parameter) for
  # Jenss-Bayley trajectories.
  else if (curveFun %in% c("Jenss-Bayley", "JB")){
    basis_fun <- function(t, param) cbind(t, exp(param * t) - 1)
    pooled_data <- .flatten_curve_data(dat_traj, dat_time)
    pooled_fit <- .select_stable_curve_fit(.grid_curve_fits(
      pooled_data$measures,
      pooled_data$times,
      .jb_param_grid(pooled_data$times),
      basis_fun
    ), rss_tol = 0.001)
    gamma_anchor <- pooled_fit$param

    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- as.numeric(ind[1:nT])
      times <- as.numeric(ind[-(1:nT)])
      person_fit <- .select_stable_curve_fit(
        .grid_curve_fits(measures, times, .jb_param_grid(times, anchor = gamma_anchor), basis_fun),
        rss_tol = 0.01
      )
      gamma_start <- .shrink_curve_param(gamma_anchor, person_fit$param, anchor_weight = 0.75)
      refit <- .fit_basis_curve(measures, times, basis_fun, gamma_start)
      if (length(refit$coef) < 3 || any(!is.finite(refit$coef[1:3]))) {
        refit <- person_fit
        gamma_start <- person_fit$param
      }
      c(eta0 = refit$coef[1], eta1 = refit$coef[2], eta2 = refit$coef[3],
        gamma = gamma_start, res_var = refit$res_var)
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2], eta2 = coef_mat[, 3],
                                gamma = coef_mat[, 4])
    attr(growth_factor, "mean_res_var") <- median(coef_mat[, ncol(coef_mat)], na.rm = TRUE)
  }
  # Calculate individual growth factors eta0 (initial values), eta1 (slopes in the 1st stage), eta2 (slopes in
  # the 2nd stage), and gamma (knot locations) for bilinear spline trajectories.
  else if (curveFun %in% c("bilinear spline", "BLS")){
    basis_fun <- function(t, param) cbind(pmin(t, param), pmax(t - param, 0))
    pooled_data <- .flatten_curve_data(dat_traj, dat_time)
    pooled_fit <- .grid_search_curve(
      pooled_data$measures,
      pooled_data$times,
      .bls_knot_grid(pooled_data$times),
      basis_fun
    )
    knot_anchor <- pooled_fit$param

    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- as.numeric(ind[1:nT])
      times <- as.numeric(ind[-(1:nT)])
      person_fit <- .grid_search_curve(measures, times, .bls_knot_grid(times, anchor = knot_anchor), basis_fun)
      c(eta0 = person_fit$coef[1], eta1 = person_fit$coef[2], eta2 = person_fit$coef[3],
        gamma = person_fit$param, res_var = person_fit$res_var)
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2], eta2 = coef_mat[, 3],
                                gamma = coef_mat[, 4])
    attr(growth_factor, "mean_res_var") <- median(coef_mat[, ncol(coef_mat)], na.rm = TRUE)
  }
  # Calculate individual growth factors eta0 (initial values) and eta1 (slopes in the 1st time interval) for
  # nonparametric trajectories.
  else if (curveFun %in% c("nonparametric", "NonP")){
    # NonP-1: Compute interval-specific mean rates using mean(dy/dt) for clarity and robustness.
    # Each slp[j] stores the mean instantaneous rate (intercept only) for interval j.
    slp <- rep(0, nT - 1)
    for (j in 1:(nT - 1)){
      dy <- dat_traj[, j + 1] - dat_traj[, j]
      dt <- dat_time[, j + 1] - dat_time[, j]
      rate <- dy / dt
      slp[j] <- mean(rate, na.rm = TRUE)
    }
    eta0 <- dat_traj[, 1]
    time_delta <- dat_time[, 2] - dat_time[, 1]
    eta1 <- (dat_traj[, 2] - dat_traj[, 1]) / time_delta
    growth_factor <- data.frame(eta0 = eta0, eta1 = eta1)
    # For nonparametric, compute residuals as deviations from piecewise-linear interpolation
    # using the interval-specific slopes
    fitted_vals <- matrix(NA, nrow = nrow(dat_traj), ncol = nT)
    fitted_vals[, 1] <- dat_traj[, 1]
    for (j in 2:nT){
      fitted_vals[, j] <- fitted_vals[, j - 1] + slp[j - 1] * (dat_time[, j] - dat_time[, j - 1])
    }
    res_var_per_person <- rowMeans((as.matrix(dat_traj) - fitted_vals)^2, na.rm = TRUE)
    attr(growth_factor, "mean_res_var") <- median(res_var_per_person, na.rm = TRUE)
  }
  return(growth_factor)
}

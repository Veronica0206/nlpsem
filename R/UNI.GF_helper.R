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
#' @importFrom stats coef residuals median
#'
.curve_span <- function(times) {
  rng <- range(times[is.finite(times)])
  span <- diff(rng)
  if (!is.finite(span) || span <= 0) 1 else span
}

# Restore the 0.3 per-person nonlinear least-squares policy. The explicitly
# recorded controls are stats::nls defaults, including the default Gauss-Newton
# algorithm. Neither this policy nor a failed NLS fit changes the SEM optimizer.
.nls_initialization_control <- function() {
  stats::nls.control(maxiter = 50L, tol = 1e-5, minFactor = 1 / 1024,
                    printEval = FALSE, warnOnly = FALSE)
}

.fit_nls_initialization <- function(formula, data, start) {
  stats::nls(formula, data = data, start = start, algorithm = "default",
             control = .nls_initialization_control())
}

# Historical seeds with bounded handling of incomplete/repeated occasions.
# Complete observations are sorted by time for seed calculation. Repeated times
# are averaged only for finite-difference seeds; NLS uses every complete pair.
.nls_curve_seed <- function(measures, times, curve) {
  unique_times <- sort(unique(times))
  y <- vapply(unique_times, function(t) mean(measures[times == t]), numeric(1))
  t <- unique_times
  n <- length(t)
  if (n < 2L || !all(is.finite(y))) return(NULL)
  span <- t[n] - t[1L]
  slope <- (y[n] - y[1L]) / span
  regression_slope <- function(x, response) {
    if (length(x) < 2L || length(unique(x)) < 2L || any(!is.finite(response))) return(NA_real_)
    unname(stats::coef(stats::lm(response ~ x))[2L])
  }
  if (curve == "EXP") {
    rate <- abs(regression_slope(t[-n], log(abs(diff(y) / diff(t)) + 1e-9)))
    # Two occasions or a constant finite-difference slope cannot identify the
    # nonlinear rate. A finite time-scale seed is reported as a heuristic.
    if (!is.finite(rate) || rate <= 0) rate <- 1 / span
    return(list(a = y[1L], b = y[n] - y[1L], c = rate))
  }
  if (curve == "JB") {
    rate <- if (n >= 3L) {
      -regression_slope(t[-c(1L, n)], diff(y, differences = 2) / (diff(t)^2)[-1L])
    } else NA_real_
    if (!is.finite(rate)) rate <- -1 / span
    return(list(a = y[1L], b = slope, c = y[n] - y[1L], d = rate))
  }
  knot <- (t[1L] + t[n]) / 2
  # Use indices in the original sorted vector, not indices in a filtered vector.
  # This corrects the historical right-segment index mapping.
  left <- which(t <= knot)
  right <- which(t >= knot)
  left_end <- utils::tail(left, 1L)
  right_start <- right[1L]
  first <- if (left_end > 1L) (y[left_end] - y[1L]) / (t[left_end] - t[1L]) else slope
  second <- if (right_start < n) (y[n] - y[right_start]) / (t[n] - t[right_start]) else slope
  list(a = y[1L], b = first, c = second, d = knot)
}

.nls_curve_growth_factors <- function(dat_traj, dat_time, nT, curve) {
  columns <- if (curve == "EXP") c("eta0", "eta1", "gamma") else c("eta0", "eta1", "eta2", "gamma")
  n <- nrow(dat_traj)
  factors <- matrix(NA_real_, nrow = n, ncol = length(columns), dimnames = list(NULL, columns))
  residual_matrix <- matrix(NA_real_, nrow = n, ncol = nT)
  residual_variances <- rep(NA_real_, n)
  audit <- data.frame(person = seq_len(n), method = rep("unavailable", n),
                      nls_converged = rep(FALSE, n), fallback_reason = rep("", n),
                      attempts = integer(n), n_complete = integer(n), warnings = rep("", n),
                      stringsAsFactors = FALSE)
  formula <- switch(EXPR = curve,
    EXP = measures ~ a + b * (1 - exp(-c * times)),
    JB = measures ~ a + b * times + c * (exp(d * times) - 1),
    BLS = measures ~ a + b * pmin(times, d) + c * pmax(times - d, 0))
  evaluate <- switch(EXPR = curve,
    EXP = function(t, p) p[1L] + p[2L] * (1 - exp(-p[3L] * t)),
    JB = function(t, p) p[1L] + p[2L] * t + p[3L] * (exp(p[4L] * t) - 1),
    BLS = function(t, p) p[1L] + p[2L] * pmin(t, p[4L]) + p[3L] * pmax(t - p[4L], 0))
  for (i in seq_len(n)) {
    measures <- as.numeric(dat_traj[i, ])
    times <- as.numeric(dat_time[i, ])
    complete <- is.finite(measures) & is.finite(times)
    audit$n_complete[i] <- sum(complete)
    order_complete <- which(complete)[order(times[complete])]
    y <- measures[order_complete]
    t <- times[order_complete]
    seed <- .nls_curve_seed(y, t, curve)
    if (is.null(seed)) {
      audit$fallback_reason[i] <- "Fewer than two distinct complete times; no finite heuristic start"
      next
    }
    coefficients <- unlist(seed, use.names = FALSE)
    warning_messages <- character()
    error_message <- ""
    fit <- NULL
    if (length(unique(t)) < length(seed)) {
      error_message <- "Too few distinct complete times for the nonlinear parameter count"
    } else {
      audit$attempts[i] <- 1L
      fit <- withCallingHandlers(tryCatch(
        .fit_nls_initialization(formula, data.frame(measures = y, times = t), seed),
        error = function(e) { error_message <<- conditionMessage(e); NULL }),
        warning = function(w) { warning_messages <<- c(warning_messages, conditionMessage(w)); invokeRestart("muffleWarning") })
    }
    if (!is.null(fit)) {
      estimate <- as.numeric(stats::coef(fit))
      fitted_values <- suppressWarnings(evaluate(t, estimate))
      converged <- isTRUE(fit$convInfo$isConv)
      if (converged && length(estimate) == length(seed) &&
          all(is.finite(estimate)) && all(is.finite(fitted_values))) {
        coefficients <- estimate
        audit$method[i] <- "nls"
        audit$nls_converged[i] <- TRUE
      } else {
        error_message <- if (!converged) "NLS did not converge" else "NLS returned nonfinite coefficients or fitted values"
      }
    }
    audit$warnings[i] <- paste(unique(warning_messages), collapse = "; ")
    if (!audit$nls_converged[i]) {
      audit$method[i] <- "heuristic_fallback"
      audit$fallback_reason[i] <- error_message
    }
    fitted_values <- suppressWarnings(evaluate(times[complete], coefficients))
    if (any(!is.finite(coefficients)) || any(!is.finite(fitted_values))) {
      audit$method[i] <- "unavailable"
      audit$fallback_reason[i] <- paste(audit$fallback_reason[i], "Heuristic coefficients or fitted values are nonfinite", sep = "; ")
      next
    }
    factors[i, ] <- coefficients
    residual_matrix[i, complete] <- measures[complete] - fitted_values
    residual_variances[i] <- mean(residual_matrix[i, complete]^2)
  }
  result <- as.data.frame(factors)
  attr(result, "mean_res_var") <- if (any(is.finite(residual_variances))) {
    stats::median(residual_variances[is.finite(residual_variances)])
  } else NA_real_
  attr(result, "residuals") <- residual_matrix
  attr(result, "initialization") <- list(method = "per_person_nls", algorithm = "default",
                                         control = .nls_initialization_control(), rows = audit)
  result
}

# Besides the growth factors, getUNI.GF() attaches two attributes: "mean_res_var" (median over persons of
# the per-person residual variance) and "residuals" (persons x occasions matrix of observed minus fitted
# values from the individual curve fits). The residual matrix feeds .residual_correlation_starts().

getUNI.GF <- function(dat_traj, dat_time, nT, curveFun){
  # Calculate individual growth factors eta0 (initial values) and eta1 (slopes) for linear trajectories.
  if (curveFun %in% c("linear", "LIN")){
    # Apply the linear model to each row
    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- ind[1:nT]
      times <- ind[-(1:nT)]
      # A person with fewer than two complete occasions cannot support a linear fit
      if (sum(is.finite(measures) & is.finite(times)) < 2){
        return(c(eta0 = NA_real_, eta1 = NA_real_, res_var = NA_real_, rep(NA_real_, nT)))
      }
      ind_model <- lm(measures ~ times)
      fitted <- unname(coef(ind_model)[1] + coef(ind_model)[2] * times)
      c(eta0 = unname(coef(ind_model)[1]), eta1 = unname(coef(ind_model)[2]),
        res_var = mean(residuals(ind_model)^2), unname(measures - fitted))
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2])
    attr(growth_factor, "mean_res_var") <- median(coef_mat[, 3], na.rm = TRUE)
    attr(growth_factor, "residuals") <- unname(coef_mat[, -(1:3), drop = FALSE])
  }
  # Calculate individual growth factors eta0 (initial values), eta1 (linear slopes), and eta2 (quadratic slopes)
  # for quadratic trajectories.
  else if (curveFun %in% c("quadratic", "QUAD")){
    # Apply the quadratic model to each row
    coef_mat <- t(apply(cbind(dat_traj, dat_time), 1, function(ind){
      measures <- ind[1:nT]
      times <- ind[-(1:nT)]
      # A person with fewer than three complete occasions cannot support a quadratic fit
      if (sum(is.finite(measures) & is.finite(times)) < 3){
        return(c(eta0 = NA_real_, eta1 = NA_real_, eta2 = NA_real_, res_var = NA_real_, rep(NA_real_, nT)))
      }
      times_sq <- times^2
      ind_model <- lm(measures ~ times + times_sq)
      fitted <- unname(coef(ind_model)[1] + coef(ind_model)[2] * times + coef(ind_model)[3] * times_sq)
      c(eta0 = unname(coef(ind_model)[1]), eta1 = unname(coef(ind_model)[2]),
        eta2 = unname(coef(ind_model)[3]), res_var = mean(residuals(ind_model)^2), unname(measures - fitted))
    }))
    growth_factor <- data.frame(eta0 = coef_mat[, 1], eta1 = coef_mat[, 2], eta2 = coef_mat[, 3])
    attr(growth_factor, "mean_res_var") <- median(coef_mat[, 4], na.rm = TRUE)
    attr(growth_factor, "residuals") <- unname(coef_mat[, -(1:4), drop = FALSE])
  }
  # Nonlinear growth factors use the same per-person NLS policy in every model
  # family. The retained audit distinguishes converged NLS fits from heuristic
  # fallbacks; no grid-search values are silently substituted.
  else if (curveFun %in% c("negative exponential", "EXP")){
    growth_factor <- .nls_curve_growth_factors(dat_traj, dat_time, nT, "EXP")
  }
  else if (curveFun %in% c("Jenss-Bayley", "JB")){
    growth_factor <- .nls_curve_growth_factors(dat_traj, dat_time, nT, "JB")
  }
  else if (curveFun %in% c("bilinear spline", "BLS")){
    growth_factor <- .nls_curve_growth_factors(dat_traj, dat_time, nT, "BLS")
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
    # Residuals of each person's own fit: an intercept and a scale for the change along the average
    # interval-specific slopes. Deviations from the average curve itself would also contain the
    # between-person differences in intercept and slope, and overstate the residual variance.
    Y <- as.matrix(dat_traj); Tm <- as.matrix(dat_time)
    fitted_vals <- matrix(NA_real_, nrow = nrow(Y), ncol = nT)
    for (i in seq_len(nrow(Y))){
      basis <- c(0, cumsum(slp * diff(Tm[i, ])))
      ok <- is.finite(Y[i, ]) & is.finite(basis)
      if (sum(ok) > 2L && stats::var(basis[ok]) > 0){
        fitted_vals[i, ok] <- Y[i, ok] - stats::lm.fit(cbind(1, basis[ok]), Y[i, ok])$residuals
      }
    }
    resid <- Y - fitted_vals
    n_fit <- rowSums(is.finite(resid))
    res_var_per_person <- ifelse(n_fit > 2L, rowSums(resid^2, na.rm = TRUE) / pmax(n_fit - 2L, 1L), NA_real_)
    attr(growth_factor, "mean_res_var") <- median(res_var_per_person, na.rm = TRUE)
    attr(growth_factor, "residuals") <- unname(resid)
  }
  return(growth_factor)
}

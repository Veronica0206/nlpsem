# Independent plotting targets from the RAM expectation, including stored fits.
.ram_curve_at <- function(model, times, outcome = "M", group = NULL, growth_tvc = FALSE){
  if (!is.null(group)){
    parent_data <- model@data
    model <- model@submodels[[group]]
    if (!length(model@data)) model <- OpenMx::mxModel(model, parent_data)
  }
  if (growth_tvc){
    labels <- names(OpenMx::omxGetParameters(model))
    model <- OpenMx::omxSetParameters(model, labels = labels[endsWith(labels, "kappa")], values = 0)
  }
  vapply(times, function(tt){
    at_time <- model
    observed <- at_time@data@observed
    columns <- grep("^T[0-9]+$", names(observed), value = TRUE)
    observed[columns] <- lapply(observed[columns], function(x) rep(tt, length(x)))
    at_time <- OpenMx::mxModel(at_time, OpenMx::mxData(observed, type = "raw"))
    as.numeric(OpenMx::mxGetExpected(at_time, "means", defvar.row = 1)[1, paste0(outcome, "1")])
  }, numeric(1))
}

.ram_parameter_jacobian <- function(model, expected){
  parameters <- OpenMx::omxGetParameters(model)
  target_length <- length(expected(model))
  matrix(vapply(seq_along(parameters), function(j){
    h <- 1e-7 * max(1, abs(parameters[j]))
    plus <- OpenMx::omxSetParameters(model, names(parameters)[j], values = parameters[j] + h)
    minus <- OpenMx::omxSetParameters(model, names(parameters)[j], values = parameters[j] - h)
    (expected(plus) - expected(minus)) / (2 * h)
  }, numeric(target_length)), nrow = target_length)
}

.expect_bls_ram_plot <- function(out, fit, outcome = "M", group = NULL, growth_tvc = FALSE){
  expected <- function(m) .ram_curve_at(m, out$time, outcome, group, growth_tvc)
  expect_equal(out$Y_hat, expected(fit), tolerance = 1e-10)
  jac <- .ram_parameter_jacobian(fit, expected)
  covariance <- stats::vcov(fit)
  expected_se <- sqrt(diag(jac %*% covariance %*% t(jac)))
  expect_equal((out$Y_hat_u - out$Y_hat_l) / (2 * qnorm(.975)), expected_se, tolerance = 3e-4)
}


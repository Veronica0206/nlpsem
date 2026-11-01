# BLS targets come from the fitted RAM construction, not a reconstructed spline.
.bls_plot_case <- function(random = TRUE, tic = TRUE, contrast = .1, knot = 1){
  covariance <- diag(.01, 4)
  starts <- list(Y_starts = list(means = c(50, 1.5, -.5, 5),
    covMatrix = covariance[lower.tri(covariance, diag = TRUE)], residuals = 1),
    TIC_starts = list(means = 2, covMatrix = 1), beta = matrix(c(0, 0, contrast, knot), ncol = 1))
  times <- c(0, 4, 5, 6, 7, 10)
  dat <- data.frame(ID = 1:8, x = seq(1, 3, length.out = 8), grp = rep(1:2, each = 4))
  for (j in seq_along(times)){
    dat[[paste0("T", j)]] <- times[j]
    dat[[paste0("Y", j)]] <- 45 + j + seq_len(8) / 10
    dat[[paste0("X", j)]] <- j + seq_len(8) / 20
  }
  model <- nlpsem:::getLGCM.mxModel(dat, "T", "Y", "BLS", random, seq_along(times), if (tic) "x" else NULL, starts)
  list(model = model, starts = starts, dat = dat, times = times)
}

.bls_tvc_case <- function(decompose, tic){
  case <- .bls_plot_case()
  st <- case$starts
  st$TVC_starts <- if (decompose == 0) list(means = 1:6, covMatrix = rep(1, 6)) else
    list(means = c(2, .3), covMatrix = c(1, 0, .1), residuals = 1, rel_rate = c(1, 2, .5, 1, 1))
  st$BL_starts <- if (decompose == 0) st$TIC_starts else
    list(means = if (tic) c(2, 2) else 2, covMatrix = if (tic) diag(2) else 1)
  if (decompose > 0) st$beta <- if (tic) cbind(c(.1, .02, .03, .1), c(0, 0, .1, 1)) else matrix(c(0, 0, .1, 1), ncol = 1)
  st$kappa <- 0; st$XY_cov <- .001
  case$model <- nlpsem:::getTVC.mxModel(case$dat, "T", "Y", "BLS", TRUE, 1:6, "LGCM", "X", decompose, if (tic) "x" else NULL, st)
  case
}

.bls_ram <- function(model, group = NULL){
  if (!is.null(group)){
    observed <- model@data
    model <- model@submodels[[group]]
    if (!length(model@data)) model <- OpenMx::mxModel(model, observed)
  }
  labels <- names(OpenMx::omxGetParameters(model))
  kappa <- labels[endsWith(labels, "kappa")]
  if (length(kappa)) model <- OpenMx::omxSetParameters(model, kappa, values = 0)
  as.numeric(OpenMx::mxGetExpected(model, "means", defvar.row = 1)[1, paste0("Y", 1:6)])
}

.expect_bls_mean <- function(model, times, prefix = "Y", family = "LGCM"){
  comp <- nlpsem:::.plot_gf_components(model, prefix, prefix, family, "BLS")
  comp$cov <- matrix(0, length(comp$est), length(comp$est))
  out <- nlpsem:::.plot_curve_delta(comp, times, "BLS")
  expect_equal(out$est, .bls_ram(model), tolerance = 1e-10)
  algebra <- OpenMx::mxModel(model, OpenMx::mxAlgebraFromString(comp$expression, name = "plot_components"))
  expect_equal(c(OpenMx::mxEvalByName("plot_components", algebra, compute = TRUE)), comp$est, tolerance = 1e-10)
  expect_gt(min(eigen(OpenMx::mxGetExpected(model, "covariance", defvar.row = 1), symmetric = TRUE)$values), 0)
  invisible(out)
}

test_that("BLS plots follow RAM through and beyond shifted knots, including controls", {
  for (random in c(FALSE, TRUE)) for (tic in c(FALSE, TRUE)){
    for (contrast in c(0, .1)) for (knot in c(-1, 0, 1)){
      case <- .bls_plot_case(random, tic, contrast, knot)
      .expect_bls_mean(case$model, case$times)
    }
  }
  case <- .bls_plot_case()
  expected <- .bls_ram(case$model)
  expect_equal(expected[c(3, 4, 6)], c(46.9997, 49.1999993500004, 53.99999995), tolerance = 1e-10)
})

test_that("BLS class curves retain their own expansion and covariate means", {
  for (container in c("multigroup", "mixture")){
    case <- .bls_plot_case()
    second <- case$starts
    second$Y_starts$means <- c(60, 2, -.3, 6)
    second$TIC_starts$means <- -1
    second$beta <- matrix(c(.2, .05, -.1, 1.5), ncol = 1)
    classes <- nlpsem:::getsub.LGCM(case$dat, 2, "T", "Y", "BLS", TRUE, 1:6, "x",
      list(case$starts, second), container, if (container == "multigroup") "grp" else NULL)
    for (k in 1:2){
      model <- OpenMx::mxModel(classes[[k]], OpenMx::mxData(case$dat, type = "raw"))
      .expect_bls_mean(model, case$times, paste0("c", k, "Y"))
    }
  }
})

test_that("BLS TVC components follow RAM in all decompositions with and without TICs", {
  for (decompose in 0:3) for (tic in c(FALSE, TRUE)){
    case <- .bls_tvc_case(decompose, tic)
    .expect_bls_mean(case$model, case$times, family = "TVC")
  }
})

test_that("BLS covariance propagation follows numerical RAM derivatives in parent models", {
  skip_on_cran()   # numerical RAM Jacobians over seven containers take about 20 s
  for (container in c("single", "multigroup", "mixture", paste0("TVC", 0:3))){
    case <- if (startsWith(container, "TVC")) .bls_tvc_case(as.integer(sub("TVC", "", container)), TRUE) else .bls_plot_case()
    grouped <- container %in% c("multigroup", "mixture")
    family <- if (startsWith(container, "TVC")) "TVC" else "LGCM"
    if (grouped){
      second <- case$starts; second$TIC_starts$means <- -1
      classes <- nlpsem:::getsub.LGCM(case$dat, 2, "T", "Y", "BLS", TRUE, 1:6, "x",
        list(case$starts, second), container, if (container == "multigroup") "grp" else NULL)
      model <- OpenMx::mxModel("Classes", classes)
      if (container == "mixture") model <- OpenMx::mxModel(model, OpenMx::mxData(case$dat, type = "raw"))
    } else model <- case$model
    parameters <- OpenMx::omxGetParameters(model)
    # Specified covariance isolates propagation, not estimated SEs or interval coverage.
    a <- diag(length(parameters))
    a[1, 2] <- -.3; a[1, nrow(a)] <- .2
    pcov <- crossprod(a) * .001
    dimnames(pcov) <- list(names(parameters), names(parameters))
    local_mocked_bindings(mxSE = function(expression, model, ...){
      OpenMx::mxSE(expression, model, cov = pcov, ...)
    }, .package = "nlpsem")
    for (k in if (grouped) 1:2 else 1){
      group <- if (grouped) k else NULL
      expected <- function(m) .bls_ram(m, group)
      jac <- .ram_parameter_jacobian(model, expected)
      eval_model <- if (grouped) model@submodels[[k]] else model
      prefix <- if (grouped) paste0("c", k, "Y") else "Y"
      se_prefix <- if (grouped) paste0("Class", k, ".", prefix) else prefix
      out <- nlpsem:::.trajectory_plot_data(model, eval_model, prefix, se_prefix, family, "BLS", "LGCM", case$times)
      expect_equal(out$Y_hat, expected(model), tolerance = 1e-10)
      expect_equal((out$Y_hat_u - out$Y_hat_l) / (2 * qnorm(.975)),
        sqrt(diag(jac %*% pcov %*% t(jac))), tolerance = 3e-4)
      # The intrinsic expression is smooth at both expansion and adjusted knots.
      expect_true(all(is.finite(out$Y_hat_u)))
      if (!length(eval_model@data)) eval_model <- OpenMx::mxModel(eval_model, model@data)
      dist <- nlpsem:::.plot_gf_distribution(model, eval_model, prefix, se_prefix, family, "BLS")
      change <- nlpsem:::.plot_curve_delta(dist, case$times, "BLS", "change", baseline_time = 5)
      target <- expected(model) - expected(model)[3]
      change_jac <- sweep(jac, 2, jac[3, ], "-")
      expect_equal(change$est, target, tolerance = 1e-10)
      expect_equal((change$upper - change$lower) / (2 * qnorm(.975)),
        sqrt(diag(change_jac %*% pcov %*% t(change_jac))), tolerance = 3e-4)
    }
  }
})

test_that("intrinsic multivariate BLS curves retain each smoothed process", {
  case <- .bls_plot_case()
  covariance <- diag(.01, 8)
  starts <- list(means = list(c(50, 1.5, -.5, 5), c(10, .3, -.1, 7)),
    covMatrix = covariance[lower.tri(covariance, diag = TRUE)], residuals = diag(2))
  model <- nlpsem:::getMGM.mxModel(case$dat, c("T", "T"), c("Y", "X"), "BLS", TRUE,
    list(1:6, 1:6), "LGCM", starts)
  expected <- OpenMx::mxGetExpected(model, "means", defvar.row = 1)
  for (prefix in c("Y", "X")){
    comp <- nlpsem:::.plot_gf_components(model, prefix, prefix, "MGM", "BLS")
    expect_true(comp$intrinsic_bls)
    comp$cov <- matrix(0, length(comp$est), length(comp$est))
    out <- nlpsem:::.plot_curve_delta(comp, case$times, "BLS")
    expect_equal(out$est, as.numeric(expected[1, paste0(prefix, 1:6)]), tolerance = 1e-10)
  }
})


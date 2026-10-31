# The RAM expectation is the independent target, including nonzero covariate
# means and effects on both amplitude and rate. No optimization is required.
.lgcm_plot_case <- function(curve, random = TRUE, tic = TRUE, shifted = TRUE){
  means <- if (curve == "EXP") c(50, 10, .2) else c(50, 2, -10, -.2)
  n <- length(means)
  covariance <- diag(.01, n)
  beta <- matrix(seq(.05, by = .03, length.out = n), ncol = 1)
  beta[n, 1] <- if (shifted) .1 else 0
  starts <- list(Y_starts = list(means = means, covMatrix = covariance[lower.tri(covariance, diag = TRUE)],
                                 residuals = 1),
                 TIC_starts = list(means = 2, covMatrix = 1), beta = beta)
  dat <- data.frame(ID = 1:8, x = seq(1, 3, length.out = 8), grp = rep(1:2, each = 4))
  times <- c(0, 5, 10, 15)
  for (j in 1:4){
    dat[[paste0("T", j)]] <- times[j]
    dat[[paste0("Y", j)]] <- 50 + j + seq_len(8) / 10
    dat[[paste0("X", j)]] <- j + seq_len(8) / 20
  }
  model <- nlpsem:::getLGCM.mxModel(dat, "T", "Y", curve, random, 1:4, if (tic) "x" else NULL, starts)
  list(model = model, starts = starts, dat = dat, times = times)
}

.lgcm_ram_mean <- function(model){
  as.numeric(OpenMx::mxGetExpected(model, "means", defvar.row = 1)[1, paste0("Y", 1:4)])
}

.expect_lgcm_plot_mean <- function(model, curve, times, prefix = "Y", family = "LGCM"){
  components <- nlpsem:::.plot_gf_components(model, prefix, prefix, family, curve)
  # Zero covariance extracts the point mean only; it is not fitted uncertainty.
  components$cov <- matrix(0, length(components$est), length(components$est))
  out <- nlpsem:::.plot_curve_delta(components, times, curve)
  expect_equal(out$est, .lgcm_ram_mean(model), tolerance = 1e-10)
  expr <- OpenMx::mxModel(model, OpenMx::mxAlgebraFromString(components$expression, name = "plot_components"))
  expect_equal(c(OpenMx::mxEvalByName("plot_components", expr, compute = TRUE)), components$est, tolerance = 1e-10)
  expect_gt(min(eigen(OpenMx::mxGetExpected(model, "covariance", defvar.row = 1), symmetric = TRUE)$values), 0)
}

test_that("LGCM plots follow RAM with shifted amplitude/rate and retain fixed-rate controls", {
  for (curve in c("EXP", "JB")) for (random in c(FALSE, TRUE)) for (tic in c(FALSE, TRUE)){
    for (shifted in c(FALSE, TRUE)){
      case <- .lgcm_plot_case(curve, random, tic, shifted)
      .expect_lgcm_plot_mean(case$model, curve, case$times)
    }
  }
})

test_that("class-specific LGCM means use their own expansion quantities and covariate shifts", {
  for (curve in c("EXP", "JB")) for (container in c("multigroup", "mixture")){
    case <- .lgcm_plot_case(curve)
    second <- case$starts
    second$TIC_starts$means <- 3
    second$Y_starts$means[2] <- second$Y_starts$means[2] * 1.5
    classes <- nlpsem:::getsub.LGCM(case$dat, 2, "T", "Y", curve, TRUE, 1:4, "x",
      list(case$starts, second), container, if (container == "multigroup") "grp" else NULL)
    for (k in 1:2){
      model <- OpenMx::mxModel(classes[[k]], OpenMx::mxData(case$dat, type = "raw"))
      .expect_lgcm_plot_mean(model, curve, case$times, paste0("c", k, "Y"))
    }
  }
})

test_that("direct and decomposed TVC growth curves retain covariate-induced rate means", {
  for (curve in c("EXP", "JB")) for (decompose in 0:3) for (tic in c(FALSE, TRUE)){
    case <- .lgcm_plot_case(curve)
    starts <- case$starts
    starts$TVC_starts <- if (decompose == 0) list(means = 1:4, covMatrix = rep(1, 4)) else
      list(means = c(1, .3), covMatrix = c(1, 0, .1), residuals = 1, rel_rate = c(1, 2, .5))
    starts$BL_starts <- list(means = if (tic) c(2, 1) else 1, covMatrix = if (tic) diag(2) else 1)
    if (decompose == 0) starts$BL_starts <- case$starts$TIC_starts
    if (decompose != 0) starts$beta <- matrix(seq(.03, by = .01,
      length.out = length(starts$Y_starts$means) * (1 + tic)), ncol = 1 + tic)
    starts$kappa <- .5; starts$XY_cov <- .001
    model <- nlpsem:::getTVC.mxModel(case$dat, "T", "Y", curve, TRUE, 1:4, "LGCM", "X", decompose,
                                    if (tic) "x" else NULL, starts)
    # Isolate the documented growth component, leaving all growth-factor paths.
    model <- OpenMx::omxSetParameters(model, labels = "kappa", values = 0)
    .expect_lgcm_plot_mean(model, curve, case$times, family = "TVC")
  }
})

test_that("LGCM plot uncertainty agrees with numerical RAM derivatives including parent covariance", {
  for (curve in c("EXP", "JB")) for (group in c(FALSE, TRUE)){
    case <- .lgcm_plot_case(curve)
    if (group){
      classes <- nlpsem:::getsub.LGCM(case$dat, 2, "T", "Y", curve, TRUE, 1:4, "x",
        rep(list(case$starts), 2), "multigroup", "grp")
      model <- OpenMx::mxModel("Groups", classes, OpenMx::mxFitFunctionMultigroup(c("Class1", "Class2")))
    } else model <- case$model
    parameters <- OpenMx::omxGetParameters(model)
    # A specified, correlated parameter covariance tests propagation. It is not
    # an estimated covariance or an interval-coverage claim for this small fixture.
    a <- diag(length(parameters)); a[1, 2] <- -.3
    pcov <- crossprod(a) * .001
    dimnames(pcov) <- list(names(parameters), names(parameters))
    local_mocked_bindings(mxSE = function(expression, model, ...){
      OpenMx::mxSE(expression, model, cov = pcov, ...)
    }, .package = "nlpsem")
    for (k in if (group) 1:2 else 1){
      eval_model <- if (group) model@submodels[[k]] else model
      pfx <- if (group) paste0("c", k, "Y") else "Y"
      se_pfx <- if (group) paste0("Class", k, ".", pfx) else pfx
      out <- nlpsem:::.trajectory_plot_data(model, eval_model, pfx, se_pfx, "LGCM", curve, "LGCM", case$times)
      expected <- function(m) .lgcm_ram_mean(if (group) m@submodels[[k]] else m)
      jac <- vapply(seq_along(parameters), function(j){
        step <- 1e-5 * max(1, abs(parameters[j]))
        plus <- OpenMx::omxSetParameters(model, labels = names(parameters)[j], values = parameters[j] + step)
        minus <- OpenMx::omxSetParameters(model, labels = names(parameters)[j], values = parameters[j] - step)
        (expected(plus) - expected(minus)) / (2 * step)
      }, numeric(4))
      expect_equal(out$Y_hat, expected(model), tolerance = 1e-10)
      actual_se <- (out$Y_hat_u - out$Y_hat_l) / (2 * qnorm(.975))
      expect_equal(actual_se, sqrt(diag(jac %*% pcov %*% t(jac))), tolerance = 3e-4)
    }
  }
})

test_that("intrinsic curve status, rate and change gradients retain all expansion quantities", {
  for (curve in c("EXP", "JB")){
    p <- if (curve == "EXP") c(50, 11, .2, 10, .2) else c(50, 2.2, -9, .2, -10, -.2)
    f <- if (curve == "EXP") function(p, t) p[1] + p[2]*(1-exp(-p[5]*t)) + p[4]*p[3]*t*exp(-p[5]*t) else
      function(p, t) p[1] + p[2]*t + p[3]*(exp(p[6]*t)-1) + p[5]*p[4]*t*exp(p[6]*t)
    a <- diag(length(p)); a[1, 3] <- -.3
    covariance <- crossprod(a) * .01
    for (quantity in c("status", "rate", "change")){
      fn <- switch(quantity, status = f, change = function(p,t) f(p,t)-f(p,2),
                   rate = function(p,t) (f(p,t+1e-4)-f(p,t-1e-4))/2e-4)
      tt <- c(1, 5, 10)
      jac <- vapply(seq_along(p), function(j){
        plus <- minus <- p; plus[j] <- p[j]+1e-4; minus[j] <- p[j]-1e-4
        (fn(plus,tt)-fn(minus,tt))/2e-4
      }, numeric(length(tt)))
      out <- nlpsem:::.plot_curve_delta(list(est=p,cov=covariance,intrinsic_rate=TRUE),tt,curve,quantity,baseline_time=2)
      expect_equal(out$est, fn(p,tt), tolerance=1e-6)
      expect_equal((out$upper-out$lower)/(2*qnorm(.975)),sqrt(diag(jac%*%covariance%*%t(jac))),tolerance=1e-6)
    }
  }
})

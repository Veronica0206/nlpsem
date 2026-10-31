# These fixtures have exact, independently specified sample moments. Per-person
# measurement errors are orthogonal to intercept/time, so person-level linear
# regression recovers the stated factors without requiring optimization.
.med_cov_fixture <- function(){
  n <- 128L
  a <- 2 * pi * seq_len(n) / n
  z <- do.call(cbind, lapply(1:5, function(k) cbind(sin(k * a), cos(k * a)))) *
    sqrt(2 * (n - 1) / n)
  X0 <- 10 * z[, 1]; X1 <- z[, 6]
  M0 <- 2 * X0 + z[, 2]; M1 <- 0.5 * X1 + z[, 3]
  Y0 <- 3 * M0 - 6 * X0 + z[, 4]
  Y1 <- 2 * M1 - X1 + z[, 5]
  make_process <- function(prefix, intercept, slope, noise){
    y <- outer(intercept, rep(1, 5)) + outer(slope, 0:4) +
      outer(noise, c(1, -2, 1, 0, 0))
    setNames(as.data.frame(y), paste0(prefix, 1:5))
  }
  dat <- cbind(data.frame(ID = seq_len(n), X = X0),
               make_process("Y", Y0, Y1, z[, 7]),
               make_process("M", M0, M1, z[, 8]),
               make_process("P", X0, X1, z[, 9]),
               setNames(as.data.frame(matrix(rep(0:4, each = n), n)), paste0("T", 1:5)))
  list(dat = dat, factors = data.frame(X0, X1, M0, M1, Y0, Y1), z = z)
}

.med_cov_args <- function(dat, x_type, curve = "LIN"){
  long <- x_type == "longitudinal"
  list(dat = dat, t_var = rep("T", if (long) 3 else 2), y_var = "Y", m_var = "M",
       x_type = x_type, x_var = if (long) "P" else "X", curveFun = curve,
       records = rep(list(1:5), if (long) 3 else 2), res_scale = NULL, res_cor = NULL)
}

test_that("opposing mediation paths give correct and admissible linear starts", {
  fixture <- .med_cov_fixture()
  for (x_type in c("baseline", "longitudinal")){
    args <- .med_cov_args(fixture$dat, x_type)
    starts <- do.call(nlpsem:::getMED.initial, args)
    expected_m <- diag(c(1, if (x_type == "baseline") 1.25 else 1))
    expected_y <- diag(c(1, if (x_type == "baseline") 1.8 else 1))
    expect_equal(unname(starts$M$covMatrix), expected_m, tolerance = 1e-10)
    expect_equal(unname(starts$Y$covMatrix), expected_y, tolerance = 1e-10)
    expect_equal(as.numeric(starts$Y$beta_x)[1], -6, tolerance = 1e-10)
    expect_equal(starts$Y$beta_m[1, 1], 3, tolerance = 1e-10)
    bargs <- args; bargs$res_scale <- NULL; bargs$starts <- starts
    model <- do.call(nlpsem:::getMED.mxModel, bargs)
    expected <- OpenMx::mxGetExpected(model, "covariance", defvar.row = 1)
    expect_gt(min(eigen(expected, symmetric = TRUE)$values), 0)
    # Evaluate likelihood only: these tests do not replace or retry the optimizer.
    evaluated <- OpenMx::mxRun(OpenMx::mxModel(model,
      OpenMx::mxComputeOnce("fitfunction", "fit")), silent = TRUE)
    expect_true(is.finite(evaluated$output$fit))
  }
})

test_that("mediation covariance uses aligned regression residuals with missing predictors", {
  fixture <- .med_cov_fixture()
  missing <- abs(fixture$factors$X0) < quantile(abs(fixture$factors$X0), 0.75)
  for (x_type in c("baseline", "longitudinal")){
    dat <- fixture$dat
    f <- fixture$factors
    f$X0[missing] <- NA_real_
    if (x_type == "baseline"){
      dat$X[missing] <- NA_real_
      m_res <- residuals(lm(cbind(M0, M1) ~ X0, data = f, na.action = na.exclude))
      y_res <- cbind(residuals(lm(Y0 ~ X0 + M0, data = f, na.action = na.exclude)),
                     residuals(lm(Y1 ~ X0 + M0 + M1, data = f, na.action = na.exclude)))
    } else{
      dat[missing, paste0("P", 1:5)] <- NA_real_
      f$X1[missing] <- NA_real_
      m_res <- cbind(residuals(lm(M0 ~ X0, data = f, na.action = na.exclude)),
                     residuals(lm(M1 ~ X0 + X1, data = f, na.action = na.exclude)))
      y_res <- cbind(residuals(lm(Y0 ~ M0 + X0, data = f, na.action = na.exclude)),
                     residuals(lm(Y1 ~ M0 + M1 + X0 + X1, data = f, na.action = na.exclude)))
    }
    starts <- do.call(nlpsem:::getMED.initial, .med_cov_args(dat, x_type))
    expect_equal(unname(starts$M$covMatrix), unname(cov(m_res, use = "complete.obs")), tolerance = 1e-10)
    expect_equal(unname(starts$Y$covMatrix), unname(cov(y_res, use = "complete.obs")), tolerance = 1e-10)
    expect_gt(min(eigen(starts$M$covMatrix, symmetric = TRUE)$values), 0)
    expect_gt(min(eigen(starts$Y$covMatrix, symmetric = TRUE)$values), 0)
  }
})

test_that("group and mixture mediation paths retain the corrected covariance starts", {
  dat <- .med_cov_fixture()$dat
  for (x_type in c("baseline", "longitudinal")){
    args <- .med_cov_args(dat, x_type)
    expected <- do.call(nlpsem:::getMED.initial, args)
    group_args <- c(args, list(nClass = 2L, grp_var = "group", sub_Model = "MED",
                              intrinsic = NULL, TVC = NULL, decompose = 0, growth_TIC = NULL))
    group_args$dat <- rbind(transform(dat, group = 1L), transform(dat, group = 2L))
    groups <- do.call(nlpsem:::getMGroup.initial, group_args)
    mix_args <- c(args, list(nClass = 1L, prop_starts = 1, sub_Model = "MED",
                            cluster_TIC = NULL, TVC = NULL, decompose = 0,
                            growth_TIC = NULL, intrinsic = NULL))
    mixture <- do.call(nlpsem:::getMIX.initial, mix_args)
    for (s in c(groups, mixture[seq_len(mix_args$nClass)])){
      expect_equal(s$M$covMatrix, expected$M$covMatrix, tolerance = 1e-10)
      expect_equal(s$Y$covMatrix, expected$Y$covMatrix, tolerance = 1e-10)
    }
  }
})

test_that("missing residual overlap has an explicit diagonal fallback", {
  outcomes <- cbind(c(1, 4, 7, NA), c(NA, NA, 3, 8))
  predictors <- cbind(0, c(NA, NA, 0, 0))
  expect_warning(covariance <- nlpsem:::.med_residual_covariance(
    outcomes, predictors, c(0, 0), matrix(0, 2, 2), list(1L, 1:2), "test residuals"),
    "Fewer than two shared residual rows")
  # Missing predictor 2 must not discard equation 1's first two residuals.
  expect_equal(unname(covariance), diag(c(9, 12.5)))
  expect_warning(fallback <- nlpsem:::.med_residual_covariance(
    cbind(1:3, NA_real_), matrix(0, 3, 1), c(0, 0), matrix(0, 2, 1),
    list(1L, 1L), "test unavailable residuals"), "using diagonal covariance starts")
  expect_equal(unname(fallback), diag(c(1, 1)))
})

test_that("BLS disturbance covariance is estimated before unchanged path damping", {
  fixture <- .med_cov_fixture()
  z <- fixture$z
  damping <- c(0.7, 0.7, 0.5)
  xm <- c(2, 0.5, 0.25)
  my <- c(3, 2, 1.5)
  for (x_type in c("baseline", "longitudinal")){
    if (x_type == "baseline"){
      x <- matrix(10 * z[, 1], ncol = 1)
      m <- x %*% t(xm) + z[, 2:4]
      y <- sweep(z[, 2:4], 2, my, "*") + z[, 5:7]
    } else{
      x <- sweep(z[, 1:3], 2, c(10, 2, 1), "*")
      m <- sweep(x, 2, xm, "*") + z[, 4:6]
      y <- sweep(z[, 4:6], 2, my, "*") + z[, 7:9]
    }
    # All structural regression disturbances have unit sample covariance;
    # damping the returned seeds must not add prediction error to that target.
    expected_m <- expected_y <- diag(3)
    # Isolate the structural regression stage. Existing NLS tests separately
    # exercise the unchanged person-level BLS estimation and knot profiling.
    raw_bls <- function(factors){
      gf <- data.frame(eta0 = factors[, 2] - 3 * factors[, 1],
                       eta1 = factors[, 1], eta2 = factors[, 3], gamma = 3)
      attr(gf, "mean_res_var") <- 1
      gf
    }
    supplied <- list(Y = raw_bls(y), M = raw_bls(m))
    if (x_type == "longitudinal") supplied$P <- raw_bls(x)
    local_mocked_bindings(getUNI.GF = function(dat_traj, ...){
      supplied[[substr(names(dat_traj)[1], 1, 1)]]
    }, .package = "nlpsem")
    starts <- do.call(nlpsem:::getMED.initial, .med_cov_args(fixture$dat, x_type, "BLS"))
    expect_equal(unname(starts$M$covMatrix), expected_m, tolerance = 1e-9)
    expect_equal(unname(starts$Y$covMatrix), expected_y, tolerance = 1e-9)
    expect_equal(diag(starts$Y$beta_m), damping * my, tolerance = 1e-9)
    if (x_type == "baseline"){
      expect_equal(as.numeric(starts$M$beta), damping * xm, tolerance = 1e-9)
    } else{
      expect_equal(diag(starts$M$beta), damping * xm, tolerance = 1e-9)
    }
    expect_gt(min(eigen(starts$M$covMatrix, symmetric = TRUE)$values), 0)
    expect_gt(min(eigen(starts$Y$covMatrix, symmetric = TRUE)$values), 0)
  }
})

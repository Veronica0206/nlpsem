.quad_truth_starts <- function(truth, prefix = "", tic = FALSE) {
  take <- function(nm) unname(truth[paste0(prefix, nm)])
  out <- list(Y_starts = list(
    means = take(paste0(if (tic) "alpha" else "mueta", 0:2)),
    covMatrix = take(c("psi00", "psi01", "psi02", "psi11", "psi12", "psi22")),
    residuals = take("residuals")))
  if (tic) {
    out$beta <- matrix(take(c("beta10", "beta11", "beta12")), 3, 1)
    out$TIC_starts <- list(means = take("mux1"), covMatrix = take("phi11"))
  }
  out
}

test_that("shifted quadratic truth matches the actual LCSM mean and covariance", {
  for (jitter in c(0, 0.25)) {
    dat <- getSimData(n = 30, curveFun = "QUAD", times = 5:14, jitter = jitter, seed = 1)
    truth <- attr(dat, "truth")
    lcsm <- attr(dat, "truth_lcsm")
    info <- attr(dat, "truth_info")
    p <- attr(dat, "settings")$params[[1]]
    expect_equal(truth[c("mueta0", "psi00")], c(mueta0 = 50, psi00 = 25))
    expect_equal(lcsm[c("mueta0", "psi00", "psi01", "psi02")],
                 c(mueta0 = 125, psi00 = 166.25, psi01 = 8.75, psi02 = 3.15))
    expect_identical(info$intercept_reference, "time_zero")
    expect_identical(info$intercept_time, 0)
    expect_equal(info$first_observation_time, 5)
    expect_identical(info$reference_records, 1:10)
    expect_true(info$lcsm$truth_available)
    expect_false(info$lcsm$same_parameter_definitions)
    model <- nlpsem:::getLCSM.mxModel(dat, "T", "Y", "QUAD", FALSE, 1:10, NULL,
                                     .quad_truth_starts(lcsm))
    # Compare the built recurrence with moments from the independent generating equation.
    for (i in c(1, 8)) {
      time <- unlist(dat[i, paste0("T", 1:10)], use.names = FALSE)
      L <- cbind(1, time, time^2)
      expected_mean <- drop(L %*% p$mean)
      expected_cov <- L %*% p$cov %*% t(L) + diag(p$residual, 10)
      expect_equal(as.numeric(OpenMx::mxGetExpected(model, "means", defvar.row = i)), expected_mean,
                   tolerance = 1e-10)
      expect_equal(unname(OpenMx::mxGetExpected(model, "covariance", defvar.row = i)), expected_cov,
                   tolerance = 1e-10)
    }
  }
})

test_that("quadratic TIC targets transform conditional means, effects and covariance together", {
  p <- list(mean = c(eta0 = 50, eta1 = 20, eta2 = -1),
            cov = matrix(c(25, 1.5, 0.45, 1.5, 1, 0.09, 0.45, 0.09, 0.09), 3), residual = 1,
            tic_mean = 2, tic_cov = matrix(1), beta = matrix(c(1, 2, 3), 3, 1))
  dat <- getSimData(n = 30, curveFun = "QUAD", times = 5:14, jitter = 0, params = p, seed = 2)
  lcsm <- attr(dat, "truth_lcsm")
  expect_equal(attr(dat, "truth")[c("alpha0", "beta10")], c(alpha0 = 48, beta10 = 1))
  expect_equal(lcsm[c("alpha0", "beta10", "mueta0", "psi00")],
               c(alpha0 = -47, beta10 = 86, mueta0 = 125, psi00 = 166.25))
  model <- nlpsem:::getLCSM.mxModel(dat, "T", "Y", "QUAD", FALSE, 1:10, "ex1",
                                   .quad_truth_starts(lcsm, tic = TRUE))
  L <- cbind(1, 5:14, (5:14)^2)
  total_cov <- p$cov + p$beta %*% p$tic_cov %*% t(p$beta)
  cov_yx <- L %*% p$beta %*% p$tic_cov
  target_cov <- rbind(cbind(L %*% total_cov %*% t(L) + diag(1, 10), cov_yx),
                      cbind(t(cov_yx), p$tic_cov))
  expect_equal(as.numeric(OpenMx::mxGetExpected(model, "means")), c(drop(L %*% p$mean), 2),
               tolerance = 1e-10)
  expect_equal(unname(OpenMx::mxGetExpected(model, "covariance")), unname(target_cov), tolerance = 1e-10)
})

test_that("quadratic class-specific targets retain class membership parameters", {
  p1 <- list(mean = c(eta0 = 50, eta1 = 20, eta2 = -1), cov = diag(c(25, 1, 0.09)), residual = 1)
  p2 <- p1; p2$mean <- c(eta0 = 60, eta1 = 10, eta2 = -0.5)
  dat <- getSimData(n = 60, curveFun = "QUAD", times = 5:14, jitter = 0,
                    params = list(p1, p2), prop = c(0.4, 0.6), seed = 3)
  lcsm <- attr(dat, "truth_lcsm")
  expect_equal(lcsm[c("c1mueta0", "c2mueta0", "p2")],
               c(c1mueta0 = 125, c2mueta0 = 97.5, p2 = log(1.5)))
  starts <- lapply(1:2, function(k) .quad_truth_starts(lcsm, paste0("c", k)))
  models <- nlpsem:::getsub.LCSM(dat, 2, "T", "Y", "QUAD", FALSE, 1:10, NULL,
                                starts, container = "multigroup", grp_var = "class")
  L <- cbind(1, 5:14, (5:14)^2)
  for (k in 1:2) {
    p <- list(p1, p2)[[k]]
    expect_equal(as.numeric(OpenMx::mxGetExpected(models[[k]], "means")), drop(L %*% p$mean),
                 tolerance = 1e-10)
    expect_equal(unname(OpenMx::mxGetExpected(models[[k]], "covariance")),
                 L %*% p$cov %*% t(L) + diag(1, 10), tolerance = 1e-10)
  }
  cluster <- list(mean = 0, cov = matrix(1), beta = rbind(c(0, 0), c(-0.5, 1)))
  clustered <- getSimData(n = 30, curveFun = "QUAD", times = 5:14, jitter = 0,
                          params = list(p1, p2), cluster = cluster, seed = 4)
  expect_identical(attr(clustered, "truth_lcsm")[c("beta20", "beta21")],
                    attr(clustered, "truth")[c("beta20", "beta21")])
})

test_that("quadratic multivariate targets retain cross-process covariance meaning", {
  dat <- getSimData(n = 50, sub_Model = "MGM", curveFun = "QUAD", times = 5:14,
                    jitter = 0, y_var = c("A", "B"), seed = 5)
  p <- attr(dat, "settings")$params[[1]]
  lcsm <- attr(dat, "truth_lcsm")
  expect_equal(lcsm[c("A_mueta0", "B_mueta0")], c(A_mueta0 = 125, B_mueta0 = 129))
  expected_baseline_cross_cov <- drop(c(1, 5, 25) %*% p$cov[1:3, 4:6] %*% c(1, 5, 25))
  expect_equal(unname(lcsm[["AB_psi00"]]), expected_baseline_cross_cov)
  starts <- nlpsem:::getMULTI.initial(dat, c("T", "T"), c("A", "B"), "QUAD", list(1:10, 1:10),
                                     res_scale = NULL, res_cor = 0.3)
  model <- nlpsem:::getMGM.mxModel(dat, c("T", "T"), c("A", "B"), "QUAD", FALSE,
                                  list(1:10, 1:10), "LCSM", starts)
  model_labels <- names(lcsm)
  model_labels[match(c("A_res", "B_res", "AB_res"), model_labels)] <- c("A_residuals", "B_residuals", "AB_RES")
  expect_setequal(model_labels, names(OpenMx::omxGetParameters(model)))
  model <- OpenMx::omxSetParameters(model, model_labels, values = unname(lcsm))
  L <- cbind(1, 5:14, (5:14)^2)
  loadings <- kronecker(diag(2), L)
  target_cov <- loadings %*% p$cov %*% t(loadings) + kronecker(p$residual, diag(10))
  expect_equal(as.numeric(OpenMx::mxGetExpected(model, "means")),
               drop(loadings %*% unlist(p$mean, use.names = FALSE)), tolerance = 1e-10)
  expect_equal(unname(OpenMx::mxGetExpected(model, "covariance")), target_cov, tolerance = 1e-10)
})

test_that("quadratic TVC targets transform both kinds of growth covariate effects", {
  for (decompose in 0:3) {
    p <- nlpsem:::.sim_default_params(list(sub_Model = "TVC", form = "QUAD", intrinsic = FALSE,
                                         nT = 10, decompose = decompose))
    p$tic_mean <- c(2, 3); p$tic_cov <- diag(2)
    p$beta <- matrix(c(1, 2, 3, -1, 1, 0.5), 3, 2)
    if (decompose != 0) p$cov_bl <- c(0, 0)
    dat <- getSimData(n = 30, sub_Model = "TVC", curveFun = "QUAD", decompose = decompose,
                      times = 5:14, jitter = 0, params = p, seed = 6)
    lcsm <- attr(dat, "truth_lcsm")
    b <- c(1, 5, 25)
    beta <- drop(b %*% p$beta)
    alpha0 <- sum(b * p$mean) - sum(beta * p$tic_mean)
    if (decompose != 0) {
      beta_tvc0 <- sum(b * p$beta_tvc)
      alpha0 <- alpha0 - beta_tvc0 * p$tvc$mean[["eta0"]]
      expect_equal(unname(lcsm[["betaTVC0"]]), beta_tvc0)
      expect_equal(lcsm[c("X_mueta0", "X_psi00")], attr(dat, "truth")[c("X_mueta0", "X_psi00")])
    }
    expect_equal(unname(lcsm[c("Y_alpha0", "betaTIC10", "betaTIC20", "Y_psi00")]),
                 c(alpha0, beta, drop(b %*% p$cov %*% b)))
    expect_equal(lcsm[["kappa"]], p$kappa)
    expect_match(attr(dat, "truth_info")$lcsm$qualification, "do not establish equivalence")
  }
})

test_that("truth metadata distinguishes supported origins and unavailable nonlinear targets", {
  zero <- getSimData(n = 10, curveFun = "QUAD", times = 0:9, seed = 7)
  expect_identical(attr(zero, "truth_lcsm"), attr(zero, "truth"))
  expect_true(attr(zero, "truth_info")$lcsm$same_parameter_definitions)
  nonp <- getSimData(n = 10, curveFun = "NonP", times = 5:14, seed = 8)
  expect_identical(attr(nonp, "truth_lcsm"), attr(nonp, "truth"))
  expect_identical(attr(nonp, "truth_info")$intercept_reference, "first_observation")
  expect_equal(attr(nonp, "truth_info")$intercept_time, 5)
  expect_true(attr(nonp, "truth_info")$lcsm$same_parameter_definitions)
  for (form in c("EXP", "JB")) for (intrinsic in c(FALSE, TRUE)) {
    dat <- getSimData(n = 10, curveFun = form, intrinsic = intrinsic, times = 5:14, seed = 9)
    expect_null(attr(dat, "truth_lcsm"))
    expect_false(attr(dat, "truth_info")$lcsm$truth_available)
    expect_false(attr(dat, "truth_info")$lcsm$same_parameter_definitions)
    expect_match(attr(dat, "truth_info")$lcsm$qualification, "No exact LCSM targets")
  }
  med <- getSimData(n = 10, sub_Model = "MED", curveFun = "BLS", x_type = "baseline", seed = 10)
  expect_identical(attr(med, "truth_info")$intercept_reference, "process_knot")
  expect_null(attr(med, "truth_info")$intercept_time)
  expect_false(attr(med, "truth_info")$lcsm$supported_curve)
  expect_null(attr(med, "truth_lcsm"))
  # Metadata is tied to the original full schedule, not a later selected/jittered baseline.
  expect_identical(attr(zero, "truth_info")$reference_records, 1:10)
  expect_match(attr(zero, "truth_info")$lcsm$qualification, "later selected baselines")
})

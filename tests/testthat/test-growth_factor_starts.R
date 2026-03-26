test_that("getUNI.GF returns finite anchored starts for EXP trajectories", {
  set.seed(123)
  times <- matrix(rep(0:4, each = 24), nrow = 24)
  eta0 <- rnorm(24, 10, 0.6)
  eta1 <- rnorm(24, 4, 0.5)
  gamma <- pmax(rnorm(24, 0.45, 0.08), 0.05)
  traj <- eta0 + eta1 * (1 - exp(-gamma * times)) + matrix(rnorm(24 * 5, sd = 0.1), nrow = 24)

  gf <- nlpsem:::getUNI.GF(as.data.frame(traj), as.data.frame(times), nT = 5, curveFun = "EXP")

  expect_true(all(is.finite(as.matrix(gf))))
  expect_gt(stats::sd(gf$gamma), 0)
})

test_that("getUNI.GF returns finite anchored starts for JB trajectories", {
  set.seed(456)
  times <- matrix(rep(seq(0, 2, length.out = 5), each = 24), nrow = 24)
  eta0 <- rnorm(24, 8, 0.5)
  eta1 <- rnorm(24, 1.1, 0.15)
  eta2 <- rnorm(24, 2.5, 0.3)
  gamma <- rnorm(24, -0.55, 0.08)
  traj <- eta0 + eta1 * times + eta2 * (exp(gamma * times) - 1) + matrix(rnorm(24 * 5, sd = 0.08), nrow = 24)

  gf <- nlpsem:::getUNI.GF(as.data.frame(traj), as.data.frame(times), nT = 5, curveFun = "JB")

  expect_true(all(is.finite(as.matrix(gf))))
  expect_gt(stats::sd(gf$gamma), 0)
})

test_that("JB starts avoid near-zero gamma collapse on package example data", {
  data("RMS_dat", package = "nlpsem")
  RMS_dat0 <- RMS_dat
  baseT <- RMS_dat0$T1
  for (j in 1:9) RMS_dat0[[paste0("T", j)]] <- RMS_dat0[[paste0("T", j)]] - baseT

  starts <- nlpsem:::getUNI.initial(
    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "JB",
    records = 1:9, growth_TIC = NULL, res_scale = NULL
  )

  expect_gt(abs(unname(starts$Y_starts$means[4])), 1e-3)
  expect_lt(abs(unname(starts$Y_starts$means[3])), 1000)
})

test_that("intrinsic nonlinear starts damp the random shape component", {
  data("RMS_dat", package = "nlpsem")
  RMS_dat0 <- RMS_dat
  baseT <- RMS_dat0$T1
  for (j in 1:9) RMS_dat0[[paste0("T", j)]] <- RMS_dat0[[paste0("T", j)]] - baseT
  RMS_datM <- RMS_dat0
  for (j in 1:9) RMS_datM[[paste0("T", j)]] <- RMS_datM[[paste0("T", j)]] / 12

  jb_fixed <- nlpsem:::getUNI.initial(
    dat = RMS_datM, t_var = "T", y_var = "M", curveFun = "JB",
    records = 1:9, growth_TIC = NULL, res_scale = NULL, intrinsic = FALSE
  )
  jb_intrinsic <- nlpsem:::getUNI.initial(
    dat = RMS_datM, t_var = "T", y_var = "M", curveFun = "JB",
    records = 1:9, growth_TIC = NULL, res_scale = NULL, intrinsic = TRUE
  )
  bgs_fixed <- nlpsem:::getUNI.initial(
    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS",
    records = 1:9, growth_TIC = NULL, res_scale = NULL, intrinsic = FALSE
  )
  bgs_intrinsic <- nlpsem:::getUNI.initial(
    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS",
    records = 1:9, growth_TIC = NULL, res_scale = NULL, intrinsic = TRUE
  )

  expect_lt(unname(jb_intrinsic$Y_starts$covMatrix[10]), unname(jb_fixed$Y_starts$covMatrix[10]))
  expect_lt(unname(bgs_intrinsic$Y_starts$covMatrix[10]), unname(bgs_fixed$Y_starts$covMatrix[10]))
  expect_equal(unname(bgs_intrinsic$Y_starts$covMatrix[4]), 0)
  expect_equal(unname(bgs_intrinsic$Y_starts$covMatrix[7]), 0)
  expect_equal(unname(bgs_intrinsic$Y_starts$covMatrix[9]), 0)
})

test_that("intrinsic TIC and TVC starts shrink covariate path starts", {
  data("RMS_dat", package = "nlpsem")
  RMS_dat0 <- RMS_dat
  baseT <- RMS_dat0$T1
  for (j in 1:9) RMS_dat0[[paste0("T", j)]] <- RMS_dat0[[paste0("T", j)]] - baseT
  RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
  BL_mean <- mean(RMS_dat0[, "R1"])
  BL_var <- var(RMS_dat0[, "R1"])
  for (j in 1:9) RMS_dat0[[paste0("Rs", j)]] <- (RMS_dat0[[paste0("R", j)]] - BL_mean)/sqrt(BL_var)

  bls_tic_fixed <- nlpsem:::getUNI.initial(
    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS",
    records = 1:9, growth_TIC = "ex1", res_scale = NULL, intrinsic = FALSE
  )
  bls_tic_intrinsic <- nlpsem:::getUNI.initial(
    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS",
    records = 1:9, growth_TIC = "ex1", res_scale = NULL, intrinsic = TRUE
  )
  bls_tvc_fixed <- nlpsem:::getTVC.initial(
    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS",
    records = 1:9, growth_TIC = "ex1", TVC = "Rs", decompose = 1,
    res_scale = NULL, intrinsic = FALSE
  )
  bls_tvc_intrinsic <- nlpsem:::getTVC.initial(
    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS",
    records = 1:9, growth_TIC = "ex1", TVC = "Rs", decompose = 1,
    res_scale = NULL, intrinsic = TRUE
  )

  expect_lt(abs(unname(bls_tic_intrinsic$beta[1, 1])), abs(unname(bls_tic_fixed$beta[1, 1])))
  expect_lt(abs(unname(bls_tvc_intrinsic$beta[1, 2])), abs(unname(bls_tvc_fixed$beta[1, 2])))
  expect_lt(abs(unname(bls_tvc_intrinsic$kappa)), abs(unname(bls_tvc_fixed$kappa)))
})

test_that("decompose 3 TVC kappa uses change-from-baseline states", {
  data("RMS_dat", package = "nlpsem")
  RMS_dat0 <- RMS_dat
  baseT <- RMS_dat0$T1
  for (j in 1:9) RMS_dat0[[paste0("T", j)]] <- RMS_dat0[[paste0("T", j)]] - baseT
  RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
  BL_mean <- mean(RMS_dat0[, "R1"])
  BL_var <- var(RMS_dat0[, "R1"])
  for (j in 1:9) RMS_dat0[[paste0("Rs", j)]] <- (RMS_dat0[[paste0("R", j)]] - BL_mean)/sqrt(BL_var)

  starts_d2 <- nlpsem:::getTVC.initial(
    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS",
    records = 1:9, growth_TIC = "ex1", TVC = "Rs", decompose = 2,
    res_scale = NULL, intrinsic = TRUE
  )
  starts_d3 <- nlpsem:::getTVC.initial(
    dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS",
    records = 1:9, growth_TIC = "ex1", TVC = "Rs", decompose = 3,
    res_scale = NULL, intrinsic = TRUE
  )

  expect_true(is.finite(starts_d3$kappa))
  expect_lt(abs(unname(starts_d3$kappa)), abs(unname(starts_d2$kappa)))
  expect_lt(abs(unname(starts_d3$kappa)), 1)
})

test_that("decompose 3 TVC starts are finite across supported functional forms", {
  data("RMS_dat", package = "nlpsem")
  RMS_dat0 <- RMS_dat
  baseT <- RMS_dat0$T1
  for (j in 1:9) RMS_dat0[[paste0("T", j)]] <- (RMS_dat0[[paste0("T", j)]] - baseT)/12
  RMS_dat0$ex1 <- as.numeric(scale(RMS_dat0$Approach_to_Learning))
  RMS_dat0$ex2 <- as.numeric(scale(RMS_dat0$Attention_focus))
  BL_mean <- mean(RMS_dat0[, "R1"])
  BL_var <- var(RMS_dat0[, "R1"])
  for (j in 1:9) RMS_dat0[[paste0("Rs", j)]] <- (RMS_dat0[[paste0("R", j)]] - BL_mean)/sqrt(BL_var)

  cases <- list(
    list(curveFun = "LIN", intrinsic = FALSE),
    list(curveFun = "QUAD", intrinsic = FALSE),
    list(curveFun = "EXP", intrinsic = FALSE),
    list(curveFun = "EXP", intrinsic = TRUE),
    list(curveFun = "JB", intrinsic = FALSE),
    list(curveFun = "JB", intrinsic = TRUE),
    list(curveFun = "BLS", intrinsic = FALSE),
    list(curveFun = "BLS", intrinsic = TRUE),
    list(curveFun = "NonP", intrinsic = FALSE)
  )

  kappas <- vapply(cases, function(case) {
    starts <- nlpsem:::getTVC.initial(
      dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = case$curveFun,
      records = 1:9, growth_TIC = c("ex1", "ex2"), TVC = "Rs", decompose = 3,
      res_scale = NULL, intrinsic = case$intrinsic
    )
    expect_true(is.finite(starts$kappa))
    unname(starts$kappa)
  }, numeric(1))

  expect_equal(length(kappas), length(cases))
})

test_that("getUNI.GF returns finite anchored starts for BLS trajectories", {
  set.seed(789)
  times <- matrix(rep(seq(0, 4, length.out = 6), each = 24), nrow = 24)
  eta0 <- rnorm(24, 12, 0.7)
  eta1 <- rnorm(24, 1.8, 0.2)
  eta2 <- rnorm(24, -0.9, 0.15)
  knot <- rnorm(24, 2.2, 0.25)
  traj <- eta0 + eta1 * pmin(times, knot) + eta2 * pmax(times - knot, 0) + matrix(rnorm(24 * 6, sd = 0.12), nrow = 24)

  gf <- nlpsem:::getUNI.GF(as.data.frame(traj), as.data.frame(times), nT = 6, curveFun = "BLS")

  expect_true(all(is.finite(as.matrix(gf))))
  expect_gt(stats::sd(gf$gamma), 0)
})

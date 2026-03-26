test_that("BLS mediation starts with longitudinal predictor are damped and PD", {
  data("RMS_dat", package = "nlpsem")
  RMS_dat0 <- RMS_dat
  baseT <- RMS_dat0$T1
  for (j in 1:9) RMS_dat0[[paste0("T", j)]] <- RMS_dat0[[paste0("T", j)]] - baseT

  reparam_bls <- function(gf){
    tmp_etar <- gf[, 1] + gf[, 2] * gf[, 4]
    tmp_eta1 <- gf[, 2]
    gf[, 1] <- tmp_eta1
    gf[, 2] <- tmp_etar
    names(gf)[1:2] <- c("eta1", "etar")
    gf
  }

  starts <- nlpsem:::getMED.initial(
    dat = RMS_dat0, t_var = rep("T", 3), y_var = "S", m_var = "M",
    x_type = "longitudinal", x_var = "R", curveFun = "BLS",
    records = list(2:9, 1:9, 1:9), res_scale = NULL, res_cor = NULL
  )

  X_growth_factor <- reparam_bls(nlpsem:::getUNI.GF(
    dat_traj = RMS_dat0[, paste0("R", 1:9)],
    dat_time = RMS_dat0[, paste0("T", 1:9)],
    nT = 9, curveFun = "BLS"
  ))
  M_growth_factor <- reparam_bls(nlpsem:::getUNI.GF(
    dat_traj = RMS_dat0[, paste0("M", 1:9)],
    dat_time = RMS_dat0[, paste0("T", 1:9)],
    nT = 9, curveFun = "BLS"
  ))
  Y_growth_factor <- reparam_bls(nlpsem:::getUNI.GF(
    dat_traj = RMS_dat0[, paste0("S", 2:9)],
    dat_time = RMS_dat0[, paste0("T", 2:9)],
    nT = 8, curveFun = "BLS"
  ))

  M_reg_1 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = M_growth_factor[, 1], X_growth_factor[, 1]),
                           na.action = na.exclude)$coefficients)
  M_reg_r <- as.numeric(lm(M_gammaV ~ ., data = data.frame(M_gammaV = M_growth_factor[, 2], X_growth_factor[, 1:2]),
                           na.action = na.exclude)$coefficients)
  M_reg_2 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = M_growth_factor[, 3], X_growth_factor[, 1:3]),
                           na.action = na.exclude)$coefficients)
  M_beta_raw <- matrix(c(M_reg_1[-1], 0, 0, M_reg_r[-1], 0, M_reg_2[-1]), byrow = TRUE, nrow = 3, ncol = 3)

  Y_reg_1 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = Y_growth_factor[, 1], M_growth_factor[, 1], X_growth_factor[, 1]),
                           na.action = na.exclude)$coefficients)
  Y_reg_r <- as.numeric(lm(Y_gammaV ~ ., data = data.frame(Y_gammaV = Y_growth_factor[, 2], M_growth_factor[, 1:2], X_growth_factor[, 1:2]),
                           na.action = na.exclude)$coefficients)
  Y_reg_2 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = Y_growth_factor[, 3], M_growth_factor[, 1:3], X_growth_factor[, 1:3]),
                           na.action = na.exclude)$coefficients)
  Y_beta_raw <- matrix(c(Y_reg_1[2], rep(0, 2), Y_reg_1[3], rep(0, 2),
                         Y_reg_r[2:3], 0, Y_reg_r[4:5], 0, Y_reg_2[-1]),
                       byrow = TRUE, nrow = 3, ncol = 6)

  expect_true(all(is.finite(unlist(starts))))
  expect_lt(abs(starts$M$beta[2, 1]), abs(M_beta_raw[2, 1]))
  expect_lt(abs(starts$Y$beta_x[2, 1]), abs(Y_beta_raw[2, 4]))
  expect_lt(abs(starts$Y$beta_m[2, 1]), abs(Y_beta_raw[2, 1]))
  expect_true(all(eigen(starts$X$covMatrix, symmetric = TRUE)$values > 0))
  expect_true(all(eigen(starts$M$covMatrix, symmetric = TRUE)$values > 0))
  expect_true(all(eigen(starts$Y$covMatrix, symmetric = TRUE)$values > 0))
})

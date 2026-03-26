test_that("getMGM preserves multi-character outcome prefixes", {
  dat <- data.frame(
    T1 = c(0, 0, 0, 0),
    T2 = c(1, 1, 1, 1),
    score1 = c(1, 2, 3, 4),
    score2 = c(2, 3, 4, 5),
    memory1 = c(2, 3, 4, 5),
    memory2 = c(3, 4, 5, 6)
  )

  captured_y_var <- NULL
  local_mocked_bindings(
    getMULTI.initial = function(dat, t_var, y_var, curveFun, records, res_scale, res_cor, intrinsic = FALSE) {
      captured_y_var <<- y_var
      list(dummy = TRUE)
    },
    getMGM.mxModel = function(dat, t_var, y_var, curveFun, intrinsic, records, y_model, starts) {
      captured_y_var <<- y_var
      OpenMx::mxModel("stub", type = "RAM")
    },
    mxRun = function(model) model,
    .package = "nlpsem"
  )

  fit <- nlpsem::getMGM(
    dat = dat,
    t_var = c("T", "T"),
    y_var = c("score", "memory"),
    curveFun = "LIN",
    intrinsic = FALSE,
    records = list(1:2, 1:2),
    y_model = "LGCM"
  )

  expect_s4_class(fit, "myMxOutput")
  expect_equal(captured_y_var, c("score", "memory"))
})

test_that("intrinsic multivariate BLS starts damp knot-related covariance", {
  data("RMS_dat", package = "nlpsem")
  RMS_dat0 <- RMS_dat
  baseT <- RMS_dat0$T1
  for (j in 1:9) RMS_dat0[[paste0("T", j)]] <- RMS_dat0[[paste0("T", j)]] - baseT

  fixed_starts <- nlpsem:::getMULTI.initial(
    dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"),
    curveFun = "BLS", records = list(1:9, 1:9),
    res_scale = NULL, res_cor = NULL, intrinsic = FALSE
  )
  intrinsic_starts <- nlpsem:::getMULTI.initial(
    dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"),
    curveFun = "BLS", records = list(1:9, 1:9),
    res_scale = NULL, res_cor = NULL, intrinsic = TRUE
  )

  expect_lt(intrinsic_starts$covMatrix[4, 4], fixed_starts$covMatrix[4, 4])
  expect_lt(intrinsic_starts$covMatrix[8, 8], fixed_starts$covMatrix[8, 8])
  expect_equal(unname(intrinsic_starts$covMatrix[4, 1]), 0)
  expect_equal(unname(intrinsic_starts$covMatrix[4, 2]), 0)
  expect_equal(unname(intrinsic_starts$covMatrix[4, 8]), 0)
})

test_that("getMGroup.initial forwards intrinsic to multivariate starts", {
  dat <- data.frame(
    grp = c(1, 1, 2, 2),
    T1 = c(0, 0, 0, 0),
    T2 = c(1, 1, 1, 1),
    R1 = c(1, 2, 3, 4),
    R2 = c(2, 3, 4, 5),
    M1 = c(2, 3, 4, 5),
    M2 = c(3, 4, 5, 6)
  )

  intrinsic_seen <- logical()
  local_mocked_bindings(
    getMULTI.initial = function(dat, t_var, y_var, curveFun, records, res_scale, res_cor, intrinsic = FALSE) {
      intrinsic_seen <<- c(intrinsic_seen, intrinsic)
      list(means = list(c(0, 0), c(0, 0)), covMatrix = diag(4), residuals = diag(2))
    },
    .package = "nlpsem"
  )

  starts <- nlpsem:::getMGroup.initial(
    dat = dat, nClass = 2, grp_var = "grp", sub_Model = "MGM",
    t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS",
    intrinsic = TRUE, records = list(1:2, 1:2), m_var = NULL, x_var = NULL,
    x_type = NULL, TVC = NULL, decompose = NULL, growth_TIC = NULL,
    res_scale = NULL, res_cor = NULL
  )

  expect_length(starts, 2)
  expect_equal(intrinsic_seen, c(TRUE, TRUE))
})

test_that("getMIX.initial forwards intrinsic to multivariate starts", {
  dat <- data.frame(
    T1 = c(0, 0, 0, 0),
    T2 = c(1, 1, 1, 1),
    R1 = c(1, 2, 8, 9),
    R2 = c(2, 3, 9, 10),
    M1 = c(2, 3, 9, 10),
    M2 = c(3, 4, 10, 11)
  )

  intrinsic_seen <- logical()
  local_mocked_bindings(
    getUNI.GF = function(dat_traj, dat_time, nT, curveFun) {
      matrix(
        c(0, 0, 0, 0,
          0, 0, 0, 0,
          10, 10, 10, 10,
          10, 10, 10, 10),
        nrow = nrow(dat_traj), byrow = TRUE
      )
    },
    getMULTI.initial = function(dat, t_var, y_var, curveFun, records, res_scale, res_cor, intrinsic = FALSE) {
      intrinsic_seen <<- c(intrinsic_seen, intrinsic)
      list(means = list(c(0, 0), c(0, 0)), covMatrix = diag(4), residuals = diag(2))
    },
    .package = "nlpsem"
  )

  starts <- nlpsem:::getMIX.initial(
    dat = dat, nClass = 2, prop_starts = c(0.5, 0.5), sub_Model = "MGM",
    cluster_TIC = NULL, t_var = c("T", "T"), records = list(1:2, 1:2),
    y_var = c("R", "M"), curveFun = "BLS", m_var = NULL, x_var = NULL,
    x_type = NULL, TVC = NULL, decompose = NULL, growth_TIC = NULL,
    res_scale = NULL, res_cor = NULL, intrinsic = TRUE
  )

  expect_length(starts, 3)
  expect_equal(intrinsic_seen, c(TRUE, TRUE))
})

test_that("getMGroup.initial damps class-specific non-intrinsic BLS MGM covariance starts", {
  dat <- data.frame(
    grp = c(1, 1, 2, 2),
    T1 = c(0, 0, 0, 0),
    T2 = c(1, 1, 1, 1),
    R1 = c(1, 2, 8, 9),
    R2 = c(2, 3, 9, 10),
    M1 = c(2, 3, 9, 10),
    M2 = c(3, 4, 10, 11)
  )
  cov_mat <- diag(8)

  local_mocked_bindings(
    getMULTI.initial = function(dat, t_var, y_var, curveFun, records, res_scale, res_cor, intrinsic = FALSE) {
      list(means = list(c(0, 0, 0, 4), c(0, 0, 0, 5)), covMatrix = cov_mat, residuals = diag(2))
    },
    .package = "nlpsem"
  )

  starts <- nlpsem:::getMGroup.initial(
    dat = dat, nClass = 2, grp_var = "grp", sub_Model = "MGM",
    t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS",
    intrinsic = FALSE, records = list(1:2, 1:2), m_var = NULL, x_var = NULL,
    x_type = NULL, TVC = NULL, decompose = NULL, growth_TIC = NULL,
    res_scale = list(c(0.1, 0.1), c(0.1, 0.1)), res_cor = list(0.3, 0.3)
  )

  keep_idx <- c(1:3, 5:7)
  expect_equal(diag(starts[[1]]$covMatrix)[keep_idx], rep(0.6, length(keep_idx)))
  expect_equal(diag(starts[[1]]$covMatrix)[c(4, 8)], c(1, 1))
})

test_that("mixture weight starts use softmax logits", {
  expect_equal(
    nlpsem:::.mixture_weight_starts(c(0.2, 0.3, 0.5)),
    c(0, log(0.3 / 0.2), log(0.5 / 0.2))
  )

  expect_equal(
    nlpsem:::.mixture_weight_starts(c(2, 3, 5)),
    c(0, log(3 / 2), log(5 / 2))
  )
})

test_that("getMIX.initial stores class starts on the softmax scale", {
  dat <- data.frame(
    T1 = c(0, 0, 0, 0),
    T2 = c(1, 1, 1, 1),
    Y1 = c(1, 2, 8, 9),
    Y2 = c(2, 3, 9, 10)
  )

  local_mocked_bindings(
    getUNI.GF = function(dat_traj, dat_time, nT, curveFun) {
      matrix(
        c(0, 0,
          0, 0,
          10, 10,
          10, 10),
        nrow = nrow(dat_traj), byrow = TRUE
      )
    },
    getUNI.initial = function(dat, t_var, y_var, curveFun, records, growth_TIC, res_scale, intrinsic = NULL) {
      list(Y_starts = list(means = c(0, 0), covMatrix = diag(2), residuals = 1))
    },
    .package = "nlpsem"
  )

  starts <- nlpsem:::getMIX.initial(
    dat = dat, nClass = 2, prop_starts = c(0.2, 0.8), sub_Model = "LGCM",
    cluster_TIC = NULL, t_var = "T", records = 1:2, y_var = "Y",
    curveFun = "LIN", m_var = NULL, x_var = NULL, x_type = NULL,
    TVC = NULL, decompose = NULL, growth_TIC = NULL,
    res_scale = NULL, res_cor = NULL, intrinsic = FALSE
  )

  expect_equal(starts[[3]], c(0, log(0.8 / 0.2)))
})

test_that("cluster TIC starts are finite and keep constant predictors at zero slope", {
  dat <- data.frame(
    label = c(1, 1, 1, 2, 2, 2),
    gx1 = c(-2, -1, 0, 0.5, 1, 2),
    gx2 = c(5, 5, 5, 5, 5, 5)
  )

  starts <- nlpsem:::.cluster_tic_weight_starts(dat, c("gx1", "gx2"), nClass = 2)

  expect_equal(dim(starts), c(2, 3))
  expect_true(all(is.finite(starts)))
  expect_equal(starts[1, ], c(0, 0, 0))
  expect_equal(unname(starts[2, 3]), 0)
})

test_that("cluster TIC starts fall back to intercept-only logits when predictors do not vary", {
  dat <- data.frame(
    label = c(1, 1, 2, 2),
    gx1 = c(3, 3, 3, 3),
    gx2 = c(7, 7, 7, 7)
  )
  starts <- NULL

  expect_warning(
    starts <- nlpsem:::.cluster_tic_weight_starts(dat, c("gx1", "gx2"), nClass = 2),
    "using intercept-only softmax starts"
  )

  expect_equal(unname(starts[, 1, drop = FALSE]), unname(nlpsem:::.cluster_tic_intercept_starts(dat$label, 2)))
  expect_equal(starts[, -1, drop = FALSE], matrix(0, nrow = 2, ncol = 2))
})

test_that("cluster TIC output names adapt to the number of predictors", {
  expect_equal(
    nlpsem:::.cluster_tic_param_names(3, c("gx1", "gx2", "gx3")),
    c("beta20", "beta21", "beta22", "beta23",
      "beta30", "beta31", "beta32", "beta33")
  )
})

test_that(".ensure_pd repairs slightly indefinite matrices", {
  mat <- matrix(c(1, 2, 2, 1), nrow = 2)
  repaired <- NULL

  expect_warning(
    repaired <- nlpsem:::.ensure_pd(mat, tol = 1e-6),
    "projecting to nearest positive-definite matrix"
  )

  expect_true(isSymmetric(repaired))
  expect_true(all(eigen(repaired, symmetric = TRUE)$values >= 1e-6 - 1e-10))
})

test_that(".append_between_output appends one pair at a time", {
  out <- list()
  out <- nlpsem:::.append_between_output(out, matrix(c(1, 2, 3, 4), nrow = 2), 0.5)
  out <- nlpsem:::.append_between_output(out, matrix(c(5, 6, 7, 8), nrow = 2), 0.9)

  expect_equal(out[[1]], c(1, 2, 3, 4, 0.5))
  expect_equal(out[[2]], c(5, 6, 7, 8, 0.9))
  expect_equal(length(unlist(out)), 10)
})

test_that(".multi_reduced_cov_starts keeps three transformed coefficients per trajectory", {
  cov_mat <- matrix(seq_len(12 * 12), nrow = 12)
  cov_mat <- cov_mat + t(cov_mat)
  keep_idx <- c(1:3, 5:7, 9:11)
  expected <- t(cov_mat[keep_idx, keep_idx])[row(t(cov_mat[keep_idx, keep_idx])) >= col(t(cov_mat[keep_idx, keep_idx]))]

  reduced <- nlpsem:::.multi_reduced_cov_starts(cov_mat, n_traj = 3)

  expect_equal(reduced, expected)
  expect_length(reduced, 45)
})

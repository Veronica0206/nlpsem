.interval_fixture <- function(records = c(1, 3, 5, 7, 9)) {
  dat <- data.frame(group = c(1, 1, 2, 2))
  for (pos in seq_along(records)) {
    time <- (pos - 1) * dat$group
    dat[[paste0("T", records[pos])]] <- time
    dat[[paste0("U", records[pos])]] <- time * 3
    dat[[paste0("Y", records[pos])]] <- 10 + 3 * time
    dat[[paste0("Z", records[pos])]] <- 20 + 2 * time
    dat[[paste0("X", records[pos])]] <- 5 + time
  }
  uni <- list(Y_starts = list(means = c(10, 3), covMatrix = c(4, 0.2, 1),
                              residuals = 1, rel_rate = rep(1, length(records) - 1)))
  multi <- list(means = list(c(10, 3), c(20, 2)), covMatrix = diag(c(4, 1, 4, 1)),
                residuals = diag(2), rel_rate = rep(list(rep(1, length(records) - 1)), 2))
  tvc <- uni
  tvc$TVC_starts <- list(means = rep(5, length(records)), covMatrix = rep(1, length(records)))
  tvc$kappa <- 0
  list(dat = dat, records = records, uni = uni, multi = multi, tvc = tvc)
}

.derived_value <- function(model, name) {
  as.numeric(OpenMx::mxEvalByName(name, model, compute = TRUE))
}

test_that("single-group LCSMs use selected wave positions and previous selected occasions", {
  for (records in list(c(1, 3, 5, 7, 9), c(2, 4, 6, 8, 10))) {
    f <- .interval_fixture(records)
    m <- nlpsem:::getLCSM.mxModel(f$dat[f$dat$group == 1, ], t_var = "T", y_var = "Y",
                                  curveFun = "NonP", intrinsic = FALSE, records = records,
                                  growth_TIC = NULL, starts = f$uni)
    expect_equal(.derived_value(m, "Ychg_inv_m"), rep(3, 4))
    expect_equal(.derived_value(m, "Ychg_inv_v"), rep(1, 4))
    expect_equal(.derived_value(m, "Ychg_bl_m"), c(3, 6, 9, 12))
    expect_equal(as.numeric(OpenMx::mxGetExpected(m, "means")), c(10, 13, 16, 19, 22))
    expect_equal(.derived_value(m, paste0("lag", records[2])), 1)
    expect_equal(.derived_value(m, paste0("mid_t", records[2])), 0.5)
  }
})

test_that("multiple-group LCSM and TVC changes use each group's observed intervals", {
  f <- .interval_fixture()
  common <- list(dat = f$dat, nClass = 2, t_var = "T", y_var = "Y", curveFun = "NonP",
                 intrinsic = FALSE, records = f$records, growth_TIC = NULL, container = "multigroup",
                 grp_var = "group")
  lcsm <- do.call(nlpsem:::getsub.LCSM, c(common, list(starts = rep(list(f$uni), 2))))
  tvc <- do.call(nlpsem:::getsub.TVC, c(common, list(y_model = "LCSM", TVC = "X", decompose = 0,
                                                  starts = rep(list(f$tvc), 2))))
  for (models in list(lcsm, tvc)) {
    for (k in 1:2) {
      m <- models[[k]]
      expect_equal(as.numeric(m$M_TIME$values), (0:4) * k)
      expect_equal(.derived_value(m, paste0("c", k, "Ychg_inv_m")), rep(3 * k, 4))
      expect_equal(.derived_value(m, paste0("c", k, "Ychg_inv_v")), rep(k^2, 4))
      expect_equal(.derived_value(m, paste0("c", k, "Ychg_bl_m")), (1:4) * 3 * k)
      expect_s4_class(OpenMx::mxRun(m, onlyFrontend = TRUE, silent = TRUE), "MxModel")
    }
  }
  # Latent classes use a common pooled grid, rather than an arbitrary provisional assignment.
  common$container <- "mixture"
  mix <- do.call(nlpsem:::getsub.LCSM, c(common, list(starts = rep(list(f$uni), 2))))
  mix_tvc <- do.call(nlpsem:::getsub.TVC, c(common, list(y_model = "LCSM", TVC = "X", decompose = 0,
                                                      starts = rep(list(f$tvc), 2))))
  for (models in list(mix, mix_tvc)) {
    expect_equal(as.numeric(models[[1]]$M_TIME$values), (0:4) * 1.5)
    expect_equal(models[[1]]$M_TIME$values, models[[2]]$M_TIME$values)
  }
})

test_that("multivariate LCSMs preserve each outcome's and group's time grid", {
  f <- .interval_fixture()
  args <- list(dat = f$dat, t_var = c("T", "U"), y_var = c("Y", "Z"), curveFun = "NonP",
               intrinsic = FALSE, records = rep(list(f$records), 2), y_model = "LCSM")
  single <- do.call(nlpsem:::getMGM.mxModel, c(args, list(starts = f$multi)))
  expect_equal(.derived_value(single, "Ychg_inv_m"), rep(4.5, 4))
  expect_equal(.derived_value(single, "Zchg_inv_m"), rep(9, 4))
  expect_s4_class(OpenMx::mxRun(single, onlyFrontend = TRUE, silent = TRUE), "MxModel")
  grouped <- do.call(nlpsem:::getsub.MGM, c(args, list(nClass = 2, starts = rep(list(f$multi), 2),
                                                    container = "multigroup", grp_var = "group")))
  for (k in 1:2) {
    m <- grouped[[k]]
    expect_equal(as.numeric(m$YM_TIME$values), (0:4) * k)
    expect_equal(as.numeric(m$ZM_TIME$values), (0:4) * 3 * k)
    expect_equal(.derived_value(m, paste0("c", k, "Ychg_inv_m")), rep(3 * k, 4))
    expect_equal(.derived_value(m, paste0("c", k, "Zchg_inv_m")), rep(6 * k, 4))
    expect_equal(.derived_value(m, paste0("c", k, "Zchg_inv_v")), rep(9 * k^2, 4))
    expect_s4_class(OpenMx::mxRun(m, onlyFrontend = TRUE, silent = TRUE), "MxModel")
  }
  mix <- do.call(nlpsem:::getsub.MGM, c(args, list(nClass = 2, starts = rep(list(f$multi), 2),
                                                container = "mixture")))
  expect_equal(mix[[1]]$YM_TIME$values, mix[[2]]$YM_TIME$values)
  expect_equal(mix[[1]]$ZM_TIME$values, mix[[2]]$ZM_TIME$values)
  expect_equal(as.numeric(mix[[1]]$ZM_TIME$values), (0:4) * 4.5)
})

test_that("decomposed TVCs map selected waves to consecutive rate parameters", {
  for (records in list(c(1, 3, 5, 7, 9), c(2, 4, 6, 8, 10))) {
    f <- .interval_fixture(records)
    starts <- f$uni
    starts$TVC_starts <- list(means = c(5, 2), covMatrix = c(1, 0, 1), residuals = 1,
                              rel_rate = c(1, 2, 3, 4))
    starts$kappa <- 0
    starts$beta <- c(0, 0)
    starts$XY_cov <- 0
    for (decompose in 1:3) {
      m <- nlpsem:::getTVC.mxModel(f$dat, t_var = "T", y_var = "Y", curveFun = "NonP",
                                  intrinsic = FALSE, records = records, y_model = "LCSM", TVC = "X",
                                  decompose = decompose, growth_TIC = NULL, starts = starts)
      expect_s4_class(OpenMx::mxRun(m, onlyFrontend = TRUE, silent = TRUE), "MxModel")
      for (pos in 2:5) {
        expect_equal(.derived_value(m, paste0("X_abs_rate", records[pos])), 2 * (pos - 1))
        expect_equal(.derived_value(m, paste0("lag", records[pos])), 1)
      }
    }
  }
})

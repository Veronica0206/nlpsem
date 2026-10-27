# Tests for the class-specific TVC Jacobian, the recoding of grouping variables, and initial
# values computed from data that contain persons with too few observations.

.rebaseline_full <- function() {
  dat <- .rebaseline_rms()
  dat$ex1 <- as.numeric(scale(dat$Approach_to_Learning))
  BL_mean <- mean(dat$R1); BL_sd <- sd(dat$R1)
  for (j in 1:9) dat[[paste0("Rs", j)]] <- (dat[[paste0("R", j)]] - BL_mean) / BL_sd
  dat
}

# Evaluate class k's Jacobian through the container: multiple-group sub-models hold their own
# data, mixture sub-models read the container's data via definition variables.
.knot_jacobian_entry <- function(model, k) {
  sub <- model@submodels[[k]]
  g <- OpenMx::mxEvalByName(paste0("Class", k, ".c", k, "grad"), model, compute = TRUE)
  mean_name <- intersect(paste0("c", k, c("Y_alpha_s", "Y_mean_s")), names(sub@algebras))[1]
  e <- c(OpenMx::mxEvalByName(paste0("Class", k, ".", mean_name), model, compute = TRUE))
  list(actual = unname(g[1, 4]), expected = unname(-(e[2] - e[3])), dims = dim(g))
}

test_that("class-specific TVC builders carry the knot-variance term in the Jacobian", {
  dat <- .rebaseline_full()
  configs <- list(list(decompose = 1, tic = "ex1"), list(decompose = 0, tic = NULL))
  for (cfg in configs) {
    st <- suppressWarnings(nlpsem:::getMGroup.initial(
      dat, nClass = 2, grp_var = "SEX", sub_Model = "TVC", t_var = "T", y_var = "M", curveFun = "BLS",
      intrinsic = TRUE, records = 1:9, m_var = NULL, x_var = NULL, x_type = NULL, TVC = "Rs",
      decompose = cfg$decompose, growth_TIC = cfg$tic, res_scale = NULL, res_cor = NULL))
    mg <- nlpsem:::getMGroup.mxModel(
      dat, grp_var = "SEX", nClass = 2, sub_Model = "TVC", t_var = "T", y_var = "M", curveFun = "BLS",
      intrinsic = TRUE, records = 1:9, y_model = "LGCM", m_var = NULL, x_var = NULL, x_type = NULL,
      TVC = "Rs", decompose = cfg$decompose, growth_TIC = cfg$tic, starts = st, res_cor = NULL)
    for (k in 1:2) {
      j <- .knot_jacobian_entry(mg, k)
      expect_equal(j$dims, c(4, 4))
      expect_equal(j$actual, j$expected, info = paste("multiple-group, decompose", cfg$decompose, "class", k))
      expect_gt(abs(j$actual), 0)
    }
  }
  stx <- suppressWarnings(nlpsem:::getMIX.initial(
    dat, nClass = 2, prop_starts = c(0.5, 0.5), sub_Model = "TVC", cluster_TIC = NULL, t_var = "T",
    records = 1:9, y_var = "M", curveFun = "BLS", m_var = NULL, x_var = NULL, x_type = NULL, TVC = "Rs",
    decompose = 1, growth_TIC = "ex1", res_scale = NULL, res_cor = NULL, intrinsic = TRUE))
  mx <- nlpsem:::getMIX.mxModel(
    dat, nClass = 2, sub_Model = "TVC", cluster_TIC = NULL, t_var = "T", y_var = "M", curveFun = "BLS",
    intrinsic = TRUE, records = 1:9, y_model = "LGCM", m_var = NULL, x_var = NULL, x_type = NULL,
    TVC = "Rs", decompose = 1, growth_TIC = "ex1", starts = stx, res_cor = NULL)
  for (k in 1:2) {
    j <- .knot_jacobian_entry(mx, k)
    expect_equal(j$actual, j$expected, info = paste("mixture class", k))
    expect_gt(abs(j$actual), 0)
  }
})

test_that(".recode_group_codes maps any coding to 1..K and reports it", {
  expect_silent(r <- nlpsem:::.recode_group_codes(c(1, 2, 1, 2), "g"))
  expect_equal(r$codes, c(1L, 2L, 1L, 2L))
  expect_equal(r$nClass, 2)
  expect_message(r <- nlpsem:::.recode_group_codes(c(0, 1, 1, 0), "g"), "0 -> 1, 1 -> 2")
  expect_equal(r$codes, c(1L, 2L, 2L, 1L))
  expect_message(r <- nlpsem:::.recode_group_codes(c("M", "F", "F"), "g"), "F -> 1, M -> 2")
  expect_equal(r$codes, c(2L, 1L, 1L))
  expect_message(r <- nlpsem:::.recode_group_codes(factor(c("b", "a", "b"), levels = c("b", "a")), "g"), "b -> 1, a -> 2")
  expect_equal(r$codes, c(1L, 2L, 1L))
  expect_warning(r <- nlpsem:::.recode_group_codes(c(1, NA, 2), "g"), "removed")
  expect_equal(r$codes, c(1L, NA, 2L))
  expect_error(nlpsem:::.recode_group_codes(c(1, 1, 1), "g"), "at least two")
})

test_that("getMGroup gives the same fit for any coding of the grouping variable", {
  skip_on_cran()
  dat <- .rebaseline_rms()
  dat$SEX01 <- dat$SEX - 1
  dat$SEXchr <- ifelse(dat$SEX == 1, "F", "M")
  fit <- function(data, g) suppressWarnings(getMGroup(data, grp_var = g, sub_Model = "LGCM", y_var = "M", t_var = "T",
                                                      records = 1:9, curveFun = "LIN", intrinsic = FALSE))
  ll <- function(f) f@mxOutput@output$Minus2LogLikelihood
  sizes <- function(f) unname(sapply(f@mxOutput@submodels, function(s) nrow(s$data$observed)))
  f1 <- suppressMessages(fit(dat, "SEX"))
  f2 <- suppressMessages(fit(dat, "SEX01"))
  msgs <- character()
  f3 <- withCallingHandlers(fit(dat, "SEXchr"),
                            message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("recoded for estimation: F -> 1, M -> 2", msgs)))
  expect_equal(ll(f2), ll(f1), tolerance = 1e-6)
  expect_equal(ll(f3), ll(f1), tolerance = 1e-6)
  expect_equal(sizes(f2), sizes(f1))
  expect_equal(sizes(f3), sizes(f1))
  dat_na <- dat
  dat_na$SEX[1:5] <- NA
  warns <- character()
  f4 <- withCallingHandlers(
    suppressMessages(getMGroup(dat_na, grp_var = "SEX", sub_Model = "LGCM", y_var = "M", t_var = "T",
                               records = 1:9, curveFun = "LIN", intrinsic = FALSE)),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("5 row\\(s\\) with a missing value", warns)))
  expect_equal(sum(sizes(f4)), nrow(dat) - 5)
})

test_that("starts tolerate persons with too few observations", {
  dat <- .rebaseline_full()
  dat[1, paste0("M", 1:9)] <- NA   # no observation at all
  dat[2, paste0("M", 2:9)] <- NA   # one observation
  dat[3, paste0("M", 3:9)] <- NA   # two observations

  gf <- nlpsem:::getUNI.GF(dat[, paste0("M", 1:9)], dat[, paste0("T", 1:9)], 9, "LIN")
  expect_true(all(is.na(gf[1:2, ])))
  expect_true(all(is.finite(as.matrix(gf[-(1:2), ]))))
  gfq <- nlpsem:::getUNI.GF(dat[, paste0("M", 1:9)], dat[, paste0("T", 1:9)], 9, "QUAD")
  expect_true(all(is.na(gfq[1:3, ])))
  expect_true(all(is.finite(as.matrix(gfq[-(1:3), ]))))

  for (cf in c("LIN", "QUAD", "EXP", "JB", "BLS", "NonP")) {
    s <- suppressWarnings(nlpsem:::getUNI.initial(dat, "T", "M", cf, 1:9, growth_TIC = NULL, res_scale = NULL))
    expect_true(all(is.finite(unlist(s$Y_starts$means))), info = cf)
    expect_true(all(is.finite(s$Y_starts$covMatrix)), info = cf)
    expect_true(is.finite(s$Y_starts$residuals), info = cf)
  }
  for (cf in c("LIN", "BLS")) {
    s <- suppressWarnings(nlpsem:::getMULTI.initial(dat, c("T", "T"), c("R", "M"), cf, list(1:9, 1:9),
                                                    res_scale = NULL, res_cor = NULL))
    expect_true(all(is.finite(unlist(s$means))), info = cf)
    expect_true(all(is.finite(s$covMatrix)), info = cf)
    expect_true(all(is.finite(s$residuals)), info = cf)
  }
  for (cf in c("LIN", "BLS")) {
    s <- suppressWarnings(nlpsem:::getMED.initial(dat, t_var = rep("T", 2), y_var = "M", m_var = "R", x_type = "baseline",
                                                  x_var = "ex1", curveFun = cf, records = list(1:9, 1:9),
                                                  res_scale = NULL, res_cor = NULL))
    expect_true(all(is.finite(unlist(s))), info = cf)
  }
  s <- suppressWarnings(nlpsem:::getTVC.initial(dat, "T", "M", "BLS", 1:9, growth_TIC = "ex1", TVC = "Rs", decompose = 1,
                                                res_scale = NULL, res_cor = NULL, intrinsic = TRUE))
  expect_true(all(is.finite(unlist(s))))
  set.seed(1)
  s <- suppressWarnings(nlpsem:::getMIX.initial(dat, nClass = 2, prop_starts = c(0.5, 0.5), sub_Model = "LGCM",
                                                cluster_TIC = NULL, t_var = "T", records = 1:9, y_var = "M", curveFun = "BLS",
                                                m_var = NULL, x_var = NULL, x_type = NULL, TVC = NULL, decompose = NULL,
                                                growth_TIC = NULL, res_scale = NULL, res_cor = NULL, intrinsic = FALSE))
  expect_length(s, 3)
  expect_true(all(is.finite(unlist(s[1:2]))))
})

test_that("provisional class labels cover every person", {
  gf <- matrix(c(0, 0, 0, 0, 10, 10, 10, 10, NA, 9.5, NA, NA), ncol = 2, byrow = TRUE)
  set.seed(2)
  lab <- nlpsem:::.provisional_class_labels(gf, 2)
  expect_length(lab, 6)
  expect_true(all(lab %in% 1:2))
  expect_equal(lab[1:4], c(1L, 1L, 2L, 2L))
  expect_equal(lab[5], 2L)
  expect_error(nlpsem:::.provisional_class_labels(gf[c(1, 5, 6), ], 2), "complete growth-factor starts")
})

test_that("provisional classes too small for their starts are formed on normal scores", {
  # Two persons far out in the second growth factor: k-means on the raw values isolates them.
  set.seed(1)
  gf <- cbind(rnorm(40, 30, 5), c(rnorm(38, 80, 10), 1200, 1900), rnorm(40, 0.5, 0.1))
  set.seed(2)
  expect_equal(min(tabulate(nlpsem:::.provisional_class_labels(gf, 2))), 2)
  set.seed(2)
  lab <- nlpsem:::.provisional_class_labels(gf, 2, min_size = nlpsem:::.min_class_size(3, 2))
  expect_gte(min(tabulate(lab, nbins = 2)), 4)
  expect_equal(nlpsem:::.min_class_size(3), 4L)
  expect_equal(nlpsem:::.min_class_size(2, 4), 6L)
  set.seed(2)
  expect_warning(small <- nlpsem:::.provisional_class_labels(gf[1:6, ], 2, min_size = 4), "each class needs at least 4")
  set.seed(2)
  expect_identical(small, nlpsem:::.provisional_class_labels(gf[1:6, ], 2))
})

test_that("covariate effects a class cannot determine start at 0", {
  expect_warning(beta <- nlpsem:::.zero_aliased_effects(matrix(c(1, NA, 2, 3), 2)), "starting them at 0")
  expect_equal(beta, matrix(c(1, 0, 2, 3), 2))
  expect_silent(nlpsem:::.zero_aliased_effects(matrix(1:4, 2)))
})

test_that("negative exponential mixture starts with two growth TICs are computed on the reference sample", {
  skip_on_cran()
  # k-means on the individual NLS growth factors forms classes of 198 and 2 persons on this sample;
  # the class of 2 cannot support an intercept and two covariate effects.
  dat <- .golden_data(n = 200)
  for (intrinsic in c(FALSE, TRUE)) {
    set.seed(20191029)
    starts <- suppressWarnings(nlpsem:::getMIX.initial(
      dat = dat, nClass = 2, prop_starts = c(0.5, 0.5), sub_Model = "LGCM", cluster_TIC = NULL, t_var = "T",
      records = 1:9, y_var = "M", curveFun = "EXP", m_var = NULL, x_var = NULL, x_type = NULL, TVC = NULL,
      decompose = NULL, growth_TIC = c("ex1", "gx1"), res_scale = NULL, res_cor = NULL, intrinsic = intrinsic))
    for (k in 1:2) {
      expect_true(all(is.finite(starts[[k]]$beta)), info = paste(intrinsic, k))
      expect_true(all(is.finite(starts[[k]]$Y_starts$covMatrix)), info = paste(intrinsic, k))
    }
  }
})

test_that("class-specific mediation models with a baseline predictor build for any predictor name", {
  dat <- .rebaseline_full()
  for (cf in c("LIN", "BLS")) {
    a <- list(t_var = rep("T", 2), y_var = "M", m_var = "R", x_type = "baseline", x_var = "ex1",
              curveFun = cf, records = list(1:9, 1:9))
    st <- suppressWarnings(do.call(nlpsem:::getMGroup.initial, c(list(
      dat = dat, nClass = 2, grp_var = "SEX", sub_Model = "MED", intrinsic = NULL, TVC = NULL,
      decompose = NULL, growth_TIC = NULL, res_scale = NULL, res_cor = NULL), a)))
    mg <- do.call(nlpsem:::getMGroup.mxModel, c(list(
      dat = dat, grp_var = "SEX", nClass = 2, sub_Model = "MED", intrinsic = NULL, y_model = NULL,
      TVC = NULL, decompose = NULL, growth_TIC = NULL, starts = st, res_cor = NULL), a))
    expect_s4_class(mg, "MxModel")
    expect_true("ex1" %in% mg@submodels[[1]]@manifestVars)
    stx <- suppressWarnings(do.call(nlpsem:::getMIX.initial, c(list(
      dat = dat, nClass = 2, prop_starts = c(0.5, 0.5), sub_Model = "MED", cluster_TIC = NULL,
      intrinsic = NULL, TVC = NULL, decompose = NULL, growth_TIC = NULL, res_scale = NULL, res_cor = NULL), a)))
    mx <- do.call(nlpsem:::getMIX.mxModel, c(list(
      dat = dat, nClass = 2, sub_Model = "MED", cluster_TIC = NULL, intrinsic = NULL, y_model = NULL,
      TVC = NULL, decompose = NULL, growth_TIC = NULL, starts = stx, res_cor = NULL), a))
    expect_s4_class(mx, "MxModel")
  }
})

test_that("per-person curve residuals are attached and residual-correlation starts ignore a shared trend", {
  set.seed(1)
  n <- 120; tt <- matrix(rep(0:5, each = n), n)
  eta0 <- rnorm(n, 10, 2); eta1 <- rnorm(n, 2, 0.5)
  y1 <- eta0 + eta1 * tt + matrix(rnorm(n * 6), n)
  y2 <- 5 + 0.8 * eta0 + 1.5 * eta1 * tt + matrix(rnorm(n * 6), n)   # same growth, independent residuals
  dat_traj <- as.data.frame(y1); dat_traj[3, 2] <- NA; dat_time <- as.data.frame(tt)
  for (cf in c("LIN", "QUAD", "EXP", "JB", "BLS", "NonP")) {
    gf <- nlpsem:::getUNI.GF(dat_traj, dat_time, nT = 6, curveFun = cf)
    res <- attr(gf, "residuals")
    expect_equal(dim(res), c(n, 6), info = cf)
    expect_true(is.na(res[3, 2]), info = cf)
    expect_true(all(is.finite(res[-3, ])), info = cf)
  }
  gf <- nlpsem:::getUNI.GF(dat_traj, dat_time, nT = 6, curveFun = "LIN")
  expect_equal(attr(gf, "residuals")[1, ], unname(y1[1, ] - (gf$eta0[1] + gf$eta1[1] * tt[1, ])))
  gfb <- nlpsem:::getUNI.GF(dat_traj, dat_time, nT = 6, curveFun = "BLS")
  k <- gfb$gamma[1]
  expect_equal(attr(gfb, "residuals")[1, ],
               unname(y1[1, ] - (gfb$eta0[1] + gfb$eta1[1] * pmin(tt[1, ], k) + gfb$eta2[1] * pmax(tt[1, ] - k, 0))))
  # Person-mean centring (the 0.4 heuristic) sees the shared growth trend; the curve residuals do not.
  r_old <- cor(c(y1 - rowMeans(y1)), c(y2 - rowMeans(y2)))
  expect_gt(r_old, 0.7)
  gf2 <- nlpsem:::getUNI.GF(as.data.frame(y2), dat_time, nT = 6, curveFun = "LIN")
  r_new <- nlpsem:::.residual_correlation_starts(list(attr(gf, "residuals"), attr(gf2, "residuals")), labels = c("y1", "y2"))
  expect_lt(abs(r_new), 0.15)
  expect_warning(fb <- nlpsem:::.residual_correlation_starts(list(matrix(NA_real_, 3, 2), matrix(NA_real_, 3, 2)),
                                                             labels = c("a", "b")), "falling back")
  expect_equal(fb, 0.3)
  # Three processes give the three pairwise correlations in upper-triangle order
  expect_length(nlpsem:::.residual_correlation_starts(list(attr(gf, "residuals"), attr(gf2, "residuals"), attr(gfb, "residuals"))), 3)
})

test_that("multivariate starts on the shipped data no longer begin from a near-perfect residual correlation", {
  dat <- .golden_data(n = 500)
  st <- nlpsem:::getMULTI.initial(dat = dat, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS",
                                  records = list(1:9, 1:9), res_scale = NULL, res_cor = NULL, intrinsic = FALSE)
  r <- st$residuals[1, 2] / sqrt(st$residuals[1, 1] * st$residuals[2, 2])
  expect_lt(abs(r), 0.6)
  expect_gt(abs(st$residuals[1, 2]), 0)
})

test_that("mixture weight rows report the standard error of the softmax weight", {
  skip_on_cran()
  dat <- .golden_data(n = 200)
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  set.seed(20191029)
  fit <- suppressMessages(getMIX(dat = dat, prop_starts = c(0.5, 0.5), sub_Model = "LGCM", t_var = "T", y_var = "M",
                                 curveFun = "LIN", intrinsic = FALSE, records = 1:9, paramOut = TRUE))
  m <- fit@mxOutput
  w <- fit@Estimates[fit@Estimates$Name == "p2", ]
  expect_equal(w$Estimate, unname(m@output$estimate["w2"]))
  expect_equal(w$SE, unname(m@output$standardErrors["w2", 1]))
  expect_false(isTRUE(all.equal(w$SE, w$Estimate)))
})

test_that("class-specific decomposed-TVC tables read the renamed absolute-rate algebras", {
  skip_on_cran()
  dat <- .golden_data(n = 200)
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  set.seed(20191029)
  fit <- suppressMessages(suppressWarnings(getMGroup(dat = dat, grp_var = "SEX", sub_Model = "TVC", t_var = "T", y_var = "M",
                                                     curveFun = "LIN", intrinsic = FALSE, records = 1:9, y_model = "LGCM",
                                                     TVC = "Rs", decompose = 1, paramOut = TRUE)))
  est <- fit@Estimates
  rate_rows <- est[grepl("X_abs_rate", est$Name), ]
  expect_equal(nrow(rate_rows), 2 * 8)
  expect_true(all(is.finite(rate_rows$Estimate)))
  for (k in 1:2) {
    direct <- unname(OpenMx::mxEvalByName(paste0("c", k, "X_abs_rate2"), fit@mxOutput@submodels[[k]]))
    expect_equal(rate_rows$Estimate[rate_rows$Name == paste0("c", k, "X_abs_rate2")], c(direct))
  }
})

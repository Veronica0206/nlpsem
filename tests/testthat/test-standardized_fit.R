# Quadratic, negative exponential and Jenss-Bayley models are fitted on the standardized scale and
# converted back exactly (R/standardized_fit.R).

.with_std <- function(value, expr) {
  old <- options(nlpsem.standardize = value)
  on.exit(options(old), add = TRUE)
  force(expr)
}

.std_unfitted <- function(fn, args, dat, starts = NULL) {
  call_args <- c(list(dat = dat), args, list(paramOut = FALSE))
  if (!is.null(starts)) call_args$starts <- starts
  suppressMessages(suppressWarnings(do.call(fn, call_args)))
}

test_that("the parameter map reproduces the likelihood exactly for every kind of label", {
  # Building without fitting: the optimizer calls return the model they are given.
  local_mocked_bindings(mxRun = function(model, ...) model, mxTryHard = function(model, ...) model, .package = "nlpsem")
  dat <- .golden_data(n = 200)
  cases <- list(
    LGCM_EXP_r_tic = list(getLGCM, "LGCM", list(t_var = "T", y_var = "M", curveFun = "EXP", intrinsic = TRUE, records = 1:9, growth_TIC = "ex1")),
    LCSM_JB_r = list(getLCSM, "LCSM", list(t_var = "T", y_var = "M", curveFun = "JB", intrinsic = TRUE, records = 1:9)),
    TVC_QUAD_d1_tic = list(getTVCmodel, "TVC", list(t_var = "T", y_var = "M", curveFun = "QUAD", intrinsic = FALSE, records = 1:9,
                                                    y_model = "LGCM", TVC = "Rs", decompose = 1, growth_TIC = "ex1")),
    TVC_LCSM_EXP_r_d0 = list(getTVCmodel, "TVC", list(t_var = "T", y_var = "M", curveFun = "EXP", intrinsic = TRUE, records = 1:9,
                                                      y_model = "LCSM", TVC = "Rs", decompose = 0)),
    MGM_JB_r = list(getMGM, "MGM", list(t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "JB", intrinsic = TRUE,
                                        records = list(1:9, 1:9), y_model = "LGCM")),
    MGroup_TVC_EXP_r_d2 = list(getMGroup, "TVC", list(grp_var = "SEX", sub_Model = "TVC", t_var = "T", y_var = "M", curveFun = "EXP",
                                                      intrinsic = TRUE, records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 2)),
    MGroup_MGM_QUAD = list(getMGroup, "MGM", list(grp_var = "SEX", sub_Model = "MGM", t_var = c("T", "T"), y_var = c("R", "M"),
                                                  curveFun = "QUAD", intrinsic = FALSE, records = list(1:9, 1:9), y_model = "LGCM")),
    MIXc_QUAD = list(getMIX, "LGCM", list(prop_starts = c(0.5, 0.5), sub_Model = "LGCM", cluster_TIC = "gx1", t_var = "T", y_var = "M",
                                          curveFun = "QUAD", intrinsic = FALSE, records = 1:9)),
    # a time column that is also a growth or class-membership covariate
    LGCM_QUAD_timeTIC = list(getLGCM, "LGCM", list(t_var = "T", y_var = "M", curveFun = "QUAD", intrinsic = FALSE, records = 1:9,
                                                   growth_TIC = c("ex1", "T9"))),
    MIXc_EXP_timeTIC = list(getMIX, "LGCM", list(prop_starts = c(0.5, 0.5), sub_Model = "LGCM", cluster_TIC = "T9", t_var = "T",
                                                 y_var = "M", curveFun = "EXP", intrinsic = FALSE, records = 1:9))
  )
  for (nm in names(cases)) {
    fn <- cases[[nm]][[1]]; sub <- cases[[nm]][[2]]; a <- cases[[nm]][[3]]
    set.seed(20191029)
    fit <- .std_unfitted(fn, a, dat)
    raw <- fit@mxOutput
    plan <- nlpsem:::.std_plan(dat, sub, a$curveFun, a$y_var, a$records, a$t_var, a$TVC, a$decompose,
                               covariates = c(a$growth_TIC, a$cluster_TIC),
                               groups = if (!is.null(a$grp_var)) dat[[a$grp_var]])
    expect_false(is.null(plan), info = nm)
    a_std <- a
    a_std$growth_TIC <- nlpsem:::.std_cov(plan, a$growth_TIC)
    a_std$cluster_TIC <- nlpsem:::.std_cov(plan, a$cluster_TIC)
    std <- .with_std(FALSE, .std_unfitted(fn, a_std, nlpsem:::.std_data(dat, plan), starts = getSpecification(fit)$fitting$starts))@mxOutput
    theta <- OpenMx::omxGetParameters(raw)
    map <- nlpsem:::.std_param_map(names(theta), plan, nlpsem:::.std_container_labels(raw))
    expect_length(map$unknown, 0)
    K <- sum(vapply(plan$processes, function(p) 2 * sum(!is.na(as.matrix(dat[, p$columns]))) * log(p$scale), numeric(1)))
    for (point in list(theta, theta * (1 + 0.003 * sin(seq_along(theta))))) {
      o <- OpenMx::mxRun(OpenMx::omxSetParameters(raw, labels = names(point), values = point), useOptimizer = FALSE, silent = TRUE)
      z <- OpenMx::mxRun(OpenMx::omxSetParameters(std, labels = names(point),
                                                 values = (point - map$shift[names(point)]) / map$factor[names(point)]),
                         useOptimizer = FALSE, silent = TRUE)
      expect_equal(o@output$Minus2LogLikelihood - z@output$Minus2LogLikelihood, K, tolerance = 1e-10, info = nm)
    }
  }
})

test_that("unit factors follow the units of each growth factor", {
  dat <- .golden_data(n = 200)
  plan <- nlpsem:::.std_plan(dat, "TVC", "JB", "M", 1:9, "T", "Rs", 1L)
  sy <- sd(dat$M1); sx <- sd(dat$Rs1); b <- mean(dat$T9)
  expect_equal(plan$processes$Y$center, mean(dat$M1))
  expect_equal(plan$time$scale, b)
  map <- nlpsem:::.std_param_map(c("Y_mueta0", "Y_mueta1", "Y_mueta2", "Y_acc_ratio", "Y_psi1g", "Y_psigg", "Y_residuals",
                                   "beta1TIC1", "betagTVC", "kappa", "Cov_XYres", "X_mueta1", "X_psi01", "X_rel_rate3"), plan)
  expect_equal(map$factor, c(Y_mueta0 = sy, Y_mueta1 = sy / b, Y_mueta2 = sy, Y_acc_ratio = 1 / b, Y_psi1g = sy / b^2,
                             Y_psigg = 1 / b^2, Y_residuals = sy^2, beta1TIC1 = sy / b, betagTVC = 1 / (b * sx),
                             kappa = sy * b / sx, Cov_XYres = sy * sx, X_mueta1 = sx / b, X_psi01 = sx^2 / b, X_rel_rate3 = 1))
  expect_equal(map$shift[["Y_mueta0"]], mean(dat$M1))
  expect_equal(unname(map$shift[-1]), rep(0, 13))
  quad <- nlpsem:::.std_param_map(c("Y_mueta2", "Y_psi12"), nlpsem:::.std_plan(dat, "LGCM", "QUAD", "M", 1:9, "T"))
  expect_equal(quad$factor, c(Y_mueta2 = sy / b^2, Y_psi12 = sy^2 / b^3))
  expect_equal(nlpsem:::.std_param_map(c("Y_mueta0", "Y_unknown"), plan)$unknown, "Y_unknown")
})

test_that("multiple-group and mixture fits scale outcomes within groups", {
  dat <- .golden_data(n = 200)
  within <- function(x, g) sqrt(sum((x - ave(x, g))^2) / (length(x) - length(unique(g))))
  plan <- nlpsem:::.std_plan(dat, "TVC", "EXP", "M", 1:9, "T", "Rs", 0, groups = dat$SEX)
  expect_equal(plan$processes$Y$scale, within(dat$M1, dat$SEX))
  expect_equal(plan$processes$X$scale, within(dat$Rs1, dat$SEX))
  expect_identical(plan$processes$Y$scale_basis, "within groups")
  expect_equal(plan$processes$Y$center, mean(dat$M1))
  # one group, no groups, or groups of the wrong length: the overall standard deviation
  for (g in list(rep(1, nrow(dat)), NULL, dat$SEX[-1])) {
    p <- nlpsem:::.std_plan(dat, "LGCM", "EXP", "M", 1:9, "T", groups = g)
    expect_equal(p$processes$Y$scale, sd(dat$M1))
    expect_identical(p$processes$Y$scale_basis, "overall")
  }
})

test_that("standardized fits converge when groups or classes are far apart", {
  skip_on_cran()
  # Two classes of negative exponential trajectories whose intercepts are 50 apart: the overall
  # first-occasion standard deviation is about five times the within-class one.
  set.seed(20261027)
  n <- 300; class <- rep(1:2, each = n / 2)
  eta0 <- c(0, 50)[class] + rnorm(n, sd = 5); eta1 <- 25 + rnorm(n, sd = 3)
  dat <- data.frame(id = seq_len(n), class = class)
  for (j in 1:10) {
    tj <- if (j %in% c(1, 10)) rep(j - 1, n) else j - 1 + runif(n, -0.25, 0.25)
    dat[[paste0("T", j)]] <- tj
    dat[[paste0("Y", j)]] <- eta0 + eta1 * (1 - exp(-0.6 * tj)) + rnorm(n)
  }
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  fit_both <- function(f) {
    std <- f(); raw <- .with_std(FALSE, f())
    list(std = std@mxOutput@output, raw = raw@mxOutput@output, spec = getSpecification(std)$fitting$standardization)
  }
  multigroup <- fit_both(function() suppressMessages(getMGroup(dat = dat, grp_var = "class", sub_Model = "LGCM", y_var = "Y",
                                                               t_var = "T", records = 1:10, curveFun = "EXP", intrinsic = FALSE)))
  mixture <- fit_both(function() { set.seed(20191029)
    suppressMessages(getMIX(dat = dat, prop_starts = c(0.5, 0.5), sub_Model = "LGCM", y_var = "Y", t_var = "T", records = 1:10,
                            curveFun = "EXP", intrinsic = FALSE)) })
  for (res in list(multigroup, mixture)) {
    expect_identical(res$spec$processes$scale_basis, "within groups")
    expect_lt(res$spec$processes$scale, sd(dat$Y1) / 3)
    expect_equal(res$std$status$code, 0)
    expect_equal(res$std$Minus2LogLikelihood, res$raw$Minus2LogLikelihood, tolerance = 1e-8)
  }
})

test_that("a covariate that is also a time column is not rescaled", {
  skip_on_cran()
  dat <- .golden_data(n = 200)
  # time in years since birth, so that the first time column is baseline age
  for (j in 1:9) dat[[paste0("T", j)]] <- RMS_dat[[paste0("T", j)]][1:200] / 12
  dat$baseline_age <- dat$T1
  plan <- nlpsem:::.std_plan(dat, "LGCM", "QUAD", "M", 1:9, "T", covariates = "T1")
  expect_identical(plan$copies, c(T1 = "nlpsemTIC_T1"))
  sdat <- nlpsem:::.std_data(dat, plan)
  expect_identical(sdat$nlpsemTIC_T1, dat$T1)
  expect_equal(sdat$T1, dat$T1 / plan$time$scale)
  expect_identical(nlpsem:::.std_cov(plan, c("ex1", "T1")), c("ex1", "nlpsemTIC_T1"))
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  fit_tic <- function(tic){
    set.seed(20191029)
    suppressMessages(getLGCM(dat = dat, t_var = "T", y_var = "M", curveFun = "QUAD", intrinsic = FALSE, records = 1:9,
                             growth_TIC = tic, paramOut = TRUE))
  }
  shared <- fit_tic("T1"); copied <- fit_tic("baseline_age")
  expect_true(getSpecification(shared)$fitting$standardization$applied)
  expect_equal(shared@mxOutput@output$Minus2LogLikelihood, copied@mxOutput@output$Minus2LogLikelihood, tolerance = 1e-10)
  expect_equal(coef(shared), coef(copied), tolerance = 1e-8)
  expect_equal(shared@mxOutput@output$standardErrors, copied@mxOutput@output$standardErrors, tolerance = 1e-8)
  expect_equal(unname(coef(shared)["mux1"]), mean(dat$T1), tolerance = 1e-5)
})

test_that("a standardized fit carries a consistent raw-scale solution and uncertainty", {
  skip_on_cran()
  dat <- .golden_data(n = 200)
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  set.seed(20191029)
  fit <- suppressMessages(getLGCM(dat = dat, t_var = "T", y_var = "M", curveFun = "EXP", intrinsic = FALSE, records = 1:9, paramOut = TRUE))
  m <- fit@mxOutput; out <- m@output
  rec <- getSpecification(fit)$fitting$standardization
  expect_true(rec$applied)
  expect_equal(rec$processes$center, mean(dat$M1)); expect_equal(rec$processes$scale, sd(dat$M1))
  expect_equal(rec$time$scale, mean(dat$T9))
  expect_lt(rec$condition_number_standardized, rec$condition_number_raw_scale)
  # the stored likelihood is the raw-data likelihood at the stored estimates
  again <- OpenMx::mxRun(m, useOptimizer = FALSE, silent = TRUE)
  expect_equal(out$Minus2LogLikelihood, again@output$Minus2LogLikelihood, tolerance = 1e-10)
  expect_equal(out$estimate, OpenMx::omxGetParameters(m)[names(out$estimate)])
  # standard errors, parameter covariance matrix and Hessian agree with each other
  expect_equal(unname(out$standardErrors[, 1]), unname(sqrt(diag(out$vcov))[rownames(out$standardErrors)]))
  expect_equal(unname(out$vcov %*% out$hessian), diag(2, nrow(out$vcov)), tolerance = 1e-8)
  expect_equal(summary(m)$parameters$Std.Error, unname(out$standardErrors[summary(m)$parameters$name, 1]))
  expect_equal(fit@Estimates$SE[fit@Estimates$Name == "mueta0"], unname(out$standardErrors["Y_mueta0", 1]))
  # the raw-scale fit reaches the same solution on this case
  set.seed(20191029)
  raw <- .with_std(FALSE, suppressMessages(getLGCM(dat = dat, t_var = "T", y_var = "M", curveFun = "EXP", intrinsic = FALSE, records = 1:9)))
  expect_null(getSpecification(raw)$fitting$standardization)
  expect_equal(raw@mxOutput@output$Minus2LogLikelihood, out$Minus2LogLikelihood, tolerance = 1e-8)
  expect_equal(OpenMx::omxGetParameters(raw@mxOutput), OpenMx::omxGetParameters(m), tolerance = 1e-3)
})

test_that("other functional forms and mediation models are fitted as built", {
  dat <- .golden_data(n = 200)
  expect_null(nlpsem:::.std_plan(dat, "LGCM", "LIN", "M", 1:9, "T"))
  expect_null(nlpsem:::.std_plan(dat, "LGCM", "BLS", "M", 1:9, "T"))
  expect_null(nlpsem:::.std_plan(dat, "LCSM", "NonP", "M", 1:9, "T"))
  expect_null(nlpsem:::.std_plan(dat, "MED", "LIN", "M", list(1:9, 1:9), c("T", "T")))
  expect_null(.with_std(FALSE, nlpsem:::.std_plan(dat, "LGCM", "EXP", "M", 1:9, "T")))
  set.seed(20191029)
  fit <- suppressMessages(getLGCM(dat = dat, t_var = "T", y_var = "M", curveFun = "LIN", intrinsic = FALSE, records = 1:9))
  expect_null(getSpecification(fit)$fitting$standardization)
})

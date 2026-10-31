# Independent-review regressions for precision, class parsing and one-class output.

.precision_output_table <- function(multiplier = 1, se = 0.00006) {
  model <- OpenMx::mxModel("precision_probe",
    OpenMx::mxMatrix("Full", 2, 1, values = c(0.00014 * multiplier, 0.5), name = "Y_mean0"),
    OpenMx::mxMatrix("Symm", 2, 2, values = matrix(c(.2, .01, .01, .3), 2), name = "Y_psi0"))
  model@output$estimate <- c(Y_residuals = .1, unrelated = .2)
  model@output$standardErrors <- matrix(c(.01, .02), 2, 1,
    dimnames = list(c("Y_residuals", "unrelated"), "SE"))
  local_mocked_bindings(mxSE = function(expression, model, ...) {
    target <- deparse(substitute(expression))
    if (target == "Y_mean0") return(c(se * multiplier, .1))
    if (target == "Y_psi0") return(matrix(c(.03, .01, .01, .03), 2))
    stop("Unexpected algebra: ", target)
  }, .package = "nlpsem")
  nlpsem:::getLGCM.output(model, growth_TIC = NULL,
    names = c("mu0", "mu1", "psi00", "psi01", "psi11", "residual"))
}

test_that("stored precision preserves inference when measurement units change", {
  original <- .precision_output_table()
  scaled <- .precision_output_table(1000)
  tiny <- .precision_output_table(se = .000000049)
  expect_identical(original$Estimate[1], .00014)
  expect_identical(original$SE[1], .00006)
  expect_identical(tiny$SE[1], .000000049)
  original_stats <- getEstimateStats(est_in = original)@wald["mu0", ]
  scaled_stats <- getEstimateStats(est_in = scaled)@wald["mu0", ]
  tiny_stats <- getEstimateStats(est_in = tiny)@wald["mu0", ]
  expect_equal(original_stats$Estimate / original_stats$SE, 7 / 3)
  expect_equal(scaled_stats$Estimate / scaled_stats$SE, original_stats$Estimate / original_stats$SE)
  expect_identical(original_stats$p.value, "0.0196")
  expect_identical(scaled_stats$p.value, original_stats$p.value)
  expect_identical(tiny_stats$p.value, "<0.0001")
  expect_equal(c(original_stats$wald_lbound, original_stats$wald_ubound),
    .00014 + c(-1, 1) * qnorm(.975) * .00006)
  expect_equal(c(scaled_stats$wald_lbound, scaled_stats$wald_ubound) / 1000,
    c(original_stats$wald_lbound, original_stats$wald_ubound))
  expect_gt(tiny_stats$wald_ubound, tiny_stats$wald_lbound)
})

test_that("numeric bootstrap and likelihood intervals retain sub-decimal precision", {
  estimates <- c(a = .00000123456789)
  raw <- data.frame(statusCode = OpenMx::as.statusCode(rep(0, 4)),
    a = c(.000000123456, .000001234567, .000002345678, .000003456789))
  boot <- nlpsem:::.bootstrap_percentile_CI(raw, OpenMx::as.statusCode(0),
    estimates, 4, .95)
  expect_identical(boot$Estimate, unname(estimates))
  expect_equal(c(boot$boot_lbound, boot$boot_ubound),
    unname(quantile(raw$a, c(.025, .975))))
  expect_gt(boot$boot_ubound, boot$boot_lbound)

  model <- OpenMx::mxModel("interval_probe",
    OpenMx::mxMatrix("Full", 1, 1, values = estimates, free = TRUE, labels = "a", name = "A"))
  intervals <- matrix(c(.000000123456, estimates, .000003456789), 1, 3,
    dimnames = list("a", c("lbound", "estimate", "ubound")))
  local_mocked_bindings(mxTryHard = function(model, ...) {
    model@output$confidenceIntervals <- intervals
    model
  }, .package = "nlpsem")
  out <- getEstimateStats(model, data.frame(Name = "a", Estimate = estimates, SE = 1e-7),
    CI_type = "likelihood")@likelihood
  expect_equal(unlist(out[1, ], use.names = FALSE), c(intervals[1, 2], intervals[1, 1], intervals[1, 3]))
})

test_that("all seven output entry points preserve algebra values and small standard errors", {
  # Real saved OpenMx models supply the estimates and algebra dimensions.
  # Only the SE provider is controlled, to test precision independently of optimization.
  tiny_se <- .000049123456
  local_mocked_bindings(mxSE = function(expression, model, forceName = FALSE, ...) {
    target <- if (forceName) expression else paste(deparse(substitute(expression)), collapse = "")
    value <- OpenMx::mxEvalByName(target, model)
    value[] <- tiny_se
    value
  }, .package = "nlpsem")

  lgcm <- .load_extdata("getLGCM_examples.RData")$Math_LGCM_TIC_BLS_f
  lcsm <- .load_extdata("getLCSM_examples.RData")$Read_LCSM_QUAD
  tvc <- .load_extdata("getTVCmodel_examples.RData")$Math_TVCslp_BLS_f
  mgm <- .load_extdata("getMGM_examples.RData")$RM_PLGCM.f
  med <- .load_extdata("getMediation_examples.RData")$Med2_LGCM_BLS
  group <- .load_extdata("getMGroup_examples.RData")$MGroup_Read_EXP_LGCM_TIC_r
  mix <- .load_extdata("getMIX_examples.RData")$RM_BLS_PLGCM3

  group_names <- nlpsem:::.auto_names_subModel("LGCM", "EXP", FALSE, c("ex1", "ex2"),
    1:9, NULL, NULL, "R", NULL, NULL, NULL)
  mix_names <- nlpsem:::.auto_names_subModel("MGM", "BLS", FALSE, NULL,
    list(1:9, 1:9), "LGCM", NULL, c("R", "M"), NULL, NULL, NULL)
  outputs <- list(
    LGCM = nlpsem:::getLGCM.output(lgcm@mxOutput, c("ex1", "ex2"), lgcm@Estimates$Name),
    LCSM = nlpsem:::getLCSM.output(lcsm@mxOutput, "QUAD", NULL, lcsm@Estimates$Name),
    TVC = nlpsem:::getTVC.output(tvc@mxOutput, "BLS", 1:9, "LGCM", 1, "ex1", tvc@Estimates$Name),
    MGM = nlpsem:::getMGM.output(mgm@mxOutput, c("R", "M"), list(1:9, 1:9), "BLS", "LGCM", mgm@Estimates$Name),
    MED = nlpsem:::getMED.output(med@mxOutput, "M", "R", "baseline", "ex1", "BLS", med@Estimates$Name),
    MGroup = nlpsem:::getMGroup.output(group@mxOutput, 2, "LGCM", "R", "EXP", NULL,
      1:9, c("ex1", "ex2"), NULL, NULL, group_names),
    MIX = nlpsem:::getMIX.output(mix@mxOutput, 3, "MGM", c("gx1", "gx2"), c("R", "M"), "BLS",
      NULL, list(1:9, 1:9), NULL, "LGCM", NULL, mix_names))
  first_values <- c(
    OpenMx::mxEvalByName("Y_alpha0", lgcm@mxOutput)[1],
    OpenMx::mxEvalByName("Y_mean0", lcsm@mxOutput)[1],
    OpenMx::mxEvalByName("Y_alpha0", tvc@mxOutput)[1],
    OpenMx::mxEvalByName("R_mean0", mgm@mxOutput)[1],
    unname(med@mxOutput@output$estimate["muX"]),
    OpenMx::mxEvalByName("Class1.c1Y_alpha0", group@mxOutput)[1],
    OpenMx::mxEvalByName("Class1.c1R_mean0", mix@mxOutput)[1])
  for (i in seq_along(outputs)) {
    expect_equal(outputs[[i]]$Estimate[1], first_values[i], tolerance = 1e-12,
      info = names(outputs)[i])
    expect_true(any(outputs[[i]]$SE == tiny_se), info = names(outputs)[i])
    expect_true(any(abs(outputs[[i]]$Estimate - round(outputs[[i]]$Estimate, 4)) > 1e-8),
      info = names(outputs)[i])
  }
  # A derived LCSM slope and mixture cluster coefficient retain their exact values.
  expect_equal(outputs$LCSM$Estimate[outputs$LCSM$Name == "slp_val_est1"],
    unname(OpenMx::mxEvalByName("Yslp_m", lcsm@mxOutput)[1]), tolerance = 1e-12)
  beta_label <- mix@mxOutput$classbeta$labels[2, 1]
  expect_equal(outputs$MIX$Estimate[outputs$MIX$Name == beta_label],
    unname(mix@mxOutput$classbeta$values[2, 1]), tolerance = 1e-12)
})

test_that("summary residuals retain complete class and group tokens", {
  local_mocked_bindings(AIC = function(...) 120, BIC = function(...) 140,
    getPosterior = function(model, nClass, ...) {
      methods::new("postOutput", membership = as.numeric(seq_len(nClass)))
    }, .package = "nlpsem")
  for (n in 9:11) {
    subs <- lapply(seq_len(n), function(k) {
      OpenMx::mxModel(paste0("Class", k),
        OpenMx::mxMatrix("Full", 1, 1, values = k, free = TRUE,
          labels = paste0("c", k, "Y_residuals"), name = "residual"),
        OpenMx::mxData(data.frame(Y = c(1, 2)), type = "raw"))
    })
    group <- OpenMx::mxModel("group_probe", subs)
    group@output <- list(Minus2LogLikelihood = 100)
    mixture <- OpenMx::mxModel(group, OpenMx::mxData(data.frame(Y = seq_len(n)), type = "raw"))
    for (model in list(group, mixture)) {
      out <- getSummary(list(model), HetModels = TRUE)
      columns <- paste0("Y_res_c", seq_len(n))
      expect_true(all(columns %in% names(out)))
      expect_equal(unname(unlist(out[1, columns])), as.numeric(seq_len(n)))
      expect_equal(sum(grepl("^Y_res_c[0-9]+$", names(out))), n)
      expect_false(any(grepl("^[0-9]+Y_res|[.]x$|[.]y$", names(out))))
    }
  }
})

.with_output_optimizer <- function(code) {
  old_optimizer <- OpenMx::mxOption(NULL, "Default optimizer")
  old_threads <- OpenMx::mxOption(NULL, "Number of Threads")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    OpenMx::mxOption(NULL, "Default optimizer", old_optimizer)
    OpenMx::mxOption(NULL, "Number of Threads", old_threads)
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  OpenMx::mxOption(NULL, "Default optimizer", "CSOLNP")
  OpenMx::mxOption(NULL, "Number of Threads", 1L)
  force(code)
}

test_that("one-class mixtures return parameters without nonexistent free weights", {
  # Fast, well-conditioned public fits exercise actual output extraction.
  .with_output_optimizer({
  set.seed(20261024)
  n <- 80L
  dat <- data.frame(ID = seq_len(n))
  intercept <- rnorm(n, 10, 2)
  slope <- rnorm(n, 1.5, .3)
  for (j in seq_len(5L)) {
    dat[[paste0("T", j)]] <- j - 1
    dat[[paste0("Y", j)]] <- intercept + slope * (j - 1) + rnorm(n, sd = .5)
  }
  set.seed(20261024)
  one <- suppressMessages(getMIX(dat, prop_starts = 1, sub_Model = "LGCM", y_var = "Y",
    t_var = "T", records = 1:5, curveFun = "LIN", intrinsic = FALSE, paramOut = TRUE))
  expect_equal(one@mxOutput@output$status$code, 0)
  expect_equal(nrow(one@Estimates), 6)
  expect_true(all(is.finite(one@Estimates$Estimate)))
  expect_true(all(is.finite(one@Estimates$SE)))
  expect_false(any(grepl("^p[0-9]+$", one@Estimates$Name)))
  expect_false(any(grepl("^w[0-9]+$", names(OpenMx::omxGetParameters(one@mxOutput)))))

  for (j in seq_len(5L)) dat[[paste0("Y", j)]][seq_len(n / 2)] <-
    dat[[paste0("Y", j)]][seq_len(n / 2)] + 20
  set.seed(20261024)
  two <- suppressMessages(getMIX(dat, prop_starts = c(.5, .5), sub_Model = "LGCM", y_var = "Y",
    t_var = "T", records = 1:5, curveFun = "LIN", intrinsic = FALSE, paramOut = TRUE))
  expect_equal(two@mxOutput@output$status$code, 0)
  weights <- two@Estimates[two@Estimates$Name == "p2", ]
  expect_equal(nrow(weights), 1L)
  expect_equal(weights$Estimate, unname(two@mxOutput@output$estimate["w2"]))
  expect_equal(weights$SE, unname(two@mxOutput@output$standardErrors["w2", 1]))
  })
})

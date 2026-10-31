# Multivariate parameter identity: between-process labels join the prefixes of two processes, so prefixes
# that make a joined pair equal a process prefix or another pair's are rejected, and estimates are taken
# from the fitted model by their exact labels rather than by partial matches.

.mgm_name_data <- function(y_var, n = 30) {
  dat <- data.frame(id = seq_len(n), g = rep(1:2, length.out = n))
  for (j in 1:10) dat[[paste0("T", j)]] <- j - 1
  for (p in y_var) for (j in 1:10) dat[[paste0(p, j)]] <- j + rnorm(n)
  dat
}

test_that("prefixes that give ambiguous multivariate labels are rejected with the conflicting meanings", {
  expect_error(nlpsem:::validate_mgm_names(c("A", "B", "AB")),
               'process "AB" and pair \\("A", "B"\\) both use the identifier "AB"')
  expect_error(nlpsem:::validate_mgm_names(c("A", "BC", "AB", "C")),
               'pairs \\("A", "BC"\\) and \\("AB", "C"\\) both use the identifier "ABC"')
  expect_error(nlpsem:::validate_mgm_names(c("Y", "y"), "LCSM"), "more than letter case")
  expect_invisible(nlpsem:::validate_mgm_names(c("Y", "y"), "LGCM"))
  expect_invisible(nlpsem:::validate_mgm_names(c("Y", "MY")))
  expect_invisible(nlpsem:::validate_mgm_names(c("Y", "XY", "Z")))
  expect_invisible(nlpsem:::validate_mgm_names(c("R", "M")))
})

test_that("every multivariate entry point rejects ambiguous prefixes before fitting or simulating", {
  y_var <- c("A", "B", "AB")
  dat <- .mgm_name_data(y_var)
  common <- list(dat = dat, t_var = rep("T", 3), y_var = y_var, curveFun = "LIN", intrinsic = FALSE,
                 records = rep(list(1:10), 3), y_model = "LGCM")
  pattern <- "Ambiguous multivariate names"
  expect_error(do.call(getMGM, common), pattern)
  expect_error(do.call(getMGM, c(common, list(paramOut = TRUE))), pattern)
  expect_error(do.call(getMGM, c(common, list(starts = list()))), pattern)
  op <- options(nlpsem.standardize = FALSE)
  expect_error(do.call(getMGM, common), pattern)
  options(op)
  expect_error(do.call(getMGroup, c(common, list(grp_var = "g", sub_Model = "MGM"))), pattern)
  expect_error(do.call(getMIX, c(common, list(prop_starts = c(0.5, 0.5), sub_Model = "MGM"))), pattern)
  expect_error(nlpsem:::getMGM.mxModel(dat = dat, t_var = rep("T", 3), y_var = y_var, curveFun = "LIN", intrinsic = FALSE,
                                       records = rep(list(1:10), 3), y_model = "LGCM", starts = NULL), pattern)
  # the simulator rejects them without drawing random numbers
  set.seed(1); before <- .Random.seed
  expect_error(getSimData(n = 20, sub_Model = "MGM", curveFun = "LIN", y_var = y_var), pattern)
  expect_error(getSimData(n = 20, sub_Model = "MGM", curveFun = "QUAD", y_var = c("A", "BC", "AB", "C")), pattern)
  expect_error(getSimData(n = 20, sub_Model = "MGM", curveFun = "NonP", y_var = c("Y", "y")), "letter case")
  expect_identical(.Random.seed, before)
})

test_that("multivariate output takes each residual from its exact label", {
  skip_on_cran()
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  check <- function(fit, pairs) {
    est <- fit@Estimates; par <- OpenMx::omxGetParameters(fit@mxOutput)
    se <- stats::setNames(fit@mxOutput@output$standardErrors[, 1], rownames(fit@mxOutput@output$standardErrors))
    for (nm in names(pairs)) {
      expect_equal(est$Estimate[est$Name == nm], unname(par[[pairs[[nm]]]]), info = nm)
      expect_equal(est$SE[est$Name == nm], unname(se[[pairs[[nm]]]]), info = nm)
    }
  }
  # "Y_residuals" is part of "MY_residuals"
  dat <- getSimData(n = 100, sub_Model = "MGM", curveFun = "LIN", y_var = c("Y", "MY"), seed = 3)
  set.seed(20191029)
  fit <- suppressMessages(getMGM(dat = dat, t_var = c("T", "T"), y_var = c("Y", "MY"), curveFun = "LIN", intrinsic = FALSE,
                                 records = list(1:10, 1:10), y_model = "LGCM", paramOut = TRUE))
  expect_equal(nrow(fit@Estimates), 17)
  check(fit, list(Y_res = "Y_residuals", MY_res = "MY_residuals", YMY_res = "YMY_RES"))
  expect_identical(names(nlpsem:::.free_by_label(fit@mxOutput, c("MY_residuals", "Y_residuals"))$est), c("MY_residuals", "Y_residuals"))
  expect_error(nlpsem:::.free_by_label(fit@mxOutput, "Z_residuals"), "not found in the fitted model: Z_residuals")
  # "YZ_RES" is part of "XYZ_RES"
  dat <- getSimData(n = 100, sub_Model = "MGM", curveFun = "LIN", y_var = c("Y", "XY", "Z"), seed = 5)
  set.seed(20191029)
  fit <- suppressMessages(getMGM(dat = dat, t_var = rep("T", 3), y_var = c("Y", "XY", "Z"), curveFun = "LIN", intrinsic = FALSE,
                                 records = rep(list(1:10), 3), y_model = "LGCM", paramOut = TRUE))
  check(fit, list(YZ_res = "YZ_RES", XYZ_res = "XYZ_RES", YXY_res = "YXY_RES", XY_res = "XY_residuals"))
})

test_that("class-specific multivariate output and relative rates use exact labels", {
  skip_on_cran()
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  dat <- getSimData(n = 200, sub_Model = "MGM", curveFun = "NonP", y_var = c("Y", "MY"), seed = 7)
  dat$g <- rep(1:2, 100)
  set.seed(20191029)
  fit <- suppressMessages(getMGroup(dat = dat, grp_var = "g", sub_Model = "MGM", t_var = c("T", "T"), y_var = c("Y", "MY"),
                                    curveFun = "NonP", intrinsic = FALSE, records = list(1:10, 1:10), y_model = "LCSM", paramOut = TRUE))
  est <- fit@Estimates; par <- OpenMx::omxGetParameters(fit@mxOutput)
  for (k in 1:2) for (p in c("Y", "MY")) {
    expect_equal(est$Estimate[est$Name == paste0("c", k, p, "_res")], unname(par[[paste0("c", k, p, "_residuals")]]))
    rates <- est[est$Name %in% paste0("c", k, p, "_rel_rate", 2:9), ]
    expect_equal(rates$Name, paste0("c", k, p, "_rel_rate", 2:9))
    expect_equal(rates$Estimate, unname(par[paste0("c", k, p, "_rel_rate", 2:9)]))
  }
  expect_equal(est$Estimate[est$Name == "c2YMY_res"], unname(par[["c2YMY_RES"]]))
})

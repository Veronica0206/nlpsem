test_that("bootstrap intervals exclude failed and non-finite replications by parameter name", {
  raw <- data.frame(seed = 1:6, fit = 1:6,
                    statusCode = OpenMx::as.statusCode(c(0, 1, 0, 6, 0, NA)),
                    b = c(10, 20, 30, -1e6, Inf, 1e6),
                    a = c(1, 2, 3, 1e6, 4, -1e6))
  expect_warning(out <- nlpsem:::.bootstrap_percentile_CI(
    raw, OpenMx::as.statusCode(c(0, 1)), c(a = 2, b = 20), 6, .95),
    "Discarded 3 of 6")
  expect_identical(rownames(out), c("a", "b"))
  expect_equal(out$Estimate, c(2, 20))
  expect_equal(out$boot_lbound, c(1.05, 10.5))
  expect_equal(out$boot_ubound, c(2.95, 29.5))
  info <- attr(out, "bootstrap_diagnostics")
  expect_equal(unlist(info[c("requested", "attempted", "successful", "discarded")]),
               c(requested = 6, attempted = 6, successful = 3, discarded = 3))
  expect_equal(sum(info$status_counts), 6)
  expect_warning(strict <- nlpsem:::.bootstrap_percentile_CI(
    raw, OpenMx::as.statusCode(0), c(a = 2, b = 20), 6, .95), "Discarded 4")
  expect_equal(strict$boot_lbound, c(1.05, 10.5))
  expect_error(nlpsem:::.bootstrap_percentile_CI(
    raw[4:6, ], OpenMx::as.statusCode(c(0, 1)), c(a = 2, b = 20), 3, .95),
    "at least two")
})

test_that("the public bootstrap API preserves diagnostics and filters bad fits", {
  stub <- OpenMx::mxModel("stub", OpenMx::mxMatrix(
    "Full", 1, 1, free = TRUE, values = 2, labels = "a", name = "A"))
  local_mocked_bindings(mxBootstrap = function(model, replications) {
    list(compute = list(OK = OpenMx::as.statusCode(c(0, 1)), output = list(raw =
      data.frame(seed = 1:4, fit = 1:4, a = c(1, 2, 3, 1000),
                 statusCode = OpenMx::as.statusCode(c(0, 1, 0, 6))))))
  }, .package = "nlpsem")
  est <- data.frame(Name = "a", Estimate = 2, SE = .1)
  expect_warning(out <- getEstimateStats(stub, est, CI_type = "bootstrap", rep = 4),
                 "Discarded 1 of 4")
  expect_equal(out@bootstrap$boot_ubound, 2.95)
  expect_equal(attr(out@bootstrap, "bootstrap_diagnostics")$successful, 3)
})

test_that("invalid interval inputs are rejected before estimation", {
  est <- data.frame(Name = "a", Estimate = 1, SE = .1)
  bad_levels <- list(0, 1, -.5, 1.5, Inf, NA_real_, NaN, numeric(), c(.9, .95), "0.95")
  for (level in bad_levels) {
    for (type in c("Wald", "likelihood", "bootstrap", "all")) {
      expect_error(getEstimateStats(est_in = est, CI_type = type, rep = 20,
                                    conf.level = level), "between 0 and 1")
    }
    expect_error(getLatentKappa(c(1, 2), c(1, 2), conf.level = level), "between 0 and 1")
  }
  for (reps in list(0, -1, 2.5, Inf, NaN, NA_real_, numeric(), c(2, 3), "20", 2^31)) {
    expect_error(getEstimateStats(est_in = est, CI_type = "bootstrap", rep = reps), "positive integer")
    expect_error(getLRT(NULL, NULL, boot = TRUE, rep = reps), "positive integer")
  }
})

test_that("agreement matches values rather than factor positions", {
  a <- factor(c("A", "A", "B", "B"), levels = c("A", "B", "unused"))
  b <- factor(c("A", "A", "B", "B"), levels = c("B", "A"))
  expect_equal(getLatentKappa(a, b)@kappa_value, "1 (1, 1)")
  expect_equal(getLatentKappa(a, b)@judgment, "Perfect Agreement")
  # Different names are different categories; arbitrary class labels are not auto-aligned.
  expect_match(getLatentKappa(c(1, 1, 2, 2), c(2, 2, 3, 3))@kappa_value,
               "^-0.3333 ")
  expect_equal(getLatentKappa(c("A", "B", "B", NA), c("A", "B", "B", "C"))@kappa_value,
               "1 (1, 1)")
  x <- c("A", "A", "B", "B", "C", "C")
  y <- c("A", "B", "B", "C", "C", "A")
  expected <- getLatentKappa(x, y)@kappa_value
  expect_identical(getLatentKappa(factor(x, levels = c("C", "B", "A")),
                                 factor(y, levels = c("B", "A", "C")))@kappa_value, expected)
})

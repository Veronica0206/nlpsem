test_that("nonlinear initialization actually calls NLS and recovers identified curves", {
  original <- nlpsem:::.fit_nls_initialization
  calls <- list()
  testthat::local_mocked_bindings(
    .fit_nls_initialization = function(formula, data, start) {
      calls[[length(calls) + 1L]] <<- list(formula = formula, data = data, start = start)
      original(formula, data, start)
    }, .package = "nlpsem")
  fixtures <- list(
    EXP = list(t = seq(0, 6, length.out = 16), truth = c(10, 6, .5),
               mean = function(t) 10 + 6 * (1 - exp(-.5 * t))),
    JB = list(t = seq(0, 2, length.out = 12), truth = c(8, .1, -3, -2),
              mean = function(t) 8 + .1 * t - 3 * (exp(-2 * t) - 1)),
    BLS = list(t = seq(0, 6, length.out = 16), truth = c(10, 2, -.7, 2.7),
               mean = function(t) 10 + 2 * pmin(t, 2.7) - .7 * pmax(t - 2.7, 0)))
  # Small deterministic perturbations avoid the zero-residual convergence
  # limitation of the historical default NLS controls.
  for (curve in names(fixtures)) {
    fixture <- fixtures[[curve]]
    noise <- if (curve == "BLS") .01 * cos(seq_along(fixture$t)) else .0001 * sin(seq_along(fixture$t))
    y <- fixture$mean(fixture$t) + noise
    gf <- nlpsem:::getUNI.GF(as.data.frame(t(y)), as.data.frame(t(fixture$t)),
                            length(y), curve)
    audit <- attr(gf, "initialization")
    expect_identical(audit$method, "per_person_nls", info = curve)
    expect_identical(audit$algorithm, "default", info = curve)
    expect_true(audit$rows$nls_converged, info = curve)
    expect_identical(audit$rows$method, "nls", info = curve)
    expect_identical(audit$rows$attempts, 1L, info = curve)
    expect_equal(unname(as.numeric(gf)), fixture$truth, tolerance = .002, info = curve)
    expect_true(all(is.finite(attr(gf, "residuals"))), info = curve)
    expect_equal(attr(gf, "mean_res_var"), mean(attr(gf, "residuals")^2), info = curve)
  }
  expect_length(calls, 3L)
  expect_identical(nlpsem:::.nls_initialization_control(),
                   stats::nls.control(maxiter = 50L, tol = 1e-5, minFactor = 1 / 1024,
                                      printEval = FALSE, warnOnly = FALSE))
})

test_that("NLS initialization preserves ordered results and original residual positions", {
  times <- seq(0, 6, length.out = 16)
  measures <- 10 + 2 * pmin(times, 2.7) - .7 * pmax(times - 2.7, 0) +
    .01 * cos(seq_along(times))
  measures[c(3, 14)] <- NA_real_
  order <- c(9:16, 1:8)
  original <- nlpsem:::getUNI.GF(as.data.frame(t(measures)), as.data.frame(t(times)), 16, "BLS")
  reordered <- nlpsem:::getUNI.GF(as.data.frame(t(measures[order])),
                                 as.data.frame(t(times[order])), 16, "BLS")
  expect_true(attr(original, "initialization")$rows$nls_converged)
  expect_true(attr(reordered, "initialization")$rows$nls_converged)
  expect_equal(as.numeric(original), as.numeric(reordered), tolerance = 1e-10)
  expect_equal(attr(original, "residuals")[, order, drop = FALSE], attr(reordered, "residuals"),
               tolerance = 1e-10)
  expect_identical(attr(original, "initialization")$rows$n_complete, 14L)
  expect_true(all(is.na(attr(original, "residuals")[1, c(3, 14)])))
})

test_that("NLS failures are explicit audited heuristic fallbacks", {
  testthat::local_mocked_bindings(
    .fit_nls_initialization = function(...) stop("controlled NLS failure"),
    .package = "nlpsem")
  times <- 0:5
  measures <- c(0, 2, 4, 7, 9, 12)
  gf <- nlpsem:::getUNI.GF(as.data.frame(t(measures)), as.data.frame(t(times)), 6, "BLS")
  # Right segment begins at original index four (time 3), not index one
  # of the full vector as in the historical filtered-index bug.
  seed <- c(0, 2, 2.5, 2.5)
  expect_equal(as.numeric(gf), seed)
  audit <- attr(gf, "initialization")$rows
  expect_identical(audit$method, "heuristic_fallback")
  expect_false(audit$nls_converged)
  expect_identical(audit$attempts, 1L)
  expect_match(audit$fallback_reason, "controlled NLS failure", fixed = TRUE)
  expected_residuals <- measures - (seed[1] + seed[2] * pmin(times, seed[4]) +
                                      seed[3] * pmax(times - seed[4], 0))
  expect_equal(as.numeric(attr(gf, "residuals")), expected_residuals)
  expect_equal(attr(gf, "mean_res_var"), mean(expected_residuals^2))
})

test_that("incomplete and degenerate occasions are not reported as successful NLS", {
  traj <- as.data.frame(rbind(c(1, NA, NA, NA), c(1, 2, NA, NA), c(1, 1, 1, 1)))
  time <- as.data.frame(rbind(0:3, 0:3, rep(0, 4)))
  gf <- nlpsem:::getUNI.GF(traj, time, 4, "EXP")
  audit <- attr(gf, "initialization")$rows
  expect_identical(audit$method, c("unavailable", "heuristic_fallback", "unavailable"))
  expect_false(any(audit$nls_converged))
  expect_identical(audit$attempts, c(0L, 0L, 0L))
  expect_identical(audit$n_complete, c(1L, 2L, 4L))
  expect_true(all(is.na(as.matrix(gf)[c(1, 3), ])))
  expect_true(all(is.na(attr(gf, "residuals")[c(1, 3), ])))
  expect_true(all(is.finite(as.matrix(gf)[2, ])))
  expect_true(all(nzchar(audit$fallback_reason)))
})

test_that("NLS initialization neither changes SEM optimizer nor consumes RNG", {
  optimizer <- OpenMx::mxOption(NULL, "Default optimizer")
  set.seed(112358)
  seed <- .Random.seed
  times <- seq(0, 6, length.out = 16)
  measures <- 10 + 6 * (1 - exp(-.5 * times)) + .0001 * sin(seq_along(times))
  nlpsem:::getUNI.GF(as.data.frame(t(measures)), as.data.frame(t(times)), 16, "EXP")
  expect_identical(OpenMx::mxOption(NULL, "Default optimizer"), optimizer)
  expect_identical(.Random.seed, seed)
})

.diagnostics_example <- local({
  fitted <- NULL
  function() {
    if (!is.null(fitted)) return(fitted)
    old_optimizer <- OpenMx::mxOption(NULL, "Default optimizer")
    old_threads <- OpenMx::mxOption(NULL, "Number of Threads")
    on.exit({
      OpenMx::mxOption(NULL, "Default optimizer", old_optimizer)
      OpenMx::mxOption(NULL, "Number of Threads", old_threads)
    })
    OpenMx::mxOption(NULL, "Default optimizer", "CSOLNP")
    OpenMx::mxOption(NULL, "Number of Threads", 1L)
    # Deterministic Gaussian quantiles provide a well-conditioned small real fit.
    dat <- data.frame(Y = stats::qnorm((seq_len(80) - .5) / 80))
    model <- OpenMx::mxModel("diagnosticExample", type = "RAM", manifestVars = "Y", latentVars = "F",
      OpenMx::mxPath("F", "Y", arrows = 1, free = FALSE, values = 1),
      OpenMx::mxPath("F", arrows = 2, free = TRUE, values = .8, labels = "variance"),
      OpenMx::mxPath("Y", arrows = 2, free = FALSE, values = .1),
      OpenMx::mxPath("one", "F", free = TRUE, values = 0, labels = "mean"),
      OpenMx::mxData(dat, type = "raw"))
    fitted <<- suppressMessages(OpenMx::mxRun(model, silent = TRUE))
    fitted
  }
})

test_that("diagnostics distinguish numerical checks from fit acceptance without mutation", {
  model <- .diagnostics_example()
  before <- serialize(model, NULL)
  diagnostics <- getDiagnostics(model)
  expect_true(diagnostics$optimizer$status_zero)
  expect_true(diagnostics$optimizer$finite_likelihood)
  expect_true(all(diagnostics$covariance$status == "positive_definite"))
  expect_identical(diagnostics$hessian$status, "positive_definite")
  expect_identical(diagnostics$standard_errors$status, "no_issues_detected")
  expect_length(diagnostics$issues, 0L)
  expect_identical(diagnostics$repeated_starts, "not_assessed")
  expect_identical(diagnostics$component_sizes$status, "not_applicable")
  expect_identical(serialize(model, NULL), before)
  wrapped <- methods::new("myMxOutput", mxOutput = model, Estimates = data.frame())
  expect_identical(getDiagnostics(wrapped), diagnostics)
})

test_that("status zero does not hide inadmissible shipped latent covariances", {
  examples <- .load_extdata("getTVCmodel_examples.RData")
  output <- getDiagnostics(examples$Math_TVCchgBL_BLS_r)
  expect_true(output$optimizer$status_zero)
  latent <- subset(output$covariance, block == "Latent disturbances")
  expect_true(any(latent$minimum_variance < -1))
  expect_true(any(latent$status == "negative_variance"))
  original <- subset(output$covariance, matrix == "Y_psi_r")
  expect_equal(original$minimum_variance, -1.395537, tolerance = 1e-5)
  expect_equal(original$minimum_eigenvalue, -5.197491, tolerance = 1e-5)
  expect_true(any(grepl("Covariance problem", output$issues, fixed = TRUE)))
})

test_that("covariance diagnostics are scale-aware and allow deterministic zero disturbances", {
  valid <- matrix(c(4, 1, 1, 9), 2)
  scaling <- diag(c(.00001, 10000))
  a <- nlpsem:::.diagnostic_matrix(valid, 1e-8)
  b <- nlpsem:::.diagnostic_matrix(scaling %*% valid %*% scaling, 1e-8)
  expect_identical(a$status, b$status)
  expect_equal(a$scaled_minimum_eigenvalue, b$scaled_minimum_eigenvalue)
  invalid <- matrix(c(1, 2, 2, 1), 2)
  expect_identical(nlpsem:::.diagnostic_matrix(invalid, 1e-8)$status, "indefinite")
  expect_identical(nlpsem:::.diagnostic_matrix(scaling %*% invalid %*% scaling, 1e-8)$status, "indefinite")
  zeros <- nlpsem:::.diagnostic_matrix(diag(c(1, 0)), 1e-8)
  expect_identical(zeros$status, "near_boundary")
  expect_identical(zeros$negative_diagonals, "")
  expect_identical(nlpsem:::.diagnostic_matrix(matrix(0, 2, 2), 1e-8)$status, "near_boundary")
  expect_identical(nlpsem:::.diagnostic_matrix(matrix(c(1, 1, 1, 0), 2), 1e-8)$status, "indefinite")
  expect_identical(nlpsem:::.diagnostic_matrix(matrix(c(1, 0, 2, 1), 2), 1e-8)$status, "not_symmetric")
})

test_that("missing or invalid information remains explicit", {
  model <- .diagnostics_example()
  model@output <- list()
  output <- getDiagnostics(model)
  expect_true(is.na(output$optimizer$status))
  expect_true(all(output$covariance$status == "unavailable"))
  expect_identical(output$hessian$status, "unavailable")
  expect_identical(output$standard_errors$status, "unavailable")
  model <- .diagnostics_example()
  model@output$status$code <- 6L
  model@output$standardErrors[1, 1] <- NA_real_
  model@output$calculatedHessian <- diag(c(1, -1))
  model@output$infoDefinite <- FALSE
  output <- getDiagnostics(model)
  expect_false(output$optimizer$status_zero)
  expect_identical(output$standard_errors$status, "issues_detected")
  expect_false(output$hessian$openmx_info_definite)
  expect_true(any(grepl("Hessian", output$issues)))
  expect_true(any(grepl("Standard errors", output$issues)))
  unsupported <- OpenMx::mxModel("unsupported")
  expect_identical(getDiagnostics(unsupported)$covariance$status, "unavailable")
  expect_error(getDiagnostics(NULL), "myMxOutput or MxModel")
  for (tol in list(NA, NaN, Inf, 0, -1, 1, c(.1, .2), "tiny")) {
    expect_error(getDiagnostics(model, tolerance = tol), "tolerance")
  }
})

test_that("row-dependent covariances are not assessed from cached values", {
  model <- .diagnostics_example()
  model@matrices$S@labels[1, 1] <- "data.rowVariance"
  output <- getDiagnostics(model)
  expect_identical(output$covariance$status, "unavailable")
  expect_match(output$covariance$detail, "row-specific")
  model <- .diagnostics_example()
  model <- OpenMx::mxModel(model,
    OpenMx::mxMatrix("Full", 1, 1, labels = "data.rowVariance", name = "rowCovariance"),
    OpenMx::mxAlgebra(rowCovariance, name = "dependent"),
    OpenMx::mxAlgebra(dependent, name = "Y_psi0"))
  # Test transitive algebra dependency detection independently of refitting.
  expect_true(nlpsem:::.diagnostic_has_definition_variables(model, "Y_psi0"))
  expect_false(nlpsem:::.diagnostic_has_definition_variables(model, "S"))
})

test_that("jointly fitted children and mixture component sizes are inspected", {
  examples <- .load_extdata("getMIX_examples.RData")
  output <- getDiagnostics(examples$Math_BLS_LGCM3)
  expect_equal(length(unique(output$covariance$model)), 3L)
  expect_true(any(grepl("Class3", output$covariance$model)))
  sizes <- output$component_sizes$table
  expect_equal(nrow(sizes), 3L)
  expect_equal(sum(sizes$assigned_n), 500)
  expect_equal(sum(sizes$posterior_n), 500, tolerance = 1e-10)
  expect_equal(sum(sizes$assigned_proportion), 1)
  covariate <- getDiagnostics(examples$RM_BLS_PLGCM3)
  expect_identical(covariate$component_sizes$status, "available")
  expect_equal(sum(covariate$component_sizes$table$posterior_n), 500, tolerance = 1e-10)
  group <- OpenMx::mxModel("joint", OpenMx::mxRename(.diagnostics_example(), "First"),
                           OpenMx::mxRename(.diagnostics_example(), "Second"))
  group@output <- .diagnostics_example()@output
  group@submodels$Second@matrices$S@values["F", "F"] <- -1
  output <- getDiagnostics(group)
  expect_equal(length(unique(output$covariance$model)), 2L)
  expect_true(any(grepl("Second", output$covariance$model[output$covariance$status == "negative_variance"])))
  expect_identical(output$component_sizes$status, "not_applicable")
})

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
  # The shipped multiple-group random-knot example reaches status 0 with positive variances but an
  # indefinite growth-factor covariance matrix in group 1.
  examples <- .load_extdata("getMGroup_examples.RData")
  output <- getDiagnostics(examples$MGroup_Math_BLS_LGCM_f)
  expect_true(output$optimizer$status_zero)
  latent <- subset(output$covariance, block == "Latent disturbances")
  expect_true(all(latent$minimum_variance > 0))
  expect_identical(latent$status[grepl("Class1", latent$model)], "indefinite")
  expect_identical(latent$status[grepl("Class2", latent$model)], "positive_definite")
  original <- subset(output$covariance, matrix == "c1Y_psi0")
  expect_equal(original$minimum_eigenvalue, -0.005265433, tolerance = 1e-5)
  expect_true(any(grepl("Covariance problem", output$issues, fixed = TRUE)))
  # The change-from-baseline TVC example, inadmissible in its 1.0.0.9000 fit, now passes.
  expect_identical(getDiagnostics(.load_extdata("getTVCmodel_examples.RData")$Math_TVCchgBL_BLS_r)$issues, character(0))
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

# A stored-results fixture with a known joint covariance. It does not fit or
# silently recompute cached algebras; the independent joint matrix is explicit.
.diagnostics_cross_fixture <- function(prefix = "") {
  variables <- c("Reading0", "Reading1", "Math0", "Math1")
  left <- matrix(paste0(prefix, c("read00", "read01", "read01", "read11")), 2L)
  right <- matrix(paste0(prefix, c("math00", "math01", "math01", "math11")), 2L)
  cross_labels <- matrix(paste0(prefix, c("cross00", "cross10", "cross01", "cross11")), 2L)
  labels <- rbind(cbind(left, cross_labels), cbind(t(cross_labels), right))
  cross <- matrix(c(-.2, .3, .1, -.1), 2L)
  joint <- rbind(cbind(diag(2), cross), cbind(t(cross), diag(2)))
  block_algebra <- function(grid, name) {
    text <- paste0("rbind(", paste(apply(grid, 1L, function(row) {
      paste0("cbind(", paste(row, collapse = ","), ")")
    }), collapse = ","), ")")
    OpenMx::mxAlgebraFromString(text, name = name)
  }
  cross_name <- paste0(prefix, "ReadingEfficiencyMathReasoning_psi_s")
  left_name <- paste0(prefix, "ReadingEfficiency_psi_s")
  model <- OpenMx::mxModel("crossBlockExample", type = "RAM",
    manifestVars = "Y", latentVars = variables,
    OpenMx::mxPath(from = variables, to = variables, connect = "unique.pairs",
      arrows = 2, free = FALSE, values = joint[lower.tri(joint, diag = TRUE)],
      labels = labels[lower.tri(labels, diag = TRUE)]),
    OpenMx::mxPath("Y", arrows = 2, free = FALSE, values = 1),
    block_algebra(cross_labels, cross_name), block_algebra(left, left_name))
  model@algebras[[cross_name]]@result <- cross
  model@algebras[[left_name]]@result <- diag(2)
  model@output <- list(Minus2LogLikelihood = 1, status = list(code = 0L))
  list(model = model, cross_name = cross_name, left_name = left_name,
       variables = variables, joint = joint)
}

test_that("cross blocks use their joint covariance regardless of names or class prefixes", {
  for (prefix in c("", "c12")) {
    fixture <- .diagnostics_cross_fixture(prefix)
    expect_gt(min(eigen(fixture$joint, symmetric = TRUE)$values), 0)
    output <- getDiagnostics(fixture$model)
    expect_true(all(output$covariance$status == "positive_definite"))
    expect_false(fixture$cross_name %in% output$covariance$matrix)
    expect_identical(output$cross_covariances$matrix, fixture$cross_name)
    expect_identical(output$cross_covariances$status, "joint_covariance_checked")
    expect_identical(output$cross_covariances$row_variables, "Reading0, Reading1")
    expect_identical(output$cross_covariances$column_variables, "Math0, Math1")
    expect_length(output$issues, 0L)
    parent <- OpenMx::mxModel("joint", OpenMx::mxRename(fixture$model, "Class12"))
    parent@output <- fixture$model@output
    child <- getDiagnostics(parent)
    expect_identical(child$cross_covariances$model, "joint/Class12")
    expect_length(child$issues, 0L)
  }
})

test_that("valid marginal blocks cannot conceal an inadmissible joint covariance", {
  fixture <- .diagnostics_cross_fixture("c123")
  model <- fixture$model
  cross <- matrix(c(-1.2, .3, .1, -1.1), 2L)
  joint <- rbind(cbind(diag(2), cross), cbind(t(cross), diag(2)))
  expect_lt(min(eigen(joint, symmetric = TRUE)$values), 0)
  model@matrices$S@values[fixture$variables, fixture$variables] <- joint
  model@algebras[[fixture$cross_name]]@result <- cross
  output <- getDiagnostics(model)
  expect_true(all(subset(output$covariance, matrix == "S" & block != "Manifest residuals")$status == "indefinite"))
  expect_identical(subset(output$covariance, matrix == fixture$left_name)$status, "positive_definite")
  expect_true(any(grepl("Covariance problem", output$issues, fixed = TRUE)))
  expect_identical(output$cross_covariances$status, "joint_covariance_checked")

  # Numerical asymmetry alone never exempts a genuine covariance algebra.
  model <- fixture$model
  model@algebras[[fixture$left_name]]@result[1L, 2L] <- .2
  output <- getDiagnostics(model)
  expect_identical(subset(output$covariance, matrix == fixture$left_name)$status, "not_symmetric")
  expect_true(any(grepl(fixture$left_name, output$issues, fixed = TRUE)))
})

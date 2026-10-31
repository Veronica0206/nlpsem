test_that("kappa printing uses the requested confidence level without changing S4 slots", {
  x <- c(1, 1, 1, 2, 2, 2)
  y <- c(1, 1, 2, 1, 2, 2)
  for (level in c(.80, .95, .99)) {
    result <- getLatentKappa(x, y, conf.level = level)
    text <- capture.output(returned <- withVisible(printTable(result)))
    expect_true(any(grepl(paste0(100 * level, "% CI"), text, fixed = TRUE)))
    expect_identical(returned$value, result)
    expect_false(returned$visible)
    expect_setequal(methods::slotNames(result), c("kappa_value", "judgment"))
  }
  legacy <- result
  attr(legacy, "conf.level") <- NULL
  text <- capture.output(printTable(legacy))
  expect_true(any(grepl("level not recorded", text, fixed = TRUE)))
  expect_false(any(grepl("95%", text, fixed = TRUE)))
})

test_that("table methods return their objects invisibly", {
  objects <- list(
    getEstimateStats(est_in = data.frame(Name = "a", Estimate = 1, SE = .1)),
    methods::new("FSOutput", scores_est = matrix(1), scores_se = matrix(.1)),
    methods::new("postOutput", prob = matrix(1), membership = 1, entropy = NA_real_),
    .load_extdata("getLGCM_examples.RData")$Math_LGCM_BLS_r
  )
  for (object in objects) {
    invisible(capture.output(result <- withVisible(printTable(object))))
    expect_identical(result$value, object)
    expect_false(result$visible)
  }
  figure <- methods::new("figOutput", figures = list())
  invisible(capture.output(result <- withVisible(methods::show(figure))))
  expect_identical(result$value, figure)
  expect_false(result$visible)
})

test_that("factor scoring rejects unsupported model containers before calling OpenMx", {
  expect_error(getIndFS(NULL), "fitted OpenMx model")
  for (value in list(NA_character_, character(), c("ML", "Regression"), "wrong", 1)) {
    expect_error(getIndFS(NULL, value), "FS_type must be")
  }
  mixture <- .load_extdata("getMIX_examples.RData")$Math_BLS_LGCM3@mxOutput
  expect_error(getIndFS(mixture), "multiple-group and mixture parent models")
  groups <- .load_extdata("getMGroup_examples.RData")
  fit <- Filter(function(x) methods::is(x, "myMxOutput"), as.list(groups))[[1]]
  expect_error(getIndFS(fit@mxOutput), "multiple-group and mixture parent models")
})

test_that("single-factor scoring preserves matrix dimensions", {
  model <- .load_extdata("getLGCM_examples.RData")$Math_LGCM_BLS_r@mxOutput
  # The return-shape contract is independent of a costly factor-scoring fit.
  local_mocked_bindings(mxFactorScores = function(model, type) array(
    c(2, .2), c(1, 1, 2), dimnames = list("person", "factor", c("estimate", "SE"))),
    .package = "nlpsem")
  result <- getIndFS(model)
  expect_identical(dim(result@scores_est), c(1L, 1L))
  expect_identical(dim(result@scores_se), c(1L, 1L))
  expect_equal(unname(result@scores_est), matrix(2))
  expect_equal(unname(result@scores_se), matrix(.2))
})

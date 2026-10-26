# Review regressions: user-visible posterior and summary contracts.

.review_probability_model <- function(likelihood, logits, covariates = FALSE) {
  n <- nrow(likelihood)
  k <- ncol(likelihood)
  children <- lapply(seq_len(k), function(j) OpenMx::mxModel(paste0("Class", j),
    OpenMx::mxMatrix("Full", n, 1, values = likelihood[, j], name = "component")))
  dat <- data.frame(Y = seq_len(n), z = seq_len(n), z2 = seq_len(n) / 10)
  if (covariates) {
    OpenMx::mxModel("probability_probe", children, OpenMx::mxData(dat, type = "raw"),
      OpenMx::mxMatrix("Full", k, 3, values = logits, name = "classbeta"),
      OpenMx::mxMatrix("Full", 3, 1, labels = c(NA, "data.z", "data.z2"), values = 1, name = "weightsV"))
  } else {
    OpenMx::mxModel("probability_probe", children, OpenMx::mxData(dat, type = "raw"),
      OpenMx::mxMatrix("Full", k, 1, values = logits, name = "weights"))
  }
}

.review_posterior_mock <- function() {
  local_mocked_bindings(mxEval = function(object, model, ...) model$component$values,
    .package = "nlpsem", .env = parent.frame())
}

test_that("posterior normalization is stable and follows Bayes' rule", {
  .review_posterior_mock()
  likelihood <- matrix(c(.2, .4, .6, .8, .1, .1), 3, 2)
  model <- .review_probability_model(likelihood, log(c(.25, .75)))
  expected <- sweep(likelihood, 2, c(.25, .75), "*")
  expected <- expected / rowSums(expected)
  out <- getPosterior(model, 2)
  expect_equal(out@prob, expected, tolerance = 1e-14)
  expect_equal(out@membership, c(2, 1, 1))
  for (shift in c(-1000, 1000)) {
    shifted <- model
    shifted$weights$values <- shifted$weights$values + shift
    expect_equal(getPosterior(shifted, 2)@prob, expected, tolerance = 1e-12)
  }
  extreme <- .review_probability_model(matrix(c(0, 1), 1, 2), c(0, -800))
  expect_equal(getPosterior(extreme, 2)@prob, matrix(c(0, 1), 1, 2))
  extreme <- .review_probability_model(matrix(c(.1, .2, .3), 1, 3), c(0, 800, 801))
  expected <- c(0, .2 * exp(-1), .3)
  expected <- matrix(expected / sum(expected), 1, 3)
  expect_equal(getPosterior(extreme, 3)@prob, expected, tolerance = 1e-14)
})

test_that("cluster covariates are recovered and checked against fitted order", {
  .review_posterior_mock()
  likelihood <- matrix(c(.2, .4, .6, .8, .1, .1), 3, 2)
  beta <- rbind(c(800, 0, 0), c(801, .2, -.1))
  model <- .review_probability_model(likelihood, beta, covariates = TRUE)
  logits <- cbind(1, model$data$observed$z, model$data$observed$z2) %*% t(beta)
  weights <- exp(sweep(logits, 1, apply(logits, 1, max), "-"))
  expected <- likelihood * weights
  expected <- expected / rowSums(expected)
  expect_equal(getPosterior(model, 2)@prob, expected, tolerance = 1e-12)
  expect_equal(getPosterior(model, 2, cluster_TIC = c("z", "z2"))@prob, expected, tolerance = 1e-12)
  for (cols in list("z", c("z2", "z"), c("z", "absent"))) {
    expect_error(getPosterior(model, 2, cluster_TIC = cols), "original order")
  }
  model$data$observed$z[1] <- NA_real_
  expect_error(getPosterior(model, 2), "numeric and finite")
})

test_that("posterior rejects malformed counts, flags and component likelihoods", {
  .review_posterior_mock()
  model <- .review_probability_model(matrix(.5, 2, 2), c(0, 0))
  for (k in list(NA_real_, Inf, 0, -1, 1.5, c(2, 2), "2")) {
    expect_error(getPosterior(model, k), "positive integer")
  }
  for (k in c(1, 3)) expect_error(getPosterior(model, k), "number of fitted")
  expect_error(getPosterior(NULL, 2), "mxOutput")
  for (flag in list(NA, 1, logical(), c(TRUE, FALSE))) {
    expect_error(getPosterior(model, 2, label = flag), "TRUE or FALSE")
  }
  for (cols in list(NA_character_, character(), c("z", "z"), 1, "")) {
    expect_error(getPosterior(model, 2, cluster_TIC = cols), "distinct nonmissing")
  }
  expect_error(getPosterior(model, 2, cluster_TIC = "z"), "original order")
  for (bad in c(-1, NA_real_, Inf)) {
    modified <- model
    modified@submodels[[1]]$component$values[1] <- bad
    expect_error(getPosterior(modified, 2), "invalid component likelihoods")
  }
  zeros <- .review_probability_model(matrix(0, 2, 2), c(0, 0))
  expect_error(getPosterior(zeros, 2), "all class likelihoods are zero")
  model$weights$values[1] <- Inf
  expect_error(getPosterior(model, 2), "finite log weight")
})

test_that("shipped mixture probabilities retain dimensions and coherent summaries", {
  objects <- .load_extdata("getMIX_examples.RData")
  model <- objects$Math_BLS_LGCM3@mxOutput
  out <- getPosterior(model, 3)
  expect_equal(dim(out@prob), c(500L, 3L))
  expect_true(all(is.finite(out@prob)))
  expect_equal(rowSums(out@prob), rep(1, 500), tolerance = 1e-12)
  expect_equal(out@entropy, .6504447, tolerance = 1e-6)
  expect_error(getPosterior(model, 2), "number of fitted")
  expect_error(getPosterior(model, 4), "number of fitted")
  multivariate <- objects$RM_BLS_PLGCM3@mxOutput
  expect_equal(getPosterior(multivariate, 3)@prob,
    getPosterior(multivariate, 3, cluster_TIC = c("gx1", "gx2"))@prob)
  out <- getSummary(list(multivariate), HetModels = TRUE)
  expect_equal(grep("^%Class", names(out), value = TRUE), paste0("%Class", 1:3))
  expect_false(any(grepl("[.]x$|[.]y$", names(out))))
  expect_equal(unlist(out[paste0("%Class", 1:3)], use.names = FALSE), c("22.2%", "52.6%", "25.2%"))
  mixed <- getSummary(list(model, multivariate), HetModels = TRUE)
  expect_equal(nrow(mixed), 2L)
  expect_identical(mixed$Model, c("Model1", "Model2"))
  expect_true(is.na(mixed$R_res_c1[1]))
})

test_that("three-process group summaries have one proportion per group", {
  local_mocked_bindings(AIC = function(...) 120, BIC = function(...) 140,
    .package = "nlpsem")
  children <- lapply(1:3, function(k) OpenMx::mxModel(paste0("Class", k),
    OpenMx::mxMatrix("Full", 3, 1, values = c(k, 2*k, 3*k), free = TRUE,
      labels = paste0("c", k, c("Y_residuals", "M_residuals", "X_residuals")), name = "residuals"),
    OpenMx::mxData(data.frame(Y = seq_len(k)), type = "raw")))
  model <- OpenMx::mxModel("three_process_group", children)
  model@output <- list(Minus2LogLikelihood = 100)
  out <- getSummary(list(model), HetModels = TRUE)
  expect_equal(sum(grepl("^%Class", names(out))), 3L)
  expect_equal(sum(grepl("_res_c", names(out))), 9L)
  expect_equal(out$X_res_c3, 9)
  expect_equal(as.numeric(sub("%", "", unlist(out[paste0("%Class", 1:3)]))), 100*c(1, 2, 3)/6)
  expect_error(getSummary(list()), "nonempty list")
  expect_error(getSummary(list(NULL)), "nonempty list")
  expect_error(getSummary(list(model), HetModels = NA), "TRUE or FALSE")
})

test_that("a fitted one-class mixture supports posterior summary and plotting", {
  old_optimizer <- OpenMx::mxOption(NULL, "Default optimizer")
  old_threads <- OpenMx::mxOption(NULL, "Number of Threads")
  on.exit({
    OpenMx::mxOption(NULL, "Default optimizer", old_optimizer)
    OpenMx::mxOption(NULL, "Number of Threads", old_threads)
  }, add = TRUE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(20261024)
  OpenMx::mxOption(NULL, "Default optimizer", "CSOLNP")
  OpenMx::mxOption(NULL, "Number of Threads", 1L)
  n <- 80L
  dat <- data.frame(ID = seq_len(n))
  intercept <- rnorm(n, 10, 2)
  slope <- rnorm(n, 1.5, .3)
  for (j in 1:5) {
    dat[[paste0("T", j)]] <- j - 1
    dat[[paste0("Y", j)]] <- intercept + slope * (j - 1) + rnorm(n, sd = .5)
  }
  fit <- suppressMessages(getMIX(dat, prop_starts = 1, sub_Model = "LGCM", y_var = "Y",
    t_var = "T", records = 1:5, curveFun = "LIN", intrinsic = FALSE, paramOut = TRUE))
  expect_equal(fit@mxOutput@output$status$code, 0)
  posterior <- getPosterior(fit@mxOutput, 1)
  expect_equal(posterior@prob, matrix(1, n, 1))
  expect_identical(posterior@membership, rep(1, n))
  expect_identical(posterior@entropy, NA_real_)
  summary <- getSummary(list(fit@mxOutput), HetModels = TRUE)
  expect_identical(summary[["%Class1"]], "100%")
  figure <- plot(fit, draw = FALSE)
  expect_s4_class(figure, "figOutput")
  expect_true(length(figure@figures) > 0)
})

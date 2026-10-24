# Regression tests for post-processing functions, using synthetic inputs and the
# fitted example objects shipped in inst/extdata. The extdata objects were fitted
# under an earlier OpenMx version, so loading them and touching their slots can
# emit version warnings; those are suppressed where expected.

.load_extdata <- function(file) {
  path <- system.file("extdata", file, package = "nlpsem")
  skip_if(path == "", paste(file, "not available"))
  e <- new.env()
  suppressWarnings(load(path, envir = e))
  e
}

test_that("getEstimateStats computes Wald CIs and formatted p-values", {
  est <- data.frame(Name = c("a", "b", "c"),
                    Estimate = c(1, 0.5, 0),
                    SE = c(0.1, 0.25, 0.2))
  out <- getEstimateStats(est_in = est, CI_type = "Wald")
  expect_s4_class(out, "StatsOutput")
  wald <- out@wald
  expect_equal(wald$wald_lbound, est$Estimate - qnorm(0.975) * est$SE)
  expect_equal(wald$wald_ubound, est$Estimate + qnorm(0.975) * est$SE)
  expect_equal(wald["a", "p.value"], "<0.0001")
  expect_equal(wald["b", "p.value"], "0.0455")
  expect_equal(wald["c", "p.value"], ">0.9999")
})

test_that("getEstimateStats rejects an unrecognized CI_type informatively", {
  est <- data.frame(Name = "a", Estimate = 1, SE = 0.1)
  expect_error(getEstimateStats(est_in = est, CI_type = "wald"), "should be one of")
})

test_that("getEstimateStats requires a positive integer rep for bootstrap CIs", {
  est <- data.frame(Name = "a", Estimate = 1, SE = 0.1)
  stub_model <- OpenMx::mxModel("stub")
  expect_error(getEstimateStats(model = stub_model, est_in = est, CI_type = "bootstrap"),
               "positive integer")
  expect_error(getEstimateStats(model = stub_model, est_in = est, CI_type = "all", rep = -5),
               "positive integer")
})

test_that("getEstimateStats requires a model for likelihood CIs", {
  est <- data.frame(Name = "a", Estimate = 1, SE = 0.1)
  expect_error(getEstimateStats(est_in = est, CI_type = "likelihood"), "must be provided")
})

test_that("getSummary reports fit statistics for fitted LGCM examples", {
  e <- .load_extdata("getLGCM_examples.RData")
  models <- list(e$Math_LGCM_BLS_f@mxOutput, e$Math_LGCM_BLS_r@mxOutput)
  s <- suppressWarnings(getSummary(model_list = models, HetModels = FALSE))
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 2)
  expect_true(all(c("Model", "No_Params", "-2ll", "AIC", "BIC") %in% names(s)))
  expect_equal(s[1, "-2ll"],
               round(suppressWarnings(models[[1]]@output$Minus2LogLikelihood), 2))
  expect_true(all(s$AIC >= s[["-2ll"]]))
  expect_true(all(s$BIC >= s$AIC))
})

test_that("getPosterior returns coherent probabilities for a 3-class mixture", {
  e <- .load_extdata("getMIX_examples.RData")
  mx <- e$Math_BLS_LGCM3@mxOutput
  post <- suppressWarnings(getPosterior(model = mx, nClass = 3))
  expect_s4_class(post, "postOutput")
  n_obs <- suppressWarnings(nrow(mx$data$observed))
  expect_equal(dim(post@prob), c(n_obs, 3))
  expect_equal(unname(rowSums(post@prob)), rep(1, n_obs), tolerance = 1e-8)
  expect_true(all(post@membership %in% 1:3))
  expect_true(post@entropy <= 1)
})

test_that("getLatentKappa handles perfect agreement and input errors", {
  labels <- c(1, 1, 2, 2, 1, 2)
  kap <- getLatentKappa(label1 = labels, label2 = labels)
  expect_s4_class(kap, "KappaOutput")
  expect_equal(kap@judgment, "Perfect Agreement")
  expect_match(kap@kappa_value, "^1 ")
  expect_error(getLatentKappa(label1 = c(1, 2), label2 = c(1, 2, 1)), "Lengths")
  expect_error(getLatentKappa(label1 = c(1, 1, 1), label2 = c(1, 1, 1)), "2 levels")
})

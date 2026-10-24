# A small real fit verifies the public methods without relying only on old fixtures.
.fit_methods_example <- local({
  fitted <- NULL
  function(){
    if (!is.null(fitted)) return(fitted)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    old_optimizer <- OpenMx::mxOption(NULL, "Default optimizer")
    old_threads <- OpenMx::mxOption(NULL, "Number of Threads")
    on.exit({
      if (had_seed) assign(".Random.seed", seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
      OpenMx::mxOption(NULL, "Default optimizer", old_optimizer)
      OpenMx::mxOption(NULL, "Number of Threads", old_threads)
    })
    OpenMx::mxOption(NULL, "Default optimizer", "CSOLNP")
    OpenMx::mxOption(NULL, "Number of Threads", 1L)
    set.seed(5374)
    n <- 120L
    dat <- data.frame(ID = seq_len(n))
    intercept <- rnorm(n, 10, 2)
    slope <- rnorm(n, 1, .3)
    for (j in 1:4){
      dat[[paste0("Time", j)]] <- (j - 1) + runif(n, -.08, .08)
      dat[[paste0("Outcome", j)]] <- intercept + slope * dat[[paste0("Time", j)]] + rnorm(n, 0, .5)
    }
    fitted <<- suppressMessages(getLGCM(dat, t_var = "Time", y_var = "Outcome", curveFun = "LIN",
                                        intrinsic = FALSE, records = 1:4, paramOut = TRUE))
    fitted
  }
})

test_that("standard print and summary methods work on a real fitted object", {
  fit <- .fit_methods_example()
  expect_equal(fit@mxOutput@output$status$code, 0)
  expect_setequal(methods::slotNames(fit), c("mxOutput", "Estimates"))
  text <- capture.output(returned <- withVisible(print(fit)))
  expect_true(any(grepl("nlpsem fit:", text, fixed = TRUE)))
  expect_true(any(grepl("Observations: 120", text, fixed = TRUE)))
  expect_true(any(grepl("Optimizer status: 0", text, fixed = TRUE)))
  expect_false(returned$visible)
  expect_identical(returned$value, fit)
  shown <- capture.output(returned_show <- withVisible(methods::show(fit)))
  expect_identical(shown, text)
  expect_false(returned_show$visible)
  expect_identical(returned_show$value, fit)
  expect_equal(summary(fit), summary(fit@mxOutput))
  expect_equal(summary(fit, verbose = TRUE), summary(fit@mxOutput, verbose = TRUE))
  expect_error(print(fit, digits = NA_real_), "digits")
})

test_that("coef returns named unrounded free estimates independently of the output table", {
  fit <- .fit_methods_example()
  expected <- OpenMx::omxGetParameters(fit@mxOutput, free = TRUE)
  expect_identical(coef(fit), expected)
  expect_true(is.numeric(coef(fit)))
  expect_true(length(names(coef(fit))) == length(coef(fit)))
  expect_true(any(abs(coef(fit) - round(coef(fit), 4)) > 1e-8))
  expect_equal(length(coef(fit)), length(fit@mxOutput@output$estimate))
  without_table <- fit
  without_table@Estimates <- data.frame()
  expect_identical(coef(without_table), expected)
})

test_that("logLik, AIC and BIC preserve OpenMx counts for every container", {
  single <- .fit_methods_example()
  grouped <- .load_extdata("getMGroup_examples.RData")$MGroup_Math_BLS_LGCM_f
  mixture <- .load_extdata("getMIX_examples.RData")$Math_BLS_LGCM3
  for (fit in list(single, grouped, mixture)){
    ll <- logLik(fit)
    expected <- logLik(fit@mxOutput)
    expect_s3_class(ll, "logLik")
    expect_identical(ll, expected)
    expect_equal(attr(ll, "df"), length(coef(fit)))
    n <- if (length(fit@mxOutput@data$observed)) nrow(fit@mxOutput@data$observed) else
      sum(vapply(fit@mxOutput@submodels, function(x) nrow(x@data$observed), integer(1)))
    expect_equal(attr(ll, "nobs"), n)
    expect_equal(AIC(fit), AIC(fit@mxOutput))
    expect_equal(BIC(fit), BIC(fit@mxOutput))
    expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * length(coef(fit)))
    expect_equal(BIC(fit), -2 * as.numeric(ll) + log(n) * length(coef(fit)))
  }
  # Standard multi-model comparisons use a single logLik object per fit.
  same_fit <- single
  expect_equal(AIC(single, same_fit)$AIC, rep(AIC(single), 2))
  expect_equal(BIC(single, same_fit)$BIC, rep(BIC(single), 2))
})

test_that("plot uses saved model structure and only permits display changes", {
  fit <- .fit_methods_example()
  out <- withVisible(plot(fit, draw = FALSE))
  expect_false(out$visible)
  expect_s4_class(out$value, "figOutput")
  reference <- getFigure(fit@mxOutput, sub_Model = "LGCM", y_var = "Outcome", curveFun = "LIN",
                         t_var = "Time", records = 1:4, xstarts = 0, outcome = "Outcome")
  expect_equal(out$value@figures[[1]][[1]]$data, reference@figures[[1]][[1]]$data)
  shifted <- plot(fit, draw = FALSE, xstarts = 60, xlab = "Month", outcome = "Mathematics")
  expect_equal(shifted@figures[[1]][[1]]$data$time, out$value@figures[[1]][[1]]$data$time + 60)
  expect_equal(shifted@figures[[1]][[1]]$data$Y_hat, out$value@figures[[1]][[1]]$data$Y_hat)
  expect_equal(shifted@figures[[1]][[1]]$labels$x, "Month")
  expect_equal(shifted@figures[[1]][[1]]$labels$y, "Measurement of Mathematics")
  expect_error(plot(fit, draw = FALSE, curveFun = "QUAD"), "saved fit specification")
  expect_error(plot(fit, draw = FALSE, records = 1:3), "saved fit specification")
  expect_error(plot(fit, draw = FALSE, y_var = "Other"), "saved fit specification")
  expect_error(plot(fit, draw = NA), "draw")
  expect_error(plot(fit, draw = FALSE, xstarts = NA_real_), "xstarts")
  expect_error(plot(fit, 2, draw = FALSE), "by name")
})

test_that("plot draws by default and returns the figure invisibly", {
  fit <- .fit_methods_example()
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  output <- capture.output(result <- withVisible(plot(fit)))
  expect_true(any(grepl("figOutput", output, fixed = TRUE)))
  expect_false(result$visible)
  expect_s4_class(result$value, "figOutput")
})

test_that("legacy serialized fits retain inspection and explicit plotting support", {
  fit <- .load_extdata("getLGCM_examples.RData")$Math_LGCM_BLS_r
  expect_null(getSpecification(fit))
  expect_output(print(fit), "nlpsem fit")
  expect_equal(summary(fit), summary(fit@mxOutput))
  expect_identical(coef(fit), OpenMx::omxGetParameters(fit@mxOutput, free = TRUE))
  expect_identical(logLik(fit), logLik(fit@mxOutput))
  expect_error(plot(fit, draw = FALSE), "legacy fit has no saved specification")
  expect_error(plot(fit, sub_Model = "LGCM", draw = FALSE), "y_var, curveFun, t_var, records")
  expect_error(plot(fit, model = fit@mxOutput, draw = FALSE), "Unknown plotting argument")
  fig <- plot(fit, sub_Model = "LGCM", y_var = "M", curveFun = "BLS", t_var = "T",
              records = 1:9, xstarts = 0, draw = FALSE)
  expect_s4_class(fig, "figOutput")
  expect_equal(nrow(fig@figures[[1]][[1]]$layers[[4]]$data), nrow(fit@mxOutput@data$observed) * 9)
})

test_that("fit methods and saved plotting specifications survive saveRDS reload", {
  fit <- .fit_methods_example()
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  recovered <- readRDS(path)
  unlink(path)
  expect_identical(getSpecification(recovered), getSpecification(fit))
  expect_identical(coef(recovered), coef(fit))
  expect_identical(logLik(recovered), logLik(fit))
  expect_equal(summary(recovered), summary(fit))
  expect_equal(plot(recovered, draw = FALSE)@figures[[1]][[1]]$data,
                plot(fit, draw = FALSE)@figures[[1]][[1]]$data)
})

test_that("real multigroup fits print original labels and plot without repeated structure", {
  dat <- as.data.frame(.fit_methods_example()@mxOutput@data$observed)
  dat$arm <- rep(c("Control", "Program"), length.out = nrow(dat))
  grouped <- suppressMessages(getMGroup(dat, grp_var = "arm", sub_Model = "LGCM",
    t_var = "Time", y_var = "Outcome", curveFun = "LIN", intrinsic = FALSE,
    records = 1:4, paramOut = TRUE))
  expect_equal(grouped@mxOutput@output$status$code, 0)
  output <- capture.output(print(grouped))
  expect_true(any(grepl("1 = Control; 2 = Program", output, fixed = TRUE)))
  figure <- plot(grouped, draw = FALSE)
  expect_s4_class(figure, "figOutput")
  expect_setequal(unique(figure@figures[[1]][[1]]$data$Class), 1:2)
  expect_equal(attr(logLik(grouped), "nobs"), nrow(dat))
  expect_equal(AIC(grouped), AIC(grouped@mxOutput))
  expect_equal(BIC(grouped), BIC(grouped@mxOutput))
})

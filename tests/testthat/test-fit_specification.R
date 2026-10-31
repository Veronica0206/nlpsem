# Phase C fit provenance: preserve the existing S4 layout and attach portable metadata.
.fit_specification_data <- function(n = 12L) {
  records <- c(2L, 4L, 6L, 8L)
  dat <- data.frame(baseline = seq_len(n) / n, covariate = cos(seq_len(n)),
                    group = rep(c("control", "treated"), length.out = n))
  for (j in seq_along(records)) {
    dat[[paste0("Time", records[j])]] <- rep(j - 1L, n)
    dat[[paste0("Score", records[j])]] <- seq_len(n) + j + sin(seq_len(n) * j)
    dat[[paste0("Mediator", records[j])]] <- seq_len(n) / 2 + j + cos(seq_len(n) * j)
  }
  rownames(dat) <- paste0("person-", seq_len(n))
  dat
}

.fit_specification_mock_fitters <- function(.env = parent.frame()) {
  builder <- function(...) OpenMx::mxModel("specification_stub", type = "RAM")
  # the stubs are never fitted, so the fit-time diagnostics warning is not exercised here
  testthat::local_mocked_bindings(
    getLGCM.mxModel = builder, getLCSM.mxModel = builder,
    getTVC.mxModel = builder, getMGM.mxModel = builder,
    getMED.mxModel = builder, getMGroup.mxModel = builder,
    getMIX.mxModel = builder, mxRun = function(model) model,
    mxTryHard = function(model, ...) model, .warn_fit_issues = function(object) object,
    .package = "nlpsem", .env = .env
  )
}

.fit_specification_cases <- function() {
  uni <- list(t_var = "Time", y_var = "Score", curveFun = "LIN", intrinsic = FALSE,
              records = c(2L, 4L, 6L, 8L))
  list(
    getLGCM = uni,
    getLCSM = modifyList(uni, list(curveFun = "NonP")),
    getTVCmodel = c(uni, list(y_model = "LGCM", TVC = "Mediator", decompose = 0L,
                             growth_TIC = "covariate")),
    getMGM = list(t_var = c("Time", "Time"), y_var = c("Score", "Mediator"),
                  curveFun = "LIN", intrinsic = FALSE,
                  records = rep(list(c(2L, 4L, 6L, 8L)), 2), y_model = "LGCM"),
    getMediation = list(t_var = c("Time", "Time"), y_var = "Score", m_var = "Mediator",
                        x_type = "baseline", x_var = "baseline", curveFun = "LIN",
                        records = rep(list(c(2L, 4L, 6L, 8L)), 2)),
    getMGroup = c(list(grp_var = "group", sub_Model = "LGCM"), uni),
    getMIX = c(list(prop_starts = c(0.4, 0.6), sub_Model = "LGCM", cluster_TIC = "covariate"), uni)
  )
}

.fit_specification_contains <- function(x, predicate) {
  if (predicate(x)) return(TRUE)
  if (is.call(x) || is.pairlist(x)) x <- as.list(x)
  if (is.list(x)) return(any(vapply(x, .fit_specification_contains, logical(1), predicate = predicate)))
  FALSE
}

.fit_specification_preserve_rng <- function(expr, fit_control = FALSE) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  if (fit_control) {
    old_optimizer <- OpenMx::mxOption(NULL, "Default optimizer")
    old_threads <- OpenMx::mxOption(NULL, "Number of Threads")
    on.exit({
      OpenMx::mxOption(NULL, "Default optimizer", old_optimizer)
      OpenMx::mxOption(NULL, "Number of Threads", old_threads)
    }, add = TRUE)
    OpenMx::mxOption(NULL, "Default optimizer", "CSOLNP")
    OpenMx::mxOption(NULL, "Number of Threads", 1L)
  }
  eval.parent(substitute(expr))
}

test_that("all seven public fitters attach a portable canonical specification", {
  .fit_specification_mock_fitters()
  dat <- .fit_specification_data()
  cases <- .fit_specification_cases()
  family <- c(getLGCM = "LGCM", getLCSM = "LCSM", getTVCmodel = "TVC",
              getMGM = "MGM", getMediation = "MED", getMGroup = "LGCM", getMIX = "LGCM")
  required <- c("schema_version", "call", "function_name", "specification", "fitting",
                "data", "versions", "diagnostics")
  supplied_starts <- list(marker = 42)
  for (fn in names(cases)) {
    fit <- suppressMessages(do.call(get(fn, envir = asNamespace("nlpsem")),
                   c(list(dat = dat, starts = supplied_starts, tries = 2L), cases[[fn]])))
    info <- getSpecification(fit)
    expect_s4_class(fit, "myMxOutput")
    expect_setequal(methods::slotNames(fit), c("mxOutput", "Estimates"))
    expect_identical(info, attr(fit, "nlpsem", exact = TRUE))
    expect_setequal(names(info), required)
    expect_identical(info$schema_version, 1L)
    expect_identical(info$function_name, fn)
    expect_true(is.call(info$call))
    expect_identical(info$specification$sub_Model, unname(family[[fn]]))
    expect_identical(info$specification$nClass, if (fn %in% c("getMGroup", "getMIX")) 2L else NULL)
    for (arg in setdiff(names(cases[[fn]]), "prop_starts")) {
      expect_true(arg %in% names(info$specification), info = paste(fn, arg))
      expect_equal(info$specification[[arg]], cases[[fn]][[arg]], info = paste(fn, arg))
    }
    expect_true(all(c("t_var", "y_var", "curveFun", "records", "sub_Model", "y_model", "nClass") %in%
                      names(info$specification)))
    expect_identical(info$fitting$starts_source, "user")
    expect_identical(info$fitting$starts, supplied_starts)
    expect_identical(info$fitting$tries, 2L)
    expect_identical(info$fitting$optimizer, OpenMx::mxOption(NULL, "Default optimizer"))
    expect_identical(info$data$original_n, nrow(dat))
    expect_identical(info$data$retained_rows, seq_len(nrow(dat)))
    expect_length(info$data$dropped_rows, 0L)
    expect_identical(info$call$dat, quote(.data))
    expect_identical(info$call$starts, quote(.starts))
    expect_identical(info$call[[1]], as.name(fn))
    expect_true(is.na(info$diagnostics$status))
    expect_true(is.na(info$diagnostics$minus2LL))
    expect_identical(info$versions$OpenMx, as.character(packageVersion("OpenMx")))
    expect_identical(info$versions$nlpsem, as.character(packageVersion("nlpsem")))
    # R.version.string, which on development builds carries no version number
    expect_identical(info$versions$R, R.version.string)
    expect_false(.fit_specification_contains(info, is.environment))
    expect_false(.fit_specification_contains(info, is.function))
    expect_false(.fit_specification_contains(info, function(x) identical(x, dat)))
    expect_identical(getSpecification(unserialize(serialize(fit, NULL))), info)
  }
})

test_that("legacy fit objects remain readable without fabricated provenance", {
  legacy <- methods::new("myMxOutput", mxOutput = OpenMx::mxModel("legacy", type = "RAM"))
  expect_null(getSpecification(legacy))
  expect_null(attr(legacy, "nlpsem", exact = TRUE))
  expect_setequal(methods::slotNames(legacy), c("mxOutput", "Estimates"))
})

test_that("metadata does not create RNG state when no RNG was used", {
  .fit_specification_mock_fitters()
  .fit_specification_preserve_rng({
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    fit <- getLGCM(.fit_specification_data(), "Time", "Score", "LIN", intrinsic = FALSE,
                   records = c(2L, 4L, 6L, 8L), starts = list(marker = 42))
    rng_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    expect_false(rng_exists)
    expect_null(getSpecification(fit)$fitting$RNG_state)
    expect_identical(getSpecification(fit)$fitting$RNG_kind, RNGkind())
  })
})

test_that("metadata captures pre-initialization RNG state without consuming random numbers", {
  .fit_specification_mock_fitters()
  local_mocked_bindings(getUNI.initial = function(...) list(marker = runif(1)), .package = "nlpsem")
  .fit_specification_preserve_rng({
    set.seed(3184)
    before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    expected_start <- runif(1)
    expected_after <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    assign(".Random.seed", before, envir = .GlobalEnv)
    fit <- getLGCM(.fit_specification_data(), "Time", "Score", "linear", intrinsic = FALSE,
                   records = c(2L, 4L, 6L, 8L))
    after <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    info <- getSpecification(fit)
    expect_identical(info$fitting$RNG_state, before)
    expect_identical(after, expected_after)
    expect_identical(info$fitting$starts_source, "automatic")
    expect_identical(info$fitting$starts, list(marker = expected_start))
    expect_identical(info$specification$curveFun, "LIN")
  })
})

test_that("group provenance retains original labels and original row positions", {
  .fit_specification_mock_fitters()
  dat <- .fit_specification_data()
  # Factor level order is intentional and differs from alphabetical label order.
  dat$group <- factor(dat$group, levels = c("treated", "control", "unused"))
  dat$group[c(2L, 9L)] <- NA
  fit <- suppressWarnings(suppressMessages(getMGroup(dat, grp_var = "group", sub_Model = "LGCM",
    t_var = "Time", y_var = "Score", curveFun = "LIN", intrinsic = FALSE,
    records = c(2L, 4L, 6L, 8L), starts = list(marker = 42))))
  info <- getSpecification(fit)
  keep <- which(!is.na(dat$group))
  expect_identical(info$data$original_n, nrow(dat))
  expect_identical(info$data$retained_rows, keep)
  expect_identical(info$data$dropped_rows, c(2L, 9L))
  expect_identical(info$data$row_names, rownames(dat))
  expect_identical(info$data$group_mapping$code, 1:2)
  expect_identical(as.character(info$data$group_mapping$label), c("treated", "control"))
  expect_equal(unname(info$data$group_rows), lapply(c("treated", "control"), function(label) which(dat$group == label)))
  expect_identical(info$specification$grp_var, "group")
  expect_identical(info$specification$nClass, 2L)
  expect_identical(info$specification$container, "multigroup")
  expect_identical(getSpecification(unserialize(serialize(fit, NULL))), info)
})

test_that("real fitted metadata survives serialization outside the original data environment", {
  .fit_specification_preserve_rng({
    # Only serialized objects leave this local environment; provenance cannot rely
    # on lazy promises or symbols such as local_dat or record_indices remaining live.
    bytes <- local({
      set.seed(744)
      local_dat <- .fit_specification_data(180L)
      record_indices <- c(2L, 4L, 6L, 8L)
      intercept <- rnorm(nrow(local_dat), mean = 4, sd = 1)
      slope <- rnorm(nrow(local_dat), mean = 1.2, sd = 0.3)
      for (j in seq_along(record_indices)) {
        local_dat[[paste0("Score", record_indices[j])]] <- intercept + slope * (j - 1) + rnorm(nrow(local_dat), sd = 0.4)
      }
      local_dat$group[c(6L, 17L)] <- NA
      single <- suppressMessages(getLGCM(local_dat, "Time", "Score", "LIN", intrinsic = FALSE,
        records = record_indices, paramOut = TRUE))
      grouped <- suppressWarnings(suppressMessages(getMGroup(local_dat, grp_var = "group", sub_Model = "LGCM",
        t_var = "Time", y_var = "Score", curveFun = "LIN", intrinsic = FALSE,
        records = record_indices, paramOut = TRUE)))
      serialize(list(single = single, grouped = grouped), NULL)
    })
    fits <- unserialize(bytes)
    for (fit in fits) {
      info <- getSpecification(fit)
      expect_equal(fit@mxOutput@output$status$code, 0)
      expect_equal(info$diagnostics$status, fit@mxOutput@output$status$code)
      expect_equal(info$diagnostics$minus2LL, unname(fit@mxOutput@output$Minus2LogLikelihood))
      expect_equal(info$diagnostics$n_parameters, length(OpenMx::omxGetParameters(fit@mxOutput)))
      expect_identical(info$specification$records, c(2L, 4L, 6L, 8L))
      expect_identical(info$specification$t_var, "Time")
      expect_identical(info$specification$y_var, "Score")
      expect_false(.fit_specification_contains(info, is.environment))
      expect_false(.fit_specification_contains(info, is.function))
      expect_gt(nrow(fit@Estimates), 0L)
      expect_identical(getSpecification(unserialize(serialize(fit, NULL))), info)
    }
    expect_equal(getSpecification(fits$single)$diagnostics$n_observations, 180L)
    expect_equal(getSpecification(fits$grouped)$diagnostics$n_observations, 178L)
    expect_identical(getSpecification(fits$grouped)$data$dropped_rows, c(6L, 17L))
    expect_identical(getSpecification(fits$single)$data$group_mapping, NULL)
  }, fit_control = TRUE)
})


test_that("provenance reports nonconverged optimizer status without converting it to success", {
  .fit_specification_mock_fitters()
  local_mocked_bindings(mxRun = function(model) {
    model@output$status <- list(code = 6L)
    model@output$Minus2LogLikelihood <- 123.4
    model
  }, .package = "nlpsem")
  fit <- getLGCM(.fit_specification_data(), "Time", "Score", "LIN", intrinsic = FALSE,
                 records = c(2L, 4L, 6L, 8L), starts = list(marker = 42))
  expect_identical(getSpecification(fit)$diagnostics$status, 6L)
  expect_identical(getSpecification(fit)$diagnostics$minus2LL, 123.4)
})

test_that("a new fit with diagnostic issues is returned with a classed warning", {
  builder <- function(...) OpenMx::mxModel("specification_stub", type = "RAM")
  local_mocked_bindings(getLGCM.mxModel = builder, mxRun = function(model) {
    model@output$status <- list(code = 6L)
    model@output$Minus2LogLikelihood <- 123.4
    model
  }, .package = "nlpsem")
  expect_warning(fit <- getLGCM(.fit_specification_data(), "Time", "Score", "LIN", intrinsic = FALSE,
                                records = c(2L, 4L, 6L, 8L), starts = list(marker = 42)),
                 class = "nlpsem_diagnostics_warning", regexp = "Optimizer status is 6")
  expect_s4_class(fit, "myMxOutput")
  expect_identical(getSpecification(fit)$diagnostics$status, 6L)
})

# Metadata must observe ordinary R argument evaluation, not change it.
.lazy_metadata_data <- function(){
  data.frame(Time1 = 0, Time2 = 1, Time3 = 2, Time4 = 3,
              Outcome1 = 1:5, Outcome2 = 2:6, Outcome3 = 3:7, Outcome4 = 4:8)
}

.lazy_metadata_preserve_rng <- function(expr){
  present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved <- if (present) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (present) assign(".Random.seed", saved, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  eval.parent(substitute(expr))
}

.lazy_metadata_stub <- function(.env = parent.frame()){
  testthat::local_mocked_bindings(
    getLGCM.mxModel = function(...) OpenMx::mxModel("lazy_metadata_stub", type = "RAM"),
    mxRun = function(model) model,
    .warn_fit_issues = function(object) object, .package = "nlpsem", .env = .env
  )
}

test_that("metadata never forces unused retry controls", {
  .lazy_metadata_stub()
  .lazy_metadata_preserve_rng({
    set.seed(475)
    before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    fit <- getLGCM(.lazy_metadata_data(), "Time", "Outcome", "LIN", intrinsic = FALSE,
      records = 1:4, starts = list(marker = 42), tries = NULL,
      OKStatus = stop("unused OKStatus evaluated"),
      jitterD = {runif(1); "runif"},
      loc = stop("unused loc evaluated"), scale = stop("unused scale evaluated"))
    expect_s4_class(fit, "myMxOutput")
    expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), before)
    expect_identical(getSpecification(fit)$fitting$RNG_state, before)
    expect_identical(getSpecification(fit)$fitting$starts, list(marker = 42))
  })
})

test_that("unused random control expressions do not create RNG state", {
  .lazy_metadata_stub()
  .lazy_metadata_preserve_rng({
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    fit <- getLGCM(.lazy_metadata_data(), "Time", "Outcome", "LIN", intrinsic = FALSE,
      records = 1:4, starts = list(marker = 42), tries = NULL,
      jitterD = {runif(1); "runif"}, loc = runif(1), scale = runif(1))
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    expect_null(getSpecification(fit)$fitting$RNG_state)
  })
})

test_that("user-supplied starts retain their evaluation order before retries", {
  events <- character()
  local_mocked_bindings(
    getLGCM.mxModel = function(...){
      events <<- c(events, "builder")
      OpenMx::mxModel("lazy_metadata_stub", type = "RAM")
    },
    mxTryHard = function(model, extraTries, ...){
      force(extraTries)
      events <<- c(events, "retry")
      model
    }, .warn_fit_issues = function(object) object, .package = "nlpsem"
  )
  .lazy_metadata_preserve_rng({
    set.seed(303)
    before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    expected_start <- runif(1)
    expected_tries <- sample(1:2, 1)
    expected_after <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    assign(".Random.seed", before, envir = .GlobalEnv)
    fit <- getLGCM(.lazy_metadata_data(), "Time", "Outcome", "LIN", intrinsic = FALSE,
      records = 1:4,
      starts = {events <- c(events, "starts"); list(marker = runif(1))},
      tries = {events <- c(events, "tries"); sample(1:2, 1)})
    expect_identical(events, c("starts", "builder", "tries", "retry"))
    expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), expected_after)
    info <- getSpecification(fit)
    expect_identical(info$fitting$RNG_state, before)
    expect_identical(info$fitting$starts_source, "user")
    expect_identical(info$fitting$starts, list(marker = expected_start))
    expect_identical(info$fitting$tries, expected_tries)
  })
})

test_that("automatic-start random draws precede lazy retry arguments", {
  events <- character()
  local_mocked_bindings(
    getUNI.initial = function(...){
      events <<- c(events, "automatic starts")
      list(marker = runif(1))
    },
    getLGCM.mxModel = function(...){
      events <<- c(events, "builder")
      OpenMx::mxModel("lazy_metadata_stub", type = "RAM")
    },
    mxTryHard = function(model, extraTries, ...){
      force(extraTries)
      events <<- c(events, "retry")
      model
    }, .warn_fit_issues = function(object) object, .package = "nlpsem"
  )
  .lazy_metadata_preserve_rng({
    set.seed(517)
    before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    expected_start <- runif(1)
    expected_tries <- sample(1:2, 1)
    expected_after <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    assign(".Random.seed", before, envir = .GlobalEnv)
    fit <- getLGCM(.lazy_metadata_data(), "Time", "Outcome", "LIN", intrinsic = FALSE,
      records = 1:4, tries = {events <- c(events, "tries"); sample(1:2, 1)})
    expect_identical(events, c("automatic starts", "builder", "tries", "retry"))
    expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), expected_after)
    info <- getSpecification(fit)
    expect_identical(info$fitting$RNG_state, before)
    expect_identical(info$fitting$starts_source, "automatic")
    expect_identical(info$fitting$starts, list(marker = expected_start))
    expect_identical(info$fitting$tries, expected_tries)
  })
})

test_that("canonical curve names use factor labels rather than integer codes", {
  .lazy_metadata_stub()
  fit <- getLGCM(.lazy_metadata_data(), "Time", "Outcome", factor("EXP"), intrinsic = FALSE,
                 records = 1:4, starts = list(marker = 42))
  expect_identical(getSpecification(fit)$specification$curveFun, "EXP")
})

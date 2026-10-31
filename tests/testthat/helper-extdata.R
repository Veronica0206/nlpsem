# Shared helpers for tests that use the fitted example objects shipped in inst/extdata.

.load_extdata <- function(file) {
  path <- system.file("extdata", file, package = "nlpsem")
  testthat::skip_if(path == "", paste(file, "not available"))
  e <- new.env()
  suppressWarnings(load(path, envir = e))
  e
}

.rebaseline_rms <- function() {
  env <- new.env()
  data("RMS_dat", package = "nlpsem", envir = env)
  dat <- env$RMS_dat
  baseT <- dat$T1
  for (j in 1:9) dat[[paste0("T", j)]] <- dat[[paste0("T", j)]] - baseT
  dat
}

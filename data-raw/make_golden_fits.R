# Generate tests/testthat/golden/fit_golden.rds: reference values for the fitting-layer tests.
# Run from the package root after any intentional change to the estimation code:
#   Rscript data-raw/make_golden_fits.R
# Review the printed status codes, diagnostic issues and timings; cases that do not converge (status
# other than 0) are kept in the file but flagged, and the test file skips them. A case that converges
# to an inadmissible solution (issues > 0) is not a positive reference: list it in
# .known_inadmissible_golden in tests/testthat/test-fitting_golden.R.
pkgload::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(OpenMx))
source("tests/testthat/helper-golden.R")
dat <- .golden_data()
cases <- .golden_cases()
golden <- list()
for (nm in names(cases)) {
  t0 <- Sys.time()
  res <- tryCatch(.fit_golden_case(cases[[nm]], .golden_case_data(cases[[nm]], dat)),
                  error = function(e) list(error = conditionMessage(e)))
  res$n <- if (is.null(cases[[nm]]$n)) nrow(dat) else cases[[nm]]$n
  res$elapsed <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)
  golden[[nm]] <- res
  if (is.null(res$error)) {
    message(sprintf("%-28s status=%s  -2LL=%.3f  params=%d  issues=%d  %.1fs", nm, res$status, res$minus2LL,
                    length(res$estimates), length(res$issues), res$elapsed))
  } else {
    message(sprintf("%-28s ERROR: %s", nm, res$error))
  }
}
attr(golden, "session") <- list(date = format(Sys.time(), "%Y-%m-%d"), R = R.version.string,
                                OpenMx = as.character(packageVersion("OpenMx")),
                                nlpsem = as.character(packageVersion("nlpsem")), n = nrow(dat))
dir.create("tests/testthat/golden", showWarnings = FALSE, recursive = TRUE)
saveRDS(golden, "tests/testthat/golden/fit_golden.rds", version = 2)
message("written tests/testthat/golden/fit_golden.rds; total ", round(sum(sapply(golden, `[[`, "elapsed"))), " s")

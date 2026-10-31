# Invoke with Rscript -e 'source("data-raw/run_validation_tests.R")' <source> <filter> <output-prefix>.
# The source() call parses this complete driver before any long fit starts.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
Sys.setenv(NOT_CRAN = "true")
OpenMx::mxOption(NULL, "Number of Threads", 1L)
if (args[2] == "all-except-golden") {
  files <- list.files(file.path(args[1], "tests/testthat"), pattern = "^test-.*[.]R$")
  ids <- sub("^test-", "", sub("[.]R$", "", files))
  args[2] <- paste0("^(", paste(setdiff(ids, "fitting_golden"), collapse = "|"), ")$")
}
started <- Sys.time()
results <- testthat::test_local(args[1], filter = args[2], reporter = "summary", stop_on_failure = FALSE)
tab <- as.data.frame(results)
saveRDS(tab, paste0(args[3], ".rds"))
utils::write.csv(tab[, setdiff(names(tab), "result")], paste0(args[3], ".csv"), row.names = FALSE)
cat("\nTotals:\n")
print(colSums(tab[, c("failed", "skipped", "error", "warning", "passed")]))
cat("Elapsed seconds:", as.numeric(difftime(Sys.time(), started, units = "secs")), "\n")
print(sessionInfo())
if (any(tab$failed > 0 | tab$error | tab$warning > 0)) quit(status = 1L)

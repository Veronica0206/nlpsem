# Run from the package directory:
# Rscript data-raw/run_portable_qualification.R <output-dir> [replicates=50] [seed=20261025] [source-dir=.]
# Requires pkgload in addition to package dependencies. New independent synthetic
# evidence; does not replace historical references or establish full V1 qualification.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 4L) stop("Usage: output-dir [replicates=50] [seed=20261025] [source-dir=.]")
out_dir <- args[1]
replicates <- if (length(args) >= 2L) suppressWarnings(as.numeric(args[2])) else 50
seed <- if (length(args) >= 3L) suppressWarnings(as.numeric(args[3])) else 20261025
source_dir <- normalizePath(if (length(args) >= 4L) args[4] else ".", mustWork = TRUE)
if (!is.finite(replicates) || replicates != floor(replicates) || !(replicates == 0 || (replicates >= 20 && replicates <= 1000))) stop("replicates must be zero (smoke only) or an integer 20 through 1000")
if (!is.finite(seed) || seed < 1 || seed != floor(seed) || seed > .Machine$integer.max - 20000) stop("Invalid seed")
if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload to load the reviewed source")
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Use a new or empty output directory to preserve earlier qualification evidence")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)
pkgload::load_all(source_dir, quiet = TRUE)
OpenMx::mxOption(NULL, "Default optimizer", "CSOLNP")
OpenMx::mxOption(NULL, "Number of Threads", 1L)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
started <- Sys.time()
# Criteria are recorded before fitting. With 50 replicates, nominal 95% coverage
# has MCSE sqrt(.95*.05/50)=.031: this is a coarse screen, not precise calibration.
protocol <- list(scope = "Seven basic fitting-family smoke cases plus limited linear LGCM mean recovery only",
  generator = "Independent Gaussian effects and explicit measurement equations; no package generation helpers",
  seed = seed, replicates = replicates, optimizer = "CSOLNP", threads = 1,
  smoke_n = 300L, recovery_n = 200L, waves = 6L,
  smoke_acceptance = "status 0; finite likelihood/parameters; no warnings or diagnostic issues; available covariance, positive-definite Hessian and complete positive finite SEs",
  recovery_minimum_accepted_fraction = .9,
  recovery_bias_screen = "abs(bias) <= max(practical tolerance, 3 * empirical SD / sqrt(accepted replicates))",
  practical_bias_tolerance = c(mueta0 = .10, mueta1 = .03),
  coverage_screen = "two-sided exact binomial test of coverage=.95, p>=.01 among accepted fits; failures also counted as uncovered in all-attempt coverage",
  intervals = "Normal 95% Wald intervals for true means 10 and 1.5",
  unassessed = c("nonlinear and random-knot models", "five unavailable historical references", "global optimum and multiple starts", "missingness mechanisms", "mediation/mixture inference", "nonlinear LRT/bootstrap", "all curve/container combinations"))
dput(protocol, file.path(out_dir, "protocol.R"))
files <- sort(c(file.path(source_dir, c("DESCRIPTION", "NAMESPACE")), list.files(file.path(source_dir, "R"), pattern = "[.]R$", full.names = TRUE), file.path(source_dir, "data-raw/run_portable_qualification.R")))
manifest <- data.frame(file = substring(files, nchar(source_dir) + 2L), md5 = unname(tools::md5sum(files)))
write.csv(manifest, file.path(out_dir, "source-manifest.csv"), row.names = FALSE)
writeLines(c(paste("Package", as.character(packageVersion("nlpsem"))), paste("Source manifest MD5", unname(tools::md5sum(file.path(out_dir, "source-manifest.csv")))), paste("UTC start", format(started, tz = "UTC", usetz = TRUE)), capture.output(sessionInfo())), file.path(out_dir, "environment.txt"))

# Known linear truth: E(b)=(10,1.5), var(b0)=4, var(b1)=.16,
# cov(b0,b1)=.16 and independent residual variance=.5625.
make_linear <- function(n, data_seed) {
  set.seed(data_seed)
  z <- matrix(rnorm(n * 2L), n, 2)
  b0 <- 10 + 2*z[, 1]
  b1 <- 1.5 + .4*(.2*z[, 1] + sqrt(1-.2^2)*z[, 2])
  dat <- data.frame(ID = seq_len(n))
  for (j in 1:6) {
    dat[[paste0("T", j)]] <- j-1 + runif(n, -.08, .08)
    dat[[paste0("Y", j)]] <- b0 + b1*dat[[paste0("T", j)]] + rnorm(n, sd = .75)
  }
  dat
}
make_smoke <- function(family, data_seed, n = 300L) {
  dat <- make_linear(n, data_seed)
  set.seed(data_seed + 500L)
  if (family == "TVC") {
    for (j in 1:6) {
      dat[[paste0("U", j)]] <- rnorm(n)
      dat[[paste0("Y", j)]] <- dat[[paste0("Y", j)]] + .4*dat[[paste0("U", j)]]
    }
  } else if (family == "MGM") {
    m0 <- 5+rnorm(n, sd = 1.5); m1 <- .8+rnorm(n, sd = .3)
    for (j in 1:6) dat[[paste0("M", j)]] <- m0 + m1*dat[[paste0("T", j)]] + rnorm(n, sd = .6)
  } else if (family == "MED") {
    x <- rnorm(n); m0 <- 4+.5*x+rnorm(n, sd = 1.5); m1 <- .8+.15*x+rnorm(n, sd = .3)
    y0 <- 10+.2*x+.6*m0+rnorm(n, sd = 2)
    y1 <- 1.5+.1*x+.2*m0+.3*m1+rnorm(n, sd = .4)
    dat$X <- x
    for (j in 1:6) {
      time <- dat[[paste0("T", j)]]
      dat[[paste0("M", j)]] <- m0+m1*time+rnorm(n, sd = .6)
      dat[[paste0("Y", j)]] <- y0+y1*time+rnorm(n, sd = .75)
    }
  } else if (family %in% c("MGroup", "MIX")) {
    dat$group <- rep(1:2, length.out = n)
    for (j in 1:6) dat[[paste0("Y", j)]] <- dat[[paste0("Y", j)]] + (dat$group-1)*(if (family == "MIX") 12 else 4)
  }
  dat
}
fit_family <- function(family, dat) {
  basic <- list(dat = dat, t_var = "T", y_var = "Y", records = 1:6, curveFun = "LIN", intrinsic = FALSE, paramOut = TRUE)
  if (family == "LGCM") return(do.call(getLGCM, basic))
  if (family == "LCSM") { basic$curveFun <- "NonP"; return(do.call(getLCSM, basic)) }
  if (family == "TVC") return(do.call(getTVCmodel, c(basic, list(y_model = "LGCM", TVC = "U", decompose = 0))))
  if (family == "MGM") {
    basic$t_var <- c("T", "T"); basic$y_var <- c("Y", "M"); basic$records <- list(1:6, 1:6)
    return(do.call(getMGM, c(basic, list(y_model = "LGCM"))))
  }
  if (family == "MED") {
    basic$t_var <- c("T", "T"); basic$records <- list(1:6, 1:6); basic$intrinsic <- NULL
    return(do.call(getMediation, c(basic, list(m_var = "M", x_type = "baseline", x_var = "X"))))
  }
  if (family == "MGroup") return(do.call(getMGroup, c(basic, list(grp_var = "group", sub_Model = "LGCM"))))
  if (family == "MIX") return(do.call(getMIX, c(basic, list(prop_starts = c(.5, .5), sub_Model = "LGCM"))))
  stop("Unknown family")
}
run_fit <- function(family, dat, fit_seed) {
  set.seed(fit_seed); warnings <- character(); elapsed <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(suppressMessages(fit_family(family, dat)), warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - elapsed
  if (inherits(fit, "error")) return(list(fit = NULL, diagnostics = NULL, result = data.frame(family = family, accepted = FALSE, status = NA_integer_, minus2LL = NA_real_, finite_parameters = FALSE, covariance_available = FALSE, hessian_status = "unavailable", standard_errors = "unavailable", diagnostic_issues = "fit error", warnings = paste(unique(warnings), collapse = " | "), error = conditionMessage(fit), elapsed_seconds = elapsed)))
  diagnostics <- tryCatch(getDiagnostics(fit), error = function(e) e)
  if (inherits(diagnostics, "error")) {
    return(list(fit = fit, diagnostics = NULL, result = data.frame(
      family = family, accepted = FALSE, status = fit@mxOutput@output$status$code,
      minus2LL = fit@mxOutput@output$Minus2LogLikelihood, finite_parameters = all(is.finite(coef(fit))),
      covariance_available = FALSE, hessian_status = "unavailable", standard_errors = "unavailable",
      diagnostic_issues = "diagnostic error", warnings = paste(unique(warnings), collapse = " | "),
      error = conditionMessage(diagnostics), elapsed_seconds = elapsed)))
  }
  covariance_available <- any(diagnostics$covariance$status != "unavailable")
  finite_parameters <- length(coef(fit)) > 0L && all(is.finite(coef(fit)))
  accepted <- isTRUE(diagnostics$optimizer$status_zero) && isTRUE(diagnostics$optimizer$finite_likelihood) && finite_parameters && covariance_available && !length(diagnostics$issues) && !length(warnings) && identical(diagnostics$hessian$status, "positive_definite") && identical(diagnostics$standard_errors$status, "no_issues_detected")
  list(fit = fit, diagnostics = diagnostics, result = data.frame(family = family, accepted = accepted, status = diagnostics$optimizer$status, minus2LL = fit@mxOutput@output$Minus2LogLikelihood, finite_parameters = finite_parameters, covariance_available = covariance_available, hessian_status = diagnostics$hessian$status, standard_errors = diagnostics$standard_errors$status, diagnostic_issues = paste(diagnostics$issues, collapse = " | "), warnings = paste(unique(warnings), collapse = " | "), error = "", elapsed_seconds = elapsed))
}
smoke <- list()
for (i in seq_along(families <- c("LGCM", "LCSM", "TVC", "MGM", "MED", "MGroup", "MIX"))) {
  family <- families[i]; cat("Smoke", family, "\n")
  result <- run_fit(family, make_smoke(family, seed+i), seed+100L+i)
  result$result$data_seed <- seed+i; result$result$fit_seed <- seed+100L+i
  smoke[[family]] <- result$result
  saveRDS(result[c("fit", "diagnostics")], file.path(out_dir, paste0("smoke-", family, ".rds")))
  write.csv(do.call(rbind, smoke), file.path(out_dir, "smoke-results.csv"), row.names = FALSE)
  print(result$result)
}
recovery <- list()
for (i in seq_len(replicates)) {
  result <- run_fit("LGCM", make_linear(200L, seed+1000L+i), seed+10000L+i)
  row <- result$result; row$replicate <- i; row$data_seed <- seed+1000L+i; row$fit_seed <- seed+10000L+i
  for (name in c("mueta0", "mueta1")) {
    estimate <- if (is.null(result$fit)) NULL else result$fit@Estimates[result$fit@Estimates$Name == name, ]
    row[[paste0(name, "_estimate")]] <- if (is.null(estimate) || nrow(estimate) != 1L) NA_real_ else estimate$Estimate
    row[[paste0(name, "_se")]] <- if (is.null(estimate) || nrow(estimate) != 1L) NA_real_ else estimate$SE
  }
  recovery[[i]] <- row
  write.csv(do.call(rbind, recovery), file.path(out_dir, "recovery-replicates.csv"), row.names = FALSE)
  if (i %% 10L == 0L) cat("Recovery attempted", i, "of", replicates, "\n")
}
smoke_table <- do.call(rbind, smoke)
passed <- all(smoke_table$accepted)
if (replicates) {
  tab <- do.call(rbind, recovery)
  summaries <- lapply(c("mueta0", "mueta1"), function(name) {
    truth <- c(mueta0 = 10, mueta1 = 1.5)[name]
    available <- tab$accepted & is.finite(tab[[paste0(name, "_estimate")]]) & is.finite(tab[[paste0(name, "_se")]])
    estimate <- tab[[paste0(name, "_estimate")]][available]; se <- tab[[paste0(name, "_se")]][available]; n <- length(estimate)
    bias <- if (n) mean(estimate)-truth else NA_real_
    empirical_sd <- if (n >= 2L) sd(estimate) else NA_real_
    bias_mcse <- empirical_sd/sqrt(n)
    covered <- if (n) sum(abs(estimate-truth) <= qnorm(.975)*se) else 0L
    interval <- if (n) binom.test(covered, n)$conf.int else c(NA_real_, NA_real_)
    coverage_p <- if (n) binom.test(covered, n, p = .95)$p.value else NA_real_
    tolerance <- max(protocol$practical_bias_tolerance[name], 3*bias_mcse)
    screen <- n/replicates >= protocol$recovery_minimum_accepted_fraction && is.finite(bias) && is.finite(tolerance) && abs(bias) <= tolerance && is.finite(coverage_p) && coverage_p >= .01
    data.frame(parameter = name, truth = truth, attempted = replicates, accepted_with_estimate = n, excluded = replicates-n, bias = bias, bias_mcse = bias_mcse, bias_tolerance = tolerance, empirical_sd = empirical_sd, mean_standard_error = if (n) mean(se) else NA_real_, se_to_empirical_sd = if (n >= 2L) mean(se)/empirical_sd else NA_real_, covered = covered, coverage_accepted = if (n) covered/n else NA_real_, coverage_all_attempts = covered/replicates, nominal_coverage_mcse = if (n) sqrt(.95*.05/n) else NA_real_, coverage_binomial_lower = interval[1], coverage_binomial_upper = interval[2], exact_coverage_test_p = coverage_p, limited_screen_pass = screen)
  })
  summary_table <- do.call(rbind, summaries)
  write.csv(summary_table, file.path(out_dir, "recovery-summary.csv"), row.names = FALSE)
  print(summary_table); passed <- passed && all(summary_table$limited_screen_pass)
}
writeLines(c(paste("Basic smoke accepted:", sum(smoke_table$accepted), "of", nrow(smoke_table)), paste("Limited combined screen passed:", passed), paste("Elapsed seconds:", as.numeric(difftime(Sys.time(), started, units = "secs"))), "Scope: basic family smoke and one limited linear recovery screen only.", "No claim of nonlinear qualification, global optimality, cross-platform success, or comprehensive V1 readiness."), file.path(out_dir, "SUMMARY.txt"))
if (!passed) quit(status = 1L)

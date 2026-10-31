# Matched initialization / optimizer study. This does not replace golden references.
# Rscript data-raw/run_matched_nls_validation.R prepare SOURCE NEW_INPUT_DIR [CASE1,CASE2,...]
# Rscript data-raw/run_matched_nls_validation.R fit SOURCE INPUT_DIR CASE NEW_RESULT_DIR
# Rscript data-raw/run_matched_nls_validation.R compare-preparations INPUT_A INPUT_B NEW_OUTPUT
# Rscript data-raw/run_matched_nls_validation.R compare-fits INPUT_DIR RESULTS_A RESULTS_B NEW_OUTPUT
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Specify prepare, fit, compare-preparations, or compare-fits")
mode <- args[1]
default_cases <- c("LGCM_JB_r", "LCSM_JB_f", "LCSM_JB_r",
                   "TVC_LGCM_BLS_r_d1", "MGM_LGCM_BLS_r", "MED_BLS_longitudinal", "LGCM_EXP_r")
for (package in c("digest", "testthat")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Required validation dependency: ", package)
}
sha <- function(path) digest::digest(file = path, algo = "sha256")
new_dir <- function(path) {
  if (file.exists(path)) stop("Output already exists; use a new path: ", path)
  if (!dir.create(path, recursive = TRUE)) stop("Cannot create ", path)
  normalizePath(path, mustWork = TRUE)
}
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
one <- function(x, default = NA_real_) if (length(x) == 1L) unname(x) else default
same_reference <- function(a, b, tolerance) {
  # Match the package's declared test edition regardless of the caller's cwd.
  testthat::local_edition(3)
  tryCatch({testthat::expect_equal(a, b, tolerance = tolerance); TRUE},
           error = function(e) FALSE)
}
set_protocol <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  OpenMx::mxOption(NULL, "Default optimizer", "CSOLNP")
  OpenMx::mxOption(NULL, "Number of Threads", 1L)
  set.seed(20191029)
}
load_source <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required")
  pkgload::load_all(path, quiet = TRUE)
  set_protocol()
  path
}
source_manifest <- function(src) {
  paths <- sort(c("DESCRIPTION", "NAMESPACE",
    paste0("R/", list.files(file.path(src, "R"), pattern = "[.]R$")),
    "tests/testthat/helper-golden.R", "tests/testthat/golden/fit_golden.rds",
    "data-raw/run_matched_nls_validation.R"))
  data.frame(file = paths, sha256 = vapply(file.path(src, paths), sha, character(1)),
             stringsAsFactors = FALSE, row.names = NULL)
}
environment_record <- function() list(
  R = R.version.string, platform = R.version$platform,
  OpenMx = as.character(utils::packageVersion("OpenMx")),
  nlpsem = as.character(utils::packageVersion("nlpsem")),
  stats = as.character(utils::packageVersion("stats")),
  nls_control = stats::nls.control(),
  nls_algorithm = "default (Gauss-Newton)",
  nls_function_sha256 = digest::digest(paste(deparse(body(stats::nls)), collapse = "\n"),
                                     algo = "sha256", serialize = FALSE),
  optimizer = OpenMx::mxOption(NULL, "Default optimizer"),
  threads = OpenMx::mxOption(NULL, "Number of Threads"),
  rng_kind = RNGkind(), seed = 20191029L,
  session = capture.output(sessionInfo()),
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
write_environment <- function(out) {
  record <- environment_record()
  saveRDS(record, file.path(out, "environment.rds"), version = 2)
  dput(record, file = file.path(out, "environment.R"))
}
bind_rows <- function(rows) {
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    for (nm in setdiff(columns, names(x))) x[[nm]] <- NA
    x[, columns, drop = FALSE]
  })
  do.call(rbind, rows)
}

if (mode == "prepare") {
  if (!length(args) %in% 3:4) stop("prepare SOURCE NEW_INPUT_DIR [CASE1,CASE2,...]")
  src <- load_source(args[2]); out <- new_dir(args[3])
  cases_requested <- if (length(args) == 4L) strsplit(args[4], ",", fixed = TRUE)[[1]] else default_cases
  env <- new.env(parent = globalenv())
  sys.source(file.path(src, "tests/testthat/helper-golden.R"), envir = env)
  cases <- env$.golden_cases()
  stopifnot(length(cases_requested) > 0L, !anyDuplicated(cases_requested),
            all(cases_requested %in% names(cases)))
  dat <- env$.golden_data()
  saveRDS(dat, file.path(out, "data.rds"), version = 2, compress = "xz")
  write_csv(source_manifest(src), file.path(out, "source-manifest.csv"))
  write_environment(out)
  writeLines(c(
    "Prepare once; use exactly these model RDS files on both fitting platforms.",
    "Native preparation on a second platform is a separate initialization-drift audit.",
    "All final SEM fits use CSOLNP, one thread, seed 20191029 and identical case retry budgets.",
    "NLS initialization uses per-person default Gauss-Newton; failures and heuristics remain visible.",
    "Fit acceptance requires status zero, finite likelihood/parameters, available admissible",
    "stored covariances, positive-definite Hessian, valid SEs and no diagnostic issues.",
    "Historical reference comparisons are reported separately and never establish acceptance alone.",
    "This bounded matched-input study does not establish global optimality or full V1 qualification."
  ), file.path(out, "PROTOCOL.txt"))
  base_config <- list(schema_version = 1L, seed = 20191029L,
    optimizer = "CSOLNP", threads = 1L, OKstatuscodes = 0L,
    jitterDistrib = "runif", loc = 1, scale = 0.25,
    cases = cases_requested, data_sha256 = sha(file.path(out, "data.rds")),
    reference_tolerances = c(likelihood = 1e-6, parameters = 1e-4),
    nls_algorithm = "default (Gauss-Newton)", nls_control = stats::nls.control(),
    preparation_environment = environment_record())
  saveRDS(base_config, file.path(out, "protocol.rds"), version = 2)
  dput(base_config, file.path(out, "protocol.R"))
  rows <- list()
  for (id in cases_requested) {
    case_out <- new_dir(file.path(out, id))
    set_protocol()
    case <- cases[[id]]
    warning_messages <- character()
    captured <- new.env(parent = emptyenv())
    captured$calls <- list(); captured$optimizer_calls <- character()
    original <- get("getUNI.GF", envir = asNamespace("nlpsem"))
    recorded_gf <- function(dat_traj, dat_time, nT, curveFun) {
      result <- original(dat_traj, dat_time, nT, curveFun)
      captured$calls[[length(captured$calls) + 1L]] <- list(
        curve = curveFun, trajectories = colnames(dat_traj), times = colnames(dat_time),
        nT = nT, growth_factors = result,
        initialization = attr(result, "initialization", exact = TRUE))
      result
    }
    build_case <- function() {
      testthat::local_mocked_bindings(
        getUNI.GF = recorded_gf,
        mxRun = function(model, ...) {
          captured$optimizer_calls <- c(captured$optimizer_calls, "mxRun")
          model
        },
        mxTryHard = function(model, ...) {
          captured$optimizer_calls <- c(captured$optimizer_calls, "mxTryHard")
          model
        }, .package = "nlpsem")
      a <- c(list(dat = dat), case$args, list(paramOut = FALSE))
      if (!is.null(case$tries)) a$tries <- case$tries
      do.call(get(case$fn, envir = asNamespace("nlpsem")), a)
    }
    start <- proc.time()[["elapsed"]]
    object <- tryCatch(withCallingHandlers(build_case(), warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    saveRDS(captured$calls, file.path(case_out, "initialization-calls.rds"), version = 2)
    audit_rows <- list()
    for (i in seq_along(captured$calls)) {
      call <- captured$calls[[i]]
      audit <- call$initialization
      if (is.null(audit) || !is.data.frame(audit$rows)) next
      audit_rows[[length(audit_rows) + 1L]] <- cbind(
        data.frame(case = id, call = i, curve = call$curve,
                   trajectory_columns = paste(call$trajectories, collapse = ";")),
        audit$rows, stringsAsFactors = FALSE)
    }
    audit <- bind_rows(audit_rows)
    write_csv(audit, file.path(case_out, "nls-audit.csv"))
    fallback_n <- if (nrow(audit)) sum(audit$method == "heuristic_fallback") else 0L
    converged_n <- if (nrow(audit)) sum(audit$nls_converged %in% TRUE) else 0L
    row <- data.frame(case = id, run_state = "prepared",
      retry_budget = if (is.null(case$tries)) 0L else case$tries,
      preparation_seconds = proc.time()[["elapsed"]] - start,
      nls_rows = nrow(audit), nls_converged = converged_n, nls_fallback = fallback_n,
      model_sha256 = "", initial_parameters_sha256 = "",
      warnings = paste(unique(warning_messages), collapse = " | "), error = "",
      stringsAsFactors = FALSE)
    if (inherits(object, "error")) {
      row$run_state <- "preparation_error"; row$error <- conditionMessage(object)
    } else if (!length(captured$optimizer_calls) || !nrow(audit)) {
      row$run_state <- "preparation_error"
      row$error <- "No mocked fitting call or per-person NLS audit was captured"
    } else {
      model <- object@mxOutput
      stopifnot(!length(model@output$Minus2LogLikelihood))
      parameters <- OpenMx::omxGetParameters(model)
      stopifnot(length(parameters) > 0L, !anyDuplicated(names(parameters)), all(is.finite(parameters)))
      saveRDS(model, file.path(case_out, "model.rds"), version = 2, compress = "xz")
      saveRDS(parameters, file.path(case_out, "initial-parameters.rds"), version = 2)
      write_csv(data.frame(parameter = names(parameters), initial_value = as.numeric(parameters)),
                file.path(case_out, "initial-parameters.csv"))
      specification <- nlpsem::getSpecification(object)
      saveRDS(specification, file.path(case_out, "specification.rds"), version = 2)
      config <- c(list(case = id, definition = case,
                       fitting_call = captured$optimizer_calls,
                       tries = case$tries, initial_parameters = parameters), base_config)
      config$reference <- readRDS(file.path(src, "tests/testthat/golden/fit_golden.rds"))[[id]]
      config$model_sha256 <- sha(file.path(case_out, "model.rds"))
      config$initial_parameters_sha256 <- sha(file.path(case_out, "initial-parameters.rds"))
      saveRDS(config, file.path(case_out, "config.rds"), version = 2)
      dput(config[setdiff(names(config), c("initial_parameters", "reference"))],
           file.path(case_out, "config.R"))
      row$model_sha256 <- config$model_sha256
      row$initial_parameters_sha256 <- config$initial_parameters_sha256
    }
    write_csv(row, file.path(case_out, "preparation.csv"))
    rows[[id]] <- row
    write_csv(bind_rows(rows), file.path(out, "preparation.csv"))
    cat(id, row$run_state, "NLS", converged_n, "/", nrow(audit), "fallback", fallback_n, "\n")
  }
  quit(status = if (all(vapply(rows, function(r) r$run_state == "prepared", logical(1)))) 0L else 1L)
}

if (mode == "fit") {
  if (length(args) != 5L) stop("fit SOURCE INPUT_DIR CASE NEW_RESULT_DIR")
  src <- load_source(args[2]); input <- normalizePath(args[3], mustWork = TRUE)
  id <- args[4]; out <- new_dir(args[5])
  write_environment(out)
  recorded_manifest <- utils::read.csv(file.path(input, "source-manifest.csv"), stringsAsFactors = FALSE)
  current_manifest <- source_manifest(src)
  stopifnot(identical(recorded_manifest, current_manifest))
  case_in <- file.path(input, id)
  config <- readRDS(file.path(case_in, "config.rds"))
  stopifnot(identical(config$model_sha256, sha(file.path(case_in, "model.rds"))),
            identical(config$initial_parameters_sha256, sha(file.path(case_in, "initial-parameters.rds"))),
            identical(config$data_sha256, sha(file.path(input, "data.rds"))),
            identical(config$preparation_environment$OpenMx, as.character(utils::packageVersion("OpenMx"))))
  model <- readRDS(file.path(case_in, "model.rds"))
  initial <- OpenMx::omxGetParameters(model)
  stopifnot(identical(initial, config$initial_parameters),
            identical(initial, readRDS(file.path(case_in, "initial-parameters.rds"))))
  write_csv(data.frame(parameter = names(initial), initial_value = as.numeric(initial)),
            file.path(out, "verified-initial-parameters.csv"))
  saveRDS(config, file.path(out, "input-config.rds"), version = 2)
  set_protocol()
  row <- data.frame(case = id, run_state = "started", accepted = FALSE,
    status = NA_integer_, minus2LL = NA_real_, reference_minus2LL = config$reference$minus2LL,
    reference_difference = NA_real_, max_parameter_difference = NA_real_,
    reference_likelihood_match = FALSE, reference_parameters_match = FALSE,
    reference_names_match = FALSE, finite_parameters = FALSE,
    covariance_available = FALSE, covariance_admissible = FALSE,
    hessian_status = "unavailable", standard_errors = "unavailable",
    diagnostic_issues = "", warnings = "", error = "", elapsed_seconds = NA_real_,
    optimizer = "CSOLNP", threads = 1L, seed = 20191029L,
    retry_budget = if (is.null(config$tries)) 0L else config$tries,
    model_sha256 = config$model_sha256, data_sha256 = config$data_sha256,
    initial_parameters_identical = TRUE, stringsAsFactors = FALSE)
  write_csv(row, file.path(out, "started.csv"))
  warnings <- character(); started <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers({
    if (is.null(config$tries)) OpenMx::mxRun(model)
    else OpenMx::mxTryHard(model, extraTries = config$tries,
      OKstatuscodes = config$OKstatuscodes, jitterDistrib = config$jitterDistrib,
      loc = config$loc, scale = config$scale)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  row$elapsed_seconds <- proc.time()[["elapsed"]] - started
  row$warnings <- paste(unique(warnings), collapse = " | ")
  writeLines(warnings, file.path(out, "warnings.txt"))
  if (inherits(fit, "error")) {
    row$run_state <- "fit_error"; row$error <- conditionMessage(fit)
  } else {
    row$run_state <- "completed"
    saveRDS(fit, file.path(out, "fit.rds"), version = 2, compress = "xz")
    estimates <- OpenMx::omxGetParameters(fit)
    saveRDS(estimates, file.path(out, "estimates.rds"), version = 2)
    write_csv(data.frame(parameter = names(estimates), estimate = as.numeric(estimates)),
              file.path(out, "estimates.csv"))
    row$status <- one(fit@output$status$code, NA_integer_)
    row$minus2LL <- one(fit@output$Minus2LogLikelihood)
    row$finite_parameters <- length(estimates) > 0L && all(is.finite(estimates))
    row$reference_difference <- row$minus2LL - config$reference$minus2LL
    row$reference_names_match <- identical(names(estimates), names(config$reference$estimates))
    if (row$reference_names_match) {
      row$max_parameter_difference <- max(abs(estimates - config$reference$estimates))
      row$reference_parameters_match <- same_reference(estimates, config$reference$estimates, 1e-4)
    }
    row$reference_likelihood_match <- same_reference(row$minus2LL, config$reference$minus2LL, 1e-6)
    diagnostics <- tryCatch(nlpsem::getDiagnostics(fit), error = function(e) e)
    if (inherits(diagnostics, "error")) {
      row$error <- paste("Diagnostics:", conditionMessage(diagnostics))
    } else {
      saveRDS(diagnostics, file.path(out, "diagnostics.rds"), version = 2)
      capture.output(diagnostics, file = file.path(out, "diagnostics.txt"))
      write_csv(diagnostics$covariance, file.path(out, "covariance.csv"))
      row$covariance_available <- nrow(diagnostics$covariance) > 0L &&
        all(diagnostics$covariance$status != "unavailable")
      row$covariance_admissible <- row$covariance_available &&
        all(diagnostics$covariance$status %in% c("positive_definite", "near_boundary"))
      row$hessian_status <- diagnostics$hessian$status
      row$standard_errors <- diagnostics$standard_errors$status
      row$diagnostic_issues <- paste(diagnostics$issues, collapse = " | ")
      row$accepted <- isTRUE(row$status == 0L) && is.finite(row$minus2LL) &&
        row$finite_parameters && row$covariance_admissible &&
        identical(row$hessian_status, "positive_definite") &&
        identical(row$standard_errors, "no_issues_detected") &&
        !length(diagnostics$issues)
    }
  }
  write_csv(row, file.path(out, "result.csv"))
  print(row[, c("case", "run_state", "accepted", "status", "minus2LL",
                "reference_difference", "hessian_status", "diagnostic_issues")], row.names = FALSE)
  quit(status = 0L) # Completion is distinct from scientific acceptance in result.csv.
}

if (mode == "compare-preparations") {
  if (length(args) != 4L) stop("compare-preparations INPUT_A INPUT_B NEW_OUTPUT")
  a <- normalizePath(args[2], mustWork = TRUE); b <- normalizePath(args[3], mustWork = TRUE)
  out <- new_dir(args[4])
  pa <- utils::read.csv(file.path(a, "preparation.csv"), stringsAsFactors = FALSE)
  pb <- utils::read.csv(file.path(b, "preparation.csv"), stringsAsFactors = FALSE)
  data_identical <- identical(readRDS(file.path(a, "data.rds")), readRDS(file.path(b, "data.rds")))
  rows <- list()
  for (id in intersect(pa$case, pb$case)) {
    aa <- pa[pa$case == id, , drop = FALSE]; bb <- pb[pb$case == id, , drop = FALSE]
    row <- data.frame(case = id, data_identical = data_identical,
      prepared_a = aa$run_state == "prepared", prepared_b = bb$run_state == "prepared",
      parameter_names_identical = FALSE, parameters_identical = FALSE,
      max_initial_parameter_difference = NA_real_, nls_converged_a = aa$nls_converged,
      nls_converged_b = bb$nls_converged, nls_fallback_a = aa$nls_fallback,
      nls_fallback_b = bb$nls_fallback)
    if (row$prepared_a && row$prepared_b) {
      x <- readRDS(file.path(a, id, "initial-parameters.rds"))
      y <- readRDS(file.path(b, id, "initial-parameters.rds"))
      row$parameter_names_identical <- identical(names(x), names(y))
      row$parameters_identical <- identical(x, y)
      if (row$parameter_names_identical) row$max_initial_parameter_difference <- max(abs(x-y))
    }
    rows[[id]] <- row
  }
  write_csv(bind_rows(rows), file.path(out, "native-initialization-comparison.csv"))
  writeLines("Native NLS initialization drift only. The matched fitting study must use one shared input set.", file.path(out, "SCOPE.txt"))
  quit(status = 0L)
}

if (mode == "compare-fits") {
  if (length(args) != 5L) stop("compare-fits INPUT_DIR RESULTS_A RESULTS_B NEW_OUTPUT")
  input <- normalizePath(args[2], mustWork = TRUE)
  a <- normalizePath(args[3], mustWork = TRUE); b <- normalizePath(args[4], mustWork = TRUE)
  out <- new_dir(args[5])
  protocol <- readRDS(file.path(input, "protocol.rds"))
  rows <- list()
  for (id in protocol$cases) {
    paths <- c(file.path(a, id, "result.csv"), file.path(b, id, "result.csv"))
    row <- data.frame(case = id, completed_a = file.exists(paths[1]), completed_b = file.exists(paths[2]),
      initial_model_identical = FALSE, data_identical = FALSE, accepted_a = FALSE, accepted_b = FALSE,
      status_a = NA_integer_, status_b = NA_integer_, likelihood_difference_b_minus_a = NA_real_,
      max_parameter_difference = NA_real_, likelihood_matches = FALSE, parameters_match = FALSE,
      paired_gate_pass = FALSE)
    if (all(file.exists(paths))) {
      aa <- utils::read.csv(paths[1], stringsAsFactors = FALSE)
      bb <- utils::read.csv(paths[2], stringsAsFactors = FALSE)
      row$completed_a <- aa$run_state == "completed"; row$completed_b <- bb$run_state == "completed"
      row$initial_model_identical <- identical(aa$model_sha256, bb$model_sha256) &&
        isTRUE(aa$initial_parameters_identical) && isTRUE(bb$initial_parameters_identical)
      valid_sha <- function(x) is.character(x) && length(x) == 1L &&
        !is.na(x) && grepl("^[0-9a-f]{64}$", x)
      row$data_identical <- valid_sha(aa$data_sha256) && valid_sha(bb$data_sha256) &&
        identical(aa$data_sha256, bb$data_sha256)
      row$accepted_a <- aa$accepted; row$accepted_b <- bb$accepted
      row$status_a <- one(aa$status, NA_integer_)
      row$status_b <- one(bb$status, NA_integer_)
      likelihood_a <- one(aa$minus2LL); likelihood_b <- one(bb$minus2LL)
      row$likelihood_difference_b_minus_a <- likelihood_b - likelihood_a
      row$likelihood_matches <- is.finite(likelihood_a) && is.finite(likelihood_b) &&
        same_reference(likelihood_a, likelihood_b, 1e-6)
      ep <- c(file.path(a, id, "estimates.rds"), file.path(b, id, "estimates.rds"))
      if (all(file.exists(ep))) {
        x <- readRDS(ep[1]); y <- readRDS(ep[2])
        if (identical(names(x), names(y))) {
          row$max_parameter_difference <- max(abs(x-y))
          row$parameters_match <- same_reference(x, y, 1e-4)
        }
      }
      row$paired_gate_pass <- all(row$completed_a, row$completed_b, row$initial_model_identical,
        row$data_identical, row$accepted_a, row$accepted_b, row$likelihood_matches, row$parameters_match)
    }
    rows[[id]] <- row
  }
  write_csv(bind_rows(rows), file.path(out, "matched-platform-comparison.csv"))
  writeLines("Pairwise gate combines identical input models/data, both fit acceptance screens, and unchanged historical numerical tolerances. Not global optimum or full scientific qualification.", file.path(out, "SCOPE.txt"))
  quit(status = 0L)
}
stop("Unknown mode: ", mode)

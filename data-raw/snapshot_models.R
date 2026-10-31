# Build (without fitting) a 688-specification grid through the public wrappers.
# Usage: Rscript snapshot_models.R <pkg> <out.rds> [family regex]
#        [auto|capture|replay] [inputs.rds] [historical-baseline.rds]
# capture requires an exact match to the historical baseline before saving inputs.
# replay supplies frozen data/specifications/resolved starts, so construction is
# tested independently of platform-dependent automatic initialization.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Supply package path and output snapshot path")
W <- args[1]; out <- args[2]; fam <- if (length(args) >= 3) args[3] else "."
mode <- match.arg(if (length(args) >= 4L) args[4] else "auto", c("auto", "capture", "replay"))
if (mode != "auto" && length(args) < 5L) stop("Supply an inputs.rds path for capture or replay")
if (mode == "capture" && length(args) < 6L) stop("Capture requires a historical baseline path")
inputs <- if (mode == "replay") readRDS(args[5]) else NULL
if (mode == "replay") {
  stopifnot(identical(inputs$schema_version, 1L), is.data.frame(inputs$data),
            length(inputs$specs) == 688L, identical(names(inputs$specs), names(inputs$starts)),
            all(vapply(inputs$starts, function(x) !is.null(x), logical(1))))
}
suppressPackageStartupMessages({pkgload::load_all(W, quiet = TRUE); library(OpenMx)})
imp <- parent.env(asNamespace("nlpsem"))
assign("mxRun", function(model, ...) model, envir = imp)
assign("mxTryHard", function(model, ...) model, envir = imp)
source(file.path(W, "tests/testthat/helper-golden.R"))
dat <- if (mode == "replay") inputs$data else .golden_data(200)

slot_canon <- function(obj) {
  if (is.null(obj)) return(NULL)
  sn <- methods::slotNames(obj); res <- list(class = class(obj)[1])
  for (s in sn) {
    v <- tryCatch(methods::slot(obj, s), error = function(e) NULL)
    if (is.null(v) || is.environment(v) || is.function(v)) next
    if (is.character(v) || is.numeric(v) || is.logical(v)) res[[s]] <- v
    else if (methods::is(v, "MxMatrix")) res[[s]] <- list(values = v@values, free = v@free, labels = v@labels)
  }
  res
}
.sorted <- function(x) if (length(x)) x[order(names(x))] else x
canon <- function(m) {
  list(
    name = m@name, class = class(m)[1], manifestVars = m@manifestVars, latentVars = m@latentVars,
    data = if (is.null(m@data)) NULL else {
      obs <- m@data@observed; num <- vapply(obs, is.numeric, logical(1))
      list(dim = dim(obs), cols = colnames(obs), type = m@data@type,
           checksum = sum(as.matrix(obs[, num, drop = FALSE]), na.rm = TRUE), n_na = sum(is.na(obs)))
    },
    matrices = lapply(.sorted(m@matrices), function(x) list(class = class(x)[1], values = x@values, free = x@free, labels = x@labels,
                                                   lbound = x@lbound, ubound = x@ubound, dimnames = dimnames(x@values))),
    algebras = lapply(.sorted(m@algebras), function(a) list(formula = paste(deparse(a@formula, width.cutoff = 500), collapse = " "))),
    expectation = slot_canon(m@expectation), fitfunction = slot_canon(m@fitfunction),
    intervals = sort(names(m@intervals)), constraints = lapply(m@constraints, function(cc) paste(deparse(cc@formula), collapse = " ")),
    submodels = lapply(m@submodels, canon)
  )
}

# ---- specification grid ----
lgcm_forms <- list(c("LIN", FALSE), c("QUAD", FALSE), c("EXP", FALSE), c("EXP", TRUE), c("JB", FALSE), c("JB", TRUE), c("BLS", FALSE), c("BLS", TRUE))
lcsm_forms <- list(c("NonP", FALSE), c("QUAD", FALSE), c("EXP", FALSE), c("EXP", TRUE), c("JB", FALSE), c("JB", TRUE))
tics <- list(NULL, "ex1", c("ex1", "gx1"))
tag <- function(f) paste0(f[1], if (as.logical(f[2])) "_r" else "_f")
specs <- list()
add <- function(id, fn, args) specs[[id]] <<- list(fn = fn, args = args)
for (f in lgcm_forms) for (i in seq_along(tics)) add(sprintf("LGCM_%s_tic%d", tag(f), i - 1), "getLGCM",
  list(t_var = "T", y_var = "M", curveFun = f[1], intrinsic = as.logical(f[2]), records = 1:9, growth_TIC = tics[[i]]))
for (f in lcsm_forms) for (i in seq_along(tics)) add(sprintf("LCSM_%s_tic%d", tag(f), i - 1), "getLCSM",
  list(t_var = "T", y_var = "M", curveFun = f[1], intrinsic = as.logical(f[2]), records = 1:9, growth_TIC = tics[[i]]))
for (ym in c("LGCM", "LCSM")) for (f in if (ym == "LGCM") lgcm_forms else lcsm_forms) for (d in 0:3) for (i in 1:2)
  add(sprintf("TVC_%s_%s_d%d_tic%d", ym, tag(f), d, i - 1), "getTVCmodel",
      list(t_var = "T", y_var = "M", curveFun = f[1], intrinsic = as.logical(f[2]), records = 1:9, y_model = ym, TVC = "Rs",
           decompose = d, growth_TIC = tics[[i]]))
for (ym in c("LGCM", "LCSM")) for (f in if (ym == "LGCM") lgcm_forms else lcsm_forms)
  add(sprintf("MGM_%s_%s", ym, tag(f)), "getMGM",
      list(t_var = c("T", "T"), y_var = c("R", "M"), curveFun = f[1], intrinsic = as.logical(f[2]), records = list(1:9, 1:9), y_model = ym))
for (cf in c("LIN", "BLS")) {
  add(sprintf("MED_%s_baseline", cf), "getMediation", list(t_var = rep("T", 2), y_var = "M", m_var = "R", x_type = "baseline", x_var = "ex1", curveFun = cf, records = list(1:9, 1:9)))
  add(sprintf("MED_%s_longitudinal", cf), "getMediation", list(t_var = rep("T", 3), y_var = "S", m_var = "M", x_type = "longitudinal", x_var = "R", curveFun = cf, records = list(2:9, 1:9, 1:9)))
}
# class-specific containers reuse the single-group specs
sub_of <- c(getLGCM = "LGCM", getLCSM = "LCSM", getTVCmodel = "TVC", getMGM = "MGM", getMediation = "MED")
single <- specs
for (id in names(single)) {
  sp <- single[[id]]; a <- sp$args
  specs[[paste0("MGroup_", id)]] <- list(fn = "getMGroup", args = c(list(grp_var = "SEX", sub_Model = sub_of[[sp$fn]]), a))
  specs[[paste0("MIX_", id)]] <- list(fn = "getMIX", args = c(list(prop_starts = c(0.5, 0.5), sub_Model = sub_of[[sp$fn]]), a))
  specs[[paste0("MIXc_", id)]] <- list(fn = "getMIX", args = c(list(prop_starts = c(0.5, 0.5), sub_Model = sub_of[[sp$fn]], cluster_TIC = "gx1"), a))
}
if (mode == "replay") {
  if (!identical(specs, inputs$specs)) stop("Replay specifications differ from the historical grid")
  specs <- inputs$specs
}
specs <- specs[grepl(fam, names(specs))]
if (mode == "capture" && length(specs) != 688L) stop("Capture requires the complete 688-model grid")
cat("specs:", length(specs), "; input mode:", mode, "\n")
captured_starts <- stats::setNames(vector("list", length(specs)), names(specs))
snap <- list(); t0 <- Sys.time()
for (id in names(specs)) {
  sp <- specs[[id]]; fn <- get(sp$fn, envir = asNamespace("nlpsem"))
  set.seed(20191029)
  calls <- NULL
  res <- tryCatch(withCallingHandlers({
    call_args <- c(list(dat = dat), sp$args, list(paramOut = FALSE))
    if (mode == "replay") call_args$starts <- inputs$starts[[id]]
    fit <- suppressMessages(suppressWarnings(do.call(fn, call_args)))
    if (mode == "capture") {
      metadata <- getSpecification(fit)
      stopifnot(identical(metadata$fitting$starts_source, "automatic"), !is.null(metadata$fitting$starts))
      captured_starts[id] <- list(metadata$fitting$starts)
    }
    if (mode == "replay") {
      metadata <- getSpecification(fit)
      stopifnot(identical(metadata$fitting$starts_source, "user"),
                identical(metadata$fitting$starts, inputs$starts[[id]]))
    }
    canon(fit@mxOutput)
  }, error = function(e) calls <<- sys.calls()), error = function(e) list(error = conditionMessage(e),
       trace = paste(vapply(calls, function(cl) substr(paste(deparse(cl), collapse = " "), 1, 160), ""), collapse = "\n")))
  snap[[id]] <- res
}
cat("built", length(snap), "specs in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min;",
    sum(vapply(snap, function(x) !is.null(x$error), logical(1))), "errors\n")
errs <- names(snap)[vapply(snap, function(x) !is.null(x$error), logical(1))]
for (e in errs) cat("  ERROR", e, ":", substr(snap[[e]]$error, 1, 120), "\n")
saveRDS(snap, out); cat("written", out, "\n")

if (length(errs)) quit(status = 1L)

if (mode == "capture") {
  historical <- readRDS(args[6])
  if (!identical(snap, historical)) {
    stop("Captured structures differ from the historical baseline; inputs were not saved. ",
         paste(all.equal(snap, historical, tolerance = 0), collapse = "; "))
  }
  revision <- system2("git", c("-C", shQuote(normalizePath(W)), "rev-parse", "HEAD"), stdout = TRUE)
  bundle <- list(schema_version = 1L, source_commit = revision,
                 source_version = as.character(utils::packageVersion("nlpsem")),
                 session = list(R = as.character(getRversion()), platform = R.version$platform,
                                OpenMx = as.character(utils::packageVersion("OpenMx")),
                                pkgload = as.character(utils::packageVersion("pkgload")),
                                nlpsem = as.character(utils::packageVersion("nlpsem"))),
                 data = dat, specs = specs, starts = captured_starts)
  saveRDS(bundle, args[5], compress = "xz")
  cat("Captured resolved inputs for", length(specs), "models; exact historical equality verified.\n")
  cat("Inputs written to", args[5], "\n")
}

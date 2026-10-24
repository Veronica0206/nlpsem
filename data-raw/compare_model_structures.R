# Compare the historical Phase C grid with a newly generated snapshot.
# Usage: Rscript data-raw/compare_model_structures.R baseline.rds current.rds
# Only the F5 split of an identical ADD1 constant into process-specific constants
# is normalized. A 1e-10 * max(1, abs(old), abs(new)) elementwise allowance covers
# cross-platform roundoff in numeric matrix values/bounds and data checksums only.
# ADD1 constants, types,
# dimensions, free masks, labels, and formulas remain exact.
compare_model_structures <- function(baseline, current, expected_count = 688L,
                                     numeric_tolerance = 1e-10) {
  if (!is.numeric(numeric_tolerance) || length(numeric_tolerance) != 1L ||
      !is.finite(numeric_tolerance) || numeric_tolerance < 0 || numeric_tolerance > 1e-10) {
    stop("numeric_tolerance must be between 0 and 1e-10")
  }
  if (length(baseline) != expected_count || length(current) != expected_count) {
    stop("Expected ", expected_count, " model structures in both snapshots")
  }
  if (is.null(names(baseline)) || anyDuplicated(names(baseline)) ||
      !identical(names(baseline), names(current))) {
    stop("Model names differ or are not unique")
  }
  renamed <- character()
  normalize <- function(old, new, path) {
    if (!is.null(old$error) || !is.null(new$error)) stop("Construction error at ", path)
    if (!is.null(old$matrices$ADD1) && is.null(new$matrices$ADD1)) {
      old_exprs <- lapply(old$algebras, function(a) str2lang(a$formula))
      affected <- names(old_exprs)[vapply(old_exprs, function(expr) {
        "ADD1" %in% all.names(expr, functions = FALSE)
      }, logical(1))]
      if (!length(affected) || !all(grepl("r_loads$", affected))) {
        stop("Unrecognized original ADD1 references at ", path)
      }
      # Every outcome must receive its own constant, including class prefixes.
      # Checking only equal contents could accept two outcomes sharing one new
      # constant, reintroducing the unequal-wave defect that F5 fixes.
      expected_keys <- sub("r_loads$", "ADD1", affected)
      keys <- grep("ADD1$", names(new$matrices), value = TRUE)
      if (!identical(sort(keys), sort(expected_keys))) {
        stop("Outcome-specific ADD1 constant names differ at ", path)
      }
      if (!all(vapply(new$matrices[keys], identical, logical(1), old$matrices$ADD1))) {
        stop("ADD1 constants differ in content at ", path)
      }
      replace_symbol <- function(expr, replacement) {
        if (identical(expr, as.name("ADD1"))) return(as.name(replacement))
        if (is.call(expr)) {
          return(as.call(lapply(as.list(expr), replace_symbol, replacement = replacement)))
        }
        expr
      }
      for (i in seq_along(affected)) {
        algebra <- affected[i]
        expected_expr <- replace_symbol(old_exprs[[algebra]], expected_keys[i])
        if (is.null(new$algebras[[algebra]]) ||
            !identical(str2lang(new$algebras[[algebra]]$formula), expected_expr)) {
          stop("Incorrect ADD1 outcome mapping or formula at ", path, "/", algebra)
        }
        # Normalize only the formula just proved to be the approved rename.
        # All other algebra fields and all unaffected formulas remain untouched.
        new$algebras[[algebra]]$formula <- old$algebras[[algebra]]$formula
      }
      new$matrices[keys] <- NULL
      new$matrices$ADD1 <- old$matrices$ADD1
      new$matrices <- new$matrices[order(names(new$matrices))]
      renamed <<- c(renamed, path)
    }
    if (!identical(names(old$submodels), names(new$submodels))) {
      stop("Submodel names differ at ", path)
    }
    for (i in seq_along(old$submodels)) {
      new$submodels[[i]] <- normalize(old$submodels[[i]], new$submodels[[i]],
                                     paste0(path, "/", i))
    }
    new
  }
  normalized <- current
  had_renames <- setNames(rep(FALSE, length(baseline)), names(baseline))
  for (id in names(baseline)) {
    before <- length(renamed)
    normalized[[id]] <- normalize(baseline[[id]], current[[id]], id)
    had_renames[id] <- length(renamed) > before
  }
  exact <- vapply(names(baseline), function(id) identical(baseline[[id]], current[[id]]), logical(1))
  same_after_rename <- vapply(names(baseline), function(id) {
    identical(baseline[[id]], normalized[[id]])
  }, logical(1))
  roundoff_fields <- character()
  max_roundoff <- 0
  max_scaled_roundoff <- 0
  normalize_number <- function(old, new, path) {
    if (identical(old, new)) return(new)
    if (!is.numeric(old) || !is.numeric(new) || !identical(typeof(old), typeof(new)) ||
        !identical(attributes(old), attributes(new)) || length(old) != length(new)) return(new)
    finite <- is.finite(old)
    if (!identical(finite, is.finite(new)) || !identical(old[!finite], new[!finite])) return(new)
    differences <- abs(old[finite] - new[finite])
    scales <- pmax(1, abs(old[finite]), abs(new[finite]))
    scaled_differences <- differences / scales
    if (!all(scaled_differences <= numeric_tolerance)) return(new)
    roundoff_fields <<- c(roundoff_fields, path)
    max_roundoff <<- max(c(max_roundoff, differences))
    max_scaled_roundoff <<- max(c(max_scaled_roundoff, scaled_differences))
    old
  }
  normalize_roundoff <- function(old, new, path) {
    # This is an explicit allowlist, not a tolerant recursive comparison of all
    # numeric fields. In particular, constraints and parameter labels stay exact.
    for (key in intersect(names(old$matrices), names(new$matrices))) {
      if (grepl("ADD1$", key)) next
      for (field in c("values", "lbound", "ubound")) {
        value <- normalize_number(old$matrices[[key]][[field]], new$matrices[[key]][[field]],
                                  paste0(path, "/matrices/", key, "/", field))
        if (!is.null(new$matrices[[key]][[field]])) new$matrices[[key]][[field]] <- value
      }
    }
    if (!is.null(old$data$checksum) && !is.null(new$data$checksum)) {
      new$data$checksum <- normalize_number(old$data$checksum, new$data$checksum,
                                            paste0(path, "/data/checksum"))
    }
    for (i in seq_along(old$submodels)) {
      new$submodels[[i]] <- normalize_roundoff(old$submodels[[i]], new$submodels[[i]],
                                               paste0(path, "/", i))
    }
    new
  }
  for (id in names(baseline)) {
    normalized[[id]] <- normalize_roundoff(baseline[[id]], normalized[[id]], id)
  }
  same <- vapply(names(baseline), function(id) identical(baseline[[id]], normalized[[id]]), logical(1))
  roundoff <- same & !same_after_rename
  rename_only <- same_after_rename & !exact
  roundoff_only <- roundoff & !had_renames
  both <- roundoff & had_renames
  cat("Models:", length(baseline), "; exact:", sum(exact),
      "; approved ADD1 renames only:", sum(rename_only),
      "; roundoff only:", sum(roundoff_only),
      "; approved renames plus roundoff:", sum(both),
      "; unexplained:", sum(!same), "\n")
  if (length(renamed)) cat("Renamed constant scopes:", paste(renamed, collapse = ", "), "\n")
  cat("Scaled numeric tolerance:", numeric_tolerance,
      "; maximum scaled roundoff:", max_scaled_roundoff,
      "; maximum absolute roundoff:", max_roundoff, "\n")
  if (any(!same)) {
    for (id in names(baseline)[!same]) {
      cat(id, ":", paste(all.equal(baseline[[id]], normalized[[id]], tolerance = 0),
                         collapse = "; "), "\n")
    }
    stop("Unexpected structural differences")
  }
  invisible(list(models = length(baseline), exact = sum(exact),
                 approved_renames = sum(same & had_renames), renamed_scopes = renamed,
                 roundoff_equivalent = sum(roundoff), renames_only = sum(rename_only),
                 roundoff_only = sum(roundoff_only), renames_and_roundoff = sum(both),
                 numeric_tolerance = numeric_tolerance, max_abs_roundoff = max_roundoff,
                 max_scaled_roundoff = max_scaled_roundoff,
                 roundoff_fields = roundoff_fields))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Supply baseline.rds and current.rds")
  compare_model_structures(readRDS(args[1]), readRDS(args[2]))
}

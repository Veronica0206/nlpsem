# Compare the historical Phase C grid with a newly generated snapshot.
# Usage: Rscript data-raw/compare_model_structures.R baseline.rds current.rds
# Approved changes are the F5 split of identical ADD1 constants and the
# reviewed longitudinal MED M-X residual starting covariance correction, and
# two exact centered random-rate mean formulas in derived LCSM output with
# covariates. No fitting equations, free masks, labels or dimensions are
# normalized. A 1e-10 * max(1, abs(old), abs(new)) elementwise allowance covers
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
  mediation_starts <- character()
  lcsm_means <- character()
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
    # The immutable 688-case grid uses S (outcome), M (mediator), R
    # (predictor), and the default residual correlation 0.3. The old builder
    # used var(S) instead of var(M) only for M-R covariance starts. Approve
    # only these named cases and cells, after checking both formulas against
    # unchanged residual diagonal starts. This is not a generic tolerance.
    scope <- strsplit(path, "/", fixed = TRUE)[[1L]]
    approved_id <- grepl("^(MGroup_|MIX_|MIXc_)?MED_(LIN|BLS)_longitudinal$", scope[1L])
    class_scope <- length(scope) == 2L && grepl("^[1-9][0-9]*$", scope[2L])
    single_scope <- length(scope) == 1L && grepl("^MED_", scope[1L])
    if (approved_id && (single_scope || class_scope) &&
        !is.null(old$matrices$S) && !is.null(new$matrices$S) &&
        !identical(old$matrices$S$values, new$matrices$S$values)) {
      old_s <- old$matrices$S
      new_s <- new$matrices$S
      prefix <- if (class_scope) paste0("c", scope[2L]) else ""
      residual_variance <- function(role) {
        indices <- which(diag(old_s$labels) == paste0(prefix, role, "_residuals"))
        values <- unique(diag(old_s$values)[indices])
        if (length(values) != 1L || !is.finite(values) || values <= 0) {
          stop("Unrecognized MED residual diagonal at ", path)
        }
        values
      }
      vy <- residual_variance("S")
      vm <- residual_variance("M")
      vx <- residual_variance("R")
      cells <- which(!is.na(old_s$labels) & old_s$labels == paste0(prefix, "MR_RES"))
      close <- function(values, expected) {
        length(values) && all(is.finite(values)) &&
          all(abs(values - expected) <= numeric_tolerance * pmax(1, abs(values), abs(expected)))
      }
      if (!length(cells) || !identical(dim(old_s$values), dim(new_s$values)) ||
          !identical(typeof(old_s$values), typeof(new_s$values)) ||
          !close(old_s$values[cells], .3 * sqrt(vy * vx)) ||
          !close(new_s$values[cells], .3 * sqrt(vm * vx))) {
        stop("Incorrect MED M-X residual starting covariance at ", path)
      }
      # All untouched cells and all structure must still pass exact/roundoff
      # comparison below. A changed residual variance cannot justify a change.
      new$matrices$S$values[cells] <- old_s$values[cells]
      mediation_starts <<- c(mediation_starts, path)
    }
    # Only these covariate-bearing random-rate models can have a nonzero
    # mean rate deviation. Approve the two exact reporting-formula changes;
    # loadings, fitting equations, parameter structure and other output remain
    # subject to the unchanged comparison below.
    lcsm_id <- grepl(paste0("^(MGroup_|MIX_|MIXc_)?(",
      "LCSM_(EXP|JB)_r_tic[12]|",
      "TVC_LCSM_(EXP|JB)_r_(d0_tic1|d[123]_tic[01]))$"), scope[1L])
    lcsm_single <- length(scope) == 1L && !grepl("^(MGroup_|MIX_|MIXc_)", scope[1L])
    if (lcsm_id && (lcsm_single || class_scope)) {
      prefix <- if (class_scope) paste0("c", scope[2L]) else ""
      keys <- paste0(prefix, c("Yslp_m", "Ychg_bl_m"))
      changed <- vapply(keys, function(key) {
        !identical(old$algebras[[key]], new$algebras[[key]])
      }, logical(1))
      if (any(changed)) {
        exp_curve <- grepl("_EXP_r_", scope[1L])
        slice <- if (exp_curve) "2" else "2:3"
        last <- if (exp_curve) "3" else "4"
        for (j in seq_along(keys)) {
          loading <- paste0(prefix, if (j == 1L) "r_loads" else "s_loads[-1, ]")
          head <- paste0(loading, " %*% rbind(", prefix, "Y_mean0[", slice, ", ], ")
          old_formula <- paste0(head, "0)")
          new_formula <- paste0(head, prefix, "Y_mean0[", last, ", ] - ", prefix, "Y_mug)")
          a <- old$algebras[[keys[j]]]; b <- new$algebras[[keys[j]]]
          if (is.null(a$formula) || is.null(b$formula) ||
              !identical(str2lang(a$formula), str2lang(old_formula)) ||
              !identical(str2lang(b$formula), str2lang(new_formula))) {
            stop("Incorrect LCSM centered-rate reporting formula at ", path, "/", keys[j])
          }
          new$algebras[[keys[j]]]$formula <- a$formula
        }
        lcsm_means <<- c(lcsm_means, path)
      }
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
  had_mediation_starts <- had_renames
  had_lcsm_means <- had_renames
  for (id in names(baseline)) {
    before <- length(renamed)
    before_mediation <- length(mediation_starts)
    before_lcsm <- length(lcsm_means)
    normalized[[id]] <- normalize(baseline[[id]], current[[id]], id)
    had_renames[id] <- length(renamed) > before
    had_mediation_starts[id] <- length(mediation_starts) > before_mediation
    had_lcsm_means[id] <- length(lcsm_means) > before_lcsm
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
  rename_only <- same_after_rename & !exact & had_renames & !had_mediation_starts
  mediation_only <- same_after_rename & !exact & had_mediation_starts & !had_renames
  roundoff_only <- roundoff & !had_renames & !had_mediation_starts & !had_lcsm_means
  lcsm_only <- same_after_rename & !exact & had_lcsm_means & !had_renames & !had_mediation_starts
  both <- roundoff & had_renames
  cat("Models:", length(baseline), "; exact:", sum(exact),
      "; approved ADD1 renames only:", sum(rename_only),
      "; approved MED starts only:", sum(mediation_only),
      "; approved LCSM reporting means only:", sum(lcsm_only),
      "; roundoff only:", sum(roundoff_only),
      "; approved renames plus roundoff:", sum(both),
      "; approved MED starts plus roundoff:", sum(roundoff & had_mediation_starts),
      "; approved LCSM means plus roundoff:", sum(roundoff & had_lcsm_means),
      "; unexplained:", sum(!same), "\n")
  if (length(mediation_starts)) cat("Corrected MED start scopes:", paste(mediation_starts, collapse = ", "), "\n")
  if (length(lcsm_means)) cat("Corrected LCSM reporting scopes:", paste(lcsm_means, collapse = ", "), "\n")
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
                 approved_mediation_starts = sum(same & had_mediation_starts),
                 mediation_start_scopes = mediation_starts, mediation_starts_only = sum(mediation_only),
                 mediation_starts_and_roundoff = sum(roundoff & had_mediation_starts),
                 approved_lcsm_means = sum(same & had_lcsm_means),
                 lcsm_mean_scopes = lcsm_means, lcsm_means_only = sum(lcsm_only),
                 lcsm_means_and_roundoff = sum(roundoff & had_lcsm_means),
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

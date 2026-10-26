#' Inspect Numerical Diagnostics of a Stored Fit
#'
#' @description Separates optimizer completion from covariance admissibility,
#' Hessian quality, standard-error availability and mixture component sizes.
#' The function inspects stored results without refitting, changing constraints,
#' or treating an optimizer status of zero as statistical acceptance.
#'
#' @param object A fitted \code{myMxOutput} object or an OpenMx \code{MxModel}.
#' @param tolerance Positive finite tolerance below one for symmetry and
#' eigenvalue checks after diagonal scaling. Default is \code{1e-8}.
#'
#' @return A list with \code{optimizer}, a \code{covariance} data frame,
#' a \code{cross_covariances} data frame recording recognized cross blocks,
#' \code{hessian}, \code{standard_errors}, \code{component_sizes},
#' \code{issues}, \code{repeated_starts}, and \code{limitations}.
#' Each covariance row identifies its model, matrix and block. The raw and
#' scaled minimum eigenvalues are reported separately. Empty or unsupported
#' checks are explicitly unavailable. \code{issues} lists detected problems;
#' an empty list of issues does not certify statistical validity.
#'
#' @details
#' A fitted result requires a finite stored \code{Minus2LogLikelihood} and an
#' optimizer status. Other fit objectives, including weighted least squares,
#' are outside this diagnostic's current scope. For fitted RAM models, checks
#' cover the stored S matrix and its latent
#' disturbance and manifest residual blocks, including jointly fitted children.
#' Cached nlpsem covariance algebras ending in \code{_psi0}, \code{_psi_r} or
#' \code{_psi_s} are also inspected. Scalar-label block algebras are matched to
#' their positions in RAM S. A block joining distinct variable sets is reported
#' separately in \code{cross_covariances}: its entries may be negative and it need not be
#' symmetric. Admissibility is checked on the assembled joint RAM S, not on
#' the cross block alone. This structural matching does not infer roles from
#' outcome-name lengths or from the numerical symmetry of a fitted block.
#' Covariance matrices depending on definition
#' variables are not evaluated for individual rows. Other expectation types,
#' missing cached results and unfitted models are marked unavailable.
#'
#' Eigenvalues are calculated after a diagonal congruence transformation using
#' the square roots of absolute diagonal entries; a zero diagonal uses its row's
#' largest absolute entry (or the smallest positive double for an all-zero row).
#' This permits structurally zero disturbance variances. A negative diagonal is
#' always reported. Scaled eigenvalues below minus \code{tolerance} indicate an
#' indefinite matrix; near-zero eigenvalues are reported as \code{near_boundary}
#' and require substantive interpretation rather than automatic rejection.
#'
#' Hessian diagnostics use the stored calculated Hessian, falling back to the
#' stored Hessian, and separately report OpenMx's \code{infoDefinite} flag.
#' Positive semidefiniteness is insufficient for an invertible Hessian, so a
#' singular or near-singular Hessian is flagged. Standard errors are checked for
#' completeness, finiteness and positivity for the estimated free parameters.
#' These checks do not establish identification or calibrated uncertainty.
#'
#' For recognized nlpsem mixture models, component sizes include hard assignment
#' counts and sums of posterior probabilities. Empty hard-assignment classes
#' are flagged; there is no universal acceptable class-size threshold.
#' Multiple starts, likelihood optimality, row-specific implied covariances,
#' interval coverage and statistical model assumptions remain unassessed.
#' Results describe the stored fit; manually changing it can leave cached
#' outputs stale. Refit a modified model before interpreting these diagnostics.
#' Covariance checks assume references resolve within each RAM model; custom
#' cross-model covariance expressions require additional independent checks.
#'
#' @examples
#' load(system.file("extdata", "getTVCmodel_examples.RData", package = "nlpsem"))
#' diagnostics <- getDiagnostics(Math_TVCchgBL_BLS_r)
#' diagnostics$optimizer
#' diagnostics$covariance
#' diagnostics$issues
#'
#' @seealso \code{\link{getSpecification}}, \code{\link{getPosterior}}
#' @export
getDiagnostics <- function(object, tolerance = 1e-8) {
  if (methods::is(object, "myMxOutput")) object <- object@mxOutput
  if (!methods::is(object, "MxModel")) {
    stop("'object' must be a myMxOutput or MxModel object.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1) {
    stop("'tolerance' must be one finite number between zero and one.", call. = FALSE)
  }
  code <- object@output$status$code
  status <- if (length(code) == 1L && is.finite(code)) as.integer(code) else NA_integer_
  fit <- object@output$Minus2LogLikelihood
  fitted <- length(fit) == 1L && is.finite(fit) && !is.na(status)
  issues <- character()
  if (!is.na(status) && status != 0L) issues <- c(issues, paste("Optimizer status is", status))
  if (!fitted) issues <- c(issues, "A finite fitted likelihood and optimizer status are not both available")

  rows <- cross_rows <- list()
  add <- function(model, matrix, block, values = NULL, unavailable = NULL) {
    check <- .diagnostic_matrix(values, tolerance)
    if (!is.null(unavailable)) {
      check$status <- "unavailable"
      check$detail <- unavailable
    }
    rows[[length(rows) + 1L]] <<- data.frame(
      model = model, matrix = matrix, block = block, dimension = check$dimension,
      status = check$status, minimum_variance = check$minimum_diagonal,
      minimum_eigenvalue = check$minimum_eigenvalue,
      scaled_minimum_eigenvalue = check$scaled_minimum_eigenvalue,
      negative_variances = check$negative_diagonals, detail = check$detail,
      stringsAsFactors = FALSE)
  }
  walk <- function(model, path, jointly_fitted) {
    ram <- methods::is(model@expectation, "MxExpectationRAM")
    if (ram) {
      ref <- model@expectation@S
      S <- if (length(ref) == 1L) model@matrices[[ref]] else NULL
      if (!jointly_fitted) {
        add(path, if (length(ref) == 1L) ref else "S", "RAM disturbances",
            unavailable = "Model has no fitted results; starts are not checked as estimates")
      } else if (is.null(S)) {
        add(path, if (length(ref) == 1L) ref else "S", "RAM disturbances",
            unavailable = "RAM S is not a local stored matrix")
      } else if (.diagnostic_has_definition_variables(model, ref)) {
        add(path, ref, "RAM disturbances", unavailable = "Definition-variable covariance requires row-specific evaluation")
      } else {
        values <- S@values
        add(path, ref, "RAM disturbances", values)
        labels <- rownames(values)
        if (!is.null(labels)) {
          latent <- intersect(labels, model@latentVars)
          manifest <- intersect(labels, model@manifestVars)
          if (length(latent)) add(path, ref, "Latent disturbances", values[latent, latent, drop = FALSE])
          if (length(manifest)) add(path, ref, "Manifest residuals", values[manifest, manifest, drop = FALSE])
        }
      }
      covariances <- grep("_psi(0|_r|_s)$", names(model@algebras), value = TRUE)
      for (name in covariances) {
        unavailable <- if (!jointly_fitted) "Model has no fitted results" else if (
          .diagnostic_has_definition_variables(model, name)) "Definition-variable covariance requires row-specific evaluation" else NULL
        values <- model@algebras[[name]]@result
        cross <- if (!is.null(S)) {
          .diagnostic_cross_covariance(model@algebras[[name]]@formula, S@labels)
        } else NULL
        if (is.null(cross)) {
          add(path, name, "Stored covariance algebra", values, unavailable)
        } else {
          stored <- is.matrix(values) && length(values) &&
            identical(dim(values), c(length(cross$rows), length(cross$columns)))
          finite <- stored && all(is.finite(values))
          joint_available <- jointly_fitted && !.diagnostic_has_definition_variables(model, ref)
          cross_status <- if (!is.null(unavailable) || !stored) "unavailable" else if (!finite) "nonfinite" else if (
            joint_available) "joint_covariance_checked" else "joint_covariance_unavailable"
          cross_rows[[length(cross_rows) + 1L]] <<- data.frame(model = path, matrix = name,
            rows = length(cross$rows), columns = length(cross$columns),
            row_variables = paste(cross$rows, collapse = ", "),
            column_variables = paste(cross$columns, collapse = ", "),
            status = cross_status,
            detail = if (cross_status == "joint_covariance_checked") {
              "Cross-covariance of distinct RAM variable sets; admissibility is checked in the full RAM S matrix"
            } else "Cross-covariance of distinct RAM variable sets; inspect its stored-result status and RAM S availability",
            stringsAsFactors = FALSE)
          if (cross_status == "nonfinite") issues <<- c(issues, paste0("Cross-covariance problem: ", path, "/", name, " (nonfinite)"))
        }
      }
    } else if (!length(model@submodels)) {
      add(path, NA_character_, "Covariance", unavailable = "Only stored RAM covariance matrices are supported")
    }
    for (name in names(model@submodels)) {
      walk(model@submodels[[name]], paste(path, name, sep = "/"), jointly_fitted)
    }
  }
  walk(object, object@name, fitted)
  covariance <- if (length(rows)) do.call(rbind, rows) else data.frame(
    model = character(), matrix = character(), block = character(), dimension = integer(),
    status = character(), minimum_variance = numeric(), minimum_eigenvalue = numeric(),
    scaled_minimum_eigenvalue = numeric(), negative_variances = character(), detail = character())
  cross_covariances <- if (length(cross_rows)) do.call(rbind, cross_rows) else data.frame(
    model = character(), matrix = character(), rows = integer(), columns = integer(),
    row_variables = character(), column_variables = character(), status = character(), detail = character())
  bad_covariance <- covariance$status %in% c("indefinite", "negative_variance", "nonfinite", "not_symmetric")
  if (any(bad_covariance)) issues <- c(issues, paste0("Covariance problem: ",
    covariance$model[bad_covariance], "/", covariance$matrix[bad_covariance], " (",
    covariance$block[bad_covariance], ": ", covariance$status[bad_covariance], ")"))

  H <- object@output$calculatedHessian
  hessian_source <- "calculatedHessian"
  if (!is.matrix(H) || !length(H)) {
    H <- object@output$hessian
    hessian_source <- "hessian"
  }
  hessian <- .diagnostic_matrix(if (fitted) H else NULL, tolerance)
  if (hessian$status == "negative_variance") hessian$status <- "negative_diagonal"
  hessian$source <- if (is.matrix(H) && length(H)) hessian_source else NA_character_
  flag <- object@output$infoDefinite
  hessian$openmx_info_definite <- if (is.logical(flag) && length(flag) == 1L) flag else NA
  if (!hessian$status %in% c("positive_definite", "unavailable")) {
    issues <- c(issues, paste("Hessian is", hessian$status))
  }
  if (identical(hessian$openmx_info_definite, FALSE)) issues <- c(issues, "OpenMx reports infoDefinite = FALSE")
  n_parameters <- length(OpenMx::omxGetParameters(object, free = TRUE))
  se <- object@output$standardErrors
  se <- if (fitted && is.numeric(se)) as.numeric(se) else numeric()
  standard_errors <- list(status = if (!length(se)) "unavailable" else if (
    length(se) != n_parameters || any(!is.finite(se)) || any(se <= 0, na.rm = TRUE)) "issues_detected" else "no_issues_detected",
    expected = n_parameters, available = length(se),
    nonfinite = sum(!is.finite(se)), nonpositive = sum(is.finite(se) & se <= 0))
  if (standard_errors$status == "issues_detected") issues <- c(issues, "Standard errors are incomplete, nonfinite or nonpositive")

  component_sizes <- list(status = "not_applicable", table = data.frame(), detail = "Not a recognized nlpsem mixture")
  if (methods::is(object@expectation, "MxExpectationMixture")) {
    component_sizes$status <- "unavailable"
    component_sizes$detail <- "Model has no fitted results"
    if (fitted) {
      covariates <- if (is.null(object@matrices$weightsV)) NULL else {
        labels <- as.vector(object@matrices$weightsV@labels)
        sub("^data[.]", "", labels[!is.na(labels) & grepl("^data[.]", labels)])
      }
      if (!length(covariates)) covariates <- NULL
      posterior <- tryCatch(getPosterior(object, nClass = length(object@submodels),
                                        cluster_TIC = covariates), error = function(e) e)
      if (inherits(posterior, "error")) {
        component_sizes$detail <- conditionMessage(posterior)
        issues <- c(issues, "Mixture posterior sizes could not be calculated")
      } else {
        k <- ncol(posterior@prob)
        hard <- as.integer(table(factor(posterior@membership, levels = seq_len(k))))
        expected <- colSums(posterior@prob)
        component_sizes$table <- data.frame(component = names(object@submodels),
          assigned_n = hard, posterior_n = expected,
          assigned_proportion = hard / nrow(posterior@prob), stringsAsFactors = FALSE)
        component_sizes$status <- if (any(hard == 0L)) "empty_assignment_class" else "available"
        component_sizes$detail <- "Posterior counts do not establish class identification or solution stability"
        if (any(hard == 0L)) issues <- c(issues, "At least one mixture class has no maximum-posterior assignments")
      }
    }
  }
  list(optimizer = list(status = status, status_zero = if (is.na(status)) NA else status == 0L,
                         finite_likelihood = length(fit) == 1L && is.finite(fit)),
       covariance = covariance, cross_covariances = cross_covariances,
       hessian = hessian, standard_errors = standard_errors,
       component_sizes = component_sizes, issues = unique(issues),
       repeated_starts = "not_assessed",
       limitations = c("No refitting or statistical acceptance is performed",
         "Definition-variable covariances and other unsupported expectations are not fully checked",
         "Identification, global optimality, repeated-start stability, and interval calibration are not established",
         "Stored outputs can be stale if a fitted model was manually modified"))
}

# Recognize scalar-label blocks by their placement in S, not by their names or
# fitted values. A covariance block is a principal submatrix; a cross block
# occupies two disjoint variable sets and is assessed within the joint S.
.diagnostic_cross_covariance <- function(formula, labels) {
  scalar_grid <- function(expr) {
    if (is.symbol(expr)) return(matrix(as.character(expr), 1L, 1L))
    if (!is.call(expr) || !as.character(expr[[1L]]) %in% c("rbind", "cbind")) return(NULL)
    parts <- lapply(as.list(expr)[-1L], scalar_grid)
    if (!length(parts) || any(vapply(parts, is.null, logical(1)))) return(NULL)
    bind <- as.character(expr[[1L]])
    sizes <- vapply(parts, if (bind == "rbind") ncol else nrow, integer(1))
    if (length(unique(sizes)) != 1L) return(NULL)
    do.call(bind, parts)
  }
  grid <- scalar_grid(formula)
  if (is.null(grid) || !is.matrix(labels)) return(NULL)
  first <- which(!is.na(labels) & labels == grid[1L, 1L], arr.ind = TRUE)
  if (!nrow(first)) return(NULL)
  found <- list()
  for (i in seq_len(nrow(first))) {
    rows <- lapply(grid[, 1L], function(label) which(!is.na(labels[, first[i, 2L]]) & labels[, first[i, 2L]] == label))
    columns <- lapply(grid[1L, ], function(label) which(!is.na(labels[first[i, 1L], ]) & labels[first[i, 1L], ] == label))
    # Equality-constrained or otherwise ambiguous mappings are not guessed.
    if (any(lengths(rows) != 1L) || any(lengths(columns) != 1L)) next
    rows <- unlist(rows, use.names = FALSE)
    columns <- unlist(columns, use.names = FALSE)
    if (anyDuplicated(rows) || anyDuplicated(columns) ||
        !identical(unname(labels[rows, columns, drop = FALSE]), unname(grid))) next
    # A principal (or overlapping) block must not be exempted from checking.
    if (length(intersect(rows, columns))) return(NULL)
    found[[length(found) + 1L]] <- list(rows = rows, columns = columns)
  }
  if (!length(found)) return(NULL)
  result <- found[[1L]]
  result$rows <- if (is.null(rownames(labels))) as.character(result$rows) else rownames(labels)[result$rows]
  result$columns <- if (is.null(colnames(labels))) as.character(result$columns) else colnames(labels)[result$columns]
  result
}

# Inspect stored expressions, following references so indirect definition
# variables cannot make a cached covariance look constant. Unknown references
# are left to the unavailable-matrix handling above.
.diagnostic_has_definition_variables <- function(model, name, seen = character()) {
  if (name %in% seen) return(FALSE)
  seen <- c(seen, name)
  matrix <- model@matrices[[name]]
  if (!is.null(matrix)) {
    labels <- as.vector(matrix@labels)
    if (any(grepl("(^|[.])data[.]", labels), na.rm = TRUE)) return(TRUE)
    references <- sub("\\[.*$", "", labels[!is.na(labels) & grepl("\\[", labels)])
  } else {
    algebra <- model@algebras[[name]]
    if (is.null(algebra)) return(FALSE)
    references <- all.names(algebra@formula, functions = FALSE)
    if (any(grepl("(^|[.])data[.]", references))) return(TRUE)
  }
  references <- intersect(references, c(names(model@matrices), names(model@algebras)))
  any(vapply(references, function(ref) .diagnostic_has_definition_variables(model, ref, seen), logical(1)))
}

.diagnostic_matrix <- function(x, tolerance) {
  out <- list(dimension = 0L, status = "unavailable", minimum_diagonal = NA_real_,
              minimum_eigenvalue = NA_real_, scaled_minimum_eigenvalue = NA_real_,
              negative_diagonals = "", detail = "No nonempty square stored matrix")
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || nrow(x) != ncol(x)) return(out)
  out$dimension <- nrow(x)
  if (any(!is.finite(x))) {
    out$status <- "nonfinite"; out$detail <- "Matrix contains nonfinite values"; return(out)
  }
  diagonal <- diag(x)
  out$minimum_diagonal <- min(diagonal)
  names <- rownames(x)
  if (is.null(names)) names <- as.character(seq_len(nrow(x)))
  out$negative_diagonals <- paste(names[diagonal < 0], collapse = ", ")
  scale <- sqrt(abs(diagonal))
  zero <- scale == 0
  if (any(zero)) scale[zero] <- sqrt(pmax(apply(abs(x[zero, , drop = FALSE]), 1L, max), .Machine$double.xmin))
  normalized <- sweep(sweep(x, 1L, scale, "/"), 2L, scale, "/")
  if (any(!is.finite(normalized))) {
    out$status <- "nonfinite"; out$detail <- "Diagonal scaling overflowed"; return(out)
  }
  if (max(abs(normalized - t(normalized))) > tolerance * max(1, max(abs(normalized)))) {
    out$status <- "not_symmetric"; out$detail <- "Scaled matrix is not symmetric within tolerance"; return(out)
  }
  out$minimum_eigenvalue <- tryCatch(min(eigen(x, symmetric = TRUE, only.values = TRUE)$values), error = function(e) NA_real_)
  out$scaled_minimum_eigenvalue <- min(eigen(normalized, symmetric = TRUE, only.values = TRUE)$values)
  out$status <- if (any(diagonal < 0)) "negative_variance" else if (out$scaled_minimum_eigenvalue < -tolerance) {
    "indefinite"
  } else if (out$scaled_minimum_eigenvalue <= tolerance) "near_boundary" else "positive_definite"
  out$detail <- "Stored matrix; diagonal-scaled eigenvalue classification"
  out
}

#' Estimates and standard errors of free parameters selected by their exact labels
#'
#' @param model A fitted MxModel.
#' @param labels Character vector of raw parameter labels, in the order wanted.
#' @return A list with the named vectors \code{est} and \code{se}; a standard error that is not available is
#'   \code{NA}.
#' @keywords internal
#' @noRd
.free_by_label <- function(model, labels) {
  est <- model@output$estimate
  missing <- setdiff(labels, names(est))
  if (length(missing)) {
    stop("Parameter label(s) not found in the fitted model: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  se_all <- model@output$standardErrors
  se <- stats::setNames(rep(NA_real_, length(labels)), labels)
  if (length(se_all)) {
    se_all <- stats::setNames(se_all[, 1], rownames(se_all))
    have <- intersect(labels, names(se_all))
    se[have] <- se_all[have]
  }
  list(est = est[labels], se = se)
}

#' Append a single between-process output block
#'
#' @param output_list Existing list of output blocks.
#' @param psi A single between-process psi estimate or standard-error object.
#' @param res A single between-process residual estimate or standard-error object.
#' @keywords internal
#' @noRd
.append_between_output <- function(output_list, psi, res) {
  output_list[[length(output_list) + 1]] <- c(unlist(psi), unlist(res))
  output_list
}

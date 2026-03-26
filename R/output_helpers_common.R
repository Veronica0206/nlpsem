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

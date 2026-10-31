#' Reduce a block-structured covariance matrix to the free growth-factor block
#'
#' Keeps the first `n_keep` transformed growth factors from each trajectory block of
#' size `stride`, then returns the lower-triangular free values in the ordering used
#' by `mxPath(connect = "unique.pairs")`.
#'
#' @param cov_mat Symmetric covariance matrix.
#' @param n_traj Number of trajectories/processes.
#' @param n_keep Number of transformed growth factors to keep per trajectory.
#' @param stride Size of each trajectory block in the full covariance matrix.
#'
#' @return Numeric vector of lower-triangular covariance starts.
#' @keywords internal
#' @noRd
.multi_reduced_cov_starts <- function(cov_mat, n_traj, n_keep = 3L, stride = 4L){
  keep_idx <- unlist(lapply(seq_len(n_traj), function(traj){
    start <- (traj - 1L) * stride + 1L
    seq.int(start, length.out = n_keep)
  }))
  reduced <- cov_mat[keep_idx, keep_idx, drop = FALSE]
  t(reduced)[row(t(reduced)) >= col(t(reduced))]
}

# A unique index for the outcome pair (i, j), with i < j, in lexicographic
# order: (1, 2), (1, 3), ..., (2, 3), ... . Adding the two indices collides
# once there are four outcomes and silently overwrites covariance components.
#' @keywords internal
#' @noRd
.outcome_pair_index <- function(i, j, n_traj){
  (i - 1L) * n_traj - (i - 1L) * i / 2L + j - i
}

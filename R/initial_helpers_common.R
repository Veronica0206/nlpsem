#' @keywords internal
#' @noRd

## Helper: project a symmetric matrix to the nearest positive-definite matrix.
## Used to ensure residual covariance matrices (computed as Cov(GF) - B*Cov(X)*B')
## are valid when sampling variability causes small negative eigenvalues.
.ensure_pd <- function(mat, tol = 1e-6, label = ""){
  eig <- eigen(mat, symmetric = TRUE)
  if (any(eig$values < tol)){
    warning("Residual covariance matrix", if (nzchar(label)) paste0(" (", label, ")") else "",
            " had negative eigenvalues (min = ", round(min(eig$values), 8),
            "); projecting to nearest positive-definite matrix.", call. = FALSE)
    eig$values <- pmax(eig$values, tol)
    mat <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
    mat <- (mat + t(mat)) / 2
  }
  mat
}

## Helper: convert class proportions to softmax/logit starts relative to class 1.
.mixture_weight_starts <- function(prop_starts){
  prop_starts <- as.numeric(prop_starts)
  if (!length(prop_starts)) return(prop_starts)

  if (any(!is.finite(prop_starts)) || any(prop_starts <= 0)){
    warning("Mixture proportion starts must be positive and finite; using equal softmax starts.",
            call. = FALSE)
    prop_starts <- rep(1, length(prop_starts))
  }

  prop_starts <- prop_starts / sum(prop_starts)
  c(0, log(prop_starts[-1] / prop_starts[1]))
}

## Helper: intercept-only softmax starts from provisional class labels.
.cluster_tic_intercept_starts <- function(labels, nClass){
  counts <- tabulate(as.integer(labels), nbins = nClass) + 0.5
  logits <- .mixture_weight_starts(counts)
  cbind(logits)
}

## Helper: robust multinomial-logit starts for cluster-forming TICs.
.cluster_tic_weight_starts <- function(dat, cluster_TIC, nClass, shrink = 0.7){
  dat_nnet <- dat[dat$label %in% seq_len(nClass), c("label", cluster_TIC), drop = FALSE]
  dat_nnet <- dat_nnet[stats::complete.cases(dat_nnet), , drop = FALSE]

  intercept_only <- .cluster_tic_intercept_starts(dat_nnet$label, nClass)
  n_coef <- length(cluster_TIC) + 1
  fallback <- matrix(0, nrow = nClass, ncol = n_coef)
  fallback[, 1] <- intercept_only[, 1]

  if (!nrow(dat_nnet) || length(unique(dat_nnet$label)) < 2){
    warning("Cluster TIC starts had insufficient labeled observations; using intercept-only softmax starts.",
            call. = FALSE)
    return(fallback)
  }

  x_raw <- dat_nnet[, cluster_TIC, drop = FALSE]
  x_means <- vapply(x_raw, function(col) mean(col, na.rm = TRUE), numeric(1))
  x_sds <- vapply(x_raw, function(col) stats::sd(col, na.rm = TRUE), numeric(1))
  active <- is.finite(x_sds) & x_sds > 0

  if (!any(active)){
    warning("Cluster TIC starts had no varying predictors; using intercept-only softmax starts.",
            call. = FALSE)
    return(fallback)
  }

  x_scaled <- sweep(x_raw[, active, drop = FALSE], 2, x_means[active], FUN = "-")
  x_scaled <- sweep(x_scaled, 2, x_sds[active], FUN = "/")
  fit_dat <- data.frame(label = factor(dat_nnet$label, levels = seq_len(nClass)), x_scaled)

  mod_nnet <- tryCatch(
    suppressWarnings(nnet::multinom(label ~ ., data = fit_dat, trace = FALSE)),
    error = function(e) NULL
  )

  if (is.null(mod_nnet)){
    warning("Cluster TIC multinomial start model failed; using intercept-only softmax starts.",
            call. = FALSE)
    return(fallback)
  }

  coef_scaled <- summary(mod_nnet)$coefficient
  if (is.null(dim(coef_scaled))){
    coef_scaled <- matrix(coef_scaled, nrow = 1)
  } else{
    coef_scaled <- as.matrix(coef_scaled)
  }
  if (!nrow(coef_scaled)){
    return(fallback)
  }

  coef_unscaled <- matrix(0, nrow = nClass - 1, ncol = n_coef)
  coef_unscaled[, 1] <- coef_scaled[, 1, drop = TRUE]
  slope_scaled <- coef_scaled[, -1, drop = FALSE]
  slope_unscaled <- sweep(slope_scaled, 2, x_sds[active], FUN = "/")
  coef_unscaled[, 1] <- coef_unscaled[, 1] - slope_scaled %*% (x_means[active] / x_sds[active])
  coef_unscaled[, 1 + which(active)] <- slope_unscaled

  coef_unscaled[, 1] <- intercept_only[-1, 1] + shrink * (coef_unscaled[, 1] - intercept_only[-1, 1])
  coef_unscaled[, -1] <- coef_unscaled[, -1, drop = FALSE] * shrink

  out <- matrix(0, nrow = nClass, ncol = n_coef)
  out[-1, ] <- coef_unscaled
  out
}

## Helper: residual-correlation starts between longitudinal processes, from the per-person curve residuals
## that getUNI.GF() attaches to its growth factors ("residuals" attribute). Residuals are paired occasion by
## occasion over the records the two processes share. Centring the raw scores at the person mean (the 0.4
## heuristic) leaves the growth trend in the deviations, so two processes that grow together looked almost
## perfectly correlated (r > 0.9 on the shipped example data although the fitted residual correlation is 0.24);
## that start made the class-specific multivariate models unfittable. Falls back to 0.3 with a warning when
## a pair has too few jointly observed residuals.
.residual_correlation_starts <- function(res_list, labels = seq_along(res_list), fallback = 0.3, clamp = 0.9,
                                         records = NULL){
  n <- length(res_list)
  out <- numeric(0)
  if (n < 2) return(out)
  for (i in seq_len(n - 1)){
    for (j in (i + 1):n){
      ri <- res_list[[i]]; rj <- res_list[[j]]
      r <- NA_real_
      if (is.matrix(ri) && is.matrix(rj) && nrow(ri) == nrow(rj)){
        # Multivariate outcomes may have different record sets. Match the actual
        # occasions, just as the residual covariance paths in the model do.
        if (is.null(records)){
          common <- seq_len(min(ncol(ri), ncol(rj)))
          ii <- jj <- common
        } else{
          common <- intersect(records[[i]], records[[j]])
          ii <- match(common, records[[i]])
          jj <- match(common, records[[j]])
        }
        ri <- ri[, ii, drop = FALSE]; rj <- rj[, jj, drop = FALSE]
        ok <- is.finite(ri) & is.finite(rj)
        if (sum(ok) > 2) r <- suppressWarnings(stats::cor(ri[ok], rj[ok]))
      }
      if (!is.finite(r)){
        warning("Residual correlation start for ", labels[i], " and ", labels[j],
                " could not be estimated from the per-person curve residuals; falling back to ", fallback, ".",
                call. = FALSE)
        r <- fallback
      }
      out <- c(out, max(min(r, clamp), -clamp))
    }
  }
  out
}

## Helper: names for cluster-TIC multinomial coefficients in output.
.cluster_tic_param_names <- function(nClass, cluster_TIC){
  coef_labels <- c("0", seq_along(cluster_TIC))
  paste0("beta", rep(2:nClass, each = length(coef_labels)), rep(coef_labels, times = nClass - 1))
}

## Helper: damp non-intrinsic BLS multivariate covariance starts for class-specific fits.
## Single-group MGM can tolerate the data-driven transformed covariance, but mixture and
## multigroup submodels are much more sensitive because each class/group uses a smaller sample.
.stabilize_class_specific_mgm_starts <- function(starts, curveFun, intrinsic, shrink = 0.6){
  if (isTRUE(intrinsic) || !curveFun %in% c("bilinear spline", "BLS")) return(starts)
  if (is.null(starts$means) || is.null(starts$covMatrix) || !is.matrix(starts$covMatrix)) return(starts)

  n_traj <- length(starts$means)
  if (!n_traj) return(starts)

  keep_idx <- unlist(lapply(seq_len(n_traj), function(traj) ((traj - 1L) * 4L + 1L):((traj - 1L) * 4L + 3L)))
  if (max(keep_idx) > nrow(starts$covMatrix)) return(starts)

  reduced_cov <- starts$covMatrix[keep_idx, keep_idx, drop = FALSE] * shrink
  starts$covMatrix[keep_idx, keep_idx] <- .ensure_pd(reduced_cov, label = "class-specific fixed-knot multivariate starts")
  starts
}

## Helper: provisional class labels for mixture starts. k-means runs on the persons whose
## growth-factor starts are complete; the remaining persons are assigned to the nearest
## class centre over their finite coordinates (or to the largest class when none is finite).
## Classes are ordered by the mean of their centres so labels are reproducible across runs.
.provisional_class_labels <- function(growth_factor, nClass){
  gf <- as.matrix(growth_factor)
  complete <- stats::complete.cases(gf)
  if (sum(complete) < nClass){
    stop("Only ", sum(complete), " person(s) have complete growth-factor starts; at least ", nClass,
         " are needed to form the provisional classes of a mixture model.")
  }
  km <- kmeans(x = gf[complete, , drop = FALSE], centers = nClass, iter.max = 100)
  order_indices <- order(apply(km$centers, 1, mean))
  labels <- rep(NA_integer_, nrow(gf))
  labels[complete] <- match(km$cluster, order_indices)
  if (any(!complete)){
    centers <- km$centers[order_indices, , drop = FALSE]
    largest <- which.max(tabulate(labels[complete], nbins = nClass))
    for (i in which(!complete)){
      ok <- is.finite(gf[i, ])
      if (!any(ok)){
        labels[i] <- largest
      } else{
        d2 <- colSums((t(centers[, ok, drop = FALSE]) - gf[i, ok])^2)
        labels[i] <- which.min(d2)
      }
    }
  }
  labels
}

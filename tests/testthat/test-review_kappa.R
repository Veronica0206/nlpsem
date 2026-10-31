# Independent reference checks for review finding F7.
.review_labels <- function(counts) {
  cells <- which(counts > 0, arr.ind = TRUE)
  list(first = rep(cells[, 1], counts[cells]),
       second = rep(cells[, 2], counts[cells]))
}

# Numerical derivatives of kappa on the probability simplex avoid reusing the
# analytical gradient implemented by the package.
.review_kappa_reference <- function(counts, level) {
  p <- as.vector(counts / sum(counts))
  nc <- nrow(counts)
  statistic <- function(x) {
    m <- matrix(x / sum(x), nc)
    po <- sum(diag(m))
    pe <- sum(rowSums(m) * colSums(m))
    (po - pe) / (1 - pe)
  }
  step <- 1e-6
  g <- vapply(seq_along(p), function(i) {
    delta <- rep(0, length(p)); delta[i] <- step
    (statistic(p + delta) - statistic(p - delta)) / (2 * step)
  }, numeric(1))
  sampling_cov <- (diag(p) - tcrossprod(p)) / sum(counts)
  se <- sqrt(drop(t(g) %*% sampling_cov %*% g))
  kappa <- statistic(p)
  bounds <- kappa + c(-1, 1) * qnorm((1 + level)/2) * se
  paste0(round(kappa, 4), " (", round(bounds[1], 4), ", ",
         round(bounds[2], 4), ")")
}

test_that("kappa includes empirical marginal uncertainty for imbalanced labels", {
  counts <- matrix(c(95, 5, 15, 5), 2, byrow = TRUE)
  labels <- .review_labels(counts)
  expect_identical(getLatentKappa(labels$first, labels$second)@kappa_value,
                   "0.25 (0.0219, 0.4781)")
  expect_identical(getLatentKappa(labels$second, labels$first)@kappa_value,
                   "0.25 (0.0219, 0.4781)")
})

test_that("kappa intervals agree with independent numerical delta references", {
  cases <- list(matrix(c(95, 5, 15, 5), 2, byrow = TRUE),
                matrix(c(20, 5, 2, 8, 30, 4, 1, 6, 24), 3, byrow = TRUE),
                matrix(c(1, 19, 9, 1), 2, byrow = TRUE),
                matrix(c(0, 10, 0, 0, 5, 15, 0, 0, 20), 3, byrow = TRUE))
  for (counts in cases) {
    labels <- .review_labels(counts)
    for (level in c(.90, .95, .99)) {
      expected <- .review_kappa_reference(counts, level)
      expect_identical(getLatentKappa(labels$first, labels$second,
                                      conf.level = level)@kappa_value, expected)
    }
  }
})

#' @title Compute Posterior Probabilities, Cluster Assignments, and Model Entropy for a Longitudinal Mixture Model
#'
#' @description This function computes posterior probabilities, cluster assignments, and model entropy for a given mixture model
#' with a predefined number of classes. If the true labels are available, it can also compute the model accuracy.
#'
#' @param model A fitted mxModel object. Specifically, this should be the \code{mxOutput} slot from the output of \code{getMIX()}.
#' @param nClass A positive integer equal to the number of latent classes in the fitted model.
#' @param label A logical value indicating whether the data contains true labels, which are often available in a simulated data set.
#' Default is FALSE.
#' @param cluster_TIC A string or character vector representing the column name(s) for time-invariant covariate(s)
#' indicating cluster formations. Default is \code{NULL}, which recovers the fitted
#' cluster covariates from the model when present.
#' An explicit vector must match the fitted covariates in their original order.
#'
#' @return An object of class \code{postOutput}. Depending on the \code{label} argument, the object may contain the following slots:
#' \itemize{
#'   \item \code{prob}: A matrix of posterior probabilities.
#'   \item \code{membership}: A vector indicating class membership based on maximum posterior probability.
#'   \item \code{entropy}: The entropy of the model, a measure of uncertainty in class assignment.
#'   \item \code{accuracy} (optional): If \code{label = TRUE}, the model's accuracy based on true labels.
#' }
#' The content of these slots can be printed using the \code{printTable()} method for S4 objects.
#'
#' @details Probabilities are normalized in log space to avoid overflow or underflow
#' from class weights. Nonfinite or negative component likelihoods, and observations
#' with zero likelihood in every class, produce an informative error. Entropy is
#' numeric NA for a one-class model. With \code{label = TRUE}, the observed
#' \code{class} labels must already use the fitted class numbering; labels are not
#' permuted to maximize accuracy.
#'
#' @references
#' \itemize{
#'   \item {Peugh, J., & Fan, X. (2015). Enumeration Index Performance in Generalized Growth Mixture Models:
#' A Monte Carlo Test of Muthén's (2003) Hypothesis. Structural Equation Modeling: A Multidisciplinary Journal, 22(1), 115-131.
#' Routledge. \doi{10.1080/10705511.2014.919823}}
#'   \item {Lubke, G., & Muthén, B.O. (2007). Performance of Factor Mixture Models as a Function of Model Size,
#' Covariate Effects, and Class-Specific Parameters. Structural Equation Modeling: A Multidisciplinary Journal, 14(1), 26-47.
#' Routledge. \doi{10.1080/10705510709336735}}
#' }
#'
#' @seealso \code{\link{getMIX}}, \code{\link{getLatentKappa}}
#'
#' @export
#'
#' @examples
#' # Fitted mixture models shipped with the package
#' load(system.file("extdata", "getMIX_examples.RData", package = "nlpsem"))
#' # Three-class bilinear spline mixture without cluster-forming covariates
#' post1 <- getPosterior(model = Math_BLS_LGCM3@mxOutput, nClass = 3)
#' head(post1@prob)
#' post1@entropy
#' # Three-class bivariate mixture with cluster-forming covariates
#' post2 <- getPosterior(
#'   model = RM_BLS_PLGCM3@mxOutput, nClass = 3, cluster_TIC = c("gx1", "gx2")
#' )
#' table(post2@membership)
#'
#' @importFrom OpenMx mxEval
#' @importFrom methods new
#'
getPosterior <- function(model, nClass, label = FALSE, cluster_TIC = NULL){
  if (!methods::is(model, "MxModel")) {
    stop("'model' must be the fitted mxOutput from getMIX().", call. = FALSE)
  }
  if (!is.numeric(nClass) || is.complex(nClass) || length(nClass) != 1L || !is.finite(nClass) ||
      nClass < 1 || nClass != floor(nClass)) {
    stop("'nClass' must be one positive integer.", call. = FALSE)
  }
  if (!is.logical(label) || length(label) != 1L || is.na(label)) {
    stop("'label' must be TRUE or FALSE.", call. = FALSE)
  }
  subs <- model@submodels
  if (!length(subs) || nClass != length(subs)) {
    stop("'nClass' must equal the number of fitted mixture classes (", length(subs), ").", call. = FALSE)
  }
  observed <- model$data$observed
  if (is.null(observed) || !(is.data.frame(observed) || is.matrix(observed)) || !nrow(observed)) {
    stop("The mixture model must contain nonempty observed raw data at the parent level.", call. = FALSE)
  }
  observed <- as.data.frame(observed)
  n <- nrow(observed)
  has_covariates <- !is.null(model$classbeta)
  expected_covariates <- if (has_covariates) {
    definition_labels <- as.vector(model$weightsV$labels)[-1L]
    if (!length(definition_labels) || anyNA(definition_labels) ||
        any(!startsWith(definition_labels, "data."))) {
      stop("The model has malformed cluster-covariate definitions.", call. = FALSE)
    }
    sub("^data[.]", "", definition_labels)
  } else character()
  if (!is.null(cluster_TIC) &&
      (!is.character(cluster_TIC) || !length(cluster_TIC) || anyNA(cluster_TIC) ||
       any(!nzchar(cluster_TIC)) || anyDuplicated(cluster_TIC))) {
    stop("'cluster_TIC' must be NULL or distinct nonmissing column names.", call. = FALSE)
  }
  cluster_TIC <- unname(cluster_TIC)
  if (!is.null(cluster_TIC) && !identical(cluster_TIC, expected_covariates)) {
    stop("'cluster_TIC' must match the fitted cluster covariates in their original order.", call. = FALSE)
  }
  cluster_TIC <- expected_covariates
  if (has_covariates) {
    if (any(!cluster_TIC %in% names(observed)) ||
        any(!vapply(observed[cluster_TIC], is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(observed[cluster_TIC])))) {
      stop("Fitted cluster covariates must be present, numeric and finite in the model data.", call. = FALSE)
    }
    beta <- model$classbeta$values
    if (!is.numeric(beta) || !identical(dim(beta), c(as.integer(nClass), length(cluster_TIC) + 1L)) ||
        any(!is.finite(beta))) {
      stop("The fitted class-coefficient matrix has invalid dimensions or values.", call. = FALSE)
    }
    logits <- cbind(1, as.matrix(observed[cluster_TIC])) %*% t(beta)
  } else {
    weights <- model$weights$values
    if (!is.numeric(weights) || length(weights) != nClass || any(!is.finite(weights))) {
      stop("The model must contain one finite log weight per mixture class.", call. = FALSE)
    }
    logits <- matrix(rep(as.numeric(weights), each = n), nrow = n, ncol = nClass)
  }
  if (any(!is.finite(logits))) {
    stop("The fitted class logits are nonfinite; check the covariates and coefficients.", call. = FALSE)
  }
  likelihood <- vapply(seq_along(subs), function(k) {
    values <- as.numeric(mxEval(objective, subs[[k]]))
    if (length(values) != n || any(!is.finite(values)) || any(values < 0)) {
      stop("Class ", k, " has invalid component likelihoods; expected one finite nonnegative value per observation.",
           call. = FALSE)
    }
    values
  }, numeric(n))
  likelihood <- matrix(likelihood, nrow = n, ncol = nClass)
  # The prior softmax denominator cancels in Bayes' rule. Work directly with
  # log likelihood + centered logits, then normalize once at the posterior level.
  log_mass <- log(likelihood) + sweep(logits, 1L, apply(logits, 1L, max), "-")
  row_max <- apply(log_mass, 1L, max)
  if (any(!is.finite(row_max))) {
    bad <- which(!is.finite(row_max))
    stop("Posterior probabilities are undefined for observation(s) ",
         paste(utils::head(bad, 10L), collapse = ", "),
         ": all class likelihoods are zero or their weighted likelihoods are not representable. Check the fit and scaling.",
         call. = FALSE)
  }
  mass <- exp(sweep(log_mass, 1L, row_max, "-"))
  prob <- mass / rowSums(mass)
  membership <- as.numeric(max.col(prob, ties.method = "first"))
  entropy <- if (nClass == 1L) NA_real_ else {
    positive <- prob > 0
    1 + sum(prob[positive] * log(prob[positive])) / (n * log(nClass))
  }
  accuracy <- numeric()
  if (label) {
    if (is.null(observed$class)) {
      stop("Column 'class' not found in the dataset. Set label = FALSE or ensure the data contains a 'class' column.",
           call. = FALSE)
    }
    if (anyNA(observed$class)) stop("True class labels must not contain missing values.", call. = FALSE)
    accuracy <- mean(observed$class == membership)
  }
  new("postOutput", prob = prob, membership = membership, entropy = entropy, accuracy = accuracy)
}

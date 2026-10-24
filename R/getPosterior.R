#' @title Compute Posterior Probabilities, Cluster Assignments, and Model Entropy for a Longitudinal Mixture Model
#'
#' @description This function computes posterior probabilities, cluster assignments, and model entropy for a given mixture model
#' with a predefined number of classes. If the true labels are available, it can also compute the model accuracy.
#'
#' @param model A fitted mxModel object. Specifically, this should be the \code{mxOutput} slot from the output of \code{getMIX()}.
#' @param nClass An integer representing the predefined number of latent classes in the model.
#' @param label A logical value indicating whether the data contains true labels, which are often available in a simulated data set.
#' Default is FALSE.
#' @param cluster_TIC A string or character vector representing the column name(s) for time-invariant covariate(s)
#' indicating cluster formations. Default is \code{NULL}, indicating that no such time-invariant covariates are present
#' in the model.
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
  # Extract submodels from the input model
  subs <- model@submodels
  # Define a function to get the likelihood for a submodel
  getLIK <- function(num){
    return(mxEval(objective, subs[[num]]))
  }
  # Compute the likelihood for each submodel
  likelihood <- sapply(seq_along(names(subs)), getLIK)
  # Initialize a probability matrix
  prob <- matrix(0, nrow = nrow(model$data$observed), ncol = nClass)
  # Calculate class probabilities based on model weights or time-invariant covariates that inform class formation
  if (is.null(cluster_TIC)){
    # Compute class probabilities based on model weights
    log_cp <- model$weights$values
    classProbs <- exp(log_cp)/sum(exp(log_cp))
    # Compute posterior probabilities for each observation (vectorized)
    unnorm <- sweep(likelihood, 2, classProbs, `*`)
    prob <- unnorm / rowSums(unnorm)
  }
  else {
    # Compute class probabilities based on time-invariant covariates that inform class formation
    oneX <- as.matrix(cbind(matrix(1, nrow = nrow(model$data$observed), ncol = 1),
                            model$data$observed[, cluster_TIC]))
    prop <- exp(oneX %*% t(model$classbeta$values))
    classProbs <- prop/apply(prop, 1, sum)
    # Compute posterior probabilities for each observation (vectorized)
    unnorm <- classProbs * likelihood
    prob <- unnorm / rowSums(unnorm)
  }
  # Assign cluster membership based on maximum posterior probability
  membership <- apply(prob, 1, which.max)
  # Calculate entropy
  if (nClass == 1) {
    entropy <- NA
  } else {
    prob_safe <- pmax(prob, .Machine$double.xmin)
    entropy <- 1 - sum(-prob_safe * log(prob_safe))/(nrow(model$data$observed) * log(nClass))
  }
  if (label){
    if (is.null(model$data$observed$class)) {
      stop("Column 'class' not found in the dataset. Set label = FALSE or ensure the data contains a 'class' column.")
    }
    true_labels <- model$data$observed$class
    # Calculate accuracy
    accuracy <- sum(true_labels == membership)/length(true_labels)
    postOut <- new("postOutput", prob = prob, membership = membership, entropy = entropy, accuracy = accuracy)
  }
  else {
    postOut <- new("postOutput", prob = prob, membership = membership, entropy = entropy)
  }
  return(postOut)
}

#' @title Compute Latent Kappa Coefficient for Agreement between Two Latent Label Sets
#'
#' @description This function calculates the latent kappa, a measure of agreement between two sets of latent categorical
#' labels. It also computes the confidence interval and provides a qualitative interpretation of the agreement level.
#' Categories are matched by their label values, independently of factor-level order.
#' Latent class numbers have no intrinsic correspondence across models: align their meanings
#' explicitly before calling this function. The function does not permute classes to maximize agreement.
#'
#' @param label1 A factor vector representing the first set of latent categorical labels.
#' @param label2 A factor vector representing the second set of latent categorical labels.
#' @param conf.level A numeric value representing the confidence level for the confidence interval of the kappa statistic.
#' The default value is \code{0.95}.
#'
#' @return An object of class \code{KappaOutput} with the following slots:
#' \itemize{
#'   \item \code{kappa_value}: A string representing the kappa statistic along with its confidence interval.
#'   \item \code{judgment}: A string describing the level of agreement, such as "Perfect Agreement", "Slight Agreement", etc.
#' }
#' The content of these slots can be printed using the \code{printTable()} method for S4 objects.
#'
#' @references
#' \itemize{
#'   \item {Dumenci, L. (2011). The Psychometric Latent Agreement Model (PLAM) for Discrete Latent Variables Measured by Multiple
#'   Items. Organizational Research Methods, 14(1), 91-115. SAGE Publications.
#'   \doi{10.1177/1094428110374649}}
#'   \item {Landis, J., & Koch, G. (1977). The Measurement of Observer Agreement for Categorical Data. Biometrics, 33(1), 159-174.
#'   \doi{10.2307/2529310}}
#'   \item {Agresti, A. (2012). Models for Matched Pairs. In Categorical Data Analysis (pp. 413-454). Wiley.}
#' }
#'
#' @seealso \code{\link{getPosterior}}, \code{\link{getMIX}}
#'
#' @export
#'
#' @examples
#' # Class memberships from two fitted three-class mixtures shipped with the package
#' load(system.file("extdata", "getMIX_examples.RData", package = "nlpsem"))
#' label1 <- getPosterior(model = Math_BLS_LGCM3@mxOutput, nClass = 3)@membership
#' label2 <- getPosterior(
#'   model = RM_BLS_PLGCM3@mxOutput, nClass = 3, cluster_TIC = c("gx1", "gx2")
#' )@membership
#' kappa <- getLatentKappa(label1 = label1, label2 = label2)
#' printTable(kappa)
#'
#' @importFrom stats complete.cases qnorm pnorm
#' @importFrom methods new
#'
getLatentKappa <- function(label1, label2, conf.level = 0.95){
  .ci_probs(conf.level)
  if (length(label1) != length(label2)){
    stop("Lengths of 'label1' and 'label2' differ!")
  }
  # Remove both if any missing in label1 or label2
  complete <- complete.cases(label1, label2)
  label1 <- as.character(label1[complete])
  label2 <- as.character(label2[complete])
  if (length(unique(label1)) < 2 || length(unique(label2)) < 2){
    stop("Both label1 and label2 must have at least 2 levels!")
  }
  categories <- union(label1, label2)
  tbl <- table(factor(label1, levels = categories), factor(label2, levels = categories))
  N <- sum(tbl)
  Po <- sum(diag(tbl))/N
  Pe <- sum(rowSums(tbl) * colSums(tbl)/N)/N
  kappa <- (Po - Pe)/(1 - Pe)
  if (kappa >= 1){
    judge <- "Perfect Agreement"
  }
  else if (kappa <= 0){
    judge <- "No Agreement"
  }
  else {
    JUDGEMENT <- c("Slight Agreement", "Fair Agreement",
                   "Moderate Agreement", "Substantial Agreement", "Almost Perfect Agreement")
    judge <- JUDGEMENT[min(which(kappa < seq(0.2, 1, 0.2)))]
  }
  seK0 <- sqrt(Pe/(N * (1 - Pe)))
  seK <- sqrt(Po * (1 - Po)/(N * (1 - Pe)^2))
  norm.pp <- qnorm(1 - (1 - conf.level)/2)
  Z <- kappa/seK0
  p.v <- 1 - pnorm(Z)
  kappaL <- round(kappa - norm.pp * seK, 4)
  kappaU <- round(kappa + norm.pp * seK, 4)
  KappaOut <- new("KappaOutput",
                  kappa_value = paste0(round(kappa, 4), " (", kappaL, ", ", kappaU, ")"),
                  judgment = judge)
  return(KappaOut)
}

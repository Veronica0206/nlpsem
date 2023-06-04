#' @title Compute Latent Kappa Coefficient for Agreement between Two Latent Label Sets
#'
#' @description This function calculates the latent kappa, a measure of agreement between two sets of latent categorical
#' labels. It also computes the confidence interval and provides a qualitative interpretation of the agreement level.
#'
#' @param label1 A factor vector representing the first set of latent categorical labels.
#' @param label2 A factor vector representing the second set of latent categorical labels.
#' @param conf.level A numeric value representing the confidence level for the confidence interval of the kappa statistic.
#' The default value is \code{0.95}.
#'
#' @return A list with two elements: The first is a string that provides the kappa statistic along with its confidence
#' interval. The second is a string describing the level of agreement (such as "Perfect Agreement", "Slight Agreement", etc.).
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
#' @export
#'
#' @examples
#' OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
#' data("RMS_dat")
#' RMS_dat0 <- RMS_dat
#' baseT <- RMS_dat0$T1
#' RMS_dat0$T1 <- RMS_dat0$T1 - baseT
#' RMS_dat0$T2 <- RMS_dat0$T2 - baseT
#' RMS_dat0$T3 <- RMS_dat0$T3 - baseT
#' RMS_dat0$T4 <- RMS_dat0$T4 - baseT
#' RMS_dat0$T5 <- RMS_dat0$T5 - baseT
#' RMS_dat0$T6 <- RMS_dat0$T6 - baseT
#' RMS_dat0$T7 <- RMS_dat0$T7 - baseT
#' RMS_dat0$T8 <- RMS_dat0$T8 - baseT
#' RMS_dat0$T9 <- RMS_dat0$T9 - baseT
#' RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
#' RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
#' RMS_dat0$gx1 <- scale(RMS_dat0$INCOME)
#' RMS_dat0$gx2 <- scale(RMS_dat0$EDU)
#' \donttest{
#' ## Fit a growth mixture model with no TICs
#' set.seed(20191029)
#' MIX_BLS_LGCM_r <-  getMIX(
#'   dat = RMS_dat0, prop_starts = c(0.33, 0.34, 0.33), sub_Model = "LGCM",
#'   cluster_TIC = NULL, y_var = "M", t_var = "T", records = 1:9,
#'   curveFun = "BLS", intrinsic = FALSE, res_scale = list(0.3, 0.3, 0.3),
#'   growth_TIC = NULL, tries = 10
#' )
#' ## Membership of each individual from growth mixture model with no TICs
#' label1 <- getPosterior(
#'   model = MIX_BLS_LGCM_r, nClass = 3, label = FALSE, cluster_TIC = NULL
#' )
#' set.seed(20191029)
#' ## Fit a growth mixture model with growth TICs and cluster TICs
#' MIX_BLS_LGCM.TIC_r <-  getMIX(
#'   dat = RMS_dat0, prop_starts = c(0.33, 0.34, 0.33), sub_Model = "LGCM",
#'   cluster_TIC = c("gx1", "gx2"), y_var = "M", t_var = "T", records = 1:9,
#'   curveFun = "BLS", intrinsic = FALSE, res_scale = list(0.3, 0.3, 0.3),
#'   growth_TIC = c("ex1", "ex2"), tries = 10
#' )
#' ## Membership of each individual from growth mixture model with growth TICs and cluster TICs
#' label2 <- getPosterior(
#'   model = MIX_BLS_LGCM.TIC_r, nClass = 3, label = FALSE,
#'   cluster_TIC = c("gx1", "gx2")
#' )
#' ## Calcualte membership between two sets of labels
#' getLatentKappa(label1 = label1$membership, label2 = label2$membership)
#' }
#'
#' @importFrom stats complete.cases qnorm pnorm
#'
getLatentKappa <- function(label1, label2, conf.level = 0.95){
  if (length(label1) != length(label2)){
    stop("Lengths of 'label1' and 'label2' differ!")
  }
  # Remove both if any missing in label1 or label2
  complete <- complete.cases(label1, label2)
  label1 <- factor(label1[complete])
  label2 <- factor(label2[complete])
  if (nlevels(label1) < 2 || nlevels(label2) < 2){
    stop("Both label1 and label2 must have at least 2 levels!")
  }
  tbl <- table(label1, label2)
  if (nrow(tbl) != ncol(tbl)){
    stop("Levels for 2 dimensions are different!")
  }
  N <- sum(tbl)
  Po <- sum(diag(tbl))/N
  Pe <- sum(rowSums(tbl) * colSums(tbl)/N)/N
  kappa <- (Po - Pe)/(1 - Pe)
  if (kappa == 1){
    judge <- "Perfect Agreement"
  }
  else {
    JUDGEMENT <- c("No Agreement", "Slight Agreement", "Fair Agreement",
                   "Moderate Agreement", "Substantial Agreement", "Almost Perfect Agreement")
    judge <- JUDGEMENT[min(which(kappa < seq(0, 1, 0.2)))]
  }
  seK0 <- sqrt(Pe/(N * (1 - Pe)))
  seK <- sqrt(Po * (1 - Po)/(N * (1 - Pe)^2))
  norm.pp <- qnorm(1 - (1 - conf.level)/2)
  Z <- kappa/seK0
  p.v <- 1 - pnorm(Z)
  kappaL <- round(kappa - norm.pp * seK, 4)
  kappaU <- round(kappa + norm.pp * seK, 4)
  kappa_out <- list(paste0(kappa, " (", kappaL, ", ", kappaU, ")"), judge)
  return(kappa_out)
}

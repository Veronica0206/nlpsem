% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSummary.R
\name{getSummary}
\alias{getSummary}
\title{Summarize Model Fit Statistics for Fitted Models}
\usage{
getSummary(model_list, HetModels = FALSE)
}
\arguments{
\item{model_list}{A list of fitted mxModel objects. Specifically, each element of the list should be the \code{mxOutput}
slot from the result returned by one of the estimation functions provided by this package.}

\item{HetModels}{A logical flag indicating whether a mixture model or a multiple group model is included in the list.
If set to \code{TRUE}, the function can also be used for the enumeration process, allowing the determination of the
optimal number of latent classes based on model fit statistics such as BIC. The default value is \code{FALSE}.}
}
\value{
A data frame summarizing model fit statistics (number of parameters, estimated likelihood, AIC, BIC, etc.)
for each model.
}
\description{
This function summarizes the model fit statistics for a list of fitted models. The summary includes the
number of parameters, estimated likelihood (-2ll), AIC, BIC, and other relevant statistics.
}
\examples{
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
# Load ECLS-K (2011) data
data("RMS_dat")
RMS_dat0 <- RMS_dat
# Re-baseline the data so that the estimated initial status is for the starting point of the study
baseT <- RMS_dat0$T1
RMS_dat0$T1 <- RMS_dat0$T1 - baseT
RMS_dat0$T2 <- RMS_dat0$T2 - baseT
RMS_dat0$T3 <- RMS_dat0$T3 - baseT
RMS_dat0$T4 <- RMS_dat0$T4 - baseT
RMS_dat0$T5 <- RMS_dat0$T5 - baseT
RMS_dat0$T6 <- RMS_dat0$T6 - baseT
RMS_dat0$T7 <- RMS_dat0$T7 - baseT
RMS_dat0$T8 <- RMS_dat0$T8 - baseT
RMS_dat0$T9 <- RMS_dat0$T9 - baseT
\donttest{
# Fit bilinear spline growth model with fix knot
## Single group model
BLS_LGCM1 <- getLGCM(
  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
  records = 1:9, res_scale = 0.1
  )
getSummary(model_list = list(BLS_LGCM1@mxOutput), HetModels = FALSE)
## Mixture model with two latent classes
set.seed(20191029)
BLS_LGCM2 <-  getMIX(
  dat = RMS_dat0, prop_starts = c(0.45, 0.55), sub_Model = "LGCM", cluster_TIC = NULL,
  y_var = "M", t_var = "T", records = 1:9, curveFun = "BLS", intrinsic = FALSE,
  res_scale = list(0.3, 0.3), growth_TIC = NULL, tries = 10
  )
## Mixture model with three latent classes
set.seed(20191029)
BLS_LGCM3 <-  getMIX(
  dat = RMS_dat0, prop_starts = c(0.33, 0.34, 0.33), sub_Model = "LGCM", cluster_TIC = NULL,
  y_var = "M", t_var = "T", records = 1:9, curveFun = "BLS", intrinsic = FALSE,
  res_scale = list(0.3, 0.3, 0.3), growth_TIC = NULL, tries = 10
  )

getSummary(model_list = list(BLS_LGCM1@mxOutput, BLS_LGCM2@mxOutput, BLS_LGCM3@mxOutput),
  HetModels = TRUE)
}

}

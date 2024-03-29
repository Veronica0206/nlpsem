% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getLCSM.R
\name{getLCSM}
\alias{getLCSM}
\title{Fit a Latent Change Score Model with a Time-invariant Covariate (If Any)}
\usage{
getLCSM(
  dat,
  t_var,
  y_var,
  curveFun,
  intrinsic = TRUE,
  records,
  growth_TIC = NULL,
  starts = NULL,
  res_scale = NULL,
  tries = NULL,
  OKStatus = 0,
  jitterD = "runif",
  loc = 1,
  scale = 0.25,
  paramOut = FALSE,
  names = NULL
)
}
\arguments{
\item{dat}{A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
repeated measurements and occasions, and time-invariant covariates (TICs) if any.}

\item{t_var}{A string specifying the prefix of the column names corresponding to the time variable at each study wave.}

\item{y_var}{A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.}

\item{curveFun}{A string specifying the functional form of the growth curve. Supported options for latent change score
models include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}), \code{"Jenss-Bayley"}
(or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}).}

\item{intrinsic}{A logical flag indicating whether to build an intrinsically nonlinear longitudinal model. Default is
\code{TRUE}.}

\item{records}{A numeric vector specifying indices of the study waves.}

\item{growth_TIC}{A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing to the
variability of growth factors if any. Default is \code{NULL}, indicating no growth TICs are included in the model.}

\item{starts}{A list containing initial values for the parameters. Default is \code{NULL}, indicating no user-specified
initial values.}

\item{res_scale}{A numeric value representing the scaling factor for the initial calculation of the residual variance. This
value should be between \code{0} and \code{1}, exclusive. By default, this is \code{NULL}, as it is unnecessary when the
user specifies the initial values using the \code{starts} argument.}

\item{tries}{An integer specifying the number of additional optimization attempts. Default is \code{NULL}.}

\item{OKStatus}{An integer (vector) specifying acceptable status codes for convergence. Default is \code{0}.}

\item{jitterD}{A string specifying the distribution for jitter. Supported values are: \code{"runif"} (uniform
distribution), \code{"rnorm"} (normal distribution), and \code{"rcauchy"} (Cauchy distribution). Default is \code{"runif"}.}

\item{loc}{A numeric value representing the location parameter of the jitter distribution. Default is \code{1}.}

\item{scale}{A numeric value representing the scale parameter of the jitter distribution. Default is \code{0.25}.}

\item{paramOut}{A logical flag indicating whether to output the parameter estimates and standard errors. Default is \code{FALSE}.}

\item{names}{A character vector specifying parameter names. Default is \code{NULL}.}
}
\value{
An object of class \code{myMxOutput}. Depending on the \code{paramOut} argument, the object may contain the following slots:
\itemize{
  \item \code{mxOutput}: This slot contains the fitted latent change score model. A summary of this model can be obtained using the
  \code{ModelSummary()} function.
  \item \code{Estimates} (optional): If \code{paramOut = TRUE}, a data frame with parameter estimates and standard errors. The content
  of this slot can be printed using the \code{printTable()} method for S4 objects.
}
}
\description{
This function fits a latent change score model with or without time-invariant covariates to the provided data.
It manages model setup, optimization, and if requested, outputs parameter estimates and standard errors.
}
\examples{
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
# Load ECLS-K (2011) data
data("RMS_dat")
RMS_dat0 <- RMS_dat
# Re-baseline the data so that the estimated initial status is for the starting point of the study
baseT <- RMS_dat0$T1
RMS_dat0$T1 <- (RMS_dat0$T1 - baseT)/12
RMS_dat0$T2 <- (RMS_dat0$T2 - baseT)/12
RMS_dat0$T3 <- (RMS_dat0$T3 - baseT)/12
RMS_dat0$T4 <- (RMS_dat0$T4 - baseT)/12
RMS_dat0$T5 <- (RMS_dat0$T5 - baseT)/12
RMS_dat0$T6 <- (RMS_dat0$T6 - baseT)/12
RMS_dat0$T7 <- (RMS_dat0$T7 - baseT)/12
RMS_dat0$T8 <- (RMS_dat0$T8 - baseT)/12
RMS_dat0$T9 <- (RMS_dat0$T9 - baseT)/12
# Standardized time-invariant covariates
RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)

\donttest{
# Fit nonparametric change score model for reading development
## Fit model
NonP_LCSM <- getLCSM(
  dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "nonparametric",
  intrinsic = FALSE, records = 1:9, res_scale = 0.1
  )
}

}
\references{
\itemize{
  \item {Liu, J., & Perera, R. A. (2023). Estimating Rate of Change for Nonlinear Trajectories in the Framework of Individual Measurement
  Occasions: A New Perspective on Growth Curves. Behavior Research Methods. \doi{10.3758/s13428-023-02097-2}}
  \item {Liu, J. (2022). "Jenss–Bayley Latent Change Score Model With Individual Ratio of the Growth Acceleration in the Framework
  of Individual Measurement Occasions." Journal of Educational and Behavioral Statistics, 47(5), 507–543.
  \doi{10.3102/10769986221099919}}
  \item {Grimm, K. J., Zhang, Z., Hamagami, F., & Mazzocco, M. (2013). "Modeling Nonlinear Change via Latent Change and Latent
  Acceleration Frameworks: Examining Velocity and Acceleration of Growth Trajectories." Multivariate Behavioral Research, 48(1),
  117-143. \doi{10.1080/00273171.2012.755111}}
}
}

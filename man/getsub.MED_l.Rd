% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIX.SUBMODEL.MED_helper.R
\name{getsub.MED_l}
\alias{getsub.MED_l}
\title{Define Longitudinal Mediation Models as Class-specific Models (Submodels) for a Longitudinal Mixture Model}
\usage{
getsub.MED_l(
  dat,
  nClass,
  t_var,
  records,
  y_var,
  curveFun,
  m_var,
  x_var,
  x_type,
  starts,
  res_cor
)
}
\arguments{
\item{dat}{A wide-format data frame, with each row corresponding to a unique ID. It contains the observed variables with
repeated measurements and occasions for each longitudinal process, and time-invariant covariates (TICs) if any. It takes
the value passed from \code{getMIX()}.}

\item{nClass}{An integer specifying the number of latent classes for the mixture model. It takes the value passed
from \code{getMIX()}.}

\item{t_var}{A vector of strings, with each element representing the prefix for column names related to the time variable for the
corresponding longitudinal variable at each study wave. It takes the value passed from \code{getMIX()}.}

\item{records}{A list of numeric vectors, with each vector specifying the indices of the observed study waves for
the corresponding longitudinal variable. It takes the value passed from \code{getMIX()}.}

\item{y_var}{A string specifying the prefix of the column names corresponding to the outcome variable at each study wave. It takes the value
passed from \code{getMIX()}.}

\item{curveFun}{A string specifying the functional form of the growth curve. Supported options include: "linear" (or "LIN"),
and "bilinear spline" (or "BLS"). It takes the value passed from \code{getMIX()}.}

\item{m_var}{A string specifying the prefix of the column names corresponding to the mediator variable at each study wave.
It takes the value passed from \code{getMIX()}.}

\item{x_var}{A string specifying the baseline predictor if \code{x_type = "baseline"}, or the prefix of the column names
corresponding to the predictor variable at each study wave if \code{x_type = "longitudinal"}. It takes the value passed from \code{getMIX()}.}

\item{x_type}{A string indicating the type of predictor variable used in the model. Supported values are \code{"baseline"}
and \code{"longitudinal"}. It takes the value passed from \code{getMIX()}.}

\item{starts}{A list of initial values for the parameters, either takes the value passed from \code{getMIX()} or derived by the
helper function \code{getMIX.initial()}.}

\item{res_cor}{A numeric value or vector for user-specified residual correlation between any two longitudinal processes to calculate
the corresponding initial value. It takes the value passed from \code{getMIX()}.}
}
\value{
A list of manifest and latent variables and paths for an mxModel object.
}
\description{
This function defines longitudinal mediation models as class-specific models (submodels) for a longitudinal mixture model.
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TVC.info_helper.R
\name{getTVC.info}
\alias{getTVC.info}
\title{Get the Time-Varying Covariate (TVC) Information for a One-group Longitudinal Model with Time-varying Covariate}
\usage{
getTVC.info(y_var, records, growth_TIC, TVC, decompose, starts)
}
\arguments{
\item{y_var}{A string specifying the prefix of the column names corresponding to the outcome variable at each study wave.
It takes the value passed from \code{getTVCmodel()}.}

\item{records}{A numeric vector specifying the indices of the observed study waves. It takes the value passed from \code{getTVCmodel()}.}

\item{growth_TIC}{A string or character vector specifying the column name(s) of time-invariant covariate(s) that account for the
variability of growth factors, if any. Default is \code{NULL}, indicating no growth TICs present in the model. It takes the value
passed from \code{getTVCmodel()}.}

\item{TVC}{A string specifying the prefix of the column names corresponding to the time-varying covariate at each study wave. It
takes the value passed from \code{getTVCmodel()}.}

\item{decompose}{An integer specifying the decomposition option for temporal states. Supported values include \code{0} (no
decomposition), \code{1} (decomposition with interval-specific slopes as temporal states), \code{2} (decomposition with interval-
specific changes as temporal states), and \code{3} (decomposition with change-from-baseline as temporal states). It takes the
value passed from \code{getTVCmodel()}.}

\item{starts}{A list of initial values for the parameters, either takes the value passed from \code{getTVCmodel()} or
derived by the helper function \code{getTVC.initial()}.}
}
\value{
A list containing two elements: X_PARAM and KAPPA. X_PARAM is a list of OpenMx paths and parameters for the
TVC, and KAPPA is an OpenMx path for the temporal effect of the TVC on the corresponding longitudinal outcome.
}
\description{
This function constructs the OpenMx model paths and parameters for a TVC and its relationship with the parameters
related to growth factors of a longitudinal outcome.
}
\keyword{internal}

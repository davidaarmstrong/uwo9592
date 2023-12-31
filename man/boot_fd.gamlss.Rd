% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uwo9592.R
\name{boot_fd.gamlss}
\alias{boot_fd.gamlss}
\title{Bootstraps (percentile) confidence intervals for first differences from a GAMLSS model.}
\usage{
boot_fd.gamlss(
  obj,
  data,
  R = 1000,
  what = c("mu", "sigma", "nu", "tau"),
  terms,
  progress = TRUE,
  ...
)
}
\arguments{
\item{obj}{A fitted `gamlss` object.}

\item{data}{The original data used to fit `obj`.}

\item{R}{Number of bootstrap iterations}

\item{what}{Which moment of the distribution to predict - only one allowed at a time.}

\item{terms}{Named list where the name corresponds to a variable in the model.  The
values of the list should be either a vector of two numbers or a function.  See details below.}

\item{progress}{Logical indicating whether a progress bar should be printed.}

\item{...}{Currently not implemented.}
}
\description{
Bootstraps (percentile) confidence intervals for first differences from a GAMLSS model.
}
\details{
If `terms` is two values, then those two values will be used for all observations
and the average first difference will be calculated.  If the value is a character string,
that should refer to a function whose only required argument is `x` and that `x` should be the
variable named in the list.  In that case, the variable - the evaluated function and the variable
plust the evaluated function will be used.  For example if `terms = list(x = "sd")`, then
predictions will be made for `x-sd(x)` and `x+sd(x)`.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uwo9592.R
\name{bwn_eff}
\alias{bwn_eff}
\title{Calculates Effects for Between-transformed Variables}
\usage{
bwn_eff(obj, vbl, idvar, data, nvals = 10, ...)
}
\arguments{
\item{obj}{A `lmerMod` or `glmerMod` object - estimated with `lmer` or `glmer` from 
the `lme4` package.}

\item{vbl}{A character string giving the name of a variable whose effects are to be calculated.}

\item{idvar}{A character string giving the name of the grouping variable.}

\item{data}{The data frame used to fit the original model.}

\item{nvals}{The number of values to use when varying a numeric variable.  The 
values used in the prediction will be a sequence of `nvals` evenly spaced from the 
minimum to the maximum of `vbl`.}

\item{...}{Currently not implemented.}
}
\description{
Calculates predicted values for between-transformed variables holding
all other variables constant at mean values.  All within-transformed variables 
will be held constant at 0, indicating an observation at the group mean. All variables must be either factors or 
numeric, the function will fail if character variables are used.  
Currently this only works with two-levels models that have a single random effect. The function returns
the predictions on the link scale.
}

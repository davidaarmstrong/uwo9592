% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uwo9592.R
\name{gg_crplot}
\alias{gg_crplot}
\title{Creates a CR Plot using ggplot2}
\usage{
gg_crplot(
  obj,
  term,
  plot = TRUE,
  res_type = c("std_deviance", "std_pearson", "deviance", "pearson", "response",
    "working", "quantile"),
  ...
)
}
\arguments{
\item{obj}{An object of class `glm`}

\item{term}{The model term to be plotted}

\item{plot}{Logical indicating whether a plot should be returned (`TRUE`) or just the data (`FALSE`)}

\item{res_type}{The type of residuals to be used.}

\item{...}{Currently not implemented}
}
\description{
Creates a CR Plot using ggplot2
}

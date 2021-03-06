% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_ci.R
\name{boot_ci_bca}
\alias{boot_ci_bca}
\title{BCa Method}
\usage{
boot_ci_bca(bt_resamples, alpha = 0.05, stat, func, ...)
}
\arguments{
\item{bt_resamples}{An `rsplit` object created by the `bootstraps` function}

\item{alpha}{level of significance}

\item{stat}{A statistic of interest}

\item{func}{A function which when applied to data returns a vector containing the statistics of interest.}

\item{...}{Optional extra arguments to pass to `func`}
}
\description{
Calculate bootstrap confidence intervals for a statistic of interest.
}

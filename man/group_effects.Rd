% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_effects.R
\name{group_effects}
\alias{group_effects}
\title{Extract baggr study effects}
\usage{
group_effects(bg, summary = FALSE, interval = 0.95)
}
\arguments{
\item{bg}{baggr object}

\item{summary}{logical; if TRUE returns summary statistics as explained below.}

\item{interval}{uncertainty interval width (numeric between 0 and 1), if summarising}
}
\value{
Either a matrix with MCMC samples (if summary = FALSE)
or a summary of these samples (if summary = TRUE).
}
\description{
Given a baggr object, returns the raw MCMC draws of the posterior for
each group's effect, or a summary of these draws. This is an internal
function currently used as a helper for plotting and printing of results.
}
\details{
If summary = TRUE, the returned object contains for each study
or group: the posterior medians, the lower and upper bounds of the
uncertainty intervals using the central posterior credible interval
of width specified in the argument "interval", the posterior mean, and
the posterior standard deviation.
}
\examples{
fit1 <- baggr(schools)
study_effects(fit1, summary = TRUE, interval = 0.5)
}

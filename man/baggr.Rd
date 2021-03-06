% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baggr.R
\name{baggr}
\alias{baggr}
\title{Bayesian aggregate treatment effects model}
\usage{
baggr(data, model = NULL, prior = NULL, pooling = "partial",
  joint_prior = TRUE, test_data = NULL, quantiles = seq(0.05, 0.95,
  0.1), outcome = "outcome", group = "group",
  treatment = "treatment", warn = TRUE, ...)
}
\arguments{
\item{data}{data frame with summary or individual level data to meta-analyse}

\item{model}{if \code{NULL}, detected automatically from input data
otherwise choose from
\code{"rubin"}, \code{"mutau"}, \code{"individual"}, \code{"quantiles"}}

\item{prior}{list of prior arguments passed directly to each model (see Details)}

\item{pooling}{choose from \code{"none"}, \code{"partial"} (default) and \code{"full"}}

\item{joint_prior}{If \code{TRUE}, \code{mu} and \code{tau} will have joint distribution.
If \code{FALSE}, they have independent priors. Ignored if no control
(\code{mu}) data exists.}

\item{test_data}{data for cross-validation; NULL for no validation, otherwise a data frame
with the same columns as \code{data} argument}

\item{quantiles}{if \code{model = "quantiles"}, a vector indicating which quantiles of data to use
(with values between 0 and 1)}

\item{outcome}{character; column name in (individual-level)
\code{data} with outcome variable values}

\item{group}{character; column name in \code{data} with grouping factor;
it's necessary for individual-level data, for summarised data
it will be used as labels for groups when displaying results}

\item{treatment}{character; column name in (individual-level) \code{data} with treatment factor;}

\item{warn}{print an additional warning if Rhat exceeds 1.05}

\item{...}{extra options passed to Stan function, e.g. \code{control = list(adapt_delta = 0.99)},
number of iterations etc.}
}
\value{
\code{baggr} class structure: a list including Stan model fit
alongside input data, pooling metrics, various model properties.
If test data is used, mean value of -2*lpd is reported as \code{mean_lpd}
}
\description{
Bayesian inference on parameters of an average treatment effects model
that's appropriate to the supplied
individual- or group-level data, using Hamiltonian Monte Carlo in Stan.
(For overall package help file see \link{baggr_package})
}
\details{
Running \code{baggr} requires 1/ data preparation, 2/ choice of model, 3/ choice of priors.
All three are discussed in depth in the package vignette (\code{vignette("baggr")}).

\strong{Data.} For aggregate data models you need a data frame with columns
\code{tau} and \code{se} or \code{tau}, \code{mu}, \code{se.tau}, \code{se.mu}.
An additional column can be used to provide labels for each group
(by default column \code{group} is used if available, but this can be
customised -- see the example below).
For individual level data three columns are needed: outcome, treatment, group. These
are identified by using the \code{outcome}, \code{treatment} and \code{group} arguments.

When working with individual-level data,
many data preparation steps (summarising, standardisation etc.)
can be done through a helper function \link{prepare_ma}.
Using it will also automatically format data inputs to be
work with \code{baggr()}.

\strong{Models.} Available models are:
\itemize{
\item for the means: \code{"rubin"} model for average treatment effect, \code{"mutau"} version which takes
into account means of control groups, `"full"`` model which reduces to "mu and tau"
(if no covariates are used)
\item "quantiles" model is also available (see Meager, 2019 in references)
}

If no model is specified, the function tries to infer the appropriate model automatically.

\strong{Priors.} It is optional to specify priors yourself,
as the package will try propose an appropriate
prior for the input data if \code{prior=NULL}.
To set the priors yourself, please refer to the list in the \code{vignette("baggr")}
}
\examples{
df_pooled <- data.frame("tau" = c(1, -1, .5, -.5, .7, -.7, 1.3, -1.3),
                        "se" = rep(1, 8),
                        "state" = datasets::state.name[1:8])
baggr(df_pooled) #baggr automatically detects the input data
# correct labels, different pooling & passing some options to Stan
baggr(df_pooled, group = "state", pooling = "full", iter = 500)

}
\author{
Witold Wiecek, Rachael Meager
}

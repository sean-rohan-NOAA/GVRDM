% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hansen_broken_stick.R
\name{hansen_broken_stick}
\alias{hansen_broken_stick}
\title{Broken stick reaction distance model}
\usage{
hansen_broken_stick(
  rr,
  Eb,
  cc,
  vv,
  ww,
  yy,
  Eb.break,
  cc.break,
  rr.max,
  sigma.rr = NULL,
  fit = TRUE
)
}
\arguments{
\item{rr}{Observed reaction distance}

\item{Eb}{Light intensity}

\item{cc}{Effective attenuation coefficient, beam attenuation coefficient, or NTU}

\item{vv}{Light function intercept}

\item{ww}{Light function slope}

\item{yy}{Turbidity function slope}

\item{Eb.break}{Saturation intensity threshold for light}

\item{cc.break}{Minimum turbidity threshold}

\item{rr.max}{Maximum reaction distance}

\item{sigma.rr}{Standard deviation of reaction distance}

\item{fit}{Use the model in fitting mode?}
}
\value{
If fit == TRUE, returns the negative log likelihood. Otherwise, returns a list containing model parameters, negative log-likelihood, and fitted reaction distances.
}
\description{
Continuous version modified from Hansen et al. (2013, 2015)
}

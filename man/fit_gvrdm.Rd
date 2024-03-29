% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_gvrdm.R
\name{fit_gvrdm}
\alias{fit_gvrdm}
\title{Function for fitting the GVRDM and Aksnes and Utne model}
\usage{
fit_gvrdm(
  rr,
  cc,
  Ap = NA,
  C0 = NA,
  Eprime = NA,
  prey.types = NA,
  prey = NA,
  prey.size = NA,
  tt = NA,
  tt1 = NA,
  tt2 = NA,
  tt3 = NA,
  tt4 = NA,
  Ke,
  Ke.trans = FALSE,
  Eb,
  angle = NA,
  kd = NA,
  kd.mult = NA,
  beta = NA,
  kk = NA,
  delta = NA,
  alpha = 1,
  sigma = NA,
  rr.log = TRUE,
  NVrd = NA,
  NVthreshold = 0,
  NVsigma = NA,
  fit.model = TRUE,
  fit.obs = TRUE,
  ccoffset = NA,
  only.estimate = FALSE,
  silent = F,
  ...
)
}
\arguments{
\item{rr}{Reaction distance in meters}

\item{cc}{Effective attenuation coefficient or beam attenuation coefficient}

\item{Ap}{Prey area. If Ap and C0 are not provided, tt will be estimated.}

\item{C0}{Prey inherent contrast. If Ap and C0 are not provided, tt will be estimated.}

\item{Eprime}{Composite saturation parameter.}

\item{prey.types}{A vector specifying unique prey name categories that are passed to 'prey.' Used if parameters are simultaneously estimated for multiple prey types.}

\item{prey}{Vector of prey names for each observation.}

\item{prey.size}{Vector of prey size}

\item{tt}{T parameter for prey type 0. Default = NA}

\item{tt1}{T parameter for prey type 1. Default = NA}

\item{tt2}{T parameter for prey type 2. Default = NA}

\item{tt3}{T parameter for prey type 3. Default = NA}

\item{tt4}{T parameter for prey type 4. Default = NA}

\item{Ke}{Half-saturation constant}

\item{Ke.trans}{Should Ke be transformed in the model?}

\item{Eb}{# Light level}

\item{angle}{Nadir viewing angle in degrees}

\item{kd}{Diffuse attenuation coefficient of downwelling irradiance}

\item{kd.mult}{Multiplier to convert beam attenuation}

\item{beta}{Dynamic scaling function intercept parameter}

\item{delta}{Dynamic scaling function shape parameter}

\item{alpha}{Naka-Rushton exponent. Default = 1}

\item{sigma}{Standard deviation of visual reaction distance}

\item{rr.log}{Are visual errors lognormal distribution? Default = TRUE}

\item{NVrd}{Non-visual reaction distance. Default = NA (only fit visual component of the model)}

\item{NVthreshold}{Light threshold for non-visual reaction. Default = NA(only fit visual component of the model)}

\item{NVsigma}{Standard deviation of non-visual reaction distance. Default = NA(only fit visual component of the model)}

\item{fit.model}{Logical.Should the model be fitted? If FALSE, produces model diagnostics and predictions based on provided parameters.}

\item{fit.obs}{Logical. Should observations be used to calculate a likelihood?}

\item{ccoffset}{Shift kappa or beam attenuation? Default = NA results in no shift.}

\item{only.estimate}{Logical. Should only the estimated reaction distance be returned?}

\item{silent}{Logical. Should diagnostic plots be produced?}

\item{hh}{Dynamic scaling function rate parameter}

\item{return}{If fit.model is TRUE, returns the negative log likelihood. If fit is FALSE, returns model diagnostics and predictions using parameters that are passed to the model.}
}
\description{
Generalized visual reaction distance model and Aksnes and Utne model.
}
\details{
The order of composite tt parameters should match the order of prey.types in the model. By default, the model is set up to estimate five prey parameters (tt, tt1, tt2, tt3, tt4). Additional parameters can be added by modifying the function with additional numbered parameters that have names beginning with tt. Regardless of how many tt parameters are added, only parameters that are assigned values and have a corresponding prey category will be estimated by the model.
}

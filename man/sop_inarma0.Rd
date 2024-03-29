% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stationary_properties.R
\name{sop_inarma0}
\alias{sop_inarma0}
\title{Helper to compute the stationary second-order properties of an INARMA(1, 1) model.}
\usage{
sop_inarma0(mu_I, sigma2_I, phi, kappa)
}
\arguments{
\item{mu_I, sigma2_I}{mean and variance of the immigration distribution.}

\item{phi, kappa}{the remaining model parameters.}
}
\value{
a named list containing the following elements:
\describe{
\item{mu}{the stationary mean.}
\item{sigma2}{the stationary variance}
\item{acf1}{the autocorrelation \eqn{\rho_X(1)} at lag 1.}
\item{mu}{the decay factor of the autocorrelation function, i.e., \eqn{\rho_X(2)/\rho_X(1)}.}
}
}
\description{
This function computes the second-order properties in a generic fashion based on the
autoregressive parameters and first two moments of the immigration distribution. It
is used mainly inside the wrapper `sop_inarma`.
}

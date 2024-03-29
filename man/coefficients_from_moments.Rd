% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_models_moments.R
\name{coefficients_from_moments}
\alias{coefficients_from_moments}
\title{Computing moment estimates from stationary moments and autocorelations}
\usage{
coefficients_from_moments(mu, sigma2, rho, xi, family)
}
\arguments{
\item{mu}{the empirical marginal mean.}

\item{sigma2}{the empirical marginal variance.}

\item{rho}{the empirical first-order autocorrelation}

\item{xi}{the decay of the autocorrelation function, i.e., \eqn{\rho_X(2)/\rho_X(1)}.}
}
\value{
a named vector of parameter estimates
}
\description{
This requires solving a cubic equation, see Bracher and Sobolova (2024).
}

#' Helper to compute the stationary second-order properties of an INARMA(1, 1) model.
#'
#' This function computes the second-order properties in a generic fashion based on the
#' autoregressive parameters and first two moments of the immigration distribution. It
#' is used mainly inside the wrapper `sop_inarma`.
#' @param mu_I,sigma2_I mean and variance of the immigration distribution.
#' @param phi,kappa the remaining model parameters.
#' @return a named list containing the following elements:
#' \describe{
#' \item{mu}{the stationary mean.}
#' \item{sigma2}{the stationary variance}
#' \item{acf1}{the autocorrelation \eqn{\rho_X(1)} at lag 1.}
#' \item{mu}{the decay factor of the autocorrelation function, i.e., \eqn{\rho_X(2)/\rho_X(1)}.}
#' }
sop_inarma0 <- function(mu_I, sigma2_I, phi, kappa){
  xi <- 1 - phi + phi*kappa
  mu <- mu_I/(1 - kappa)
  sigma2 <- (1 - kappa*(1 + xi - phi*kappa)/(1 + xi))*sigma2_I/(1 - kappa) +
    kappa*(1 + xi - phi*kappa)/(1 + xi)*mu_I/(1 - kappa)
  acf1 <- phi*kappa*(1 + (kappa*(1 - phi)*(sigma2_I - mu_I))/
                       ((1 + xi - phi*kappa)*((1 - kappa)*sigma2_I + kappa*mu_I) + phi*kappa*sigma2_I))
  decay_acf <- 1 - phi + phi*kappa
  return(list(mu = mu,
              sigma2 = sigma2,
              acf1 = acf1,
              decay_acf = decay_acf))
}

#' Compute the stationary second-order properties of an INARMA(1, 1) model.
#'
#' Compute the stationary second-order properties of an INARMA(1, 1) model based on the model parameters.
#'
#' @export
#' @param tau,psi,phi,kappa the model parameters
#' @param family family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`
#' @return a named list containing the following elements:
#' \describe{
#' \item{mu}{the stationary mean.}
#' \item{sigma2}{the stationary variance}
#' \item{acf1}{the autocorrelation \eqn{\rho_X(1)} at lag 1.}
#' \item{mu}{the decay factor of the autocorrelation function, i.e., \eqn{\rho_X(2)/\rho_X(1)}.}
#' }
sop_inarma <- function(tau, psi = NULL, phi = NULL, beta = NULL, kappa,
                       family = c("Poisson", "Hermite", "NegBin")){
  # check arguments:
  if(is.null(phi) == is.null(beta)) stop("Exactly one of phi and beta needs to be specified.")
  if(tau <= 0) stop("tau needs to be positive.")
  if(!is.null(psi)){
    if(psi <= 0) stop("psi needs to be positive.")
  }
  if(!is.null(phi)){
    if(0 > phi | 1 < phi) stop("phi needs to be from [0, 1].")
  }
  if(!is.null(beta)){
    if(0 > beta | 1 < beta) stop("beta needs to be from [0, 1].")
  }

  # internal codes use phi parameterization:
  if(!is.null(beta)){
    phi <- 1 - beta
  }

  mu_I <- tau

  if(family == "Poisson"){
    sigma2_I <- tau
    if(!is.null(psi)) warning("Parameter psi will be ignored (there is no overdispersion parameter in the Poisson case).")
  }

  if(family == "Hermite"){
    sigma2_I <- tau*(1 + psi)
  }

  if(family == "NegBin"){
    sigma2_I <- tau + psi*tau^2
  }

  sop_inarma0(mu_I = mu_I, sigma2_I = sigma2_I, phi = phi, kappa = kappa)
}

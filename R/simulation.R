#' Simulate from an INARMA(1,1) model
#'
#' Simulation from an INARMA(1,1) model. For initialization an (approximation)
#' of the stationary distribution of \eqn{E[1]}.
#'
#' @export
#' @param tau,psi,phi,beta,kappa the model parameters. Note that `phi` and `beta`
#' characterize the same underlying parameter, with `phi = 1 - beta`. Only one
#' of them should be specified, but we keep both to ensure compatibility for
#' different settings.
#' @param lgt the length of the simulated time series.
#' @param family family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`.
#' @return A named list with the following elements:
#' \describe{
#' \item{X}{The process \eqn{X_t}}
#' \item{E}{The hidden process \eqn{E_t}}
#' \item{I}{The innovation process \eqn{I_t}}
#' }
sim_inarma <- function(tau, psi = NULL, phi = NULL, beta = NULL, kappa, lgt,
                       family = c("Poisson", "Hermite", "NegBin")){

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

  # collect arguments in list:
  args <- list(tau = tau, phi = phi, kappa = kappa, lgt = lgt)
  # add psi only if needed
  if(family %in% c("Hermite", "NegBin")) args$psi <- psi
  # call relevant simulation function
  sim_fct <- switch (family,
                     "Poisson" = sim_inarma_poisson,
                     "Hermite" = sim_inarma_hermite,
                     "NegBin" = sim_inarma_negbin
  )
  do.call(sim_fct, args)
}

#' Simulate from a Poisson INARMA(1,1) model
#'
#' Simulation from a Poisson INARMA(1,1) model. For initialization an (approximation)
#' of the stationary distribution of \eqn{E[1]}. Typically used within the wrapper `sim_inarma`.
#'
#' @param tau,phi,kappa the model parameters.
#' @param lgt the length of the simulated time series.
#' @param E1 an initial value for \eqn{E_1}.
#' @param distr_E1 an initial distribution for \eqn{E_1}.
#' @return A named list with the following elements:
#' \describe{
#' \item{X}{The process \eqn{X_t}}
#' \item{E}{The hidden process \eqn{E_t}}
#' \item{I}{The innovation process \eqn{I_t}}
#' }
sim_inarma_poisson <- function(tau, phi, kappa, E1 = NULL, distr_E1 = NULL, lgt = NULL){
  xi <- 1 - phi + phi*kappa

  # handle case where tau is scalar:
  if(length(tau) == 1){
    if(is.null(lgt)) stop("If tau is scalar lgt needs to be specified.")
    tau <- rep(tau, lgt)
    if(is.null(E1) & is.null(distr_E1)) E1 <- rpois(1, kappa*tau/(1 - xi))
  }else{
    if(!is.null(lgt)){
      if(length(tau) != lgt){
        stop("Conflicting arguments tau and lgt")
      }
    }
    lgt <- length(tau)
    if(is.null(E1) & is.null(distr_E1)) stop("If tau is not scalar E1 or distr_E1 needs to be provided.")
  }

  # sample E1 if distribution is provided:
  if(!is.null(E1) & !is.null(distr_E1)) stop("Provide either E1 OR distr_E1.")
  if(!is.null(distr_E1)){
    E1 <- sample(seq_along(distr_E1) - 1, 1, prob = distr_E1)
  }

  I <- rpois(n = lgt, lambda = tau)
  E <- rep(NA, lgt); E[1] <- E1
  X <- rep(NA, lgt)
  X[1] <- I[1] + rbinom(1, E[1], phi)

  for(t in 2:lgt){
    # version that works:
    E[t] <- E[t - 1] + I[t - 1] - rbinom(1, X[t - 1], 1 - kappa)
    X[t] <- I[t] + rbinom(1, E[t], phi)
  }
  return(list(X = tail(X, lgt), E = tail(E, lgt), I = tail(I, lgt)))
}

#' Simulate from a Hermite INARMA(1,1) model
#'
#' Simulation from a Hermite INARMA(1,1) model. For initialization an (approximation)
#' of the stationary distribution of \eqn{E[1]}. Typically used within the wrapper `sim_inarma`.
#'
#' @param tau,psi,phi,kappa the model parameters.
#' @param a1,a2 alternative parameterization of the Hermite, to be used instead of `tau, psi`.
#' @param lgt the length of the simulated time series.
#' @param E1 an initial value for \eqn{E_1}.
#' @param distr_E1 an initial distribution for \eqn{E_1}.
#' @return A named list with the following elements:
#' \describe{
#' \item{X}{The process \eqn{X_t}}
#' \item{E}{The hidden process \eqn{E_t}}
#' \item{I}{The innovation process \eqn{I_t}}
#' }
sim_inarma_hermite <- function(tau, psi, phi, kappa, a1 = NULL, a2 = NULL,
                               E1 = NULL, distr_E1 = NULL, lgt = NULL){
  # move to parameterization a1, a2 (and bring to same length if necessary):
  pars <- to_a1a2_herm(a1 = a1, a2 = a2, mu = tau, psi = psi)
  a1 <- pars$a1
  a2 <- pars$a2
  xi <- 1 - phi + phi*kappa

  # handle case where a1, a2 are scalar:
  if(length(a1) == 1){
    if(is.null(lgt)) stop("If a1 and a2 or tau and psi are scalar lgt needs to be specified.")
    a1 <- rep(a1, lgt)
    a2 <- rep(a2, lgt)
    if(is.null(E1) & is.null(distr_E1)){
      # sample E1 from stationary distribution if no starting value supplied:
      E1 <-  rherm(1,
                   a1 = (kappa*(1 + xi)*a1 + 2*kappa*(1 + xi - kappa)*a2)/
                     (1 - xi^2),
                   a2 = kappa^2*a2/(1 - xi^2))
    }
  }else{
    if(!is.null(lgt)){
      if(length(tau) != lgt){
        stop("Conflicting arguments tau and lgt")
      }
    }
    lgt <- length(tau)
    if(is.null(E1) & is.null(distr_E1)) stop("If a1 and a2 or tau and psi are not scalar E1 or distr_E1 needs to be provided.")
  }

  # sample E1 if distribution is provided:
  if(!is.null(E1) & !is.null(distr_E1)) stop("Provide either E1 OR distr_E1.")
  if(!is.null(distr_E1)){
    E1 <- sample(seq_along(distr_E1) - 1, 1, prob = distr_E1)
  }

  I <- rherm(n = lgt, a1 = a1, a2 = a2)
  E <- rep(NA, lgt); E[1] <- E1
  X <- rep(NA, lgt)
  X[1] <- I[1] + rbinom(1, E[1], phi)

  for(t in 2:lgt){
    # version that works:
    E[t] <- E[t - 1] + I[t - 1] - rbinom(1, X[t - 1], 1 - kappa)
    X[t] <- I[t] + rbinom(1, E[t], phi)
  }
  return(list(X = tail(X, lgt), E = tail(E, lgt), I = tail(I, lgt)))
}


#' Simulate from a negative binomial INARMA(1,1) model
#'
#' Simulation from a negative binomial INARMA(1,1) model. For initialization an (approximation)
#' of the stationary distribution of \eqn{E[1]}. Typically used within the wrapper `sim_inarma`.
#'
#' @param tau,psi,phi,kappa the model parameters.
#' @param lgt the length of the simulated time series.
#' @param E1 an initial value for \eqn{E_1}.
#' @param distr_E1 an initial distribution for \eqn{E_1}.
#' @return A named list with the following elements:
#' \describe{
#' \item{X}{The process \eqn{X_t}}
#' \item{E}{The hidden process \eqn{E_t}}
#' \item{I}{The innovation process \eqn{I_t}}
#' }
sim_inarma_negbin <- function(tau, psi, phi, kappa,
                               E1 = NULL, distr_E1 = NULL, lgt = NULL){

  xi <- 1 - phi + phi*kappa

  # handle case where tau, psi are scalar:
  if(length(tau) == 1 & length(psi)){
    if(is.null(lgt)) stop("If tau and psi are scalar lgt needs to be specified.")
    tau <- rep(tau, lgt)
    psi <- rep(psi, lgt)
    if(is.null(E1) & is.null(distr_E1)){
      # sample E1 from approximated stationary distribution if no starting value supplied:
      sigma2_I <- tau + psi*tau^2
      mu_E <- kappa*tau/(1 - xi)
      sigma2_E <- (kappa^2*sigma2_I + kappa*(1 - kappa + xi)*tau)/(1 - xi^2)
      size_E <- mu_E^2/(sigma2_E - mu_E)
      E1 <-  rnbinom(1, size = size_E, mu = mu_E)
    }
  }else{
    if(!is.null(lgt)){
      if(length(tau) != lgt){
        stop("Conflicting arguments tau and lgt")
      }
    }
    lgt <- length(tau)
    if(is.null(E1) & is.null(distr_E1)) stop("If tau and psi are not scalar E1 needs to be provided.")
  }

  # sample E1 if distribution is provided:
  if(!is.null(E1) & !is.null(distr_E1)) stop("Provide either E1 OR distr_E1.")
  if(!is.null(distr_E1)){
    E1 <- sample(seq_along(distr_E1) - 1, 1, prob = distr_E1)
  }

  I <- rnbinom(n = lgt, size = 1/psi, mu = tau)
  E <- rep(NA, lgt); E[1] <- E1
  X <- rep(NA, lgt)
  X[1] <- I[1] + rbinom(1, E[1], phi)

  for(t in 2:lgt){
    # version that works:
    E[t] <- E[t - 1] + I[t - 1] - rbinom(1, X[t - 1], 1 - kappa)
    X[t] <- I[t] + rbinom(1, E[t], phi)
  }
  return(list(X = tail(X, lgt), E = tail(E, lgt), I = tail(I, lgt)))
}



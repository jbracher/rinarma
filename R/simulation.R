#' Simulate from an INARMA(1,1) model
#'
#' Simulation from an INARMA(1,1) model. For initialization an (approximation)
#' of the stationary distribution of \eqn{E[1]}.
#'
#' @export
#' @param tau,psi,phi,beta,kappa,zeta the model parameters. Note that `phi` and `beta`
#' characterize the same underlying parameter, with `phi = 1 - beta`. Only one
#' of them should be specified, but we keep both to ensure compatibility for
#' different settings.
#' @param lgt the length of the simulated time series.
#' @param offspring type of the thinning \eqn{\kappa \bullet X_t}; one of `"binomial"`, or `"binomial-Poisson"`.
#' @param family distribution of the imports/innovations; one of `"Poisson"`, `"Hermite"` or `"NegBin"`.
#' @return A named list with the following elements:
#' \describe{
#' \item{X}{The process \eqn{X_t}}
#' \item{E}{The hidden process \eqn{E_t}}
#' \item{I}{The innovation process \eqn{I_t}}
#' }
sim_inarma <- function(tau, psi = NULL, phi = NULL, beta = NULL, kappa,
                       zeta = NULL, lgt,
                       offspring = c("binomial", "binomial-Poisson"),
                       family = c("Poisson", "Hermite", "NegBin"),
                       E1 = NULL){

  if(is.null(phi) == is.null(beta)) stop("Exactly one of phi and beta needs to be specified.")
  if(tau <= 0) stop("tau needs to be positive.")
  if(!is.null(psi)){
    if (family == "Poisson") {
      warning("Poisson INARMA does is equidispersed, the overdispersion parameter `psi` will be ignored")
      } else if(psi <= 0) stop("psi needs to be positive")
  }
  if(!is.null(phi)){
    if(any(0 > phi) | 1 < sum(phi)) stop("phi needs to be from [0, 1].")
  }
  if(!is.null(beta)){
    if(any(0 > beta) | 1 < sum(beta)) stop("beta needs to be from [0, 1].")
  }
  if(!is.null(kappa)){
    if(any(0 > kappa) | 1 < sum(kappa)) stop("kappa needs to be from [0, 1],
                                             otherwise the process is non-stationary.")
  }
  if (is.null(zeta) && offspring == "binomial-Poisson") {
    stop("For the binomial-Poisson thinning the parameter 'zeta' must
         be supplied.")
  }
  if (offspring == "binomial-Poisson" && is.null(E1)) {
    stop("For the binomial-Poisson thinning the initial value E1 must be supplied.")
  }

  # internal codes use phi parameterization:
  if(!is.null(beta)){
    phi <- 1 - beta
  }

  # collect arguments in list:
  args <- list(tau = tau, phi = phi, kappa = kappa, lgt = lgt, E1 = E1)
  # add psi only if needed
  if(family %in% c("Hermite", "NegBin")) args$psi <- psi
  if (offspring == "binomial") {
    args$zeta <- 0
  } else if (offspring == "binomial-Poisson") {
    args$zeta <- zeta
  } else {
    stop("Invalid offspring distribution.")
  }

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
#' @param tau,phi,kappa,zeta the model parameters.
#' @param lgt the length of the simulated time series.
#' @param E1 an initial value for \eqn{E_1}.
#' @param distr_E1 an initial distribution for \eqn{E_1}.
#' @return A named list with the following elements:
#' \describe{
#' \item{X}{The process \eqn{X_t}}
#' \item{E}{The hidden process \eqn{E_t}}
#' \item{I}{The innovation process \eqn{I_t}}
#' }
sim_inarma_poisson <- function(tau, phi, kappa, zeta = NULL, E1 = NULL,
                               distr_E1 = NULL, lgt = NULL){
  xi <- 1 - phi + (1 - sum(1 - phi)) * kappa

  # handle case where tau is scalar:
  if(length(tau) == 1){
    if(is.null(lgt)) stop("If tau is scalar lgt needs to be specified.")
    tau <- rep(tau, lgt)
    if(is.null(E1) & is.null(distr_E1) & length(phi) == 1) {
      E1 <- rpois(1, kappa*tau/(1 - xi))
    } else if (is.null(E1) & is.null(distr_E1) & length(phi) > 1) {
      stop("For simulating from a model with q > 1, E1 must be a vector of starting values of length q.")
    }
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
  if (max(length(phi), length(kappa)) > 1) {
    res <- sim_loop_higher_order(phi = phi, kappa = kappa, zeta = zeta,
                                 I = I, E1 = E1, lgt = lgt)
  } else if (length(phi) == 1 & length(kappa) == 1) {
    res <- sim_loop(phi = phi, kappa = kappa, zeta = zeta,
                    I = I, E1 = E1, lgt = lgt)
  } else {
    stop("One of the model parameters missing.")
  }

  return(list(X = tail(res$X, lgt), E = tail(res$E, lgt), I = tail(I, lgt)))
}

#' Simulate from a Hermite INARMA(1,1) model
#'
#' Simulation from a Hermite INARMA(1,1) model. For initialization an (approximation)
#' of the stationary distribution of \eqn{E[1]}. Typically used within the wrapper `sim_inarma`.
#'
#' @param tau,psi,phi,kappa,zeta the model parameters.
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
sim_inarma_hermite <- function(tau, psi, phi, kappa, zeta, a1 = NULL,
                               a2 = NULL, E1 = NULL, distr_E1 = NULL, lgt = NULL){
  # move to parameterization a1, a2 (and bring to same length if necessary):
  pars <- to_a1a2_herm(a1 = a1, a2 = a2, mu = tau, psi = psi)
  a1 <- pars$a1
  a2 <- pars$a2
  xi <- 1 - phi + (1 - sum(1 - phi)) * kappa

  # handle case where a1, a2 are scalar:
  if(length(a1) == 1){
    if(is.null(lgt)) stop("If a1 and a2 or tau and psi are scalar lgt needs to be specified.")
    a1 <- rep(a1, lgt)
    a2 <- rep(a2, lgt)
    if(is.null(E1) & is.null(distr_E1) & length(phi) == 1){
      # sample E1 from stationary distribution if no starting value supplied:
      E1 <-  rherm(1,
                   a1 = (kappa*(1 + xi)*a1 + 2*kappa*(1 + xi - kappa)*a2)/
                     (1 - xi^2),
                   a2 = kappa^2*a2/(1 - xi^2))
    } else if (is.null(E1) & is.null(distr_E1) & length(phi) > 1) {
      stop("For simulating from a model with q > 1, E1 must be a vector of starting values of length q.")
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
  if (max(length(phi), length(kappa)) > 1) {
    res <- sim_loop_higher_order(phi = phi, kappa = kappa, zeta = zeta,
                                 I = I, E1 = E1, lgt = lgt)
  } else if (length(phi) == 1 & length(kappa) == 1) {
    res <- sim_loop(phi = phi, kappa = kappa, zeta = zeta,
                    I = I, E1 = E1, lgt = lgt)
  } else {
    stop("One of the model parameters missing.")
  }
  return(list(X = tail(res$X, lgt), E = tail(res$E, lgt), I = tail(I, lgt)))
}


#' Simulate from a negative binomial INARMA(1,1) model
#'
#' Simulation from a negative binomial INARMA(1,1) model. For initialization an (approximation)
#' of the stationary distribution of \eqn{E[1]}. Typically used within the wrapper `sim_inarma`.
#'
#' @param tau,psi,phi,kappa,zeta the model parameters.
#' @param lgt the length of the simulated time series.
#' @param E1 an initial value for \eqn{E_1}.
#' @param distr_E1 an initial distribution for \eqn{E_1}.
#' @return A named list with the following elements:
#' \describe{
#' \item{X}{The process \eqn{X_t}}
#' \item{E}{The hidden process \eqn{E_t}}
#' \item{I}{The innovation process \eqn{I_t}}
#' }
sim_inarma_negbin <- function(tau, psi, phi, kappa, zeta,
                               E1 = NULL, distr_E1 = NULL, lgt = NULL){

  xi <- 1 - phi + (1 - sum(1 - phi)) * kappa

  # handle case where tau, psi are scalar:
  if(length(tau) == 1 & length(psi)){
    if(is.null(lgt)) stop("If tau and psi are scalar lgt needs to be specified.")
    tau <- rep(tau, lgt)
    psi <- rep(psi, lgt)
    if(is.null(E1) & is.null(distr_E1) & length(phi) == 1){
      # sample E1 from approximated stationary distribution if no starting value supplied:
      sigma2_I <- tau + psi*tau^2
      mu_E <- kappa*tau/(1 - xi)
      sigma2_E <- (kappa^2*sigma2_I + kappa*(1 - kappa + xi)*tau)/(1 - xi^2)
      size_E <- mu_E^2/(sigma2_E - mu_E)
      E1 <-  rnbinom(1, size = size_E, mu = mu_E)
    } else if (is.null(E1) & is.null(distr_E1) & length(phi) > 1) {
      stop("For simulating from a model with q > 1, E1 must be a vector of starting values of length q.")
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
  if (max(length(phi), length(kappa)) > 1) {
    res <- sim_loop_higher_order(phi = phi, kappa = kappa, zeta = zeta,
                                 I = I, E1 = E1, lgt = lgt)
  } else if (length(phi) == 1 & length(kappa) == 1) {
    res <- sim_loop(phi = phi, kappa = kappa, zeta = zeta,
                    I = I, E1 = E1, lgt = lgt)
  } else {
    stop("One of the model parameters missing.")
  }
  return(list(X = tail(res$X, lgt), E = tail(res$E, lgt), I = tail(I, lgt)))
}

#' Execute the loop of the INARMA(1,1) model simulation
#'
#' @param phi,kappa,zeta the model parameters.
#' @param I the sequence of imports/innovations with length `lgt`
#' @param E1 an initial value for \eqn{E_1}.
#' @param lgt the length of the simulated time series.
#' @return A named list with the following elements:
#' \describe{
#' \item{X}{The process \eqn{X_t}}
#' \item{E}{The hidden process \eqn{E_t}}
#' }
sim_loop <- function (phi, kappa, zeta, I, E1, lgt) {

  E <- rep(NA, lgt); E[1] <- E1
  L <- rbinom(1, E[1], 1 - phi)  # Newly added for the Pois-Binom thinning
  X <- rep(NA, lgt)
  X[1] <- I[1] + E[1] - L

  if (zeta == 0) {
    # version that was used for generating the data from the simulation study:
    for(t in 2:lgt){
        E[t] <- E[t - 1] + I[t - 1] - rbinom(1, X[t - 1], 1 - kappa)
        X[t] <- I[t] + rbinom(1, E[t], phi)
    }
  } else {
    # version that works for the Poisson-Binomial thinning too:
    for(t in 2:lgt){
        E[t] <- L + rbinom(1, X[t - 1], kappa * (1 - zeta)) +
          rpois(1, X[t - 1] * kappa * zeta)
        L <- rbinom(1, E[t], 1 - phi)
        X[t] <- I[t] + (E[t] - L)
      }
  }

  return(list(X = X, E = E))
}

#' Execute the loop of the INARMA(p,q) model simulation, for
#' at least one p, q > 1
#'
#' @param phi,kappa,zeta the model parameters.
#' @param I the sequence of imports/innovations with length `lgt`
#' @param E1 an initial values for \eqn{E_1}, ...\eqn{E_q} .
#' @param lgt the length of the simulated time series.
#' @return A named list with the following elements:
#' \describe{
#' \item{X}{The process \eqn{X_t}}
#' \item{E}{The hidden process \eqn{E_t}}
#' }
sim_loop_higher_order <- function (phi, kappa, zeta, I, E1, lgt) {

  # Grab the order of the model
  p <- length(kappa)
  q <- length(phi)
  max_order <- max(c(p, q))

  # Pad the parameter vectors so that they are of the same length
  beta_vec <- c(1 - phi, rep(0, max_order - q), 1 - sum(1 - phi))
  kappa_vec <- c(kappa, rep(0, max_order - p), 1 - sum(kappa))

  # Allocate the containers
  E <- X <- rep(0, lgt + max_order)
  L <- integer(q + 1)
  C <- integer(p + 1)

  # Initialize
  E[1:q] <- E1

  # Padding I in order not to get NA values at the end of the trajectory (which
  # is cropped anyway)

  I <- c(I, rep(0, max_order))

  for (k in 1:(lgt + max_order)) {

    # Do the MA thinning and update the processes
    L <- rmultinom(1, E[k], beta_vec)[, 1]
    E[k + 1:q] <- E[k + 1:q] + L[1:q]
    X[k] <- L[max_order + 1] + I[k]

    # Do the AR thinning and update the processes
    C <- rmultinom(1, X[k], kappa_vec)[, 1]
    E[k + 1:p] <- E[k + 1:p] + C[1:p]
  }

  ret_list <- list(E = E[1:lgt], X = X[1:lgt])

  return(ret_list)
}

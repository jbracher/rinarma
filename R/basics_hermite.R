#' Re-parameterize the Hermite distribution to a1, a2
#'
#' `a1, a2` are the parameters of the Hermite distribution in a parameterization
#' where \eqn{X \sim \text{Herm}(a_1, a_2)} is equivalent to
#' \deqn{X = Y_1 + 2Y_2, Y_1 \sim \text{Pois}(a_1), Y_2 \sim \text{Pois}(a_2).}
#'
#' @param a1,a2 parameters in format a1, a2 (in this case only some recycling is done)
#' @param mu,psi parameters in format mu, psi; see Bracher and Sobolova (2024) for details.
to_a1a2_herm <- function(a1 = NULL, a2 = NULL, mu = NULL, psi = NULL){
  # check that an admissible combination of parameters has been provided
  # and re-parameterize from mu, psi to a1, a2 if necessary:
  if(!is.null(a1) | !is.null(a2)){
    if(!is.null(mu) | !is.null(psi) | is.null(a1) | is.null(a2)){
      stop("Specify either a1 and a2 or mu and psi.")
    }
    if((length(a1) != length(a2)) & length(a1) > 1 & length(a2) > 1){
      warning("a1 and a2 are of different length, recycling.")
      # recycle if necessary:
      a1 <- rep(a1, length.out = max(length(a1), length(a2)))
      a2 <- rep(a2, length.out = max(length(a1), length(a2)))
    }
  }

  if(!is.null(mu) | !is.null(psi)){
    if(!is.null(a1) | !is.null(a2) | is.null(mu) | is.null(psi)){
      stop("Specify either a1 and a2 or mu and psi.")
    }
    if((length(mu) != length(psi)) & length(mu) > 1 & length(psi) > 1){
      warning("mu and psi are of different length, recycling.")
    }
    a1 <- mu*(1 - psi)
    a2 <- psi*mu/2
  }
  return(list(a1 = a1, a2 = a2))
}

#' Probability mass function of the Hermite distribution
#'
#' @param x the values at which to evaluate the probability mass function
#' @param a1,a2 parameters in format a1, a2 (in this case only some recycling is done)
#' @param mu,psi parameters in format mu, psi
dherm <- function(x, a1 = NULL, a2 = NULL, mu = NULL, psi = NULL, log = FALSE){

  if((length(a1) != length(a2)) & length(a1) > 1 & length(a2) > 1){
    stop("a1 and a2 are both of length larger than one, but not the same length.")
  }

  # move to parameterization a1, a2:
  pars <- to_a1a2_herm(a1, a2, mu, psi)
  a1 <- pars$a1
  a2 <- pars$a2

  # recycle
  if(length(a1) == 1){
    a1 <- rep(a1, length(x))
  }
  if(length(a2) == 1){
    a2 <- rep(a2, length(x))
  }

  # function to evaluate density for i-th value:
  dherm0 <- function(i){
    support_sum <- 0:floor(x[i]/2)
    exp(-a1[i] - a2[i]) * sum(a1[i]^(x[i] - 2*support_sum)*a2[i]^(support_sum)/
                          factorial(x[i] - 2*support_sum)/factorial(support_sum))
  }

  # apply to vector:
  lik <- sapply(seq_along(x), dherm0)

  # log-transform if needed
  if(log == TRUE) lik <- log(lik)
  return(lik)
}


#' Sampling from the Hermite distribution
#'
#' @param n the number of samples to generate
#' @param a1,a2 parameters in format a1, a2 (in this case only some recycling is done)
#' @param mu,psi parameters in format mu, psi
rherm <- function(n, a1 = NULL, a2 = NULL, mu = NULL, psi = NULL){
  # move to parameterization a1, a2:
  pars <- to_a1a2_herm(a1 = a1, a2 = a2, mu = mu, psi = psi)
  a1 <- pars$a1
  a2 <- pars$a2
  # sample:
  rpois(n, a1) + 2*rpois(n, a2)
}

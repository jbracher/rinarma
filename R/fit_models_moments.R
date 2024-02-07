#' Computing moment estimates from stationary moments and autocorelations
#'
#' This requires solving a cubic equation, see Bracher and Sobolova (2024).
#'
#' @param mu the empirical marginal mean.
#' @param sigma2 the empirical marginal variance.
#' @param rho the empirical first-order autocorrelation
#' @param xi the decay of the autocorrelation function, i.e., \eqn{\rho_X(2)/\rho_X(1)}.
#' @return a named vector of parameter estimates
#'
coefficients_from_moments <- function(mu, sigma2, rho, xi, family){
  # Poisson case: simple closed forms:
  if(family == "Poisson"){
    beta <- xi - rho
    phi <- 1 - beta
    tau <- mu*(rho - rho*xi)/(rho^2 + rho - rho*xi)
    kappa <- rho^2/(rho^2 + rho - rho*xi)

    coefficients <- c(tau = tau, phi = phi, kappa = kappa)
  }

  # overdispersed:

  if(family != "Poisson"){
    # parameters of the cubic equation:
    c3 <- (1 - xi)*(mu + sigma2) + 2*rho*sigma2 # positive
    c2 <- (1 - xi)*(-(2 + xi)*sigma2 - xi*mu) - 2*rho*sigma2*(2 + xi) # negative
    c1 <- (1 - xi)*(1 + xi)*(sigma2)  + rho*sigma2*3*(1 + xi) # positive
    c0 <- -rho*sigma2*(1 + xi) # negative

    # solve cubic equation:
    root <- polyroot(c(c0, c1, c2, c3))

    # choose correct root as kappa:
    kappa <- Re(root[Re(root) >= 0 & Re(root) <= 1 & abs(Im(root)) <= 0.0001])

    # compute other parameters:
    tau <- mu*(1 - kappa)
    beta <- (xi - kappa)/(1 - kappa)
    phi <- 1 - beta

    q <- (1 + beta)*kappa/(1 + xi)
    sigma2_tau <- (1 - kappa)*(sigma2 - q*mu)/(1 - q)
    rm(q)

    if(family == "Hermite"){
      psi <- (sigma2_tau/tau) - 1
    }

    if(family == "NegBin"){
      psi <- (sigma2_tau - tau) / tau^2
    }

    coefficients <- c(tau = tau, phi = phi, kappa = kappa, psi = psi)
  }
  return(coefficients)
}

#' Fitting an INARMA(1,1) model via method of moments.
#'
#' See `?fit_inarma` for details on the model definition.
#'
#' The function implements some heuristics to ensure model parameters fall into the allowed ranges.
#' See Bracher and Sobolova (2024) for details.
#'
#' @examples
#' data("measles")
#' X <- measles$value
#' # Note: running the fit takes a little while.
#' fit <- fit_inarma_moments(X, family = "Poisson")
#' summary(fit)
#'
#'
#' @export
#' @param observed a vector of observed count values
#' @param family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`
#' @param parameterization the function internally works with a slightly different notation as

#' @return A named list with the following elements.
#' \describe{
#' \item{coefficients}{the estimated coefficients.}
#' \item{coefficients_uncorrected}{the estimated model coefficients without any corrections
#' (i.e., coefficients can be negative).}
#' \item{nobs}{the number of observations.}
#' \item{fitting_method}{the method used to fit the model, here `"moments"`.}
#' }
#'
fit_inarma_moments <- function(observed, family = c("Poisson", "Hermite", "NegBin"), parameterization = "beta"){
  # get the four relevant moments:
  mu <- mean(observed)
  sigma2 <- var(observed)
  acor <- acf(observed, type = "correlation", plot = FALSE)$acf[-1]
  rho <- acor[1]
  xi <- acor[2]/acor[1]

  # compute model coefficients from moments:
  coefficients <- coefficients_from_moments(mu, sigma2, rho, xi, family = family)

  # function to check all coefficients are in their allowed ranges:
  check_coefficients <- function(coefficients){
    ret <- TRUE
    coefficients <- as.list(coefficients)
    if(coefficients$tau < 0){ret <- FALSE;  warning("Implausible estimation results: tau is negative")}
    if(coefficients$phi < 0 | coefficients$phi > 0.95){ret <- FALSE;   warning("phi is not in [0, 0.95) - risk of implausible results and numerical instability.")}
    if(coefficients$kappa < 0 | coefficients$kappa > 0.95){ret <- FALSE;   warning("kappa is not in [0, 0.95) - risk of implausible results and numerical instability.")}
    if(family == "Hermite"){
      if(coefficients$psi < 0 | coefficients$psi > 0.95){ret <- FALSE;   warning("psi is not in [0, 0.95) - risk of implausible results and numerical instability.")}
    }
    if(family == "NegBin"){
      if(coefficients$psi < 0){ret <- FALSE;  warning("Implausible estimation results: psi is negative.")}
    }
    return(ret)
  }

  # check if all coefficients are in their allowed ranges:
  coefficients_admissible <- check_coefficients(coefficients)

  # if not: repeat estimation using moments shifted to be compatible with estimation equations:
  coefficients_uncorrected <- NULL
  if((sigma2 < mu & family != "Poisson") | xi < rho | xi > 0.95){ # !coefficients_admissible
    coefficients_uncorrected <- coefficients # store original estimates

    # shift moments to be compatible:
    message("Repeating estimation with adapted moments...")
    sigma2 <- max(mu, sigma2)
    xi <- min(max(rho, xi), 0.95)

    # re-run estimation
    coefficients <- coefficients_from_moments(mu, sigma2, rho, xi, family = family)

    # re-check:
    coefficients_admissible <- check_coefficients(coefficients)
    if(!coefficients_admissible){
      warning("Estimation results still implausible, but returning as is.")
    }
  }

  ret <- list()
  ret$coefficients <- coefficients
  ret$coefficients_uncorrected <- coefficients_uncorrected
  ret$observed <- observed
  ret$nobs <- length(observed)
  ret$fitting_method <- "moments"

  # adapt nomenclature to paper Bracher and Sobolova (2024)
  if(parameterization == "beta"){
    names(ret$coefficients)[names(ret$coefficients) == "phi"] <- "beta"
    ret$coefficients["beta"] <- 1 - ret$coefficients["beta"]
    names(ret$coefficients)[names(ret$coefficients) == "phi"] <- "beta"
  }

  class(ret) <- "inarma"

  return(ret)
}

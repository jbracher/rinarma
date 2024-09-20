# Fitting an INGARCH(1, 1) model
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_ingarch <- function(observed, family = c("Poisson", "Hermite", "NegBin"),
                        parametrization = c("GLM-like", "thinning-based"),
                        start = NULL, return_se = TRUE, control_optim = NULL){

  # starting values
  if(is.null(start)) {
    start_val <- c(log_nu = 2, log_alpha = -1, log_beta = -1, log_lambda1 = 0.5)
  } else {
    if (parametrization == "thinning based") {
      kappa_start <- exp(start["logit_kappa"]) / exp(start["logit_kappa"] + 1)
      beta_start <- exp(start["logit_beta"]) / exp(start["logit_beta"] + 1)
      tau_start <- exp(start["log_tau"])
      start_val <- c(
        log_nu = log(tau_start * (1 - beta_start)),
        log_alpha = log(kappa_start * (1 - beta_start)),
        log_beta = log(beta_start),
        log_lambda1 = 0.5
      )
    }
  }

  # negative log-likelihood as function of parameter vector,
  # depends on the family:
  if(family == "Poisson"){
    nllik <- function (pars) {
      -llik_ingarch(observed, nu = exp(pars["log_nu"]),
                    alpha = exp(pars["log_alpha"]),
                    beta = exp(pars["log_beta"]),
                    lambda1 = exp(pars["log_lambda1"]))
    }
  }

  if(family == "Hermite"){
    if (is.null(start)) {
      start_val <- c(start_val, logit_psi = -0.5)
    } else {
      if (parametrization == "thinning based") {
        start_val <- c(start_val, start["logit_psi"])
      }
    }
    nllik <- function (pars) {
      -llik_hingarch(observed, nu = exp(pars["log_nu"]),
                     alpha = exp(pars["log_alpha"]),
                     beta = exp(pars["log_beta"]),
                     psi = exp(pars["logit_psi"]),
                     lambda1 = exp(pars["log_lambda1"]))
    }
  }

  if(family == "NegBin"){
    if (is.null(start)) {
      start_val <- c(start_val, log_psi = -0.5)
    } else {
      if (parametrization == "thinning based") {
        start_val <- c(start_val, start["log_psi"])
      }
    }
    nllik <- function (pars) {
      -llik_nbingarch2(observed, nu = exp(pars["log_nu"]),
                       alpha = exp(pars["log_alpha"]),
                       beta = exp(pars["log_beta"]),
                       psi = exp(pars["log_psi"]),
                       lambda1 = exp(pars["log_lambda1"]))
    }
  }

  # Optimize
  opt <- optim(start_val, nllik, control = control_optim)

  # structure results
  llik <- -opt$value
  coefficients_raw <- opt$par
  coefficients <- c(nu = exp(opt$par["log_nu"]),
                    alpha = exp(opt$par["log_alpha"]),
                    beta = exp(opt$par["log_beta"]),
                    lambda1 = exp(opt$par["log_lambda1"]))
  names(coefficients) <- c("nu", "alpha", "beta", "lambda1")

  # compute fitted values, depends on the family:
  if (family =="Poisson") {
    fitted_values <- do.call(
      llik_ingarch,
      c(as.list(coefficients), vect = list(observed), return_fitted = TRUE)
    )$fitted
  } else if (family =="Hermite") {
    coefficients <- c(coefficients, psi = unname(exp(opt$par["logit_psi"])))
    fitted_values <- do.call(
      llik_hingarch,
      c(as.list(coefficients), vect = list(observed), return_fitted = TRUE)
    )$fitted
  } else if (family =="NegBin") {
    coefficients <- c(coefficients, psi = unname(exp(opt$par["log_psi"])))
    fitted_values <- do.call(
      llik_nbingarch2,
      c(as.list(coefficients), vect = list(observed), return_fitted = TRUE)
    )$fitted
  }

  # compute AIC
  AIC <- 2*(-llik + length(coefficients))

  # Transform the coefficients if required for the thinning-based display
  # coefficients_raw REMAIN IN THE INGARCH PARAMETRIZATION
  if (parametrization == "thinning-based") {
    coefficients <- c(
      tau = coefficients["nu"] / (1 - coefficients["beta"]),
      kappa = coefficients["alpha"] / (1 - coefficients["beta"]),
      beta = coefficients["beta"],
      theta = get_cluster_size(coefficients["psi"], family = family)
    )
    names(coefficients) <- c("tau", "kappa", "beta", "theta")
  }

  # compute residuals
  pearson_residuals <- (observed - fitted_values)/sqrt(fitted_values)

  return(
    list(
      coefficients = coefficients,
      coefficients_raw = coefficients_raw,
      loglikelihood = llik, dim = length(coefficients), AIC = AIC, opt = opt,
      fitted_values = fitted_values,
      pearson_residuals = pearson_residuals,
      family = family,
      observed = observed
    )
  )
}

# Evaluate log-likelihood for an INGARCH(1, 1) model
# Arguments:
# vect: a vector containing the observed time series
# tau, phi, kappa, S1: parameters of the INGARCH(1) model
# return fitted: should fitted values be returned?
# Return: log-likelihood, or (if return_fitted) a list containing the log-likelihood
# and the fitted values
llik_ingarch <- function(vect, nu, alpha, beta, lambda1, return_fitted = FALSE){
  lgt <- length(vect)
  lambda <- numeric(lgt)
  lambda[1] <- lambda1
  for(i in 2:lgt){
    lambda[i] <- nu + alpha*vect[i - 1] + beta*lambda[i - 1]
  }
  llik <- sum(dpois(vect, lambda, log = TRUE))

  if(return_fitted){
    return(list(value = llik, fitted = lambda))
  }else{
    return(llik)
  }
}

# Evaluate log-likelihood for a negative binomial INGARCH(1, 1) model, version 2: size proportional to mean
#
# Arguments:
# vect: a vector containing the observed time series
# tau, phi, kappa, psi, S1: parameters of the INGARCH(1) model
# Return: log-likelihood, or (if return_fitted) a list containing the log-likelihood
# and the fitted values
llik_nbingarch2 <- function(vect, nu, alpha, beta, psi, lambda1,
                            return_fitted = FALSE){
  lgt <- length(vect)
  lambda <- numeric(lgt)
  lambda[1] <- lambda1
  for(i in 2:lgt){
    lambda[i] <- nu + alpha*vect[i - 1] + beta*lambda[i - 1]
  }
  llik <- sum(dnbinom(vect, mu = lambda, size = lambda/psi, log = TRUE))

  if(return_fitted){
    return(list(value = llik, fitted = lambda))
  }else{
    return(llik)
  }
}

# Evaluate log-likelihood for a Hermite INGARCH(1, 1) model
#
# Arguments:
# vect: a vector containing the observed time series
# tau, phi, kappa, psi, S1: parameters of the INGARCH(1) model
# Return: log-likelihood, or (if return_fitted) a list containing the log-likelihood
# and the fitted values
llik_hingarch <- function(vect, nu, alpha, beta, psi, lambda1, return_fitted = FALSE){
  lgt <- length(vect)
  lambda <- numeric(lgt)
  lambda[1] <- lambda1
  for(i in 2:lgt){
    lambda[i] <- nu + alpha*vect[i - 1] + beta*lambda[i - 1]
  }
  llik <- sum(dherm(vect, mu = lambda, psi = psi, log = TRUE))
  if(return_fitted){
    return(list(value = llik, fitted = lambda))
  }else{
    return(llik)
  }
}

# compute the clustersize for the summary table from the model parameter psi
# and the distribuiton family
get_cluster_size <- function(psi, family = c("Poisson", "Hermite", "NegBin")){
  # computing cluster size is tedious:
  if(family == "Poisson"){
    clustersize <- 1
  }
  if(family == "Hermite"){
    pi <- psi / (2 - psi)
    clustersize <- 1 + pi
  }
  if(family == "NegBin"){
    pi <- 1/(1 + psi)
    clustersize <- (1 - pi)/(-pi*log(pi))
  }
  return(clustersize)
}
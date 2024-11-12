library(tidyverse)
library(pracma)
library(tscount)

###############################################################################
## Functions for simulating an INARMA(p, q) process and calculating the moments
###############################################################################

## Function for simulation of the Poisson INARMA(p, q) process
##
## function arguments:
## tau, beta, kappa: model parameters
## E0: the initial values of process E_t, must be of the same length as beta
## burn_in: the burn in period for the process to reach the stationary distribution
## t_max: length of the series to be generated
## thinnings: indicator, whether the individual thinnings shall be returned
##
## the function returns the list of the trajectories of process X, E and if
## thinnings = TRUE, also the individual thinnings
sim_inarmapq <- function(tau, beta, kappa, E0, burn_in = 100, t_max = 1e4,
                         thinnings = FALSE) {

  # Check the validity of coefficients
  if (sum(beta) >= 1) {
    stop("sum of betas must be less than 1.")
  }
  if (sum(kappa) >= 1) {
    stop("sum of kappas must be less than 1.")
  }

  p <- length(kappa)
  q <- length(beta)
  max_order <- max(c(p, q))

  if (length(E0) != q) {
    stop("Length of 'E0' must be the length of the parameter vector 'beta'.")
  }

  # Generate the innovations
  I <- rpois(t_max + burn_in + max_order, tau)

  # Pad the parameter vectors so that they are of the same length
  beta_vec <- c(beta, rep(0, max_order - q), 1 - sum(beta))
  kappa_vec <- c(kappa, rep(0, max_order - p), 1 - sum(kappa))

  # Allocate the containers
  E <- X <- integer(t_max + burn_in + max_order)
  E_thin <- matrix(nrow = t_max + burn_in + max_order, ncol = max_order + 1)
  X_thin <- E_thin
  L <- integer(q + 1)
  C <- integer(p + 1)

  # Initialize
  E[1:q] <- E0

  for (k in 1:(t_max + burn_in + max_order)) {

    # Do the MA thinning and update the processes
    L <- rmultinom(1, E[k], beta_vec)[, 1]
    E_thin[k, ] <- L
    E[k + 1:q] <- E[k + 1:q] + L[1:q]
    X[k] <- L[max_order + 1] + I[k]

    # Do the AR thinning and update the processes
    C <- rmultinom(1, X[k], kappa_vec)[, 1]
    X_thin[k, ] <- C
    E[k + 1:p] <- E[k + 1:p] + C[1:p]
  }

  ret_list <- list(E = E[1:t_max + burn_in], X = X[1:t_max + burn_in])
  if (thinnings) {  # Add the individual thinning steps into the returns
    ret_list$E_thin <- E_thin[1:t_max + burn_in, ]
    ret_list$X_thin <- X_thin[1:t_max + burn_in, ]
  }

  return(ret_list)
}

## Function calculating the exact autocorrelation fun of the Poisson INARMA(p, q) model
## function arguments:
## max_d: the maximum lag
## beta, kappa: model parameters
##
## the function returns the vector of autocorrelation values up to lag max_d
acf_exact <- function(max_d, beta, kappa) {
  q <- length(beta)
  pi <- numeric(max_d + 1)

  pi[1] <- 1
  for (k in 1:max_d) {
    pr <- beta[1:k] * pi[k:1]
    pi[k + 1] <- sum(pr[!is.na(pr)])
  }

  p <- length(beta)
  u <- numeric(max_d + 1)
  u[1] <- 1
  for (k in 1:max_d) {
    pr <- kappa[1:k] * pi[k:1]
    u[k + 1] <- sum(pr[!is.na(pr)])
  }
  true_cor <- numeric(max_d + 1)
  true_cor[1] <- 1
  for (k in 1:max_d) {
    pr <- u[1:k + 1] * true_cor[k:1]
    true_cor[k + 1] <- (1 - sum(beta)) * (sum(pr[!is.na(pr)]))
  }

  return(true_cor)
}

## Function calculating the exact variance of the INARMA(1, 1) model
## with Poisson, Hermite, or Negative binomial innovation distribution
## function arguments:
## beta, kappa, tau, sigma2_tau: model parameters
##
## the function returns the vector of autocorrelation values up to lag max_d
var_exact_11 <- function (beta, kappa, tau, sigma2_tau) {
  xi <- beta + (1 - beta) * kappa
  mu <- tau / (1 - kappa)

  true_var <- mu * kappa * (1 + beta) / (1 + xi) +
    (1 - kappa * (1 + beta) / (1 + xi)) * sigma2_tau / (1 - kappa)
  return(true_var)
}

## Function calculating the exact autocorrelation fun of the INARMA(1, 1) model
## with Poisson, Hermite, or Negative binomial innovation distribution
## function arguments:
## max_d: the maximum lag
## beta, kappa: model parameters
## family: one of "Poisson", "Hermite", or "NegBin"
##
## the function returns the vector of autocorrelation values up to lag max_d
acf_exact_11 <- function (max_d, beta, kappa, tau, sigma2_tau) {
  xi <- beta + (1 - beta) * kappa

  # ACF from lag 1 onwards
  acf_true <- xi^(0:(max_d - 1)) * (1 - beta) * kappa *
    (1 +
       kappa * beta * (sigma2_tau - tau) /
       ((1 + beta) * ((1 - kappa) * sigma2_tau + kappa * tau) + (1 - beta) * kappa * sigma2_tau))
  acf_true <- c(1, acf_true)  # append 1
  return(acf_true)
}

## Function calculating the exact autocorrelation fun of the Poisson
## INGARCH(p, q) model borrowing the functionality from the tscount package
## function arguments:
## max_d: the maximum lag
## beta, kappa, tau: model parameters in the thinning formulation
##
## the function returns the vector of autocorrelation values up to lag max_d
acf_exact_ingarch <- function(max_d, beta, kappa, tau) {
  true_cor <- ingarch.acf(tau * (1 - sum(beta)), past_obs = kappa * (1 - sum(beta)),
                          lag.max = max_d, past_mean = beta, plot = FALSE)

  return(true_cor)
}

## Function calculating the exact mean of the Poisson INGARCH(p, q) model
## borrowing the functionality from the tscount package
## function arguments:
## beta, kappa, tau: model parameters in the thinning formulation
##
## the function returns the vector of autocorrelation values up to lag max_d
mean_exact_ingarch <- function(beta, kappa, tau) {
  true_cor <- ingarch.mean(tau * (1 - sum(beta)),
                           past_obs = kappa * (1 - sum(beta)), past_mean = beta)

  return(true_cor)
}

## Function calculating the exact variance of the Poisson INGARCH(p, q) model
## borrowing the functionality from the tscount package
## function arguments:
## beta, kappa, tau: model parameters in the thinning formulation
##
## the function returns the vector of autocorrelation values up to lag max_d
var_exact_ingarch <- function(beta, kappa, tau) {
  true_cor <- ingarch.var(tau * (1 - sum(beta)),
                           past_obs = kappa * (1 - sum(beta)), past_mean = beta)

  return(true_cor)
}

############################################################################
## Functions for mapping the parameter pairs from a triangle to a square and
## vice versa
############################################################################

# We can map the parameter pairs from the triangle to a square or the whole
# plane. For the square we would run the constrained optimization. For the
# unconstrained variant, the estimation often comes out better.

# Transform the parameter pair, for which par[1] + par[2] < 1 and
# 0 < par[1], par[2] < 1, onto the unit interval space [0, 1]^2
orig_par_to_box <- function (par) {
  u <- par[1] + par[2]
  v <- par[1] / u

  return(c(u, v))
}

# Transform the parameter pair on the unit interval space [0, 1]^2 back to the
# original space, where par[1] + par[2] < 1 and
# 0 < par[1], par[2] < 1
box_par_to_orig <- function (par) {

  orig1 <- par[1] * par[2]
  orig2 <- par[1] * (1 - par[2])

  return(c(orig1, orig2))
}

# Transform the parameter pair, for which par[1] + par[2] < 1 and
# 0 < par[1], par[2] < 1, onto the unit interval space [0, 1]^2
orig_par_to_unconstrained <- function (par) {
  u <- par[1] + par[2]
  v <- par[1] / u

  a <- log(u) - log(1 - u)
  b <- log(v) - log(1 - v)

  return(c(a, b))
}

# Transform the parameter pair on the unconstrained space R^2 back to the
# original space, where par[1] + par[2] < 1 and
# 0 < par[1], par[2] < 1
unconstrained_par_to_orig <- function (par) {
  u <- exp(par[1]) / (exp(par[1]) + 1)
  v <- exp(par[2]) / (exp(par[2]) + 1)

  orig1 <- u * v
  orig2 <- u * (1 - v)

  return(c(orig1, orig2))
}

###############################################################################
## Functions for calculating the standard errors for parameters on the original
## scale using the delta method
###############################################################################

## Function for calculating the standard errors of the parameters o the original
## scale
##
## function arguments:
## pair: a pair transformed model parameters (kappa, or beta) on the 2D plane
##  (from the function orig_par_to_unconstrained())
## cov_raw_part: block of the covariance matrix corresponding to the transformed
##  parameter pair
##
## the function returns the standard errors of the parameter pair, (kappa or
##  beta) on the original scale
delta_method_for_pairs <- function (pair, cov_raw_part) {
  # Calculate the gradient
  gradient <- matrix(
    nrow = 2,
    ncol = 2,
    byrow = TRUE,
    c(
      exp(pair[1]) * exp(pair[2]) / ((1 + exp(pair[1]))^2 * (1 + exp(pair[2]))),
      exp(pair[1]) * exp(pair[2]) / ((1 + exp(pair[1])) * (1 + exp(pair[2]))^2),
      exp(pair[1]) / ((1 + exp(pair[1]))^2 * (1 + exp(pair[2]))),
      -exp(pair[1]) * exp(pair[2]) / ((1 + exp(pair[1])) * (1 + exp(pair[2]))^2)
    )
  )
  # Do the delta method and return
  return(sqrt(diag(gradient %*% cov_raw_part %*% t(gradient))))
}

## Function for calculating the standard errors of the parameters o the original
## scale
##
## function arguments:
## pars: transformed model parameters, tau on the log-scale, kappa and beta on
##  the logit-scale, or on the 2D plane if dealing with a 2nd order INARMA
## hessian_transformed: the hessian matrix from the optimization corresponding
##  to the transforme parameters
##
## the function returns the standard errors of parameters tau, kappa and beta on
##  the original scale
get_ses_orig <- function (pars, hessian_transformed,
                          family = c("Poisson", "Hermite", "NegBin")) {

  # add small value to the diagonal to avoid numerical issues
  to_solve <- -( hessian_transformed + diag(10^-6, dim(hessian_transformed)[1]))
  cov_raw <- solve(to_solve)

  # Extract parameter values
  pars_names <- names(pars)
  log_tau <- pars[grepl("tau", pars_names)]
  where_kappa <- grepl("kappa", pars_names)
  where_beta <- grepl("beta", pars_names)
  transformed_kappa <- pars[where_kappa]
  transformed_beta <- pars[where_beta]

  # Standard error for tau
  tau_se <- as.numeric(sqrt(cov_raw["log_tau", "log_tau"]) * exp(log_tau))

  # Standard error for kappa
  if (length(transformed_kappa) == 1) {  # kappa is 1D
    kappa_se <- as.numeric(sqrt(cov_raw["logit_kappa", "logit_kappa"]) * exp(transformed_kappa) / (1 + exp(transformed_kappa))^2)
  } else {  # kappa is 2D
    kappa_se <- delta_method_for_pairs(transformed_kappa, cov_raw[where_kappa, where_kappa])
  }

  # Standard error for beta
  if (length(transformed_beta) == 1) {  # beta is 1D
    beta_se <- as.numeric(sqrt(cov_raw["logit_beta", "logit_beta"]) * exp(transformed_beta) / (1 + exp(transformed_beta))^2)
  } else {  # beta is 2D
    beta_se <- delta_method_for_pairs(transformed_beta, cov_raw[where_beta, where_beta])
  }

  ret <- c(tau_se, kappa_se, beta_se)

  # Standard error for psi if required by the distributional family
  if (family == "Hermite") {  # Hermite family
    transformed_psi <- pars[grepl("psi", pars_names)]
    psi_se <- as.numeric(sqrt(cov_raw["logit_psi", "logit_psi"]) * exp(transformed_psi) / (1 + exp(transformed_psi))^2)
    ret <- c(ret, psi_se)
  } else if (family == "NegBin") {  # beta is 2D
    transformed_psi <- pars[grepl("psi", pars_names)]
    psi_se <- as.numeric(sqrt(cov_raw["log_psi", "log_psi"]) * exp(transformed_psi))
    ret <- c(ret, psi_se)
  }

  if (anyNA(ret) || any(ret < 0)) {
    ret[which(ret < 0)] <- 0
    warning("At least one standard error will be zero.")
  }

  return(ret)
}

########################
## Functions to optimize
########################

## Function for evaluating the likelihood of an AR model approximating the
## INARMA(1, 1) model
##
## function arguments:
## pars: transformed model parameters, kappa and beta on the logit-scale, tau on the log-scale
## X: the observed count series
## lag_max: order of the AR model to be used for the approximation
##
## the function returns the likelihood of the approximating AR model
llik_ar_based_11 <- function (pars, X, lag_max = 8,
                              family = c("Poisson", "Hermite", "NegBin")) {

  # Extract parameter values
  tau <- exp(pars["log_tau"])
  kappa <- exp(pars["logit_kappa"]) / (1 + exp(pars["logit_kappa"]))
  beta <- exp(pars["logit_beta"]) / (1 + exp(pars["logit_beta"]))

  if (family == "Poisson") {
    psi <- 0
  } else if (family == "Hermite") {
    psi <- exp(pars["logit_psi"]) / (1 + exp(pars["logit_psi"]))
  } else if (family == "NegBin") {
    psi <- exp(pars["log_psi"])
  }

  sigma2_tau <- (1 + psi) * tau

  # Calculate the exact autocorrelation
  acf_inarma <- acf_exact_11(lag_max, beta, kappa, tau, sigma2_tau)

  # Solve the Yule-Walker equations
  YW_mat <- Toeplitz(acf_inarma[-length(acf_inarma)])
  ar_coeffs <- try(solve(YW_mat + diag(1e-6, lag_max), acf_inarma[-1]))

  if (class(ar_coeffs) == "try-error") {
    # Return an infinitely low value of the likelihood, if the Y-W equations
    # can not be solved
    return(-Inf)
  } else {
    # Reconstruct the noise
    lagged_obs_matrix <- sapply(0:(lag_max - 1), lag, x = X)[-(1:(lag_max - 1)), ]
    noise_reconstructed <- X[-(1:lag_max)] - lagged_obs_matrix[-nrow(lagged_obs_matrix), ] %*% ar_coeffs

    # Calculate the parameters of the white noise by matching the mean and
    # variance of the INARMA and AR models
    wn_mean <- tau * (1 - sum(ar_coeffs)) / (1 - kappa)
    wn_var <- var_exact_11(beta, kappa, tau, sigma2_tau) * (1 -  sum(ar_coeffs * acf_inarma[-1]))

    # Calculate the likelihood
    llik_ar <- sum(log(dnorm(noise_reconstructed, mean = wn_mean, sd = sqrt(wn_var))))

    return(llik_ar)
  }
}

## Function for evaluating the likelihood of an AR model approximating the
## INGARCH(p, q) model for p, q <= 2
##
## function arguments:
## pars: transformed model parameters (from the thinning formulation),
##  kappa and beta on the logit-scale, tau on the log-scale
## X: theobserved count series
## lag_max: order of the AR model to be used for the approximation
##
## the function returns the likelihood of the approximating AR model
llik_ar_based_ingarch <- function (pars, X, lag_max = 8) {

  # Extract parameter values
  pars_names <- names(pars)
  transformed_kappa <- pars[grepl("kappa", pars_names)]
  transformed_beta <- pars[grepl("beta", pars_names)]

  tau <- exp(pars[grepl("tau", pars_names)])

  # If kappa is one-dimensional, do the logit transformation, if 2-dimensional,
  # map it to the real plane (unconstrained space).
  if (length(transformed_kappa) == 1) {
    kappa <- exp(transformed_kappa) / (1 + exp(transformed_kappa))
  } else {
    kappa <- unconstrained_par_to_orig(transformed_kappa)
  }

  # If beta is one-dimensional, do the logit transformation, if 2-dimensional,
  # map it to the real plane (unconstrained space).
  if (length(transformed_beta) == 1) {
    beta <- exp(transformed_beta) / (1 + exp(transformed_beta))
  } else {
    beta <- unconstrained_par_to_orig(transformed_beta)
  }

  # Calculate the exact autocorrelation
  acf_inarma <- acf_exact_ingarch(lag_max, beta = beta, kappa = kappa, tau = tau)

  # Solve the Yule-Walker equations
  YW_mat <- Toeplitz(acf_inarma[-length(acf_inarma)])
  ar_coeffs <- try(solve(YW_mat + diag(1e-6, lag_max), acf_inarma[-1]))

  if (class(ar_coeffs) == "try-error") {
    # Return an infinitely low value of the likelihood, if the Y-W equations
    # can not be solved
    return(-Inf)
  } else {
    # Reconstruct the noise
    lagged_obs_matrix <- sapply(0:(lag_max - 1), lag, x = X)[-(1:(lag_max - 1)), ]
    noise_reconstructed <- X[-(1:lag_max)] - lagged_obs_matrix[-nrow(lagged_obs_matrix), ] %*% ar_coeffs
    ar_mean <- mean(noise_reconstructed)
    ar_var <- var(noise_reconstructed)

    # Calculate the parameters of the white noise
    wn_mean <- mean_exact_ingarch(beta, kappa, tau) * (1 - sum(ar_coeffs))
    wn_var <- var_exact_ingarch(beta, kappa, tau)  * (1 -  sum(ar_coeffs * acf_inarma[-1]))

    # Calculate the likelihood
    llik_ar <- sum(log(dnorm(noise_reconstructed, mean = wn_mean, sd = sqrt(wn_var))))

    return(llik_ar)
  }
}

## Function for evaluating the likelihood of an AR model approximating the
## INARMA(p, q) model, for p, q <= 2
##
## function arguments:
## pars: transformed model parameters, kappa and beta on the unconstrained scales,
##  tau on the log-scale
## X: the observed count series
## lag_max: order of the AR model to be used for the approximation
##
## the function returns the likelihood of the approximating AR model
llik_ar_based_higher <- function (pars, X, return_tau_hat = FALSE, lag_max = 10) {

  # Extract parameter values
  pars_names <- names(pars)
  transformed_kappa <- pars[grepl("kappa", pars_names)]
  transformed_beta <- pars[grepl("beta", pars_names)]

  tau <- exp(pars[grepl("tau", pars_names)])
  kappa <- unconstrained_par_to_orig(transformed_kappa)

  # If beta is one-dimensional, do the logit transformation, if 2-dimensional,
  # map it to the real plane (unconstrained space).
  if (length(transformed_beta) == 1) {
    beta <- exp(transformed_beta) / (1 + exp(transformed_beta))
  } else {
    beta <- unconstrained_par_to_orig(transformed_beta)
  }

  # Calculate the exact autocorrelation
  acf_inarma <- acf_exact(lag_max, beta = beta, kappa = kappa)

  # Solve the Yule-Walker equations
  YW_mat <- Toeplitz(acf_inarma[-length(acf_inarma)])
  ar_coeffs <- try(solve(YW_mat + diag(1e-6, lag_max), acf_inarma[-1]))

  ret <- list()

  if (class(ar_coeffs) == "try-error") {
    # Return an infinitely low value of the likelihood, if the Y-W equations
    # can not be solved
    llik_ar <- -Inf
  } else {
    # Reconstruct the noise
    lagged_obs_matrix <- sapply(0:(lag_max - 1), lag, x = X)[-(1:(lag_max - 1)), ]
    noise_reconstructed <- X[-(1:lag_max)] - lagged_obs_matrix[-nrow(lagged_obs_matrix), ] %*% ar_coeffs
    ar_mean <- mean(noise_reconstructed)
    ar_var <- var(noise_reconstructed)

    # Calculate the parameters of the white noise
    wn_mean <- tau * (1 - sum(ar_coeffs)) / (1 - sum(kappa))
    wn_var <- tau * (1 -  sum(ar_coeffs * acf_inarma[-1])) / (1 - sum(kappa))

    # Calculate the likelihood
    llik_ar <- sum(log(dnorm(noise_reconstructed, mean = wn_mean, sd = sqrt(wn_var))))

    return(llik_ar)
  }
}

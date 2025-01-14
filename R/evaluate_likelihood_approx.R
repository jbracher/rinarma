#' Function for evaluating the likelihood of an AR model approximating the
#' INARMA(1, 1) model
#' @param pars transformed model parameters, \eqn{\kappa} and \eqn{beta} on the
#'  logit-scale, tau on the log-scale
#' @param X the observed count series
#' @param lag_max order of the AR model to be used for the approximation
#' @param family family of the INARMA model; one of `"Poisson"`, `"Hermite"` or `"NegBin"`
#' @return the likelihood of the approximating AR model
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
  YW_mat <- pracma::Toeplitz(acf_inarma[-length(acf_inarma)])
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


#' Function for evaluating the likelihood of an AR model approximating the
#' INARMA(p, q) model, for p, q <= 2
#' @param pars transformed model parameters, \eqn{\kappa} and \eqn{\beta} on the
#'  unconstrained scales (2D plane), \eqn{\tau} on the log-scale
#' @param X the observed count series
#' @param lag_max order of the AR model to be used for the approximation
#' @return the likelihood of the approximating AR model
llik_ar_based_higher <- function (pars, X, lag_max = 10) {

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
  YW_mat <- pracma::Toeplitz(acf_inarma[-length(acf_inarma)])
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
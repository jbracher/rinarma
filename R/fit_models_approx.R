fit_inarma_approx <- function(observed, family = c("Poisson", "Hermite", "NegBin"),
                             order = c(p = 1, q = 1), lag_max = 10){

  if (all(order <= 0)) {
    stop("Invalid model order.")
  } else if(any(order > 2)) {
    stop("Fitting model of order higher than 2 is not supported.")
  }

  if (family %in% c("Poisson", "Hermite", "NegBin")) {
    if (family != "Poisson" && any(order) == 2) {
      stop("Fitting higher order Hermite, or Negative binomial INARMA models is
           not supported.")
    }
  } else {
    stop("Invalid family name.")
  }

  if (all(order == 1)) {
    ret <- ar_approx_11(observed, lag_max, family)
  } else {
    ret <- ar_approx_higher(observed, lag_max, order)
  }


  ret <- c(
    ret,
    family = family,
    observed = observed,
    lik_distr = NULL,
    fitted_values = NULL,
    fitted_variance = NULL,
    pearson_residuals = NULL,
    nobs = length(observed),
    fitting_method = "AR-approximation"
    )
  return(ret)
}

ar_approx_11 <- function (X, lag_max, family) {

  # Set 3 different starting values in case the optimization does not converge
  start_transformed <- matrix(c(1, 0, 0, 0, 0.2, 1, 1, -0.1, 0.2, -1, -1, -1,
                                rep(NA, 4)), nrow = 4, ncol = 4, byrow = TRUE)

  # Calculate the moment estimates and use it as a starting value for the
  # optmization in case the previous 3 are bad
  suppressWarnings(suppressMessages(start <- try(fit_inarma_moments(X, family = family))))

  if (class(start) == "try-error") {
    start_transformed[4, ] <- c(log_tau = 0.5, logit_kappa = 0, logit_beta = 0,
                                logit_psi = 0)
  } else {
    start <- start$coefficients
    start_transformed[4, 1:3] <- c(
      log_tau = unname(log(start["tau"])),
      logit_kappa = unname(log(start["kappa"] / (1 - start["kappa"]))),
      logit_beta = unname(log(start["beta"] / (1 - start["beta"])))
    )
    if (family == "NegBin") {
      start_transformed[4, 4] <- unname(log(start["psi"]))
      colnames(start_transformed) <- c("log_tau", "logit_kappa", "logit_beta",
                                       "log_psi")
    } else if (family == "Hermite") {
      start_transformed[4, 4] <- unname(log(start["psi"] / (1 - start["psi"])))
      colnames(start_transformed) <- c("log_tau", "logit_kappa", "logit_beta",
                                       "logit_psi")
    } else if (family == "Poisson") {
      start_transformed <- start_transformed[, 1:3]
      colnames(start_transformed) <- c("log_tau", "logit_kappa", "logit_beta")
    }
    start_transformed[4, is.nan(start_transformed[4, ]) |
                        is.infinite(start_transformed[4, ])] <- 0
  }

  # Refit until we get a non-problematic estimate
  refit_iter <- 1
  refit <- TRUE
  while (refit & refit_iter <= 4) {
    # Find the (approximate) maximum likelihood estimates
    op <- optim(
      par = start_transformed[refit_iter, ],
      fn = llik_ar_based_11,
      X = X,
      family = family,
      lag_max = lag_max,
      hessian = TRUE,
      control = list(fnscale = -1, maxit = 800)  # To change to maximization
    )

    # Extract and transform the results
    coeffs <- c(
      tau = unname(exp(op$par["log_tau"])),
      kappa = unname(exp(op$par["logit_kappa"]) / (1 + exp(op$par["logit_kappa"]))),
      beta = unname(exp(op$par["logit_beta"]) / (1 + exp(op$par["logit_beta"])))
    )
    if (family == "Hermite") {
      coeffs <- c(coeffs, psi = unname(exp(op$par["logit_psi"]) / (1 + exp(op$par["logit_psi"]))))
    } else if (family == "NegBin") {
      coeffs <- c(coeffs, psi = unname(exp(op$par["log_psi"])))
    }
    get_ses <- get_ses_orig(op$par, op$hessian, family = family)
    refit <- anyNA(get_ses$ses) | op$convergence != 0
    refit_iter <- refit_iter + 1
  }
  if (refit_iter > 4 && op$convergence != 0) {
    warning("Optimization did not converge.")
  }
  ret <- list(
    coefficients_raw = op$par,
    se_raw = sqrt(pmax(diag(get_ses$cov_raw), 0)),
    cov_raw = get_ses$cov_raw,
    coefficients = coeffs,
    se = get_ses$ses,
    dim = length(coeffs),
    loglikelihood = op$value,
    AIC = 2 * (-op$value + length(coeffs)),
    convergence = op$convergence,
    optim = op
    )
  return(ret)
}

ar_approx_higher <- function (X, lag_max, order) {

  # Set 4 different starting values
  start_transformed <- matrix(
    c(
      1, rep(0, order["p"]), rep(0, order["q"]),
      0.2, rep(1, order["p"]), rep(1, order["q"]),
      -0.1, rep(0.2, order["p"]), rep(1, order["q"]),
      -0.5, rep(0.2, order["p"]), rep(-0.2, order["q"])
    ),
    nrow = 4,
    ncol = 1 + sum(order),
    byrow = TRUE
  )
  colnames(start_transformed) <- c("log_tau", paste0("transformed_kappa", 1:order["p"]),
                paste0("transformed_beta", 1:order["q"]))

  # Refit until we get a non-problematic estimate
  refit_iter <- 1
  refit <- TRUE
  while (refit & refit_iter <= 4) {
    # Find the (approximate) maximum likelihood estimates
    op <- optim(
      par = start_transformed[refit_iter, ],
      fn = llik_ar_based_higher,
      X = X,
      lag_max = lag_max,
      hessian = TRUE,
      control = list(fnscale = -1, maxit = 800)  # To change to maximization
    )

    # Extract and transform the coefficients
    coeffs <- c(tau = unname(exp(op$par["log_tau"])))
    if (order["p"] == 2) {
      coeffs <- c(coeffs, kappa = unname(unconstrained_par_to_orig(op$par[2:3])))
    } else if (order["p"] == 1) {
      coeffs <- c(coeffs, kappa = unname(exp(op$par["transformed_kappa1"]) / (1 + exp(op$par["transformed_kappa1"]))))
    }

    if (order["q"] == 2) {
      coeffs <- c(coeffs, beta = unname(unconstrained_par_to_orig(op$par[(2 + p):(3 + p)])))
    } else if (order["q"] == 1) {
      coeffs <- c(coeffs, beta = unname(exp(op$par["transformed_beta1"]) / (1 + exp(op$par["transformed_beta1"]))))
    }

    # Extract and transform the standard errors
    get_ses <- get_ses_orig(op$par, op$hessian, family = "Poisson")
    refit <- anyNA(get_ses$ses) | op$convergence != 0
    refit_iter <- refit_iter + 1
  }
  if (refit_iter > 4 && op$convergence != 0) {
    warning("Optimization did not converge.")
  }
  ret <- list(
    coefficients_raw = op$par,
    se_raw = sqrt(pmax(diag(get_ses$cov_raw), 0)),
    cov_raw = get_ses$cov_raw,
    coefficients = coeffs,
    se = get_ses$ses,
    dim = length(coeffs),
    loglikelihood = op$value,
    AIC = 2 * (-op$value + length(coeffs)),
    convergence = op$convergence,
    optim = op
  )
  return(ret)
}
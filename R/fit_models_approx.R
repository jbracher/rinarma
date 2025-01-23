#' Fitting an INARMA model using the approximation by an AR model
#'
#' Inference for INARMA(p, q) model (p,q <= 2) using the Ar-approximation as
#' described in Bracher and Sobolova (2025).
#'
#' For the model definition see the `?fit_inarma` page of the documentation.
#'
#' The function uses a moment matching procedure to match the INARMA and AR
#' model of order `lag_max`.
#'
#' @export
#'
#' @examples
#' data("measles")
#' X <- measles$value
#' # Note: running the fit takes a little while.
#' \dontrun{
#' fit <- fit_inarma_approx(X, family = "Poisson", order = c(p = 1, q = 1))
#' summary(fit)
#' plot(fit, type = "fit")
#' }
#'
#'
#' @param observed a vector of observed count values
#' @param family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`
#' @param offspring the offspring distribution, or in other words, the type of
#' the \eqn{\kappa \bullet X_t} thinning; the AR-approximation procedure
#' supports only `"binomial"`.
#' @param order a vector with two named elements `p` and `q` denoting the order
#' of the INARMA(p, q) model to be fitted. p,q <= 2
#' @param lag_max the order of the approximate AR model.
#' @param control_optim a list of options passed to `optim`
#' @return an object of class `inarma`. This is a list with the following elements:
#' \describe{
#' \item{family}{the distribution family used.}
#' \item{coefficients_raw}{the estimated model coefficients on the internal scale.}
#' \item{se_raw}{the estimated standard errors on the internal scale.}
#' \item{cov_raw}{covariance matrix of the estimates on the internal scale.}
#' \item{coefficients}{the estimated model parameters transformed back to the natural scale.}
#' \item{se}{the estimated standard errors transformed back to the natural scale.}
#' \item{observed}{the vector of observed values provided by the user.}
#' \item{lik_distr}{`NULL` value for the AR-approximation}
#' \item{fitted_values}{fitted values calculated from the approximate AR model}
#' \item{fitted_variance}{the fitted conditional variances calculated from the
#' approximate AR model. In the case of the AR model, they are just the estimate
#'  of the white noise variance}
#' \item{pearson_residuals}{the Pearson residuals, again obtained from the AR model}
#' \item{dim}{the number of fitted parameters}
#' \item{loglikelihood}{the log-likelihood of the Approximate model}
#' \item{AIC}{the resulting (approximate) AIC}
#' \item{convergence}{indicates whether optimization converged}
#' \item{nobs}{the number of observations}
#' \item{optim}{return object of the call to `optim`}
#' \item{order}{order of the fitted INARMA model.}
#' \item{fitting_method}{the method used to fit the model, here `"AR-approximation"`.}
#' }
fit_inarma_approx <- function(observed, family = c("Poisson", "Hermite", "NegBin"),
                              offspring = "binomial",
                              order = c(p = 1, q = 1), lag_max = 10,
                              control_optim = NULL){

  if (offspring != "binomial") {
    stop("Other than binomial offspring is not supported.")
  }

  if (all(order <= 0)) {
    stop("Invalid model order.")
  } else if(any(order > 2)) {
    stop("Fitting model of order higher than 2 is not supported.")
  }

  if (all(order == 1)) {
    ret <- ar_approx_11(observed, lag_max, family,
                        control_optim = control_optim)
  } else {
    ret <- ar_approx_higher(observed, lag_max, order, family,
                            control_optim = control_optim)
  }

  ret <- c(
    ret,
    list(
      family = family,
      offspring = "binomial",
      observed = observed,
      lik_distr = NULL,
      nobs = length(observed),
      order = order,
      fitting_method = "AR-approximation"
    )
  )

  class(ret) <- "inarma"

  return(ret)
}

#' Function calculating the AR-approximation for an INARMA(1, 1) model.
#'
#' @param X a vector of observed count values
#' @param family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`
#' @param lag_max the order of the approximate AR model.
#' @param control_optim a list of options passed to `optim`
#' @return an object of class `inarma`. This is a list with the following elements:
#' \describe{
#' \item{coefficients_raw}{the estimated model coefficients on the internal scale.}
#' \item{se_raw}{the estimated standard errors on the internal scale.}
#' \item{cov_raw}{covariance matrix of the estimates on the internal scale.}
#' \item{coefficients}{the estimated model parameters transformed back to the natural scale.}
#' \item{se}{the estimated standard errors transformed back to the natural scale.}
#' \item{fitted_values}{fitted values calculated from the approximate AR model}
#' \item{fitted_variance}{the fitted conditional variances calculated from the
#' approximate AR model. In the case of the AR model, they are just the estimate
#'  of the white noise variance}
#' \item{pearson_residuals}{the Pearson residuals, again obtained from the AR model}
#' \item{dim}{the number of fitted parameters}
#' \item{loglikelihood}{the log-likelihood of the Approximate model}
#' \item{AIC}{the resulting (approximate) AIC}
#' \item{convergence}{indicates whether optimization converged}
#' \item{optim}{return object of the call to `optim`}
#' }
ar_approx_11 <- function (X, lag_max, family, control_optim = NULL) {

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

  # Check the control for the optimization procedure
  if (is.null(control_optim) || !is.list(control_optim)) {
    control_optim <- list(fnscale = -1, maxit = 800)
  } else {
    control_optim$fnscale <- -1
    control_optim$maxit <- max(control_optim$maxit, 800, na.rm = TRUE)
  }

  # Refit until we get a non-problematic estimate
  refit_iter <- 1
  refit <- TRUE
  while (refit & refit_iter <= 4) {
    # Find the (approximate) maximum likelihood estimates
    op <- optim(
      par = start_transformed[refit_iter, ],
      fn = llik_ar_based_11,
      return_fitted = FALSE,
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

  fitted_vals <- llik_ar_based_11(op$par, X, lag_max = lag_max,
                                 return_fitted = TRUE, family = family)

  ret <- list(
    coefficients_raw = op$par,
    se_raw = sqrt(pmax(diag(get_ses$cov_raw), 0)),
    cov_raw = get_ses$cov_raw,
    coefficients = coeffs,
    se = get_ses$ses,
    fitted_values = fitted_vals$fitted_values,
    fitted_variance = fitted_vals$fitted_variance,
    pearson_residuals = (X - fitted_vals$fitted_values) / sqrt(fitted_vals$fitted_variance),
    dim = length(coeffs),
    loglikelihood = op$value,
    AIC = 2 * (-op$value + length(coeffs)),
    convergence = op$convergence,
    optim = op
    )
  return(ret)
}

#' Function calculating the AR-approximation for an INARMA(p, q) model, where
#' p,q <= 2.
#'
#' @param X a vector of observed count values
#' @param family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`
#' @param lag_max the order of the approximate AR model.
#' @param order a vector with two named elements `p` and `q` denoting the order
#' of the INARMA(p, q) model to be fitted. p,q <= 2
#' @param control_optim a list of options passed to `optim`
#' @return an object of class `inarma`. This is a list with the following elements:
#' \describe{
#' \item{coefficients_raw}{the estimated model coefficients on the internal scale.}
#' \item{se_raw}{the estimated standard errors on the internal scale.}
#' \item{cov_raw}{covariance matrix of the estimates on the internal scale.}
#' \item{coefficients}{the estimated model parameters transformed back to the natural scale.}
#' \item{se}{the estimated standard errors transformed back to the natural scale.}
#' \item{fitted_values}{fitted values calculated from the approximate AR model}
#' \item{fitted_variance}{the fitted conditional variances calculated from the
#' approximate AR model. In the case of the AR model, they are just the estimate
#'  of the white noise variance}
#' \item{pearson_residuals}{the Pearson residuals, again obtained from the AR model}
#' \item{dim}{the number of fitted parameters}
#' \item{loglikelihood}{the log-likelihood of the Approximate model}
#' \item{AIC}{the resulting (approximate) AIC}
#' \item{convergence}{indicates whether optimization converged}
#' \item{optim}{return object of the call to `optim`}
#' }
ar_approx_higher <- function (X, lag_max, order, family, control_optim = NULL) {

  # Set 4 different starting values
  start_transformed <- matrix(
    c(
      1, rep(0, order["p"]), rep(0, order["q"]), 0.2,
      0.2, rep(1, order["p"]), rep(1, order["q"]), 0,
      -0.1, rep(0.2, order["p"]), rep(1, order["q"]), -0.2,
      -0.5, rep(0.2, order["p"]), rep(-0.2, order["q"]), -1
    ),
    nrow = 4,
    ncol = 2 + sum(order),
    byrow = TRUE
  )

  if (family == "NegBin") {
    colnames(start_transformed) <- c("log_tau", paste0("transformed_kappa", 1:order["p"]),
                                     paste0("transformed_beta", 1:order["q"]), "log_psi")
  } else if (family == "Hermite") {
    colnames(start_transformed) <- c("log_tau", paste0("transformed_kappa", 1:order["p"]),
                                     paste0("transformed_beta", 1:order["q"]), "logit_psi")
  } else if (family == "Poisson") {
    start_transformed <- start_transformed[, -(2 + sum(order))]
    colnames(start_transformed) <- c("log_tau", paste0("transformed_kappa", 1:order["p"]),
                                     paste0("transformed_beta", 1:order["q"]))
  }

  # Check the control for the optimization procedure
  if (is.null(control_optim) || !is.list(control_optim)) {
    control_optim <- list(fnscale = -1, maxit = 1000)
  } else {
    control_optim$fnscale <- -1
    control_optim$maxit <- max(control_optim$maxit, 1000, na.rm = TRUE)
  }

  # Check the control for the optimization procedure
  if (is.null(control_optim) || !is.list(control_optim)) {
    control_optim <- list(fnscale = -1, maxit = 800)
  } else {
    control_optim$fnscale <- -1
    control_optim$maxit <- max(control_optim$maxit, 800, na.rm = TRUE)
  }

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
      family = family,
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

    if (family == "Hermite") {
      coeffs <- c(coeffs, psi = unname(exp(op$par["logit_psi"]) / (1 + exp(op$par["logit_psi"]))))
    } else if (family == "NegBin") {
      coeffs <- c(coeffs, psi = unname(exp(op$par["log_psi"])))
    }

    # Extract and transform the standard errors
    get_ses <- get_ses_orig(op$par, op$hessian, family = family)
    refit <- anyNA(get_ses$ses) | op$convergence != 0
    refit_iter <- refit_iter + 1
  }
  if (refit_iter > 4 && op$convergence != 0) {
    warning("Optimization did not converge.")
  }

  fitted_vals <- llik_ar_based_higher(op$par, X, lag_max = lag_max,
                                      return_fitted = TRUE, family = family,
                                      return_approx_OK = TRUE)

  if (!fitted_vals$approx_OK) {
    warning("Error of the evaluation of the INARMA autocorellation is higher than 1e-16.
            Consider lowering the maximum lag used for the approximating AR model.")
  }

  ret <- list(
    coefficients_raw = op$par,
    se_raw = sqrt(pmax(diag(get_ses$cov_raw), 0)),
    cov_raw = get_ses$cov_raw,
    coefficients = coeffs,
    se = get_ses$ses,
    fitted_values = fitted_vals$fitted_values,
    fitted_variance = fitted_vals$fitted_variance,
    pearson_residuals = (X - fitted_vals$fitted_values) / sqrt(fitted_vals$fitted_variance),
    dim = length(coeffs),
    loglikelihood = op$value,
    AIC = 2 * (-op$value + length(coeffs)),
    convergence = op$convergence,
    optim = op
  )
  return(ret)
}
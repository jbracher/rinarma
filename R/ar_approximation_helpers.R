#' Transform a 2D parameter pair, for which par[1] + par[2] < 1 and
#' 0 < par[1], par[2] < 1, onto the unconstrained parametric space R^2
#' @param par 2D parameter of an INARMA model, either \eqn{(\kappa_1, \kappa_2)},
#'  or  \eqn{(\beta_1, \beta_2)}
#' @return transformed parameter (vector with two elements)
orig_par_to_unconstrained <- function (par) {
  # Transform to a unit interval
  u <- par[1] + par[2]
  v <- par[1] / u

  # Logit transformation
  a <- log(u) - log(1 - u)
  b <- log(v) - log(1 - v)

  return(c(a, b))
}

#' Back-transform a 2D parameter pairon the unconstrained space R^2 back to the
#' original space, where par[1] + par[2] < 1 and 0 < par[1], par[2] < 1
#' @param par transformed 2D parameter of an INARMA model, either \eqn{g(\kappa_1, \kappa_2)},
#'  or  \eqn{g(\beta_1, \beta_2)}
#' @return back-transformed parameter (vector with two elements)
unconstrained_par_to_orig <- function (par) {
  u <- exp(par[1]) / (exp(par[1]) + 1)
  v <- exp(par[2]) / (exp(par[2]) + 1)

  orig1 <- u * v
  orig2 <- u * (1 - v)

  return(c(orig1, orig2))
}

#' Delta method for calculating the standard errors of the 2D INARMA parameters
#'  on the original scale
#' @param pair a pair of transformed model parameters. Either
#'  \eqn{g(\kappa_1, \kappa_2)}, or \eqn{g(\beta_1, \beta_2)}
#' @param cov_raw_part block of the covariance matrix corresponding to the transformed
#'  parameter pair
#' @return the standard errors of the parameter pair on the original scale
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

#' Function calculating the standard errors of all INARMA parameters on the
#' original scale
#' @param pars transformed model parameters, \eqn{\tau} on the log-scale, \eqn{\kappa} and
#'  \eqn{\beta} on the logit-scale if fitting INARMA(1, 1), or on the 2D plane if
#'  fitting a 2nd order INARMA
#' @param hessian_transformed hessian matrix from the optimization corresponding
#'  to the transformed parameters
#' @return a list of quantities
#' \describe{
#' \item{ses}{standard errors of parameters \eqn{\tau}, \eqn{\kappa} and
#'  \eqn{\beta} on the original scale}
#' \item{cov_raw}{the covariance matrix of the parameter estimate on the transformed scale.}
#' }
get_ses_orig <- function (pars, hessian_transformed,
                          family = c("Poisson", "Hermite", "NegBin")) {

  # add small value to the diagonal to avoid numerical issues
  to_solve <- -( hessian_transformed - diag(10^-6, dim(hessian_transformed)[1]))
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
    kappa_se <- as.numeric(sqrt(cov_raw[where_kappa, where_kappa]) * exp(transformed_kappa) / (1 + exp(transformed_kappa))^2)
  } else {  # kappa is 2D
    kappa_se <- delta_method_for_pairs(transformed_kappa, cov_raw[where_kappa, where_kappa])
  }

  # Standard error for beta
  if (length(transformed_beta) == 1) {  # beta is 1D
    beta_se <- as.numeric(sqrt(cov_raw[where_beta, where_beta]) * exp(transformed_beta) / (1 + exp(transformed_beta))^2)
  } else {  # beta is 2D
    beta_se <- delta_method_for_pairs(transformed_beta, cov_raw[where_beta, where_beta])
  }

  ret <- c(tau = tau_se, kappa = kappa_se, beta = beta_se)

  # Standard error for psi if required by the distributional family
  if (family == "Hermite") {  # Hermite family
    transformed_psi <- pars[grepl("psi", pars_names)]
    psi_se <- as.numeric(sqrt(cov_raw["logit_psi", "logit_psi"]) * exp(transformed_psi) / (1 + exp(transformed_psi))^2)
    ret <- c(ret, psi = psi_se)
  } else if (family == "NegBin") {  # beta is 2D
    transformed_psi <- pars[grepl("psi", pars_names)]
    psi_se <- as.numeric(sqrt(cov_raw["log_psi", "log_psi"]) * exp(transformed_psi))
    ret <- c(ret, psi = psi_se)
  }

  if (anyNA(ret) || any(ret < 0)) {
    ret[which(ret < 0)] <- 0
    warning("At least one standard error will be zero.")
  }

  return(list(ses = ret, cov_raw = cov_raw))
}

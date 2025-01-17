
#' Function calculating the exact autocorrelation function of the Poisson
#' INARMA(p, q) model with parameters \eqn{\kappa} and \eqn{\beta}
#'
#' @param max_d maximum lag
#' @param beta,kappa parameters of the model, both `beta` and `kappa` can be
#'   vectors. Each of them must sum up to less than 1.
#' @return vector of the autocorrelation values from lag 0 to `max_d`
acf_exact_pois <- function(max_d, beta, kappa) {
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

#' Function calculating the exact (analytical) autocovariance function of the
#' INARMA(p, q) model with parameters \eqn{\kappa}, \eqn{\beta}, \eqn{tau} and
#' \eqn{sigma^2_\tau}. For the non-Poisson case (\eqn{tau} != \eqn{sigma^2_\tau}),
#' an infinite sum must be evaluated. The precision of this evaluation is typically
#' below 1e-8.
#'
#' @param max_d maximum lag
#' @param beta,kappa,tau,sigma2_tau parameters of the model, both `beta` and `kappa` can be
#'   vectors. Each of them must sum up to less than 1.
#' @return list with 2 elements
#' \describe{
#' \item{acov}{vector of the autocovariance values from lag 0 to `max_d`.}
#' \item{approx_OK}{indicator, whether the error of the infinite sum evaluation falls below 1e-8.
#' For the Poisson case it is always `TRUE`.}
#' }
acov_exact <- function(max_d, beta, kappa, tau = NULL, sigma2_tau = NULL) {
  if (is.null(sigma2_tau) || (all(!is.null(c(tau, sigma2_tau))) && sigma2_tau == tau)) {
    return(
      list(
        acov = acf_exact_pois(max_d, beta, kappa) * tau / (1 - sum(kappa)),
        approx_OK = TRUE
      )
    )
  } else {
    chunk <- max(max_d, 50)
    pois_acf <- acf_exact_pois(5 * chunk + max_d, beta, kappa)
    lagged_acf_matrix <- sapply(max_d:0, lag, x = pois_acf)[-(1:max_d), ]
    approx_sums <- t(lagged_acf_matrix[1:chunk, 1]) %*% lagged_acf_matrix[1:chunk, ]
    difference <- t(lagged_acf_matrix[chunk + 1, 1]) %*% lagged_acf_matrix[chunk + 1, ]

    iter <- 1
    while (any(difference > 1e-8) && iter < 5) {
      difference <- t(lagged_acf_matrix[(1:chunk) + iter * chunk, 1]) %*%
        lagged_acf_matrix[(1:chunk) + iter * chunk, ]
      approx_sums <- approx_sums + difference
      iter <- iter + 1
    }
    cov_true <- tau / (1 - sum(kappa)) * pois_acf[1:(max_d + 1)] + (sigma2_tau - tau) * approx_sums

    return(
      list(
        acov = as.numeric(cov_true),
        approx_OK = !(any(difference > 1e-8) && iter == 3))
    )
  }
}

#' Function calculating the exact autocorrelation function of the
#' INARMA(p, q) model with general innovation distribution and parameters
#' \eqn{\kappa}, \eqn{\beta}, \eqn{\tau} and \eqn{\sigma^2_\tau}
#'
#' @param max_d maximum lag
#' @param beta,kappa parameters of the model, both `beta` and `kappa` can be
#'   vectors. Each of them must sum up to less than 1.
#' @return vector of the autocorrelation values from lag 0 to `max_d`
acf_exact <- function(max_d, beta, kappa, tau = NULL, sigma2_tau = NULL) {
  acov <- acov_exact(max_d, beta, kappa, tau, sigma2_tau)
  return(list(acov = acov$acov / acov$acov[1], approx_OK = acov$approx_OK))
}

#' Function calculating the exact variance of the INARMA(1, 1) model with
#' general innovations.
#'
#' @param max_d maximum lag
#' @param beta,kappa,tau,sigma2_tau parameters of the model. For `sigma2_tau = tau`
#'  we get the Poisson model.
#'\describe{
#' \item{fitting_method}{the method used to fit the model, here `"maximum_likelihood"`.}
#' }
#' @return variance of the INARMA(1, 1) process
var_exact_11 <- function (beta, kappa, tau, sigma2_tau) {
  xi <- beta + (1 - beta) * kappa
  mu <- tau / (1 - kappa)

  true_var <- mu * kappa * (1 + beta) / (1 + xi) +
    (1 - kappa * (1 + beta) / (1 + xi)) * sigma2_tau / (1 - kappa)
  return(true_var)
}

#' Function calculating the exact autocorrelation function of INARMA(1, 1) with
#' general innovations.
#'
#' @param max_d maximum lag
#' @param beta,kappa,tau,sigma2_tau parameters of the model. For `sigma2_tau = tau`
#'  we get the Poisson model.
#' @return vector of the autocorrelation values from lag 0 to `max_d`
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

#' Function calculating the exact autocorrelation fun of the Poisson
#' INGARCH(p, q) model using the functionality from the tscount package
#' function arguments:
#' @param max_d maximum lag
#' @param beta,kappa,tau model parameters in the thinning formulation
#'
#' @return vector of the autocorrelation values from lag 0 to `max_d`
acf_exact_ingarch <- function(max_d, beta, kappa, tau) {
  true_cor <- tscount::ingarch.acf(tau * (1 - sum(beta)), past_obs = kappa * (1 - sum(beta)),
                          lag.max = max_d, past_mean = beta, plot = FALSE)

  return(true_cor)
}

#' Function calculating the exact mean of the Poisson INGARCH(p, q) model
#' using the functionality from the tscount package
#' function arguments:
#' @param beta,kappa,tau model parameters in the thinning formulation. Both
#' `beta` and `kappa` can be vectors. Each of them must sum up to less than 1.
#'
#' @return the mean of the INGARCH(p, q)
mean_exact_ingarch <- function(beta, kappa, tau) {
  true_mean <- tscount::ingarch.mean(tau * (1 - sum(beta)),
                           past_obs = kappa * (1 - sum(beta)), past_mean = beta)

  return(true_mean)
}

#' Function calculating the exact variance of the Poisson INGARCH(p, q) model
#' using the functionality from the tscount package
#' function arguments:
#' @param beta,kappa,tau model parameters in the thinning formulation. Both
#' `beta` and `kappa` can be vectors. Each of them must sum up to less than 1.
#'
#' @return the variance of the INGARCH(p, q)
var_exact_ingarch <- function(beta, kappa, tau) {
  true_var <- tscount::ingarch.var(tau * (1 - sum(beta)),
                          past_obs = kappa * (1 - sum(beta)), past_mean = beta)

  return(true_var)
}
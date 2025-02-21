#' Fitting an INGARCH model
#'
#' Maximum likelihood inference for IGARCH(1, 1) model as described in Bracher
#' and Sobolova (2024):
#'
#' \deqn{X_t = \theta * \bigl((1 - \beta) \circ E_t + I_t\bigr)}
#' \deqn{E_t = \beta \circ E_{t - 1} + \kappa \star X_{t - 1}}
#' where \eqn{\circ} denotes binomial thinning and \eqn{\star} denotes the
#' Poisson thinning. The type of thinning \eqn{*} depends on the distributional
#' family. See 'Details'.
#'
#' The two thinnings of \eqn{E_t} are coupled via
#' \deqn{[\beta \circ E_t, (1 - \beta) \circ E_t] \sim \text{Mult}(E_t, \beta, 1 - \beta).}
#' All other thinnings are independent of each other.
#' The immigration process \eqn{I_t} consists of independently and identically distributed random variables which
#' can be Poisson, Hermite or negative binomial. This distribution is charaterized by its mean, denoted by
#'  \eqn{\tau}, and potentially a dispersion parameter \eqn{\psi}. See Bracher and Sobolova (2024) for details.
#'
#' This display of the INGARCH model uses thinning
#' operators, which substantially differs from the common formulation of the
#' model in literature. For fitting an INGARCH model in its usual form
#' \deqn{X_t \sim \text{Pois}(\lambda_t),}
#' \deqn{\lambda_t = \nu + \alpha X_{t - 1} + \beta \lambda_{t - 1},}
#' we recommend using the `tscount` package (Liboschik T., Fokianos. K, Fried R.
#' (2017)).
#'
#' @details
#' The definition of thinning \eqn{*} depends on the distributional family.
#' For the Poisson INGARCH(1, 1) model, we have \eqn{\theta * 1 = 1}.
#'
#' For the Hermite INGARCH(1, 1) model we have
#' \deqn{\theta * 1 = \begin{cases}
#'                    1 \quad with probability 2 - \theta,\\
#'                    2 \quad with probability \theta - 1.
#'                    \end{cases}}
#' When fitting the Hermite INGARCH(1, 1) model, `fit_ingarch` function operates
#' internally with the overdispersion parameter \eqn{\psi}, from which
#' \eqn{\theta} is retrieved using formula
#' \deqn{\theta = \frac{2}{2 - \psi}}
#'
#' Finally for the Negative binomial INGARCH(1, 1) model, the thinning \eqn{ * }
#' is defined as
#' \deqn{\theta * 1 \sim \text{Log}\left(\frac{1}{1 + \psi}\right),}
#' where \eqn{\text{Log}(\pis)} denotes the logarithmic distribution with the
#' probability mass function.
#' \deqn{p(z) =  \frac{(1 - \pi)^z}{-z\log(\pi)}}
#' When fitting the Hermite INGARCH(1, 1) model, `fit_ingarch` function operates
#' internally with the overdispersion parameter \eqn{\psi}, from which
#' \eqn{\theta} is retrieved using formula
#' \deqn{\theta = \frac{\psi}{\log(1 + \psi)}.}
#'
#' @export
#'
#' @examples
#' data("measles")
#' X <- measles$value
#' \dontrun{
#' fit <- fit_ingarch(X, family = "Poisson")
#' }
#'
#'
#' @param observed a vector containing the observed count values
#' @param family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`
#' @param start initial values for the optimization routine (on the internal
#' scale; check the element `"optim"` of the return list).
#' @param return_se should standard errors be returned?
#' @param control_optim a list of options passed to `optim`
#' @return a list with the following elements:
#' \describe{
#' \item{family}{the distribution family used.}
#' \item{coefficients_raw}{the estimated model coefficients on the internal scale.}
#' \item{se_raw}{if `return_se == TRUE`: the estimated standard errors on the internal scale.}
#' \item{cov_raw}{covariance matrix of the estimates on the internal scale.}
#' \item{coefficients}{the estimated model parameters transformed back to the natural scale.}
#' \item{se}{the estimated standard errors transformed back to the natural scale.}
#' \item{observed}{the vector of observed values provided by the user.}
#' \item{fitted_values}{the fitted values as obtained from `lik_distr`}
#' \item{pearson_residuals}{the Pearson residuals}
#' \item{psi}{estimate and standard error of the parameter \eqn{\psi} from the
#'  GLM-formulation of the INARCH model. Included for making the transitions between the two formulations easier.}
#' \item{dim}{the number of fitted parameters}
#' \item{loglikelihood}{the log-likelihood of the fitted model}
#' \item{AIC}{the resulting AIC}
#' \item{convergence}{indicates whether optimization converged}
#' \item{nobs}{the number of observations}
#' \item{optim}{return object of the call to `optim`}
#' \item{fitting_method}{the method used to fit the model, here `"maximum_likelihood"`.}
#' }
fit_ingarch <- function(observed, family = c("Poisson", "Hermite", "NegBin"),
                        start = NULL, return_se = TRUE, control_optim = NULL){

  # Check the starting values
  if(is.null(start)) {
    start_val <- c(log_tau = 1, logit_beta = 0, logit_kappa = 0, log_mean_E1 = 0.5)
  } else{
    names_ok <- sum((names(start) %in% c("log_tau", "logit_kappa", "logit_beta",
                                    "log_mean_E1"))) == 4
    names_ok <- names_ok &&
      (family == "Poisson" || (family == "Hermite" && any(names(start) %in% "logit_psi")) ||
         (family == "NegBin" && any(names(start) %in% "log_psi")))
    if (names_ok) {
      start_val <- start
    } else{
      warning("Starting values must be supplied on the internal scale. Check the
              element 'optim' of the return list. Starting with the default
              starting values.")
      start_val <- c(log_tau = 1, logit_beta = 0, logit_kappa = 0, log_mean_E1 = 0.5)
    }
  }

  # negative log-likelihood as function of parameter vector,
  # depends on the family:
  if(family == "Poisson"){
    nllik <- function (pars) {
      -llik_ingarch(observed, tau = exp(pars["log_tau"]),
                    beta = exp(pars["logit_beta"]) / (1 + exp(pars["logit_beta"])),
                    kappa = exp(pars["logit_kappa"]) / (1 + exp(pars["logit_kappa"])),
                    mean_E1 = exp(pars["log_mean_E1"]))
    }
  }

  if(family == "Hermite"){
    if (is.null(start) || !names_ok) {
      start_val <- c(start_val, logit_psi = -0.5)
    }
    nllik <- function (pars) {
      -llik_hingarch(observed, tau = exp(pars["log_tau"]),
                     beta = exp(pars["logit_beta"]) / (1 + exp(pars["logit_beta"])),
                     kappa = exp(pars["logit_kappa"]) / (1 + exp(pars["logit_kappa"])),
                     psi = exp(pars["logit_psi"]) / (1 + exp(pars["logit_psi"])),
                     mean_E1 = exp(pars["log_mean_E1"]))
    }
  }

  if(family == "NegBin"){
    if (is.null(start) || !names_ok) {
      start_val <- c(start_val, log_psi = -0.5)
    }
    nllik <- function (pars) {
      -llik_nbingarch(observed, tau = exp(pars["log_tau"]),
                      beta = exp(pars["logit_beta"]) / (1 + exp(pars["logit_beta"])),
                      kappa = exp(pars["logit_kappa"]) / (1 + exp(pars["logit_kappa"])),
                      psi = exp(pars["log_psi"]),
                      mean_E1 = exp(pars["log_mean_E1"]))
    }
  }

  start_val[is.null(start_val) | is.infinite(start_val) | is.na(start_val)] <- 0

  # Optimize
  opt <- optim(start_val, nllik, control = control_optim, hessian = return_se)

  # structure results
  ret <- list()
  ret$family <- family
  ret$offspring <- "Poisson"

  # parameter estimates and standard errors:
  ret$coefficients_raw <- opt$par
  ret$se_raw <- ret$cov_raw <- NULL
  if(return_se){
    to_solve <- opt$hessian + diag(10^-6, dim(opt$hessian)[1]) # add small value to diagonal to avoid numerical issues
    ret$cov_raw <- solve(to_solve)
    ret$se_raw <- sqrt(diag(ret$cov_raw))
    if(any(diag(ret$cov_raw) < 0)){
      warning("Negative diagonal elements in inverse Fisher matrix - thresholding at zero. At least one estimated standard error will be zero.")
    }
    ret$se_raw <- sqrt(pmax(diag(ret$cov_raw), 0))
  }

  # compute standard errors on untransformed scale via delta method
  # (use as.numeric to get rid of naming of numeric vectors)
  se <- list(
    tau = as.numeric(ret$se_raw["log_tau"] * exp(ret$coefficients_raw["log_tau"])^2),
    beta = as.numeric(ret$se_raw["logit_beta"] * exp(ret$coefficients_raw["logit_beta"]) /
                        (1 + exp(ret$coefficients_raw["logit_beta"]))^2),
    kappa = as.numeric(ret$se_raw["logit_kappa"] * exp(ret$coefficients_raw["logit_kappa"])
                          / (1 + exp(ret$coefficients_raw["logit_kappa"]))^2),
    mean_E1 = as.numeric(ret$se_raw["log_mean_E1"] * exp(ret$coefficients_raw["log_mean_E1"])^2)
    )

  coefficients <- list(
    tau = as.numeric(exp(opt$par["log_tau"])),
    beta = as.numeric(exp(opt$par["logit_beta"]) / (1 + exp(opt$par["logit_beta"]))),
    kappa = as.numeric(exp(opt$par["logit_kappa"]) / (1 + exp(opt$par["logit_kappa"]))),
    mean_E1 = as.numeric(exp(ret$coefficients_raw["log_mean_E1"]))
  )

  # Add psi into the list of results and compute fitted values, depends on the family:
  if (family =="Poisson") {
    ret$fitted_values <- do.call(
      llik_ingarch,
      c(coefficients, vect = list(observed), return_fitted = TRUE)
    )$fitted

    ret$fitted_variance <- ret$fitted_values
    ret$lik_distr <- t(sapply(ret$fitted_values, FUN = dpois, x = 0:round(1.2 * max(observed))))
    ret$psi <- NA
    ret$coefficients <- unlist(coefficients)
    ret$se <- unlist(se)
  }
  if(family == "NegBin") {
    coefficients$psi <- as.numeric(exp(ret$coefficients_raw["log_psi"]))
    se$psi <- as.numeric(ret$se_raw["log_psi"] * exp(ret$coefficients_raw["log_psi"])^2)
    ret$fitted_values <- do.call(
      llik_nbingarch,
      c(coefficients, vect = list(observed), return_fitted = TRUE)
    )$fitted
    get_theta <- get_cluster_size(coefficients$psi, se$psi, family = "NegBin")
    coefficients$theta <- get_theta$theta
    se$theta <- get_theta$theta_se

    ret$lik_distr <-  t(
      sapply(
        ret$fitted_values, FUN = function (x) {dnbinom(0:round(1.2 * max(observed)), mu = x, size = 1 / coefficients$psi)}
        )
      )
    ret$fitted_variance <- ret$fitted_values * (1 + coefficients$psi)
    ret$psi <- c(coefficients["psi"], se_psi = unname(se["psi"]))
    ret$coefficients <- unlist(coefficients[c("tau", "beta", "kappa", "theta", "mean_E1")])
    ret$se <- unlist(se[c("tau", "beta", "kappa", "theta", "mean_E1")])
  }
  if(family == "Hermite") {
    coefficients$psi <- as.numeric(exp(ret$coefficients_raw["logit_psi"]) / (1 + exp(ret$coefficients_raw["logit_psi"])))
    se$psi <- as.numeric(ret$se_raw["logit_psi"] * exp(ret$coefficients_raw["logit_psi"]) /
                           (1 + exp(ret$coefficients_raw["logit_psi"]))^2)

    ret$fitted_values <- do.call(
      llik_hingarch,
      c(coefficients, vect = list(observed), return_fitted = TRUE)
    )$fitted
    get_theta <- get_cluster_size(coefficients$psi, se$psi, family = "Hermite")
    coefficients$theta <- get_theta$theta
    se$theta <- get_theta$theta_se

    ret$fitted_variance <- ret$fitted_values * (1 + coefficients$psi)
    ret$lik_distr <-  t(
      sapply(
        ret$fitted_values, FUN = function (x) {dherm(0:round(1.2 * max(observed)), mu = x, psi = coefficients$psi)}
      )
    )
    ret$psi <- c(coefficients["psi"], se_psi = unname(se["psi"]))
    ret$coefficients <- unlist(coefficients[c("tau", "beta", "kappa", "theta", "mean_E1")])
    ret$se <- unlist(se[c("tau", "beta", "kappa", "theta", "mean_E1")])
    }

  # compute residuals
  ret$pearson_residuals <- (observed - ret$fitted_values) / sqrt(ret$fitted_variance)

  # other:
  ret$observed <- observed
  ret$dim <- length(ret$coefficients)
  ret$loglikelihood <- -opt$value
  ret$AIC <- 2*(-ret$loglikelihood + ret$dim)
  ret$convergence <- opt$convergence
  ret$nobs <- length(observed)
  ret$optim <- opt
  ret$fitting_method <- "maximum_likelihood"
  ret$order <- c(p = 1, q = 1)
  class(ret) <- "inarma"

  return(ret)
}

#' Fitting an INARCH model
#'
#' Maximum likelihood inference for INARCH(1) model, which is obtained as a
#' special case of the INGARCH(1, 1) model as described in Bracher and Sobolova
#' (2024) for \eqn{\beta = 0}. For the definition of the INGARCH(1, 1) model,
#' see `?fit_ingarch`.
#'
#' @export
#'
#' @examples
#' data("measles")
#' X <- measles$value
#' \dontrun{
#' fit <- fit_inarch(X, family = "Poisson")
#' }
#'
#'
#' @param observed a vector containing the observed count values
#' @param family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`
#' @param start initial values for the optimization routine (on the internal
#' scale; check the element `"optim"` of the return list).
#' @param return_se should standard errors be returned?
#' @param control_optim a list of options passed to `optim`
#' @return a list with the following elements:
#' \describe{
#' \item{family}{the distribution family used.}
#' \item{coefficients_raw}{the estimated model coefficients on the internal scale.}
#' \item{se_raw}{if `return_se == TRUE`: the estimated standard errors on the internal scale.}
#' \item{cov_raw}{covariance matrix of the estimates on the internal scale.}
#' \item{coefficients}{the estimated model parameters transformed back to the natural scale.}
#' \item{se}{the estimated standard errors transformed back to the natural scale.}
#' \item{observed}{the vector of observed values provided by the user.}
#' \item{fitted_values}{the fitted values as obtained from `lik_distr`}
#' \item{pearson_residuals}{the Pearson residuals}
#' \item{psi}{estimate and standard error of the parameter \eqn{\psi} from the
#'  GLM-formulation of the INARCH model. Included for making the transitions between the two formulations easier.}
#' \item{dim}{the number of fitted parameters}
#' \item{loglikelihood}{the log-likelihood of the fitted model}
#' \item{AIC}{the resulting AIC}
#' \item{convergence}{indicates whether optimization converged}
#' \item{nobs}{the number of observations}
#' \item{optim}{return object of the call to `optim`}
#' \item{fitting_method}{the method used to fit the model, here `"maximum_likelihood"`.}
#' }
fit_inarch <- function(observed, family = c("Poisson", "Hermite", "NegBin"),
                        start = NULL, return_se = TRUE, control_optim = NULL){

  # Check the starting values
  if(is.null(start)) {
    start_val <- c(log_tau = 1, logit_kappa = 0, log_mean_E1 = 0.5)
  } else{
    names_ok <- sum((names(start) %in% c("log_tau", "logit_kappa",
                                         "log_mean_E1"))) == 3
    names_ok <- names_ok &&
      (family == "Hermite" && any(names(start) %in% "logit_psi") ||
         family == "NegBin" && any(names(start) %in% "log_psi"))
    if (names_ok) {
      start_val <- start
    } else{
      warning("Starting values must be supplied on the internal scale. Check the
              element 'optim' of the return list. Starting with the default
              starting values.")
      start_val <- c(log_tau = 1, logit_kappa = 0, log_mean_E1 = 0.5)
    }
  }

  # negative log-likelihood as function of parameter vector,
  # depends on the family:
  if(family == "Poisson"){
    nllik <- function (pars) {
      -llik_ingarch(observed, tau = exp(pars["log_tau"]),
                    beta = 0,
                    kappa = exp(pars["logit_kappa"]) / (1 + exp(pars["logit_kappa"])),
                    mean_E1 = exp(pars["log_mean_E1"]))
    }
  }

  if(family == "Hermite"){
    if (is.null(start) || !names_ok) {
      start_val <- c(start_val, logit_psi = -0.5)
    }
    nllik <- function (pars) {
      -llik_hingarch(observed, tau = exp(pars["log_tau"]),
                     beta = 0,
                     kappa = exp(pars["logit_kappa"]) / (1 + exp(pars["logit_kappa"])),
                     psi = exp(pars["logit_psi"]) / (1 + exp(pars["logit_psi"])),
                     mean_E1 = exp(pars["log_mean_E1"]))
    }
  }

  if(family == "NegBin"){
    if (is.null(start) || !names_ok) {
      start_val <- c(start_val, log_psi = -0.5)
    }
    nllik <- function (pars) {
      -llik_nbingarch(observed, tau = exp(pars["log_tau"]),
                      beta = 0,
                      kappa = exp(pars["logit_kappa"]) / (1 + exp(pars["logit_kappa"])),
                      psi = exp(pars["log_psi"]),
                      mean_E1 = exp(pars["log_mean_E1"]))
    }
  }

  # Optimize
  opt <- optim(start_val, nllik, control = control_optim, hessian = return_se)

  # structure results
  ret <- list()
  ret$family <- family
  ret$offspring <- "Poisson"

  # parameter estimates and standard errors:
  ret$coefficients_raw <- opt$par
  ret$se_raw <- ret$cov_raw <- NULL
  if(return_se){
    to_solve <- opt$hessian + diag(10^-6, dim(opt$hessian)[1]) # add small value to diagonal to avoid numerical issues
    ret$cov_raw <- solve(to_solve)
    ret$se_raw <- sqrt(diag(ret$cov_raw))
    if(any(diag(ret$cov_raw) < 0)){
      warning("Negative diagonal elements in inverse Fisher matrix - thresholding at zero. At least one estimated standard error will be zero.")
    }
    ret$se_raw <- sqrt(pmax(diag(ret$cov_raw), 0))
  }

  # compute standard errors on untransformed scale via delta method
  # (use as.numeric to get rid of naming of numeric vectors)
  se <- list(
    tau = as.numeric(ret$se_raw["log_tau"] * exp(ret$coefficients_raw["log_tau"])^2),
    kappa = as.numeric(ret$se_raw["logit_kappa"] * exp(ret$coefficients_raw["logit_kappa"])
                       / (1 + exp(ret$coefficients_raw["logit_kappa"]))^2),
    mean_E1 = as.numeric(ret$se_raw["log_mean_E1"] * exp(ret$coefficients_raw["log_mean_E1"])^2)
  )

  coefficients <- list(
    tau = as.numeric(exp(opt$par["log_tau"])),
    kappa = as.numeric(exp(opt$par["logit_kappa"]) / (1 + exp(opt$par["logit_kappa"]))),
    mean_E1 = as.numeric(exp(ret$coefficients_raw["log_mean_E1"]))
  )

  # Add psi into the list of results and compute fitted values, depends on the family:
  if (family =="Poisson") {
    ret$fitted_values <- do.call(
      llik_ingarch,
      c(coefficients, beta = 0, vect = list(observed), return_fitted = TRUE)
    )$fitted

    ret$lik_distr <- t(sapply(ret$fitted_values, FUN = dpois, x = 0:round(1.2 * max(observed))))
    ret$fitted_variance <- ret$fitted_values
    ret$psi <- NA
    ret$coefficients <- unlist(coefficients)
    ret$se <- unlist(se)
  }
  if(family == "NegBin") {
    coefficients$psi <- as.numeric(exp(ret$coefficients_raw["log_psi"]))
    se$psi <- as.numeric(ret$se_raw["log_psi"] * exp(ret$coefficients_raw["log_psi"])^2)
    ret$fitted_values <- do.call(
      llik_nbingarch,
      c(coefficients, beta = 0, vect = list(observed), return_fitted = TRUE)
    )$fitted

    ret$lik_distr <- dnbinom(0:1.2 * max(observed), mu = ret$fitted_values,
                             size = 1 / coefficients$psi)
    ret$fitted_variance <- ret$fitted_values * (1 + coefficients$psi)

    get_theta <- get_cluster_size(coefficients$psi, se$psi, family = "NegBin")
    coefficients$theta <- get_theta$theta
    se$theta <- get_theta$theta_se
    ret$psi <- c(coefficients["psi"], se_psi = unname(se["psi"]))
    ret$coefficients <- unlist(coefficients[c("tau", "kappa", "theta", "mean_E1")])
    ret$se <- unlist(se[c("tau", "kappa", "theta", "mean_E1")])
  }
  if(family == "Hermite") {
    coefficients$psi <- as.numeric(exp(ret$coefficients_raw["logit_psi"]) / (1 + exp(ret$coefficients_raw["logit_psi"])))
    se$psi <- as.numeric(ret$se_raw["logit_psi"] * exp(ret$coefficients_raw["logit_psi"]) /
                           (1 + exp(ret$coefficients_raw["logit_psi"]))^2)

    ret$fitted_values <- do.call(
      llik_hingarch,
      c(coefficients, beta = 0, vect = list(observed), return_fitted = TRUE)
    )$fitted
    get_theta <- get_cluster_size(coefficients$psi, se$psi, family = "Hermite")
    coefficients$theta <- get_theta$theta
    se$theta <- get_theta$theta_se

    ret$lik_distr <- dherm(0:1.2 * max(observed), mu = ret$fitted_values,
                           psi = coefficients$psi)
    ret$fitted_variance <- ret$fitted_values * (1 + coefficients$psi)
    ret$psi <- c(coefficients["psi"], se_psi = unname(se["psi"]))
    ret$coefficients <- unlist(coefficients[c("tau",  "kappa", "theta", "mean_E1")])
    ret$se <- unlist(se[c("tau", "kappa", "theta", "mean_E1")])
  }

  # compute residuals
  ret$pearson_residuals <- (observed - ret$fitted_values) / sqrt(ret$fitted_variance)

  # other:
  ret$observed <- observed
  ret$dim <- length(ret$coefficients)
  ret$loglikelihood <- -opt$value
  ret$AIC <- 2*(-ret$loglikelihood + ret$dim)
  ret$convergence <- opt$convergence
  ret$nobs <- length(observed)
  ret$optim <- opt
  ret$fitting_method <- "maximum_likelihood"
  ret$order <- c(p = 1, q = 0)
  class(ret) <- "inarma"

  return(ret)
}

#' Compute the clustersize and its standard error from the  model parameter psi
#' and the distribution family
#' @param psi the estimate of the parameter \eqn{psi} of the compound Poisson
#' distribution in the INGARCH framework
#' @param psi_se the standard error of \eqn{psi}
#' @param family family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`
#' @return a list containing the estimate of the parameter `theta` and its
#' standard error `theta_se`
get_cluster_size <- function(psi, psi_se, family = c("Poisson", "Hermite", "NegBin")){
  # computing cluster size is tedious:
  if(family == "Poisson"){
    clustersize <- 1
    clustersize_se <- 0
  }
  if(family == "Hermite"){
    pi <- psi / (2 - psi)
    clustersize <- 1 + pi
    clustersize_se <- psi_se * 2 / (2 - psi)^2
  }
  if(family == "NegBin"){
    pi <- 1/(1 + psi)
    clustersize <- (1 - pi)/(-pi*log(pi))
    clustersize_se <- psi_se * (log(1 + psi) - psi / (1 + psi)) / (log(1 + psi)^2)
  }
  return(list(theta = clustersize, theta_se = clustersize_se))
}
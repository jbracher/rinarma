#' Fitting an INARMA model
#'
#' Maximum likelihood inference for INARMA(1,1) model as described in Bracher and Sobolova (2024).
#' The model is defined as
#' \deqn{X_t = (1 - \beta) \circ E_t + I_t}
#' \deqn{E_t = \beta \circ E_{t - 1} + \kappa \circ X_{t - 1}}
#' where \eqn{\circ} denotes binomial thinning. The two thinnings of \eqn{E_t} are coupled via
#' \deqn{[\beta \circ E_t, (1 - \beta) \circ E_t] \sim \text{Mult}(E_t, \beta, 1 - \beta).}
#' The immigration process \eqn{I_t} consists of independently and identically distributed random variables which
#' can be Poisson, Hermite or negative binomial. This distribution is charaterized by its mean, denoted by
#'  \eqn{\tau}, and potentially a dispersion parameter \eqn{\psi}. See Bracher and Sobolova (2024) for details.
#'
#' The function uses a variation of the forward algorithm to fit the models. This is
#' computationally expensive, meaning that the optimization can be quite lengthy
#' already for medium-level counts. If the time series contains values above ~15
#' moment-based estimation using `fit_inarma_moments` may be more practical.
#'
#' @export
#'
#' @examples
#' data("measles")
#' X <- measles$value
#' # Note: running the fit takes a little while.
#' \dontrun{
#' fit <- fit_inarma(X, family = "Poisson")
#' summary(fit)
#' plot(fit, type = "fit")
#' }
#'
#'
#' @param observed a vector of observed count values
#' @param family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`
#' @param start initial values for the optimization routine (on the internal scale; check the element `"optim"` of the return list)
#' @param return_se should standard errors be returned?
#' @param parameterization the function internally works with a slightly different notation as
#' the paper by Bracher and Sobolova (2024). The default `parameterization = "beta"` will return results
#' formated as in the paper. `parameterization = "phi"` returns results as handled internally,
#' with `phi = 1 - beta`.
#' @param control_optim a list of options passed to `optim`
#' @return an object of class `inarma`. This is a list with the following elements:
#' \describe{
#' \item{family}{the distribution family used.}
#' \item{coefficients_raw}{the estimated model coefficients on the internal scale.}
#' \item{se_raw}{if `return_se == TRUE`: the estimated standard errors on the internal scale.}
#' \item{cov_raw}{covariance matrix of the estimates on the internal scale.}
#' \item{coefficients}{the estimated model parameters transformed back to the natural scale.}
#' \item{se}{the estimated standard errors transformed back to the natural scale.}
#' \item{observed}{the vector of observed values provided by the user.}
#' \item{lik_distr}{matrix containing for each time point the conditional probabilities for all
#' values in the support (given the past). This contains notably all likelihood contributions and
#' can serve to compute fitted values.}
#' \item{fitted_values}{the fitted values as obtained from `lik_distr`}
#' \item{fitted_variance}{the fitted conditional variances as obtained from `lik_distr`}
#' \item{pearson_residuals}{the Pearson residuals}
#' \item{dim}{the number of fitted parameters}
#' \item{loglikelihood}{the log-likelihood of the fitted model}
#' \item{AIC}{the resulting AIC}
#' \item{convergence}{indicates whether optimization converged}
#' \item{nobs}{the number of observations}
#' \item{optim}{return object of the call to `optim`}
#' \item{fitting_method}{the method used to fit the model, here `"maximum_likelihood"`.}
#' }
fit_inarma <- function(observed, family = c("Poisson", "Hermite", "NegBin"),
                       parameterization = "beta",
                       # data_tau = matrix(1, nrow = length(observed), # currently deprecated
                       #                 dimnames = list(NULL, "Intercept")),
                       start = NULL,
                       return_se = TRUE, control_optim = NULL){

  data_tau = matrix(1, nrow = length(observed), # currently deprecated
                    dimnames = list(NULL, "Intercept")) # to be removed when data is allowed again.

  # adapt colnames of data matrix:
  if(is.null(colnames(data_tau))){
    colnames(data_tau) <- c("Intercept", letters[seq_along(data_tau[1, -1])])
  }
  colnames(data_tau) <- paste0("tau.", colnames(data_tau))

  # create starting value for parameter optimization:
  tau_is_time_varying <- ncol(data_tau) > 1 | any(data_tau[, 1] != mean(data_tau))
  if(is.null(start)){
    # run moment estimation:
    suppressWarnings(
      suppressMessages(
        fit_moments <- fit_inarma_moments(observed = observed, family = family, parameterization = "phi")
      )
    )
    coefficients_moments <- fit_moments$coefficients
    # move coefficients into allowed ranges if necessary:
    coefficients_moments["tau"] <- max(0.1, coefficients_moments["tau"])
    coefficients_moments["phi"] <- min(max(0.05, coefficients_moments["phi"]), 0.95)
    coefficients_moments["kappa"] <- min(max(0.05, coefficients_moments["kappa"]), 0.95)
    if(family == "Hermite"){
      coefficients_moments["psi"] <- min(max(0.1, coefficients_moments["psi"]), 0.95)
    }
    if(family == "NegBin"){
      coefficients_moments["psi"] <- min(max(0.1, coefficients_moments["psi"]), 0.95)
    }

    # initialize vector of transformed moment estimators (scale used internally)
    start <- c("tau.Intercept" = NA,
               "logit_phi" = NA,
               "logit_kappa" = NA,
               if(family == "Hermite") "logit_psi" = NA,
               if(family == "NegBin") "log_psi" = NA,
               "log_mean_E1" = NA)

    # fill that vector:
    start["tau.Intercept"] <- log(coefficients_moments["tau"])
    start["logit_phi"] <- log(coefficients_moments["phi"]/(1 - coefficients_moments["phi"]))
    start["logit_kappa"] <- log(coefficients_moments["kappa"]/(1 - coefficients_moments["kappa"]))
    if(family == "Hermite"){
      start["logit_psi"] <- log(coefficients_moments["psi"]/(1 - coefficients_moments["psi"]))
    }
    if(family == "NegBin"){
      start["log_psi"] <- log(coefficients_moments["psi"])
    }
    start["log_mean_E1"] <- log(max(1, (observed[1] - coefficients_moments["tau"])/coefficients_moments["phi"]))


    # old starting values
    # start <- c(rep(0, ncol(data_tau)),
    #            -1, # phi
    #            -1, # kappa
    #            if(family %in% c("Hermite", "NegBin")) -1, # psi
    #            1) # mean_E1
    #
    #     names(start) <- c(colnames(data_tau),
    #                       "logit_phi", "logit_kappa",
    #                       if(family == "NegBin") "log_psi",
    #                       if(family == "Hermite") "logit_psi",
    #                       "log_mean_E1")

  }

  nllik_vect <- function(pars, return_distr){
    # print(".")
    lgt <- length(observed)

    beta_tau <- pars[colnames(data_tau)]

    tau <- if(tau_is_time_varying){
      exp(data_tau %*% beta_tau)
    }else{
      exp(beta_tau)
    }

    phi <- exp(pars["logit_phi"])/(1 + exp(pars["logit_phi"]))
    kappa <- exp(pars["logit_kappa"])/(1 + exp(pars["logit_kappa"]))

    if(family == "NegBin"){
      psi <- exp(pars["log_psi"])
    }

    if(family == "Hermite"){
      psi <- exp(pars["logit_psi"])/(1 + exp(pars["logit_psi"]))
    }

    # choose support (could likely be done more cleverly)
    support <- choose_support(observed = observed,
                              tau = tau, phi = phi,
                              kappa = kappa, psi = psi,
                              family = family)

    # initialization of E1 if necessary:
    if("log_mean_E1" %in% names(pars)){
      distr_E1 <- dpois(support, exp(pars["log_mean_E1"]))
    }else{
      distr_E1 <- NULL
    }

    if(family == "Poisson"){
      llik <- -llik_inarma_pois(vect = observed, distr_E1 = distr_E1, tau = tau,
                                phi = phi, kappa = kappa, support = support,
                                return_distr = return_distr)
      return(llik)
    }

    if(family == "Hermite"){
      llik <- -llik_inarma_herm(vect = observed, distr_E1 = distr_E1, tau = tau,
                                phi = phi, kappa = kappa, psi = psi, support = support,
                                return_distr = return_distr)
      return(llik)
    }

    if(family == "NegBin"){
      llik <- -llik_inarma_negbin(vect = observed, distr_E1 = distr_E1, tau = tau,
                                  phi = phi, kappa = kappa, psi = psi, support = support,
                                  return_distr = return_distr)
      return(llik)
    }
  }

  opt <- optim(par = start, fn = nllik_vect, return_distr = FALSE,
               hessian = return_se, control = control_optim)

  # very small overdispersion parameters indicate convergence issues.
  # try to catch these:
  if(family %in% c("Hermite", "NegBin")){
    estimate_psi <- ifelse(family == "Hermite",
                           exp(opt$par["logit_psi"])/(1 + exp(opt$par["logit_psi"])),
                           exp(opt$par["log_psi"]))

    if(estimate_psi < 0.05){
      message("Very low overdispersion parameter may indicate convergence problems, re-fitting model...")
      start_refit <- opt$par
      if(family == "Hermite") start_refit["logit_psi"] <- 1
      if(family == "NegBin") start_refit["log_psi"] <- 1

      opt <- optim(par = start_refit, fn = nllik_vect, return_distr = FALSE,
                   hessian = return_se, control = control_optim)
    }
  }
  # for(i in 1:3){
  #   opt <- optim(par = opt$par, fn = nllik_vect, control = control_optim,
  #                hessian = return_se)
  # }

  # create an informative return object inspired by hhh4:
  ret <- list()

  ret$family <- family
  # parameter estimates and standard errors on internal scale:
  ret$coefficients_raw <- opt$par
  ret$se_raw <- ret$cov_raw <- NULL
  if(return_se){
    to_solve <- opt$hessian + 10^-6 # add small value to diagonal to avoid numerical issues
    ret$cov_raw <- solve(to_solve)
    if(any(diag(ret$cov_raw) < 0)){
      warning("Negative diagonal elements in inverse Fisher matrix - thresholding at zero. At least one estimated standard error will be zero.")
    }
    ret$se_raw <- sqrt(pmax(diag(ret$cov_raw), 0))
  }

  # parameter estimates on original scale:
  # extract time-varying values of tau:
  beta_tau <- opt$par[grepl("tau.", names(opt$par))]
  tau <- if(tau_is_time_varying){
    exp(data_tau %*% beta_tau)
  }else{
    exp(beta_tau)
  }
  ret$coefficients <-
    list(tau = tau,
         phi = as.numeric(exp(ret$coefficients_raw["logit_phi"])/
                            (1 + exp(ret$coefficients_raw["logit_phi"]))),
         kappa = as.numeric(exp(ret$coefficients_raw["logit_kappa"])/
                              (1 + exp(ret$coefficients_raw["logit_kappa"]))))
  if(family == "NegBin") ret$coefficients$psi <- as.numeric(exp(ret$coefficients_raw["log_psi"]))
  if(family == "Hermite") ret$coefficients$psi <- as.numeric(exp(ret$coefficients_raw["logit_psi"])/(1 + exp(ret$coefficients_raw["logit_psi"])))
  ret$coefficients$mean_E1 <- as.numeric(exp(ret$coefficients_raw["log_mean_E1"]))
  # turn into vector:
  names_coefficients <- names(ret$coefficients)
  ret$coefficients <- unlist(ret$coefficients)
  names(ret$coefficients) <- names_coefficients

  # compute standard errors on untransformed scale via delta method
  # (use as.numeric to get rid of naming of numeric vectors)
  ret$se <- list(tau = as.numeric(ret$se_raw["tau.Intercept"]*exp(ret$coefficients_raw["tau.Intercept"])^2),
                 phi = as.numeric(ret$se_raw["logit_phi"]*exp(ret$coefficients_raw["logit_phi"])/(1 + exp(ret$coefficients_raw["logit_phi"]))^2),
                 kappa = as.numeric(ret$se_raw["logit_kappa"]*exp(ret$coefficients_raw["logit_kappa"])/(1 + exp(ret$coefficients_raw["logit_kappa"]))^2))
  if(family == "NegBin") ret$se$psi <- as.numeric(ret$se_raw["log_psi"]*exp(ret$coefficients_raw["log_psi"])^2)
  if(family == "Hermite") ret$se$psi <- as.numeric(ret$se_raw["logit_psi"]*exp(ret$coefficients_raw["logit_psi"])/(1 + exp(ret$coefficients_raw["logit_psi"]))^2)
  ret$se$mean_E1 <- as.numeric(ret$se_raw["log_mean_E1"]*exp(ret$coefficients_raw["log_mean_E1"])^2)
  ret$se <- unlist(ret$se)

  # adapt nomenclature to paper Bracher and Sobolova (2024)
  if(parameterization == "beta"){
    names(ret$coefficients)[names(ret$coefficients) == "phi"] <- "beta"
    ret$coefficients["beta"] <- 1 - ret$coefficients["beta"]
    names(ret$se)[names(ret$se) == "phi"] <- "beta"
  }

  ret$observed <- observed
  # to get fitted values:
  # obtain detailed likelihood, i.e. distributions for each t
  lik_distr <- -nllik_vect(pars = ret$coefficients_raw, return_distr = TRUE)
  # obtain expected values, variances and Pearson residuals from distribution:
  ret$fitted_values <- colSums((t(lik_distr)*(seq_along(lik_distr[1, ]) - 1)))
  ret$fitted_variance <- colSums((t(lik_distr)*(seq_along(lik_distr[1, ]) - 1)^2)) - ret$fitted_values^2
  ret$pearson_residuals <- (observed - ret$fitted_values)/sqrt(ret$fitted_variance)
  # and return full distribution as well
  ret$lik_distr <- lik_distr

  # other:
  ret$dim <- length(ret$coefficients)
  ret$loglikelihood <- -opt$value
  ret$AIC <- 2*(-ret$loglikelihood + ret$dim)
  ret$convergence <- (opt$convergence == 0)
  ret$nobs <- length(observed)
  ret$optim <- opt
  ret$fitting_method <- "maximum_likelihood"

  class(ret) <- "inarma"
  return(ret)
}

#' Choose the support for likelihood approximation depending on the observations and parameters
#' @param observed the time series of observed values
#' @param tau,phi,kappa,psi the model parameters
#' @param family the distributional family (`"Poisson"`, `"Hermite"` or `"NegBin`)
#' @return the support; a vector of consecutive integers starting at 0.
choose_support <- function(observed, tau, phi, kappa, psi = NULL, family){
  xi <- 1 - phi + phi*kappa

  mu_I <- tau
  sigma2_I <- switch(family,
                     "Poisson" = tau,
                     "Hermite" = (1 + psi)*tau,
                     "NegBin" = tau + psi*tau^2)

  mu_X <- mu_I/(1 - kappa)
  upper_X <- if(family == "Poisson"){
    qpois(0.999, lambda = mu_X)
  }else{
    sigma2_X <- sigma2_I/(1 - kappa) -
      (kappa*(1 + xi - phi*kappa)*(sigma2_I - mu_I))/
      ((1 - kappa)*(1 + xi))
    size_X <- mu_X^2/(sigma2_X - mu_X)
    qnbinom(0.999, mu = mu_X, size = size_X)
  }

  mu_E <- kappa*mu_I/(1 - xi)
  upper_E <- if(family == "Poisson"){
    qpois(0.9999, lambda = mu_E)
  }else{
    sigma2_E <- (kappa^2*sigma2_I + kappa*(1 - kappa + xi)*mu_I)/(1 - xi^2)
    size_E <- mu_E^2/(sigma2_E - mu_E)
    qnbinom(0.9999, mu = mu_E, size = size_E)
  }

  upper_observed <- ceiling(1.2*max(observed))

  return(0:max(upper_X, upper_E, upper_observed))
}


#' Fitting an INAR model
#'
#' Maximum likelihood inference for INAR(1) model. The model is defined as
#' \deqn{X_t = \kappa \circ X_{t - 1} + I_t}
#' where \eqn{\circ} denotes binomial thinning. The immigration process \eqn{I_t} consists of independently
#' and identically distributed random variables which can be Poisson, Hermite or negative binomial.
#' This distribution is charaterized by its mean, denoted by \eqn{\tau}, and potentially a dispersion parameter \eqn{\psi}.
#'
#' The function uses a variation of the forward algorithm to fit the models. It is a simplified version of
#' `fit_inarma` with the additional parameter \eqn{phi} set to 1.
#'
#'
#' @export
#'
#' @examples
#' data("measles")
#' X <- measles$value
#' # Note: running the fit takes a little while.
#' \dontrun{
#' fit <- fit_inar(X, family = "Poisson")
#' }
#'
#' @param observed a vector of observed count values.
#' @param family the distributional family; one of `"Poisson"`, `"Hermite"` or `"NegBin"`.
#' @param start initial values for the optimization routine.
#' @param return_se should standard errors be returned?
#' @param control_optim a list of options passed to `optim`.
#' @return an object of class `inarma`. This is a list with the following elements:
#' \describe{
#' \item{coefficients_raw}{the estimated model coefficients on the internal scale.}
#' \item{se_raw}{if `return_se == TRUE`: the estimated standard errors on the internal scale.}
#' \item{cov_raw}{covariance matrix of the estimates on the internal scale.}
#' \item{coefficients}{the estimated model parameters transformed back to the natural scale.}
#' \item{se}{the estimated standard errors transformed back to the natural scale.}
#' \item{observed}{the vector of observed values provided by the user.}
#' \item{lik_distr}{matrix containing for each time point the conditional probabilities for all
#' values in the support (given the past). This contains notably all likelihood contributions and
#' can serve to compute fitted values.}
#' \item{family}{the distribution family used.}
#' \item{fitted_values}{the fitted values as obtained from `lik_distr`.}
#' \item{fitted_variance}{the fitted conditional variances as obtained from `lik_distr`.}
#' \item{pearson_residuals}{the Pearson residuals.}
#' \item{dim}{the number of fitted parameters.}
#' \item{loglikelihood}{the log-likelihood of the fitted model.}
#' \item{AIC}{the resulting AIC.}
#' \item{convergence}{indicates whether optimization converged.}
#' \item{nobs}{the number of observations.}
#' \item{optim}{return object of the call to `optim`.}
#' \item{fitting_method}{the method used to fit the model, here `"maximum_likelihood"`.}
#' }
#' @export
fit_inar <- function(observed, family = c("Poisson", "Hermite", "NegBin"),
                     # data_tau = matrix(1, nrow = length(observed), # currently deprecated
                     #                 dimnames = list(NULL, "Intercept")),
                     start = NULL,
                     return_se = TRUE, control_optim = NULL){

  data_tau = matrix(1, nrow = length(observed), # currently deprecated
                    dimnames = list(NULL, "Intercept")) # to be removed when data is allowed again.

  # adapt colnames of data matrix:
  if(is.null(colnames(data_tau))){
    colnames(data_tau) <- c("Intercept", letters[seq_along(data_tau[1, -1])])
  }
  colnames(data_tau) <- paste0("tau.", colnames(data_tau))

  # create starting value for parameter optimization with reasonable names:
  tau_is_time_varying <- ncol(data_tau) > 1 | any(data_tau[, 1] != mean(data_tau))
  if(is.null(start)){
    start <- c(rep(0, ncol(data_tau)),
               -1, # kappa
               if(family %in% c("Hermite", "NegBin")) -1, # psi
               1) # mean_E1
    names(start) <- c(colnames(data_tau),
                      "logit_kappa",
                      if(family == "NegBin") "log_psi",
                      if(family == "Hermite") "logit_psi",
                      "log_mean_E1")
  }

  nllik_vect <- function(pars, return_distr = FALSE){
    lgt <- length(observed)

    beta_tau <- pars[colnames(data_tau)]

    tau <- if(tau_is_time_varying){
      exp(data_tau %*% beta_tau)
    }else{
      exp(beta_tau)
    }

    phi <- 1
    kappa <- exp(pars["logit_kappa"])/(1 + exp(pars["logit_kappa"]))

    if(family == "NegBin"){
      psi <- exp(pars["log_psi"])
    }

    if(family == "Hermite"){
      psi <- exp(pars["logit_psi"])/(1 + exp(pars["logit_psi"]))
    }

    # choose support could likely be done more cleverly)
    support <- choose_support(observed = observed,
                              tau = tau, phi = phi,
                              kappa = kappa, psi = psi,
                              family = family)

    # initialization of E1 if necessary:
    if("log_mean_E1" %in% names(pars)){
      distr_E1 <- dpois(support, exp(pars["log_mean_E1"]))
    }else{
      distr_E1 <- NULL
    }

    if(family == "Poisson"){
      llik <- -llik_inarma_pois(vect = observed, distr_E1 = distr_E1, tau = tau,
                                phi = phi, kappa = kappa, support = support,
                                return_distr = return_distr)
      return(llik)
    }

    if(family == "Hermite"){
      llik <- -llik_inarma_herm(vect = observed, distr_E1 = distr_E1, tau = tau,
                                phi = phi, kappa = kappa, psi = psi, support = support,
                                return_distr = return_distr)
      return(llik)
    }

    if(family == "NegBin"){
      llik <- -llik_inarma_negbin(vect = observed, distr_E1 = distr_E1, tau = tau,
                                  phi = phi, kappa = kappa, psi = psi, support = support,
                                  return_distr = return_distr)
      return(llik)
    }
  }

  opt <- optim(par = start, fn = nllik_vect, hessian = return_se,
               control = control_optim)

  # very small overdispersion parameters indicate convergence issues.
  # try to catch these:
  if(family %in% c("Hermite", "NegBin")){
    estimate_psi <- ifelse(family == "Hermite",
                           exp(opt$par["logit_psi"])/(1 + exp(opt$par["logit_psi"])),
                           exp(opt$par["log_psi"]))

    if(estimate_psi < 0.05){
      message("Very low overdispersion parameter may indicate convergence problems, re-fitting model...")
      start_refit <- opt$par
      if(family == "Hermite") start_refit["logit_psi"] <- 1
      if(family == "NegBin") start_refit["log_psi"] <- 1

      opt <- optim(par = start_refit, fn = nllik_vect, hessian = return_se,
                   control = control_optim)
    }
  }

  ret <- list()

  ret$family <- family

  # parameter estimates and standard errors:
  ret$coefficients_raw <- opt$par
  ret$se_raw <- ret$cov_raw <- NULL
  if(return_se){
    to_solve <- opt$hessian + 10^-6 # add small value to diagonal to avoid numerical issues
    ret$cov_raw <- solve(to_solve)
    ret$se_raw <- sqrt(diag(ret$cov_raw))
  }

  # parameter estimates on original scale:
  # extract time-varying values of tau:
  beta_tau <- opt$par[grepl("tau.", names(opt$par))]
  tau <- if(tau_is_time_varying){
    exp(data_tau %*% beta_tau)
  }else{
    exp(beta_tau)
  }
  ret$coefficients <-
    list(tau = tau,
         kappa = exp(ret$coefficients_raw["logit_kappa"])/
           (1 + exp(ret$coefficients_raw["logit_kappa"])))
  if(family == "NegBin") ret$coefficients$psi <- exp(ret$coefficients_raw["log_psi"])
  if(family == "Hermite") ret$coefficients$psi <- exp(ret$coefficients_raw["logit_psi"])/(1 + exp(ret$coefficients_raw["logit_psi"]))
  if(tau_is_time_varying) ret$coefficients$mean_E1 <- exp(ret$coefficients_raw["log_mean_E1"])
  # turn into vector:
  names_coefficients <- names(ret$coefficients)
  ret$coefficients <- unlist(ret$coefficients)
  names(ret$coefficients) <- names_coefficients

  ret$observed <- observed
  # to get fitted values:
  # obtain detailed likelihood, i.e. distributions for each t
  lik_distr <- -nllik_vect(pars = ret$coefficients_raw, return_distr = TRUE)
  # obtain expected value from distribution:
  ret$fitted_values <- colSums((t(lik_distr)*(seq_along(lik_distr[1, ]) - 1)))
  ret$fitted_variance <- colSums((t(lik_distr)*(seq_along(lik_distr[1, ]) - 1)^2)) - ret$fitted_values^2
  ret$pearson_residuals <- (observed - ret$fitted_values)/sqrt(ret$fitted_variance)
  # also return entire matrix:
  ret$lik_distr <- lik_distr

  # other:
  ret$dim <- length(ret$coefficients)
  ret$loglikelihood <- -opt$value
  ret$AIC <- 2*(-ret$loglikelihood + ret$dim)
  ret$convergence <- (opt$convergence == 0)
  ret$nobs <- length(observed)
  ret$optim <- opt
  ret$fitting_method <- "maximum_likelihood"

  class(ret) <- "inarma"
  return(ret)
}
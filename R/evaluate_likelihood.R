#' Get transition matrix from A[t] = (E[t] - L[t]) to E[t + 1]
#' @param kappa the parameter kappa (offspring mean)
#' @param X_tminus1 vector of probabilities for X[t - 1]; also implies the support
#' @param support_A the chosen support for A
#' @param support_A the chosen support for E
#' @return a matrix containing transition probabilities (values of A in rows, of E in columns)
tm_A_to_E <- function(kappa, X_tminus1, support_A, support_E){
  # the required matrix actually looks the same in all rows, just shifted
  # the vector of probabilities which needs to be re-used in each row
  vect0 <- dbinom(support_E,
                  size = X_tminus1,
                  prob = kappa)
  # fill the matrix with shifted versions of that vector
  tm_A_to_E <- matrix(rep(c(vect0, 0),
                          length.out = length(support_A)*length(support_E)),
                      nrow = length(support_A), byrow = TRUE)
  tm_A_to_E <- (1 - lower.tri(tm_A_to_E))*tm_A_to_E
  rownames(tm_A_to_E) <- paste0("A=", support_A, "X=", X_tminus1)
  colnames(tm_A_to_E) <- paste0("E=", support_E)

  return(tm_A_to_E)
}

#' Generic function for likelihood evaluation
#' @param vect the vector of observed values
#' @param distr_E1 a vector of probabilities used to initialize E1
#' @param distr_I the immigration distribution (probabilities for 0, ..., M); also implies the support for E and X
#' @param phi,kappa the model parameters
#' @param log should log-likelihood rather than likelihood be returned?
#' @param return_distr should a matrix with conditional probabilities (given the past) for each time point be returned?
#' @return the conditional (log) likelihood as a numeric or a matrix containing all relevant conditional probabilities (if return_distr == TRUE)
llik_inarma0_tv <- function(vect, distr_E1, distr_I, phi, kappa, log = TRUE, return_distr = FALSE){
  lgt <- length(vect)

  if(ncol(distr_I) != length(distr_E1)){
    stop("distr_E1 and distr_I need to be of the same length.")
  }

  # define supports:
  support <- 0:(ncol(distr_I) - 1)
  support_E <- support
  support_E_long <- rep(support_E, support_E + 1)
  support_L_long <- NULL
  for(s in support_E){
    support_L_long <- c(support_L_long, 0:s)
  }
  names(support_E_long) <- names(support_L_long) <- paste0("E", support_E_long,
                                                           "L", support_L_long)

  support_A <- support_E
  support_A_long <- support_E_long - support_L_long

  # helper functions to replace base:aggregate:
  # to get indices of entries corresponding to same values of A = E - L:
  get_inds_A <- function(i) which(support_A_long == i)
  # get these indices:
  inds_A <- sapply(support, FUN = get_inds_A)
  # function to get sums over the respective entries:
  sum_over_A <- function(j, p_EL_X_temp, inds_A) sum(p_EL_X_temp[inds_A[[j + 1]]])
  # helper function: compute p_X given the past from available quantities
  get_p_X_given_past <- function(x, distr_I_temp, p_EL_temp, support_L_long){
    p_I_shifted_temp <- 0*support_L_long
    p_I_shifted_temp[x - support_L_long >= 0] <-
      distr_I_temp[(x - support_L_long + 1)[x - support_L_long >= 0]]
    p_EL_X_unnorm_temp <- p_I_shifted_temp*p_EL_temp
    sum(p_EL_X_unnorm_temp)
  }

  if(return_distr){
    # initialize matrix
    p_X_given_past <- matrix(NA, nrow = lgt, ncol = length(support))
  }else{
    # or vector to store results:
    p_X_given_past <- numeric(lgt)
  }

  # initialize with stationary distribution of E1:
  p_E_temp <- distr_E1

  # move to stationary distribution of (E1, L1)
  p_EL_temp <- p_E_temp[support_E_long + 1]*
    dbinom(support_L_long, size = support_E_long, phi)

  # condition on X[1]
  # this is a bot tedious: shifted immigration distribution with 0 where
  # differences become negative
  p_I_shifted_temp <- 0*support_L_long
  p_I_shifted_temp[vect[1] - support_L_long >= 0] <-
    distr_I[1, (vect[1] - support_L_long + 1)[vect[1] - support_L_long >= 0]]
  p_EL_X_unnorm_temp <- p_I_shifted_temp*p_EL_temp
  p_EL_X_temp <- p_EL_X_unnorm_temp/(sum(p_EL_X_unnorm_temp))

  # if requested:
  if(return_distr){
    # store p(X[1]) for i in support:
    p_X_given_past[1, ] <- sapply(support, get_p_X_given_past,
                                  distr_I_temp = distr_I[1, ],
                                  p_EL_temp = p_EL_temp,
                                  support_L_long = support_L_long)
  }else{
    # otherwise just store store p(X[1])
    p_X_given_past[1] <- sum(p_EL_X_unnorm_temp)
  }

  # aggregate to p(E1 - L1)
  # p_A_X_temp <- aggregate(x = p_EL_X_temp, by = list(support_A_long), FUN = sum)[, "x"]
  p_A_X_temp <- sapply(support, sum_over_A,
                             p_EL_X_temp = p_EL_X_temp,
                             inds_A = inds_A)
  names(p_A_X_temp) <- paste0("A=", support_A)

  # now look over other time points:
  for(t in 2:lgt){
    # get transition matrix:
    tm_A_to_E_temp <- tm_A_to_E(kappa = kappa,
                                X_tminus1 = vect[t - 1],
                                support_A = support_A,
                                support_E = support_E)

    # get p(E[t] | X[< t])
    p_E_temp <- t(tm_A_to_E_temp)%*%p_A_X_temp

    # get p(E[t], L[t] | X[< t])
    p_EL_temp <- p_E_temp[support_E_long + 1]*
      dbinom(support_L_long, size = support_E_long, phi)

    # condition on X[t]
    # this is a bit tedious: shifted immigration distribution with 0 where
    # differences become negative
    p_I_shifted_temp <- 0*support_L_long
    p_I_shifted_temp[vect[t] - support_L_long >= 0] <-
      distr_I[t, (vect[t] - support_L_long + 1)[vect[t] - support_L_long >= 0]]
    p_EL_X_unnorm_temp <- p_I_shifted_temp*p_EL_temp
    p_EL_X_temp <- p_EL_X_unnorm_temp/(sum(p_EL_X_unnorm_temp))

    # store p(X[i]) for i in support (if requested):
    if(return_distr){
      p_X_given_past[t, ] <- sapply(support, get_p_X_given_past,
                                    distr_I_temp = distr_I[t, ],
                                    p_EL_temp = p_EL_temp,
                                    support_L_long = support_L_long)
    }else{
      # otherwise just store store p(X[t])
      p_X_given_past[t] <- sum(p_EL_X_unnorm_temp)
    }


    # aggregate to p(Et - Lt | X[< t])
    # p_A_X_temp[] <- aggregate(x = p_EL_X_temp, by = list(support_A_long), FUN = sum)[, "x"]
    p_A_X_temp <- sapply(support, sum_over_A,
                               p_EL_X_temp = p_EL_X_temp,
                               inds_A = inds_A)
  }

  if(return_distr){
    # if requested: return distributions as matrices
    return(p_X_given_past)
  }else{
    # otherwise return likelihood or log likelihood:
    if(log){
      return(sum(log(p_X_given_past)))
    }else{
      return(prod(p_X_given_past))
    }
  }
}

#' Checking the arguments of llik_inarma_pois, llik_inarma_herm, llik_inarma_nb
#' @param vect the vector of observed values
#' @param distr_E1 a vector of probabilities used to initialize E1
#' @param tau model parameter, scalar or a vector for time-varying immigration distributions
#' @param phi,kappa the model parameters, scalar
#' @param support the support for E and X
#' @param log should log-likelihood rather than likelihood be returned?
#' @param return_distr should a matrix with conditional probabilities (given the past) for each time point be returned?
#' @return the conditional (log) likelihood as a numeric or a matrix containing all relevant conditional probabilities (if return_distr == TRUE)
check_arguments_llik <- function(vect, distr_E1, tau, phi, psi, kappa, support, log, return_distr){
  if(is.null(distr_E1) & is.null(support)){
    stop("Either distr_E1 or support need to be provided.")
  }
  if(length(tau) > 1 & is.null(distr_E1)){
    stop("If length(tau) > 1 an initial distribution distr_E1 must be provided.")
  }
  if(!is.null(distr_E1) & !is.null(support)){
    if(length(distr_E1) != length(support)){
      stop("Conflicting distr_E1 and support.")
    }
  }
  if(length(vect) != length(tau) & length(tau) > 1){
    stop("vect and tau must be of the same length if length(tau) > 1.")
  }
}

#' Evaluating the likelihood in the Poisson case
#' A wrapper around llik_inarma0_tv
#' @param vect the vector of observed values
#' @param distr_E1 a vector of probabilities used to initialize E1
#' @param tau model parameter, scalar or a vector for time-varying immigration distributions
#' @param phi,kappa the model parameters, scalar
#' @param support the support for E and X
#' @param log should log-likelihood rather than likelihood be returned?
#' @param return_distr should a matrix with conditional probabilities (given the past) for each time point be returned?
#' @return the conditional (log) likelihood as a numeric or a matrix containing all relevant conditional probabilities (if return_distr == TRUE)
llik_inarma_pois <- function(vect, distr_E1 = NULL, tau, phi, kappa,
                             support = NULL, return_distr = FALSE){

  # check arguments:
  do.call(check_arguments_llik, as.list(environment()))

  if(is.null(support)) support <- seq_along(distr_E1) - 1

  xi <- 1 - phi + phi*kappa
  lgt <- length(vect)

  # use stationary distribution if nothing provided and scalar tau:
  if(length(tau) == 1 & is.null(distr_E1)){
    distr_E1 <- dpois(support, kappa*tau/(1 - xi))
  }

  # move to matrix with time-dependent distr_I:
  tau <- rep(tau, length.out = lgt)
  matr_distr_I <- t(sapply(tau, dpois, x = support))

  # do computations:
  llik_inarma0_tv(vect = vect, distr_E1 = distr_E1, distr_I = matr_distr_I,
                  phi = phi, kappa = kappa, return_distr = return_distr)
}



#' Evaluating the likelihood in the Hermite case
#' A wrapper around llik_inarma0_tv
#' @param vect the vector of observed values
#' @param distr_E1 a vector of probabilities used to initialize E1
#' @param tau model parameter, scalar or a vector for time-varying immigration distributions
#' @param phi,kappa,psi the model parameters, scalar
#' @param support the support for E and X
#' @param log should log-likelihood rather than likelihood be returned?
#' @param return_distr should a matrix with conditional probabilities (given the past) for each time point be returned?
#' @return the conditional (log) likelihood as a numeric or a matrix containing all relevant conditional probabilities (if return_distr == TRUE)
llik_inarma_herm <- function(vect, distr_E1 = NULL, tau, phi, psi, kappa,
                             support = NULL, return_distr = FALSE){

  # check arguments:
  do.call(check_arguments_llik, as.list(environment()))

  if(is.null(support)) support <- seq_along(distr_E1) - 1

  # move to parameterization a1, a2 (and bring to same length if necessary):
  pars <- to_a1a2_herm(mu = tau, psi = psi)
  a1 <- pars$a1
  a2 <- pars$a2

  xi <- 1 - phi + phi*kappa
  lgt <- length(vect)

  # use stationary distribution if nothing provided and scalar tau:
  if(length(tau) == 1 & is.null(distr_E1)){
    distr_E1 <- dherm(support,
                      a1 = (kappa*(1 + xi)*a1 + 2*kappa*(1 + xi - kappa)*a2)/
                        (1 - xi^2),
                      a2 = kappa^2*a2/(1 - xi^2))
  }

  # move to matrix with time-dependent distr_I:
  a1 <- rep(a1, length.out = lgt)
  a2 <- rep(a2, length.out = lgt)
  matr_distr_I <- matrix(nrow = lgt, ncol = length(support))
  for(i in 1:lgt){
    matr_distr_I[i, ] <- dherm(support, a1 = a1[i], a2 = a2[i])
  }

  # do computations:
  llik_inarma0_tv(vect = vect, distr_E1 = distr_E1, distr_I = matr_distr_I,
                  phi = phi, kappa = kappa, return_distr = return_distr)
}



#' Evaluating the likelihood in the negative binomial case
#' A wrapper around llik_inarma0_tv
#' @param vect the vector of observed values
#' @param distr_E1 a vector of probabilities used to initialize E1
#' @param tau model parameter, scalar or a vector for time-varying immigration distributions
#' @param phi,kappa,psi the model parameters, scalar
#' @param support the support for E and X
#' @param log should log-likelihood rather than likelihood be returned?
#' @param return_distr should a matrix with conditional probabilities (given the past) for each time point be returned?
llik_inarma_negbin <- function(vect, distr_E1 = NULL, tau, phi, psi, kappa,
                             support = NULL, return_distr = FALSE){

  # check arguments:
  do.call(check_arguments_llik, as.list(environment()))

  if(is.null(support)) support <- seq_along(distr_E1) - 1

  xi <- 1 - phi + phi*kappa
  lgt <- length(vect)

  # use stationary distribution if nothing provided and scalar tau:
  if(length(tau) == 1 & is.null(distr_E1)){
    sigma2_I <- tau + psi*tau^2
    mu_E <- kappa*tau/(1 - xi)
    sigma2_E <- (kappa^2*sigma2_I + kappa*(1 - kappa + xi)*tau)/(1 - xi^2)
    size_E <- mu_E^2/(sigma2_E - mu_E)
    distr_E1 <- dnbinom(support, size = size_E, mu = mu_E)
  }

  # move to matrix with time-dependent distr_I:
  tau <- rep(tau, length.out = lgt)
  # matr_distr_I <- t(sapply(tau, dnbinom, x = support, size = 1/psi, prob = NULL))
  matr_distr_I <- matrix(nrow = lgt, ncol = length(support))
  for(i in 1:lgt){
    matr_distr_I[i, ] <- dnbinom(support, mu = tau[i], size = 1/psi)
  }

  # do computations:
  llik_inarma0_tv(vect = vect, distr_E1 = distr_E1, distr_I = matr_distr_I,
                  phi = phi, kappa = kappa, return_distr = return_distr)
}
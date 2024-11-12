library(tidyverse)
library(pracma)
# library(inarma)
library(devtools)
load_all("~/Framework_INAR_INGARCH/Codes/Simulations_Baru/rinarma")
source("inst/AR_approximation_example/AR_approximation_functions.R")  # Load the functions

##############################################################################
## Inference for the INARMA(1, 1) model using the approximation by an AR model
##############################################################################

# Setup ------------------------------------------------------------------------

vals_tau <- c(1, 1, 1)
vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)

n_sim <- 1000 # number of simulation runs
lag_max <- 10  # The maximum lag of the approximating AR model
vals_lgt <- c(250, 500, 1000)  # lengths of simulated time series

# Run the loops ----------------------------------------------------------------

# Set a few different starting value in case the optimization does not converge
start_transformed <- matrix(c(1, 0, 0, 0.2, 1, 1, 0.2, -1, -1,
                              rep(NA, 3)), nrow = 4, ncol = 3, byrow = TRUE)
colnames(start_transformed) <- c("log_tau", "logit_kappa", "logit_beta")
for (s in 1:3) {  # Loop over the scenarios

  # Grab the parameter values
  tau <- vals_tau[s]
  beta <- vals_beta[s]
  kappa <- vals_kappa[s]

  for (lgt in vals_lgt) {  # Loop over the series lengths

    # Allocate the result containers
    res_11 <- matrix(NA, nrow = n_sim, ncol = 6)
    colnames(res_11) <- c("tau", "kappa", "beta", "tau_se", "kappa_se", "beta_se")
    conv_11 <- rep(NA, n_sim)
    sims <- matrix(NA, nrow = n_sim, ncol = lgt)

    for (k in 1:n_sim) {

      # Simulate an INARMA(1, 1) model
      set.seed(k)
      sim <- sim_inarma(tau = tau, kappa = kappa, beta = beta, lgt = lgt,
                        offspring = "binomial", family = "Poisson")

      # Calculate the moment estimates and use it as a starting value for the
      # optmization in case the previous 3 are bad
      start <- try(fit_inarma_moments(sim$X, family = "Poisson"))

      if (class(start) == "try-error") {
        start_transformed[4, ] <- c(log_tau = 0.5, logit_kappa = 0, logit_beta = 0,
                                    logit_psi = 0)
      } else {
        start <- start$coefficients
        start_transformed[4, ] <- c(
          log_tau = unname(log(start["tau"])),
          logit_kappa = unname(log(start["kappa"] / (1 - start["kappa"]))),
          logit_beta = unname(log(start["beta"] / (1 - start["beta"])))
        )
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
          X = sim$X,
          family = "Poisson",
          hessian = TRUE,
          control = list(fnscale = -1, maxit = 800)  # To change to maximization
        )

        # Extract and transform the results
        coeffs <- c(
          tau = unname(exp(op$par["log_tau"])),
          kappa = unname(exp(op$par["logit_kappa"]) / (1 + exp(op$par["logit_kappa"]))),
          beta = unname(exp(op$par["logit_beta"]) / (1 + exp(op$par["logit_beta"])))
        )
        ses <- get_ses_orig(op$par, op$hessian, family = "Poisson")
        refit <- anyNA(ses) | op$convergence != 0
        refit_iter <- refit_iter + 1
      }

      # Save the results
      res_11[k, 1:3] <- coeffs
      res_11[k, 4:6] <- ses
      conv_11[k] <- op$convergence
      sims[k, ] <- sim$X
    }

    print(paste("Finished length:", lgt))

    results_tab <- as_tibble(res_11) %>% mutate(convergence = conv_11, length = lgt)
    sim_tab <- as_tibble(sims)

    # Write the results
    write_csv(results_tab, file = paste0("inst/AR_approximation_example/Results/AR_INARMA11_s", s, "_lgt", lgt, ".csv"))
  }
}

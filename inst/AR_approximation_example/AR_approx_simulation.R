library(tidyverse)
library(pracma)
source("inst/AR_approximation_example/AR_approximation_functions.R")  # Load the functions

##############################################################################
## Inference for the INARMA(2, 1) model using the approximation by an AR model
##############################################################################

# Setup ------------------------------------------------------------------------

# Set the parameter values
# s <- 1
# kappa <- c(0.2, 0.6)
# beta <- 0.15
# tau <- 2

# Alternative set of parameters
s <- 2
kappa <- c(0.45, 0.25)
beta <- 0.35
tau <- 1.5

lag_max <- 10  # The maximum lag of the approximating AR model
vals_lgt <- c(250, 500, 1000)

# Run the loops ----------------------------------------------------------------

n_sim <- 1000  # Number of iterations

for (lgt in vals_lgt) {  # Loop over the series lengths

  # Allocate the result containers
  res_21 <- matrix(NA, nrow = n_sim, ncol = 4)
  colnames(res_21) <- c("tau", "kappa1", "kappa2", "beta")
  conv_21 <- rep(NA, n_sim)
  sims <- matrix(NA, nrow = n_sim, ncol = lgt)

  for (k in 1:n_sim) {

    # Simulate an INARMA(2, 1) model
    set.seed(k)
    sim <- sim_inarmapq(tau = tau, kappa = kappa, beta = beta, E0 = 1,
                        t_max = lgt)

    # Find the (approximate) maximum likelihood estimates
    op <- optim(
      par = c(
        log_tau = 1,
        transformed_kappa = orig_par_to_unconstrained(c(0.25, 0.25)),
        logit_beta = 0
      ),
      fn = llik_ar_based_higher,
      X = sim$X,
      control = list(fnscale = -1)  # To change to maximization
    )

    # Extract and transform the results
    coeffs <- c(
      tau = unname(exp(op$par["log_tau"])),
      kappa = unname(unconstrained_par_to_orig(op$par[2:3])),
      beta = unname(exp(op$par["logit_beta"]) / (1 + exp(op$par["logit_beta"])))
    )

    # Save the results
    res_21[k, ] <- coeffs
    conv_21[k] <- op$convergence
    sims[k, ] <- sim$X
  }

  print(paste("Finished length:", lgt))

  results_tab <- as_tibble(res_21) %>% mutate(convergence = conv_21, length = lgt)
  sim_tab <- as_tibble(sims)

  # Write the results
  write_csv(results_tab, file = paste0("inst/AR_approximation_example/Results/AR_INARMA21_s", s, "_lgt", lgt, ".csv"))
}

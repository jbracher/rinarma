library(tidyverse)

##############################################################################
## Inference for the INARMA(2, 1) model using the approximation by an AR model
##############################################################################

# Setup ------------------------------------------------------------------------

# Set the parameter values
s <- 1
kappa <- c(0.2, 0.6)
beta <- 0.15
tau <- 2

# Alternative set of parameters, must be uncommented manually to run
# s <- 2
# kappa <- c(0.45, 0.25)
# beta <- 0.35
# tau <- 1.5

lag_max <- 10  # The maximum lag of the approximating AR model
vals_lgt <- c(250, 500, 1000, 2000)

# Run the loops ----------------------------------------------------------------

n_sim <- 1000  # Number of iterations

for (lgt in vals_lgt) {  # Loop over the series lengths

  # Allocate the result containers
  res_21 <- matrix(NA, nrow = n_sim, ncol = 8)
  colnames(res_21) <- c("tau", "kappa1", "kappa2", "beta", "tau_se", "kappa1_se",
                        "kappa2_se", "beta_se")
  conv_21 <- rep(NA, n_sim)
  sims <- matrix(NA, nrow = n_sim, ncol = lgt)

  for (k in 1:n_sim) {

    # Simulate an INARMA(2, 1) model
    set.seed(k)
    sim <- sim_inarma(tau = tau, kappa = kappa, beta = beta, family = "Poisson",
                      offspring = "binomial", lgt = lgt)

    # Find the (approximate) maximum likelihood estimates
    fit <- fit_inarma_approx(sim$X, family = "Poisson", order = c(p = 2, q = 1),
                            lag_max = 10)

    # Save the results
    res_21[k, 1:4] <- fit$coefficients
    res_21[k, 5:8] <- fit$se
    conv_21[k] <- fit$convergence
    sims[k, ] <- sim$X
    if (k %% 50 == 0) print(paste("Finished iteration:", k))
  }

  print(paste("Finished length:", lgt))

  results_tab <- as_tibble(res_21) %>% mutate(convergence = conv_21, length = lgt)
  sim_tab <- as_tibble(sims)

  # Write the results
  write_csv(results_tab, file = paste0("inst/AR_approximation_example/Results/AR_INARMA21_s", s, "_lgt_", lgt, "_pois.csv"))
}

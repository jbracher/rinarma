library(tidyverse)

##############################################################################
## Inference for the negative binomial INARMA(2, 1) model using the
## approximation by an AR model
##############################################################################

# Setup ------------------------------------------------------------------------

# Set the parameter values
s <- 1
kappa <- c(0.2, 0.6)
beta <- 0.15
tau <- 2
psi <- 0.7

# Alternative set of parameters, must be uncommented manually to run
# s <- 2
# kappa <- c(0.45, 0.25)
# beta <- 0.35
# tau <- 1.5
# psi <- 0.5

lag_max <- 12  # The maximum lag of the approximating AR model
vals_lgt <- c(250, 500, 1000, 2000)

# Run the loops ----------------------------------------------------------------

n_sim <- 1000  # Number of iterations

for (lgt in vals_lgt) {  # Loop over the series lengths

  # Allocate the result containers
  res_21 <- matrix(NA, nrow = n_sim, ncol = 10)
  colnames(res_21) <- c("tau", "kappa1", "kappa2", "beta", "psi", "tau_se", "kappa1_se",
                        "kappa2_se", "beta_se", "psi_se")
  conv_21 <- rep(NA, n_sim)

  # Load the simulation trajectories
  dat <- unname(as.matrix(read_csv(
    file = paste0("inst/Simulation_trajectories/sim_INARMA21_s",
                  s, "_lgt_", lgt, "_nbin.csv"))))

  for (k in 1:n_sim) {

    # Find the (approximate) maximum likelihood estimates
    fit <- fit_inarma_approx(dat[k, ], family = "NegBin", order = c(p = 2, q = 1),
                             lag_max = lag_max)

    # Save the results
    res_21[k, 1:5] <- fit$coefficients
    res_21[k, 6:10] <- fit$se
    conv_21[k] <- fit$convergence

    if (k %% 50 == 0) print(paste0("Finished iteration, ", k, "."))
  }

  print(paste("Finished length:", lgt))

  results_tab <- as_tibble(res_21) %>% mutate(convergence = conv_21, length = lgt)

  # Write the results
  write_csv(results_tab, file = paste0("inst/AR_approximation_example/Results/AR_INARMA21_s", s, "_lgt_", lgt, "_nbin.csv"))
}

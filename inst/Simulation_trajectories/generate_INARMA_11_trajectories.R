# Script generating the trajectories of Poisson, Hermite and Negative binomial
# INARMA(1, 1) models from the simulation study
library(tidyverse)

# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)
vals_psi <- c(0.5, 0.7, 0.9)

vals_lgt <- c(250, 500, 1000)
vals_family <- c("Poisson", "Hermite", "NegBin")
vals_family_abbrev <- c("pois", "herm", "nbin")

n_sim <- 1e3  # Number of simulations

for (j in 1:3) {  # Loop over the families

  # Grab the full and abbreviated family name strings
  family <- vals_family[j]
  family_abbrev <- vals_family_abbrev[j]

  for (s in 1:3) {  # Loop over the scenarios

    # Grab the parameter values
    tau <- vals_tau[s]
    beta <- vals_beta[s]
    kappa <- vals_kappa[s]
    psi <- vals_psi[s]

    for (lgt in vals_lgt) {  # Loop over the lengths
      dat <- matrix(nrow = n_sim, ncol = lgt)

      for (k in 1:n_sim) {  # Loop over the iterations

        # Simulate an INARMA(1, 1) model
        set.seed(k)
        dat[k, ] <- sim_inarma(tau = tau, kappa = kappa, beta = beta, psi = psi,
                               family = family, offspring = "binomial",
                               lgt = lgt)$X

      }

      # Save in a file
      write_csv(
        as.data.frame(dat),
        file = paste0("inst/Simulation_trajectories/sim_INARMA11_s",
                      s, "_lgt_", lgt, "_", family_abbrev, ".csv")
      )

    }
  }
}
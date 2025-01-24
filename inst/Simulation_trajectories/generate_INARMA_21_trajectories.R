# Script generating the trajectories of Poisson, Hermite and Negative binomial
# INARMA(2, 1) models from the simulation study


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

vals_lgt <- c(250, 500, 1000, 2000)
vals_family <- c("Poisson", "Hermite", "NegBin")
vals_family_abbrev <- c("pois", "herm", "nbin")

n_sim <- 1e3  # Number of simulations

for (j in 1:3) {  # Loop over the families

  # Grab the full and abbreviated family name strings
  family <- vals_family[j]
  family_abbrev <- vals_family_abbrev[j]

  for (lgt in vals_lgt) {  # Loop over the lengths
    dat <- matrix(nrow = n_sim, ncol = lgt)

    for (k in 1:n_sim) {  # Loop over the iterations

      # Simulate an INARMA(2, 1) model
      set.seed(k)
      dat[k, ] <- sim_inarma(tau = tau, kappa = kappa, beta = beta, psi = psi,
                             family = family, offspring = "binomial",
                             lgt = lgt)$X

    }

    # Save in a file
    write_csv(
      as.data.frame(dat),
      file = paste0("inst/Simulation_trajectories/sim_INARMA21_s",
                    s, "_lgt_", lgt, "_", family_abbrev, ".csv")
      )

  }
}
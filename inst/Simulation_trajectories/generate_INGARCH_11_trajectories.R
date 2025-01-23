# Script generating the trajectories of Poisson, INGARCH(1, 1) model from the
# simulation study

# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)

vals_lgt <- c(250, 500, 1000)

n_sim <- 1e3  # Number of simulations


for (s in 1:3) {  # Loop over the scenarios

  # Grab the parameter values
  tau <- vals_tau[s]
  beta <- vals_beta[s]
  kappa <- vals_kappa[s]

  for (lgt in vals_lgt) {  # Loop over the lengths
    dat <- matrix(nrow = n_sim, ncol = lgt)

    for (k in 1:n_sim) {  # Loop over the iterations

      # Simulate an INGARCH(1, 1) model
      set.seed(k)
      dat[k, ] <- sim_inarma(tau = tau, kappa = kappa, beta = beta, zeta = 1,
                             family = "Poisson", offspring = "binomial-Poisson",
                             E1 = 2, lgt = lgt)$X

    }

    # Save in a file
    write_csv(
      as.data.frame(dat),
      file = paste0("inst/Simulation_trajectories/sim_INGARCH11_s",
                    s, "_lgt_", lgt, "_pois.csv")
    )

  }
}
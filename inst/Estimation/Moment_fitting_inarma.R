# Script fitting the INARMA model using the method of moments

vals_lgt <- c(250, 500, 1000) # lengths of simulated time series
vals_family <- c("Poisson", "Hermite", "NegBin")  # Families of the INARMA(1, 1) model
vals_family_abbrev <- c("pois", "herm", "nbin")  # Family abbreviations
names(vals_family_abbrev) <- vals_family

n_sim <- 1000 # number of simulation runs

for (family in vals_family) {  # Loop over the families
  for (s in 1:3) {  # Loop over the scenarios
    for (lgt in vals_lgt) {  # Loop over the lengths

      # Initialize the data frame to store the results from the INAR fitting
      res_moments <- data.frame(
        kappa = rep(NA, n_sim),
        beta = rep(NA, n_sim),
        tau = rep(NA, n_sim),
        psi = rep(NA, n_sim)
      )

      # run simulation:
      print(paste0("Started scenario ", s, " length ", lgt, " family ", vals_family_abbrev[family], "."))

      # Load the simulated trajectories
      dat <- unname(as.matrix(read.csv(
        paste0("inst/Simulation_trajectories/sim_INARMA11_s", s, "_lgt_", lgt,
               "_", vals_family_abbrev[family], ".csv")
      )))

      # Fit by the method of moments
      for (k in 1:n_sim) {

        fit <- fit_inarma_moments(dat[k, ], family = family)

        res_moments[k, c("tau", "beta", "kappa", "psi")] <- fit$coefficients[1:4]
      }

      # Store the results
      if (family == "Poisson") {
        res_moments <- res_moments[, c("tau", "beta", "kappa")]
      }
      write.csv(res_moments, file = paste0("inst/Estimation/Results/Estim_moments_s", s, "_", lgt, "_", vals_family_abbrev[family], ".csv"),
                row.names = FALSE)
      print(paste0("Finished scenario ", s, " length ", lgt, " family ", vals_family_abbrev[family], " check the results."))
    }
  }
}

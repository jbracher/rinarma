# Script fitting the Hermite INARMA model using the maximum likelihood

vals_lgt <- c(250, 500, 1000) # lengths of simulated time series
n_sim <- 1000 # number of simulation runs

# Keep track of scenarios and lengths that are done
# Needs to be changed manually, when selecting scenarios and lengths
progress_log <- expand.grid(scenario = 1:3, lgt = vals_lgt,
                            already_done = FALSE, continue = FALSE)

# If only a subset of the simulations is supposed to be run, change the
# `progress_log` variable manually now, before running the rest of the script

sims_to_run <- 1:9
sims_to_run <- sims_to_run[!(progress_log$already_done)]

# Set a few different starting value in case the optimization does not converge
start_transformed <- matrix(
  c(rep(NA, 5),
    1, 0, 2, 0, 0,
    0.2, -2, 0, 1, 1,
    0.2, -1, 3, -1, -1
  ), nrow = 4, ncol = 5, byrow = TRUE)
colnames(start_transformed) <- c("tau.Intercept", "logit_kappa", "logit_phi",
                                 "logit_psi", "log_mean_E1")

for (i in sims_to_run) {

  # should existing results be read (to take up after error)
  continue <- progress_log$continue[i]

  # grab the scenario number an the series length
  s <- progress_log$scenario[i]
  lgt <- progress_log$lgt[i]

  # run simulation:
  print(paste0("Started scenario ", s, " length ", lgt, "."))

  # Load the simulated trajectories
  dat <- unname(as.matrix(read.csv(
    paste0("inst/Simulation_trajectories/sim_INARMA11_s", s, "_lgt_", lgt, "_herm.csv")
  )))

  # initialize matrix to store results:
  # get existing results if desired (if simulation somehow failed halfway)
  if (continue) {
    res <- read.csv(
      paste0("inst/Estimation/Results/Estim_MLE_s", s, "_", lgt,
             "_herm.csv")
    )
    inds_to_run <- Position(res$kappa, f = is.na):n_sim
  } else {
    inds_to_run <- 1:n_sim
    res <- data.frame(
      kappa = rep(NA, n_sim),
      beta = rep(NA, n_sim),
      tau = rep(NA, n_sim),
      psi = rep(NA, n_sim),
      mean_E1 = rep(NA, n_sim),
      kappa_se = rep(NA, n_sim),
      beta_se = rep(NA, n_sim),
      tau_se = rep(NA, n_sim),
      psi_se = rep(NA, n_sim),
      mean_E1_se = rep(NA, n_sim),
      max_lik = rep(NA, n_sim),
      convergence = rep(NA, n_sim)
    )
  }

  # Simulate
  for (k in inds_to_run) {

    # Get the moment estimates
    mom_ests <- fit_inarma_moments(as.numeric(dat[k, ]), family = "Poisson")
    start_transformed[1, ] <- c(
      log(mom_ests$coefficients["tau"]),
      log(mom_ests$coefficients["kappa"] / (1 - mom_ests$coefficients["kappa"])),
      log((1 - mom_ests$coefficients["beta"]) / mom_ests$coefficients["beta"]),
      log(mom_ests$coefficients["psi"] / (1 - mom_ests$coefficients["psi"])),
      0.2
    )
    start_transformed[is.na(start_transformed) | is.infinite(start_transformed)] <- 0

    # Refit the model until we get a non-problematic estimate
    refit_iter <- 1
    refit <- TRUE
    while (refit & refit_iter <= 4) {

      # Fit INARMA
      fit <- fit_inarma(
        as.numeric(dat[k, ]),
        family = "Hermite",
        offspring = "binomial",
        return_se = TRUE,
        start = start_transformed[refit_iter, c("tau.Intercept", "logit_kappa",
                                                "logit_phi", "logit_psi", "log_mean_E1")]
      )

      # Increment the counters
      refit <- fit$optim$convergence != 0
      refit_iter <- refit_iter + 1
    }

    res$max_lik[k] <- fit$loglikelihood
    res[k, c("tau", "beta", "kappa", "psi", "mean_E1")] <- fit$coefficients[1:5]
    res[k, c("tau_se", "beta_se", "kappa_se", "psi_se", "mean_E1_se")] <- fit$se[1:5]
    res[k, "convergence"] <- fit$optim$convergence

    # Store the results every 10 iterations
    if (k %% 10 == 0) {
      print(paste0("Finished iteration: ", k, "."))
      write.csv(res, file = paste0("inst/Estimation/Results/Estim_MLE_s", s, "_", lgt, "_herm.csv"),
                row.names = FALSE)
    }
  }
  progress_log$already_done[i] <- TRUE
  print(paste0("Finished scenario ", s, " length ", lgt, " check the results."))
}

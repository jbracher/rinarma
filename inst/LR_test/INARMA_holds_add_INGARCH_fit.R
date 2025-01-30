# Script fitting the INGARCH model against the INARMA model, while
# INARMA holds, i.e.
# H0: zeta = 0 vs. H1: zeta > 0 and it happens that zeta = 1

vals_lgt <- c(250, 500, 1000)  # lengths of simulated time series
n_sim <- 1000  # number of simulation runs

# Keep track of scenarios and lengths that are done
# Needs to be changed manually, when selecting scenarios and lengths
progress_log <- expand.grid(scenario = 1:3, lgt = vals_lgt, already_done = FALSE,
                            continue = FALSE)

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
                                 "log_mean_E1", "logit_zeta")

inds_to_run <- 1:1000
for (i in sims_to_run) {

  # grab the scenario number an the series length
  s <- progress_log$scenario[i]
  lgt <- progress_log$lgt[i]

  # run simulation:
  print(paste0("Started scenario ", s, " length ", lgt, "."))

  # Load the simulated trajectories
  dat <- unname(as.matrix(read.csv(
    paste0("inst/Simulation_trajectories/sim_INARMA11_s", s, "_lgt_", lgt, "_pois.csv")
  )))

  # Load the estimates from the fitting of INARMA
  res <- read.csv(
    paste0("inst/LR_test/Results/INARMA_fits_INARMA_holds_s", s, "_", lgt,
           ".csv")
  )

  # Initialize the data frame to store the results from the INGARCH fitting
  res_ingarch <- data.frame(
    kappa = rep(NA, n_sim),
    beta = rep(NA, n_sim),
    tau = rep(NA, n_sim),
    kappa_se = rep(NA, n_sim),
    beta_se = rep(NA, n_sim),
    tau_se = rep(NA, n_sim),
    max_lik = rep(NA, n_sim),
    convergence = rep(NA, n_sim)
  )

  # Simulate
  for (k in inds_to_run) {

    # Fit the INGARCH model
    op_ingarch <- fit_ingarch(
      as.numeric(dat[k, ]),
      family = "Poisson",
      start = c(
        log_tau = log(res[k, "tau"]),
        logit_beta = unname(log(res[k, "beta"] / (1 - res[k, "beta"]))),
        logit_kappa = unname(log(res[k, "beta"] / (1 - res[k, "beta"]))),
        log_mean_E1 = -1
        ),
      return_se = TRUE,
      control_optim = list(maxit = 1000)
    )

    # Refit the INGARCH model until we get a non-problematic estimate
    refit_iter <- 1
    refit <- op_ingarch$optim$convergence != 0
    colnames(start_transformed) <- c("log_tau", "logit_kappa", "logit_beta",
                                     "log_mean_E1", "logit_zeta")
    while (refit & refit_iter <= 4) {

      # (Re)fit the INGARCH model
      op_ingarch <- fit_ingarch(
        as.numeric(dat[k, ]),
        family = "Poisson",
        start = start_transformed[refit_iter, 1:4],
        return_se = TRUE,
        control_optim = list(maxit = 1000)
      )

      # Increment the counters
      refit <- op_ingarch$optim$convergence != 0
      refit_iter <- refit_iter + 1
    }
  #  print(paste0("INGARCH refitted: ", refit_iter - 1, " times."))

    res_ingarch$max_lik[k] <- op_ingarch$loglikelihood
    res_ingarch[k, c("tau", "beta", "kappa")] <- op_ingarch$coefficients[1:3]
    res_ingarch[k, c("tau_se", "beta_se", "kappa_se")] <- op_ingarch$se[1:3]
    res_ingarch[k, "convergence"] <- op_ingarch$optim$convergence

    # Store the results every 10 iterations
    if (k %% 10 == 0) {
      print(paste0("Finished iteration: ", k, "."))
      write.csv(res_ingarch, file = paste0("inst/LR_test/Results/H0_INARMA_holds_results_INGARCH_fit_s", s, "_", lgt, ".csv"),
                row.names = FALSE)
    }
  }
  progress_log$already_done[i] <- TRUE
  print(paste0("Finished scenario ", s, " length ", lgt, " check the results."))
}

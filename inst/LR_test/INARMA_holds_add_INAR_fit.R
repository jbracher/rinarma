# Script fitting the INAR model against the INARMA model, while
# INARMA holds, i.e.
# H0: beta = 0 vs. H1: beta > 0

vals_lgt <- c(250, 500, 1000) # lengths of simulated time series
n_sim <- 1000 # number of simulation runs

# Keep track of scenarios and lengths that are done
# Needs to be changed manually, when selecting scenarios and lengths
progress_log <- expand.grid(scenario = 1:3, lgt = vals_lgt, already_done = FALSE,
                            continue = FALSE)

sims_to_run <- 1:9
sims_to_run <- sims_to_run[!(progress_log$already_done)]

# Set a few different starting value in case the optimization does not converge
start_transformed <- matrix(
  c(rep(NA, 4), 0.7,
    1, 0, 2, 0, 0.7,
    0.2, -2, 0, 1, 0,
    0.2, -1, 3, -1, 0
  ), nrow = 4, ncol = 4, byrow = TRUE)
colnames(start_transformed) <- c("tau.Intercept", "logit_kappa", "logit_phi",
                                 "log_mean_E1")

inds_to_run <- 1:1000
for (i in sims_to_run) {

  s <- progress_log$scenario[i]
  lgt <- progress_log$lgt[i]

  # run simulation:
  print(paste0("Started scenario ", s, " length ", lgt, "."))

  # Load the estimates from the fitting of INARMA
  res <- read.csv(
    paste0("inst/LR_test/Results/INARMA_fits_INARMA_holds_s", s, "_", lgt,
           ".csv")
  )

  # Load the simulated trajectories
  dat <- unname(as.matrix(read.csv(
    paste0("inst/Simulation_trajectories/sim_INARMA11_s", s, "_lgt_", lgt, "_pois.csv")
  )))

  # Initialize the data frame to store the results from the INAR fitting
  res_inar <- data.frame(
    kappa = rep(NA, n_sim),
    tau = rep(NA, n_sim),
    kappa_se = rep(NA, n_sim),
    tau_se = rep(NA, n_sim),
    max_lik = rep(NA, n_sim),
    convergence = rep(NA, n_sim)
  )

  # Simulate
  for (k in inds_to_run) {

    # Fit the INAR model
    op_inar <- fit_inar(
      as.numeric(dat[k, ]),
      family = "Poisson",
      offspring = "binomial",
      return_se = TRUE,
      start = c(
        log_tau = log(res[k, "tau"]),
        logit_kappa = unname(log(res[k, "beta"] / (1 - res[k, "beta"]))),
        log_mean_E1 = 0.7
      )
    )

    # Refit the INAR model until we get a non-problematic estimate
    refit_iter <- 1
    refit <- op_inar$optim$convergence != 0
    while (refit & refit_iter <= 4) {

      # (Re)fit the INAR model
      op_inar <- fit_inar(
        as.numeric(dat[k, ]),
        family = "Poisson",
        offspring = "binomial",
        return_se = TRUE,
        start = start_transformed[refit_iter, c("tau.Intercept", "logit_kappa")]
      )

      # Increment the counters
      refit <- op_inar$optim$convergence != 0
      refit_iter <- refit_iter + 1
    }

    res_inar$max_lik[k] <- op_inar$loglikelihood
    res_inar[k, c("tau", "kappa")] <- op_inar$coefficients[1:2]
    res_inar[k, c("tau_se", "kappa_se")] <- op_inar$se[1:2]
    res_inar[k, "convergence"] <- op_inar$optim$convergence

    # Store the results every 10 iterations
    if (k %% 10 == 0) {
      print(paste0("Finished iteration: ", k, "."))
      write.csv(res_inar, file = paste0("inst/LR_test/Results/INAR_fits_INARMA_holds_s", s, "_", lgt, ".csv"),
                row.names = FALSE)
    }
  }
  progress_log$already_done[i] <- TRUE
  print(paste0("Finished scenario ", s, " length ", lgt, " check the results."))
}

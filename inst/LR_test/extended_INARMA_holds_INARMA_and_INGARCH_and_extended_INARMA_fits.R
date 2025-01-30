# Script fitting the INARMA and INGARCH model against the extended INARMA model,
# while extended INARMA holds, i.e.
# H0: zeta = 0 vs. H1: zeta > 0 and
# H0: zeta = 1 vs. H1: zeta < 1

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
    paste0("inst/Simulation_trajectories/sim_extended_INARMA11_s", s, "_lgt_", lgt, "_pois.csv")
  )))

  # initialize matrix to store results:
  # get existing results if desired (if simulation somehow failed halfway)
  if (continue) {
    res <- read.csv(
      paste0("inst/LR_test/Results/H1_extended_INARMA_holds_results_s", s, "_", lgt,
             ".csv")
    )
    res_null <- read.csv(
      paste0("inst/LR_test/Results/H1_extended_INARMA_holds_results_null_s", s, "_", lgt,
             ".csv")
    )
    res_ingarch <- read.csv(
      paste0("inst/LR_test/Results/H1_extended_INARMA_holds_results_ingarch_s", s, "_", lgt,
             ".csv")
    )
    inds_to_run <- Position(res$kappa, f = is.na):n_sim
  } else {
    inds_to_run <- 1:n_sim
    res <- data.frame(
      kappa = rep(NA, n_sim),
      zeta = rep(NA, n_sim),
      beta = rep(NA, n_sim),
      tau = rep(NA, n_sim),
      mean_E1 = rep(NA, n_sim),
      kappa_se = rep(NA, n_sim),
      zeta_se = rep(NA, n_sim),
      beta_se = rep(NA, n_sim),
      tau_se = rep(NA, n_sim),
      mean_E1_se = rep(NA, n_sim),
      max_lik = rep(NA, n_sim),
      convergence = rep(NA, n_sim)
    )
    res_null <- data.frame(
      kappa = rep(NA, n_sim),
      beta = rep(NA, n_sim),
      tau = rep(NA, n_sim),
      mean_E1 = rep(NA, n_sim),
      kappa_se = rep(NA, n_sim),
      beta_se = rep(NA, n_sim),
      tau_se = rep(NA, n_sim),
      mean_E1_se = rep(NA, n_sim),
      max_lik = rep(NA, n_sim),
      convergence = rep(NA, n_sim)
    )
    res_ingarch <- data.frame(
      kappa = rep(NA, n_sim),
      beta = rep(NA, n_sim),
      tau = rep(NA, n_sim),
      mean_E1 = rep(NA, n_sim),
      kappa_se = rep(NA, n_sim),
      beta_se = rep(NA, n_sim),
      tau_se = rep(NA, n_sim),
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
      1.5,
      0.2
    )
    start_transformed[is.na(start_transformed) | is.infinite(start_transformed)] <- 0

    # Refit the null model until we get a non-problematic estimate
    refit_iter <- 1
    refit <- TRUE
    while (refit & refit_iter <= 4) {

      # Fit the INARMA model
      op_null <- fit_inarma(
        as.numeric(dat[k, ]),
        family = "Poisson",
        offspring = "binomial",
        return_se = TRUE,
        start = start_transformed[refit_iter, c("tau.Intercept", "logit_kappa",
                                                "logit_phi", "log_mean_E1")]
      )

      # Increment the counters
      refit <- op_null$optim$convergence != 0
      refit_iter <- refit_iter + 1
    }

  #  print(paste0("Null refitted: ", refit_iter - 2, " times."))

    # Fit the extended INARMA model
    op <- fit_inarma(as.numeric(dat[k, ]), family = "Poisson",
                     start = c(op_null$coefficients_raw, logit_zeta = 0, log_mean_E1 = 1),
                     offspring = "binomial-Poisson", return_se = TRUE,
                     control_optim = list(maxit = 1000))

    # Refit the alternative (big) model until we get a non-problematic estimate
    refit_iter <- 1
    refit <- op$optim$convergence != 0
    while (refit & refit_iter <= 4) {

      # (Re)fit the extended INARMA model
      op <- fit_inarma(as.numeric(dat[k, ]), family = "Poisson",
                       start = start_transformed[refit_iter,],
                       offspring = "binomial-Poisson", return_se = TRUE,
                       control_optim = list(maxit = 1000))

      # Increment the counters
      refit <- op$optim$convergence != 0
      refit_iter <- refit_iter + 1
    }
 #   print(paste0("Alternative refitted: ", refit_iter - 1, " times."))

    # Fit the INGARCH model
    op_ingarch <- fit_ingarch(
      as.numeric(dat[k, ]),
      family = "Poisson",
      start = c(
        log_tau = unname(log(op$coefficients["tau.Intercept"])),
        logit_beta = unname(log(op$coefficients["beta"] / (1 - op$coefficients["beta"]))),
        logit_kappa = unname(log(op$coefficients["kappa"] / (1 - op$coefficients["kappa"]))),
        log_mean_E1 = unname(log(op$coefficients["mean_E1"]))
      ),
      return_se = TRUE,
      control_optim = list(maxit = 1000)
    )

    # Refit the INGARCH model until we get a non-problematic estimate
    refit_iter <- 1
    refit <- op_ingarch$optim$convergence != 0
    while (refit & refit_iter <= 4) {

      # (Re)it the INGARCH model
      op_ingarch <- fit_ingarch(
        as.numeric(dat[k, ]),
        family = "Poisson",
        start = c(start_transformed[refit_iter, c(1:2, 4)],
                  logit_beta = unname(start_transformed[refit_iter, 3])),
        return_se = TRUE,
        control_optim = list(maxit = 1000)
      )

      # Increment the counters
      refit <- op_ingarch$optim$convergence != 0
      refit_iter <- refit_iter + 1
    }

    #  print(paste0("INGARCH refitted: ", refit_iter - 1, " times."))
    res$max_lik[k] <- op$loglikelihood
    res_null$max_lik[k] <- op_null$loglikelihood
    res_ingarch$max_lik[k] <- op_ingarch$loglikelihood
    res[k, c("tau", "beta", "kappa", "zeta", "mean_E1")] <- op$coefficients[1:5]
    res_null[k, c("tau", "beta", "kappa", "mean_E1")] <- op_null$coefficients[1:4]
    res_ingarch[k, c("tau", "beta", "kappa", "mean_E1")] <- op_ingarch$coefficients[1:4]
    res[k, c("tau_se", "beta_se", "kappa_se", "zeta_se", "mean_E1_se")] <- op$se[1:5]
    res_null[k, c("tau_se", "beta_se", "kappa_se", "mean_E1_se")] <- op_null$se[1:4]
    res_ingarch[k, c("tau_se", "beta_se", "kappa_se", "mean_E1_se")] <- op_ingarch$se[1:4]
    res[k, "convergence"] <- op$optim$convergence
    res_null[k, "convergence"] <- op_null$optim$convergence
    res_ingarch[k, "convergence"] <- op_ingarch$optim$convergence


    # Store the results every 10 iterations
    if (k %% 10 == 0) {
      print(paste0("Finished iteration: ", k, "."))
      write.csv(res, file = paste0("inst/LR_test/Results/extended_INARMA_fits_extended_INARMA_holds_s", s, "_", lgt, ".csv"),
                row.names = FALSE)
      write.csv(res_null, file = paste0("inst/LR_test/Results/INARMA_fits_extended_INARMA_holds_s", s, "_", lgt, ".csv"),
                row.names = FALSE)
      write.csv(res_null, file = paste0("inst/LR_test/Results/INGARCH_fits_extended_INARMA_holds_s", s, "_", lgt, ".csv"),
                row.names = FALSE)
    }
  }
  progress_log$already_done[i] <- TRUE
  print(paste0("Finished scenario ", s, " length ", lgt, " check the results."))
}

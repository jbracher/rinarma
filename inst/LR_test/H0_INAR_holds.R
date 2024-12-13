library(devtools)
load_all()

# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
vals_kappa <- c(0.5, 0.6, 0.8)
vals_zeta <- c(0.2, 0.5, 0.8)

vals_lgt <- c(250, 500, 1000) # lengths of simulated time series
n_sim <- 1000 # number of simulation runs

# Keep track of scenarios and lengths that are done
# Needs to be changed manually, when selecting scenarios and lengths
progress_log <- expand.grid(scenario = 1:3, lgt = vals_lgt, already_done = FALSE,
                            continue = FALSE)

sims_to_run <- 1:9
sims_to_run <- sims_to_run[!(progress_log$already_done)]

for (i in sims_to_run) {

  continue <- progress_log$continue[i] # should existing results be read (to take up after error)
  s <- progress_log$scenario[i]
  lgt <- progress_log$lgt[i]

  # run simulation:
  print(paste0("Started scenario ", s, " length ", lgt, "."))

  # grab true parameter values:
  tau <- vals_tau[s]
  kappa <- vals_kappa[s]

  # initialize matrix to store results:
  # get existing results if desired (if simulation somehow failed halfway)
  if (continue) {
    res <- read.csv(
      paste0("inst/LR_test/Results/H0_INAR_holds_results_s", s, "_", lgt,
             ".csv")
    )
    res_null <- read.csv(
      paste0("inst/LR_test/Results/H0_INAR_holds_results_null_s", s, "_", lgt,
             ".csv")
    )
    dat <- unname(read.csv(
      paste0("inst/LR_test/Results/pois_INAR_sim_", tau, "_", kappa, "_", lgt,
             ".csv"))
    )
    inds_to_run <- Position(res$kappa, f = is.na):n_sim
  } else {
    inds_to_run <- 1:n_sim
    res <- data.frame(
      kappa = rep(NA, n_sim),
      zeta = rep(NA, n_sim),
      tau = rep(NA, n_sim),
      kappa_se = rep(NA, n_sim),
      zeta_se = rep(NA, n_sim),
      tau_se = rep(NA, n_sim),
      max_lik = rep(NA, n_sim)
    )
    res_null <- data.frame(
      kappa = rep(NA, n_sim),
      tau = rep(NA, n_sim),
      kappa_se = rep(NA, n_sim),
      tau_se = rep(NA, n_sim),
      max_lik = rep(NA, n_sim)
    )
    dat <- matrix(nrow = n_sim, ncol = lgt)
  }

  # Simulate
  for (k in inds_to_run) {
    # Generate the data
    set.seed(k)
    dat[k, ] <- sim_inarma(beta = 0, kappa = kappa, tau = tau, lgt = lgt,
                           family = "Poisson", offspring = "binomial")$X
    # Fit
    op_null <- fit_inar(as.numeric(dat[k, ]), family = "Poisson",
                          offspring = "binomial", return_se = TRUE,
                        start = c(log_tau = 1, logit_kappa = 0))
    op <- fit_inar(as.numeric(dat[k, ]), family = "Poisson",
                     start = c(op_null$coefficients_raw, logit_zeta = 0),
                     offspring = "binomial-Poisson", return_se = TRUE)

    res$max_lik[k] <- op$loglikelihood
    res_null$max_lik[k] <- op_null$loglikelihood
    res[k, c("tau", "kappa", "zeta")] <- op$coefficients[1:3]
    res_null[k, c("tau", "kappa")] <- op_null$coefficients[1:2]
    res[k, c("tau_se",  "kappa_se", "zeta_se")] <- op$se[1:3]
    res_null[k, c("tau_se", "kappa_se")] <- op_null$se[1:2]

    # Store the results every 10 iterations
    if (k %% 10 == 0) {
      print(paste0("Finished iteration: ", k, "."))
      write.csv(res, file = paste0("inst/LR_test/Results/H0_INAR_holds_results_s", s, "_", lgt, ".csv"),
                row.names = FALSE)
      write.csv(res_null, file = paste0("Results/H0_INAR_holds_results_null_s", s, "_", lgt, ".csv"),
                row.names = FALSE)
      write.csv(dat, file = paste0("inst/LR_test/Results/pois_INAR_sim_", tau, "_",
                                   kappa, "_", lgt, ".csv"), row.names = FALSE)
    }
  }
  progress_log$already_done[i] <- TRUE
  print(paste0("Finished scenario ", s, " length ", lgt, " check the results."))
}

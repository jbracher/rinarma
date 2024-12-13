library(tscount)

# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
# vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)

n_sim <- 1000
vals_lgt <- c(250, 500, 1000)

#####################################################
## Run the simulation in the loop for 1000 iterations
#####################################################



continue <- FALSE # should existing results be read (to take up after error)

for (s in 1:3) {
  tau <- vals_tau[s]
  kappa <- vals_kappa[s]
  for (lgt in vals_lgt) {

    print(paste0("Started scenario ", s, " length ", lgt, "."))
    # initialize matrix to store results:
    # get existing results if required (if simulation somehow failed halfway)
    if (continue) {
      res <- read.csv(
        paste0("Results/H0_INARCH_holds_results_s", s, "_", lgt,
               ".csv")
      )
      res_null <- read.csv(
        paste0("Results/H0_INARCH_holds_results_null_s", s, "_", lgt,
               ".csv")
      )
      inds_to_run <- Position(res$kappa, f = is.na):n_sim
    } else {
      inds_to_run <- 1:n_sim
      res <- data.frame(
        kappa = rep(NA, n_sim),
        beta = rep(NA, n_sim),
        tau = rep(NA, n_sim),
        max_lik = rep(NA, n_sim)
      )
      res_null <- data.frame(
        kappa = rep(NA, n_sim),
        tau = rep(NA, n_sim),
        max_lik = rep(NA, n_sim)
      )
      dat <- matrix(nrow = n_sim, ncol = lgt)
    }

    # Simulate
    for (k in inds_to_run) {
      # Generate the data
      set.seed(k)
      dat[k, ] <- tsglm.sim(lgt, param = list(intercept = tau, past_obs = kappa),
                            model = list(past_obs = 1))$ts
      # Fit
      op <- tsglm(dat[k, ], list(past_obs = 1, past_mean = 1))
      op_null <- tsglm(dat[k, ], list(past_obs = 1))

      res$max_lik[k] <- op$logLik
      res_null$max_lik[k] <- op_null$logLik
      res[k, c("tau", "kappa", "beta")] <-
        c(op$coefficients[1:2] / (1 - op$coefficients[3]), op$coefficients[3])
      res_null[k, c("tau", "kappa")] <- op_null$coefficients[1:2]

      # Store the results every 10 iterations
      if (k %% 10 == 0) {
        print(paste0("Finished iteration: ", k, "."))
        write.csv(res, file = paste0("Results/H0_INARCH_holds_results_s", s, "_", lgt, ".csv"),
                  row.names = FALSE)
        write.csv(res_null, file = paste0("Results/H0_INARCH_holds_results_null_s", s, "_", lgt, ".csv"),
                  row.names = FALSE)
      }
    }
  }
}


library(xtable)

# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
vals_psi <- c(NA, NA, NA)
vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)

vals_lgt <- c(250, 500, 1000) # lengths of simulated time series
n_sim <- 1000 # number of simulation runs


# summarize results:
# initialize lists
results_sim <- list()
removed_est_se <- list()
n_available <- matrix(ncol = 3, nrow = 3,
                      dimnames = list(vals_lgt, paste0("sc", 1:3)))

for(s in 1:3){ # loop over scenarios
  tau <- vals_tau[s]
  psi <- vals_psi[s]
  beta <- vals_beta[s]
  kappa <- vals_kappa[s]

  scenario_temp <- matrix(NA, nrow = 3, ncol = 25,
                          dimnames = list(c(250, 500, 1000),
                                          c("tau", "psi", "beta", "kappa", "T",
                                            "mean_tau", "se_tau", "est_se_tau", "mean_moment_tau", "se_moment_tau",
                                            "mean_psi", "se_psi", "est_se_psi", "mean_moment_psi", "se_moment_psi",
                                            "mean_beta", "se_beta", "est_se_beta", "mean_moment_beta", "se_moment_beta",
                                            "mean_kappa", "se_kappa", "est_se_kappa", "mean_moment_kappa", "se_moment_kappa")))
  scenario_temp[1, c("tau", "psi", "beta", "kappa")] <- c(tau, psi, beta, kappa)
  scenario_temp[, "T"] <- c(250, 500, 1000)
  removed_est_se_temp <- matrix(NA, ncol = 3, nrow = 3,
                                dimnames = list(c(250, 500, 1000), c("tau", "beta", "kappa")))

  # loop over lengths of time series:
  for(lgt in c(250, 500, 1000)){

    # get parameter estimates:
    res <- read.csv(
      file = paste0("inst/LR_test/Results/INARMA_fits_INARMA_holds_s", s, "_", lgt, ".csv")
    )[, c("tau", "beta", "kappa", "tau_se", "beta_se", "kappa_se", "convergence")]

    pars <- res[, c("tau", "beta", "kappa")]
    ses <- res[, c("tau_se", "beta_se", "kappa_se")]

    n_available_temp <- max(which(!is.na(pars[, 1])))
    pars <- pars[1:n_available_temp, ]
    n_available[as.character(lgt), s] <- n_available_temp

    # fill in means and standard deviations:
    scenario_temp[as.character(lgt), "mean_tau"] <- mean(pars[res$convergence == 0, "tau"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_beta"] <- mean(pars[res$convergence == 0, "beta"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_kappa"] <- mean(pars[res$convergence == 0, "kappa"], na.rm = TRUE)

    scenario_temp[as.character(lgt), "se_tau"] <- sd(pars[res$convergence == 0, "tau"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "se_beta"] <- sd(pars[res$convergence == 0, "beta"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "se_kappa"] <- sd(pars[res$convergence == 0, "kappa"], na.rm = TRUE)

    # get estimated standard errors:

    ses[ses > 30] <- NA  # set numerically instable ses to NA
    ses[res$convergence != 0, ] <- NA  # set numerically instable ses to NA
    scenario_temp[as.character(lgt), c("est_se_tau", "est_se_beta", "est_se_kappa")] <-
      apply(ses, 2, mean, na.rm = TRUE)

    # how many times was se estimation not possible or obviously instable?
    removed_est_se_temp[as.character(lgt), ] <- colSums(is.na(ses))
    # get moment estimates:
    pars_moment <- read.csv(file = paste0("inst/Estimation/Results/Estim_moments_s", s, "_", lgt, "_pois.csv"))[, c("tau", "beta", "kappa")]

    scenario_temp[as.character(lgt), "mean_moment_tau"] <- mean(pars_moment[, "tau"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_moment_beta"] <- mean(pars_moment[, "beta"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_moment_kappa"] <- mean(pars_moment[, "kappa"], na.rm = TRUE)

    scenario_temp[as.character(lgt), "se_moment_tau"] <- sd(pars_moment[, "tau"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "se_moment_beta"] <- sd(pars_moment[, "beta"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "se_moment_kappa"] <- sd(pars_moment[, "kappa"], na.rm = TRUE)
  }

  # store results in list:
  results_sim[[s]] <- scenario_temp
  removed_est_se[[s]] <- removed_est_se_temp
}


# generate tex code:
for(s in 1:3){
  to_print <- matrix(as.character(format(results_sim[[s]],
                                         digits = 2,
                                         scientific = FALSE)),
                     nrow = 3)
  colnames(to_print) <- colnames(results_sim[[s]])

  to_print[, 5] <- vals_lgt
  to_print[2:3, 1:4] <- ""

  # replace NA values:
  to_print[grepl("NA", to_print)] <- "-"

  # re-order columns for table on ML estimate:
  to_print_ml <- to_print[, c("T",
                              "tau", "mean_tau", "se_tau", "est_se_tau",
                              "psi", "mean_psi", "se_psi", "est_se_psi",
                              "beta", "mean_beta", "se_beta", "est_se_beta",
                              "kappa", "mean_kappa", "se_kappa", "est_se_kappa")]


  write(
    print(xtable(to_print_ml), only.contents = TRUE,
          include.rownames = FALSE, include.colnames = FALSE,
          hline.after = NULL),
    file = paste0("inst/Estimation/Tables/sim_pois_sc", s, ".tex")
  )

  # re-order columns for table on moment estimate:
  to_print_moment <- to_print[, c("T",
                                  "tau", "mean_moment_tau", "se_moment_tau",
                                  "psi", "mean_moment_psi", "se_moment_psi",
                                  "beta", "mean_moment_beta", "se_moment_beta",
                                  "kappa", "mean_moment_kappa", "se_moment_kappa")]

  write(
    print(xtable(to_print_moment), only.contents = TRUE,
          include.rownames = FALSE, include.colnames = FALSE,
          hline.after = NULL),
    file = paste0("inst/Estimation/Tables/sim_pois_sc", s, "_moments.tex")
  )

}



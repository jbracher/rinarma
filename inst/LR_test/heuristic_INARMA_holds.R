
# function to find mean of normal for which the p-quantile equals 0:
find_mean_of_normal <- function(p, min = -2, max = 2){
  grid_mu <- seq(from = min, to = max, by = 0.001)
  p_quantiles <- qnorm(rep(p, length(grid_mu)),
                       grid_mu,
                       rep(1, length(grid_mu)))
  if(min(abs(p_quantiles)) > 0.001) warning("No good solution found, please increase search range.")
  ind_min <- which.min(abs(p_quantiles))
  mu <- grid_mu[ind_min]
  return(mu)
}

# wrapping threshold correction into function. The result is on the original chi2-scale!
correct_threshold_INARMA <- function(tau, beta, kappa, alpha = c(0.05, 0.1), n_sim = 1000){

  # initialize vector to store var / mean per path
  rel_var <- numeric(n_sim)

  # simulate (note: no need to store all runs):
  for(k in 1:n_sim) {
    # Generate the data
    set.seed(k)

    sim_temp <- sim_inarma(beta = beta, kappa = kappa, tau = tau, lgt = lgt,
                           family = "Poisson", offspring = "binomial")$X

    rel_var[k] <- var(sim_temp) / mean(sim_temp)
    if (var(sim_temp) == 0 & mean(sim_temp) == 0) rel_var[k] <- 1
  }

  # which quantile-level does 0 correspond to?
  p <- mean(rel_var < 1)

  # get mean of a corrected normal distribution with sd 1 (for square root of LR stat):
  mu_corrected <- find_mean_of_normal(p)

  # uncorrected critical value on chi2-scale:
  critical_uncorrected <- qnorm(1 - alpha, 0, 1)^2
  # get corrected critical value:
  critical_corrected <- qnorm(1 - alpha, mu_corrected, 1)^2
  names(critical_corrected) <- as.character(paste0("critical_", alpha))

  return(list(critical_corrected = critical_corrected,
              critical_uncorrected = critical_uncorrected))
}

# significance levels:
alpha <- c(0.05, 0.1, 0.2)


vals_lgt <- c(250, 500, 1000)

for (s in 1:3) {
  for (lgt in vals_lgt) {

    # get in results
    res_null <- read.csv(paste0("inst/LR_test/Results/H0_INARMA_holds_results_null_s", s, "_", lgt,".csv"))
    res <- read.csv(paste0("inst/LR_test/Results/H0_INARMA_holds_results_s", s, "_", lgt, ".csv"))

    # to store corrected thresholds:
    corrected_thresholds <- as.data.frame(
      matrix(nrow = nrow(res), ncol = length(alpha),
             dimnames = list(NULL, as.character(paste0("critical_", alpha))))
    )

    for(k in 1:nrow(res)){  # run through iterations:

      # grab fitted parameters
      tau <- res_null$tau[k]
      kappa <- res_null$kappa[k]
      beta <- res_null$beta[k]

      # compute corrected threshold:
      corrected_thresholds[k, ] <- correct_threshold_INARMA(tau, beta, kappa, n_sim = 500, alpha = alpha)$critical_corrected

      if (k %% 20 == 0) {
        print(paste0("Finished iteration ", k))
      }
    }

    # save the corrected thresholds
    write_csv(
      corrected_thresholds,
      file = paste0("inst/LR_test/Results/H0_INARMA_holds_correct_thresholds_s", s, "_", lgt,".csv")
    )
  }
}

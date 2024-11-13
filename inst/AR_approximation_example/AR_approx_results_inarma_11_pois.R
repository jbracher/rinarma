library(tidyverse)
library(xtable)

# Setup ------------------------------------------------------------------------

vals_tau <- c(1, 1, 1)
vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)

n_sim <- 1000 # number of simulation runs
lag_max <- 8  # The maximum lag of the approximating AR model
vals_lgt <- c(250, 500, 1000)  # lengths of simulated time series

lag_max <- 10  # The maximum lag of the approximating AR model
scenario_temp <- matrix(
  NA,
  nrow = 3,
  ncol = 13,
  dimnames = list(
    c(250, 500, 1000),
    c("T", "tau", "mean_tau", "se_tau", "est_se_tau",
      "beta", "mean_beta", "se_beta", "est_se_beta",
      "kappa", "mean_kappa", "se_kappa", "est_se_kappa")
  )
)

for (s in 1:3) {  # Loop over the scenarios
  for (lgt in vals_lgt) {  # Loop over the time series lengths

    # Load the results
    res <- read_csv(
      file = paste0("inst/AR_approximation_example/Results/AR_INARMA11_s", s,
                    "_lgt", lgt, "_pois.csv")
    )

    # Find problematic rows
    divergent <- which(res$convergence != 0)
    other_problems <- which(is.na(res[, 1:6]) | res[, 1:6] <= 0, arr.ind = TRUE)[, 1]
    problematic <- unique(c(divergent, other_problems))

    print(paste0("For length ", lgt, " there were ", length(problematic),
                 " problematic optimizations. Selecting only the good ones."))

    # Select only the good rows
    good <- setdiff(1:n_sim, problematic)

    # Aggregate
    scenario_temp[as.character(lgt), c("tau", "beta", "kappa")] <-
      c(vals_tau[s], vals_beta[s], vals_kappa[s])
    scenario_temp[as.character(lgt), c("mean_tau", "mean_beta", "mean_kappa")] <-
      apply(res[good, c("tau", "beta", "kappa")], 2, mean)
    scenario_temp[as.character(lgt), c("se_tau", "se_beta", "se_kappa")] <-
      apply(res[good, c("tau", "beta", "kappa")], 2, sd)
    scenario_temp[as.character(lgt), c("est_se_tau", "est_se_beta", "est_se_kappa")] <-
      apply(res[good, c("tau_se", "beta_se", "kappa_se")], 2, mean)

    scenario_temp[, "T"] <- vals_lgt
  }

  # Format the results

  to_print <- scenario_temp
  to_print[, 3:13] <- round(scenario_temp[, 3:13], digits = 3)
  to_print <- apply(to_print, 2, as.character)
  to_print[, 3:13] <- apply(to_print[, 3:13], 2, str_pad, pad = "0", width = 5, side = "right")
  to_print[1, "tau"] <- "1.000"

  to_insert <- matrix("-", nrow = 3, ncol = 4)
  colnames(to_insert) <- c("psi", "mean_psi", "se_psi", "est_se_psi")
  to_print <- cbind(
    to_print[, c("T", "tau", "mean_tau", "se_tau", "est_se_tau")],
    to_insert,
    to_print[, c("beta", "mean_beta", "se_beta", "est_se_beta", "kappa",
                 "mean_kappa", "se_kappa", "est_se_kappa")]
  )

  to_print[2:3, c("tau", "kappa", "beta", "psi")] <- ""

  write(
    print(
      xtable(to_print),
      only.contents = TRUE,
      include.rownames = FALSE, include.colnames = FALSE,
      hline.after = NULL),
    file = paste0("inst/AR_approximation_example/Tables/AR_approx_inarma11_sc", s, "_pois.tex")
  )
}

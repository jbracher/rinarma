## Script generating the result table of fitting the Hermite INARMA(2, 1) using
## the AR-approximation

library(tidyverse)
library(xtable)

# Setup ------------------------------------------------------------------------

n_sim <- 1000  # Number of simulations

# Set the parameter values
s <- 1
kappa <- c(0.2, 0.6)
beta <- 0.15
tau <- 2
psi <- 0.7

# Alternative set of parameters, must be uncommented manually to run
# s <- 2
# kappa <- c(0.45, 0.25)
# beta <- 0.35
# tau <- 1.5
# psi <- 0.5

lag_max <- 10  # The maximum lag of the approximating AR model
vals_lgt <- c(250, 500, 1000, 2000)
scenario_temp <- matrix(
  NA,
  nrow = 4,
  ncol = 21,
  dimnames = list(
    c(250, 500, 1000, 2000),
    c("T", "tau", "mean_tau", "se_tau", "est_se_tau",
      "psi", "mean_psi", "se_psi", "est_se_psi",
      "beta", "mean_beta", "se_beta", "est_se_beta",
      "kappa1", "mean_kappa1", "se_kappa1", "est_se_kappa1",
      "kappa2", "mean_kappa2", "se_kappa2", "est_se_kappa2")
  )
)

scenario_temp[, "T"] <- vals_lgt
scenario_temp[, "kappa1"] <- kappa[1]
scenario_temp[, "kappa2"] <- kappa[2]
scenario_temp[, "beta"] <- beta
scenario_temp[, "tau"] <- tau
scenario_temp[, "psi"] <- psi

for (lgt in vals_lgt) {

  # Load the results
  res <- read_csv(
    file = paste0("inst/AR_approximation_example/Results/AR_INARMA21_s", s,
                  "_lgt_", lgt, "_herm.csv")
  )

  # Find problematic rows
  divergent <- which(res$convergence != 0)
  other_problems <- which(is.na(res[, 1:10]) | res[, 1:10] <= 0, arr.ind = TRUE)[, 1]
  problematic <- unique(c(divergent, other_problems))

  print(paste0("For length ", lgt, " there were ", length(problematic),
               " problematic optimizations (divergent optimizations, zero
               standard errors). Selecting only the good ones."))

  # Select only the good rows
  good <- setdiff(1:n_sim, problematic)

  # Aggregate
  scenario_temp[as.character(lgt), c("mean_tau", "mean_psi", "mean_beta", "mean_kappa1", "mean_kappa2")] <-
    apply(res[good, c("tau", "psi", "beta", "kappa1", "kappa2")], 2, mean)
  scenario_temp[as.character(lgt), c("se_tau", "se_psi", "se_beta", "se_kappa1", "se_kappa2")] <-
    apply(res[good, c("tau", "psi", "beta", "kappa1", "kappa2")], 2, sd)
  scenario_temp[as.character(lgt), c("est_se_tau", "est_se_psi", "est_se_beta", "est_se_kappa1", "est_se_kappa2")] <-
    apply(res[good, c("tau_se", "psi_se", "beta_se", "kappa1_se", "kappa2_se")], 2, mean)
}

# Format the results
to_print <- scenario_temp
to_print[, 3:21] <- round(scenario_temp[, 3:21], digits = 3)
to_print <- apply(to_print, 2, as.character)
to_print[, 4:21] <- apply(to_print[, 4:21], 2, str_pad, pad = "0", width = 5, side = "right")
to_print[1, "tau"] <- ifelse(s == 1, "2.000", "1.500")
to_print[2:4, c("tau", "kappa1", "kappa2", "beta", "psi")] <- ""

write(
  print(
    xtable(to_print),
    only.contents = TRUE,
    include.rownames = FALSE, include.colnames = FALSE,
    hline.after = NULL),
  file = paste0("inst/AR_approximation_example/Tables/AR_approx_inarma21_sc", s, "herm.tex")
)

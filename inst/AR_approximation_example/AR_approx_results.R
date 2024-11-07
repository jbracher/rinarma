library(tidyverse)
library(xtable)

# Setup ------------------------------------------------------------------------

# Set the parameter values
s <- 1
kappa <- c(0.2, 0.6)
beta <- 0.15
tau <- 2

# Alternative set of parameters
# s <- 2
# kappa <- c(0.45, 0.25)
# beta <- 0.35
# tau <- 1.5

lag_max <- 10  # The maximum lag of the approximating AR model
vals_lgt <- c(250, 500, 1000)
scenario_temp <- matrix(
  NA,
  nrow = 3,
  ncol = 17,
  dimnames = list(
    c(250, 500, 1000),
    c("T", "tau", "mean_tau", "se_tau", "est_se_tau",
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

for (lgt in vals_lgt) {

  # Load the results
  res <- read_csv(
    file = paste0("inst/AR_approximation_example/Results/AR_INARMA21_s", s,
                  "_lgt", lgt, ".csv")
    )

  # Find problematic rows
  divergent <- which(res$convergence != 0)
  other_problems <- which(is.na(res[, 1:8]) | res[, 1:8] <= 0, arr.ind = TRUE)[, 1]
  problematic <- unique(c(divergent, other_problems))

  print(paste0("For length ", lgt, " there were ", length(problematic),
               " problematic optimizations. Selecting only the convergent ones."))

  # Select only the good rows
  good <- setdiff(1:n_sim, problematic)

  # Aggregate
  scenario_temp[as.character(lgt), c("mean_tau", "mean_beta", "mean_kappa1", "mean_kappa2")] <-
    apply(res[good, c("tau", "beta", "kappa1", "kappa2")], 2, mean)
  scenario_temp[as.character(lgt), c("se_tau", "se_beta", "se_kappa1", "se_kappa2")] <-
    apply(res[good, c("tau", "beta", "kappa1", "kappa2")], 2, sd)
  scenario_temp[as.character(lgt), c("est_se_tau", "est_se_beta", "est_se_kappa1", "est_se_kappa2")] <-
    apply(res[good, c("tau_se", "beta_se", "kappa1_se", "kappa2_se")], 2, mean)
}

write(
  print(
    xtable(scenario_temp, digits=c(0, 0, rep(3, ncol(scenario_temp) - 1))),
        only.contents = TRUE,
        include.rownames = FALSE, include.colnames = FALSE,
        hline.after = NULL),
  file = paste0("inst/AR_approximation_example/Tables/AR_approx_inarma21_sc", s, ".tex")
)

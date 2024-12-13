library(tidyverse)
library(tscount)
library(ks)

vals_tau <- c(1, 1, 1)
vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)
vals_zeta <- c(0.2, 0.5, 0.8)

#########################################################
## H0: zeta = 0 (binomial thinning only) vs. H1: zeta > 0
## H0 holds
#########################################################

# Scenario 1 is done for length 500 & 1000, for length 500 iteration 511 failed
# for some reason, I'll fix it later
# Scenario 2 is done for length 250
# Scenario 3 is done for length 250

s <- 3  # Scenario number
lgt <- 250  # Length

res <- read.csv(paste0("inst/LR_test/Results/H0_INARMA_holds_results_s", s, "_", lgt, ".csv"))
res_null <- read.csv(paste0("inst/LR_test/Results/H0_INARMA_holds_results_null_s", s, "_", lgt, ".csv"))

# Check the point estimates
mean(res$kappa, na.rm = TRUE) - vals_kappa[s]
mean(res_null$kappa, na.rm = TRUE) - vals_kappa[s]
mean(res$tau, na.rm = TRUE) - vals_tau[s]
mean(res_null$tau, na.rm = TRUE) - vals_tau[s]
mean(res$beta, na.rm = TRUE) - vals_beta[s]
mean(res_null$beta, na.rm = TRUE) - vals_beta[s]

# Calculate the LR statistic
LR_stat <- -2 * (res_null$max_lik - res$max_lik)

# Select only the convergent iterations
convergent <- which(res$convergence == 0 & res_null$convergence == 0)
length(convergent)
LR_stat_corrected <- LR_stat[convergent]

# Set the negative values of the LR statistic to zero
LR_stat_corrected2 <- ifelse(LR_stat_corrected < 0, 0, LR_stat_corrected)

# The 95% critical value of the mixture distribution
critical <- qchisq(0.9, 1)

# Display with ggplot
step <- 0.1
xx_gg <- seq(0, max(LR_stat_corrected2), by = step)  # breaks for the histogram
yy_gg <- 0.5 * diff(pchisq(c(xx_gg, xx_gg[length(xx_gg)] + step), 1)) / step
yy_gg[1] <- yy_gg[1] + 0.5 / step

ggplot() +
  geom_col(aes(x = xx_gg, y = yy_gg), alpha = 0.5, fill = "red",
           width = 0.1, just = 0) +
  geom_histogram(aes(x = LR_stat_corrected2, y = after_stat(density)), alpha = 0.5,
                 breaks = xx_gg) +
  geom_vline(xintercept = critical, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 7.5)) +
  labs(title = expression("H0:"~zeta~"= 0 (INARMA only) vs. H1:"~zeta~"> 0"),
       subtitle = bquote(Scenario == .(s)~~~Length == .(lgt)~~~"H0 holds"), x = "LR test statistic", y = "density") +
  theme_bw()

# "Zoom in"
ggplot() +
  geom_col(aes(x = xx_gg, y = yy_gg), alpha = 0.5, fill = "red",
           width = 0.1, just = 0) +
  geom_histogram(aes(x = LR_stat_corrected2, y = after_stat(density)), alpha = 0.5,
                 breaks = xx_gg) +
  geom_vline(xintercept = critical, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 1.5)) +
  labs(title = expression("H0:"~zeta~"= 0 (INARMA only) vs. H1:"~zeta~"> 0"),
       subtitle = bquote(Scenario == .(s)~~~Length == .(lgt)~~~"H0 holds"), x = "LR test statistic", y = "density") +
  theme_bw()

# Calculate the p-values
p_val <- 1 - 0.5 * (ifelse(LR_stat_corrected2 > 0, 1, 0) + pchisq(LR_stat_corrected2, 1))

# Display the histogram of p-values - THERE IS A SPIKE BEFORE 0.5
hist(p_val, breaks = seq(0, 1, by = 0.05), ylim = c(0, 510))
abline(h = 500, col = "red")

# We have too many p-values between 0.49 and 0.5
hist(p_val, breaks = seq(0, 1, by = 0.01), ylim = c(0, 510))
abline(h = 500, col = "red")

# Very small values of the LR-statistic
summary(LR_stat_corrected2[p_val >= 0.49 & p_val < 0.5])

# A lot of close to zero estimates of the additional parameter zeta
summary(res$zeta[p_val >= 0.49 & p_val < 0.5])

# Calculate the error of the first kind
alpha <- sum(p_val < 0.05) / length(LR_stat_corrected2)
# Alternative calculation (same value)
# alpha <- sum(LR_stat_corrected2 > critical) / length(LR_stat_corrected2)
print(paste("Type I error:", alpha * 100, "%"))

########################################################
## H0: zeta = 0 (INAR) vs. H1: zeta > 0
## H0 holds
########################################################

s <- 3  # Scenario number
lgt <- 500  # Length

res <- read.csv(paste0("inst/LR_test/Results/H0_INAR_holds_results_s", s, "_", lgt, ".csv"))
res_null <- read.csv(paste0("inst/LR_test/Results/H0_INAR_holds_results_null_s", s, "_", lgt, ".csv"))

# Check the point estimates
mean(res$kappa) - vals_kappa[s]
mean(res_null$kappa) - vals_kappa[s]
mean(res$tau) - vals_tau[s]
mean(res_null$tau) - vals_tau[s]

# Calculate the LR statistic
LR_stat <- -2 * (res_null$max_lik - res$max_lik)
LR_stat_corrected <- ifelse(LR_stat < 0, 0, LR_stat)

# The 95% critical value of the mixture distribution
critical <- qchisq(0.9, 1)

# # Display with ggplot
# step <- 0.1
# xx_gg <- seq(0, max(LR_stat_corrected), by = step)  # breaks for the histogram
# yy_gg <- 0.5 * diff(pchisq(c(xx_gg, xx_gg[length(xx_gg)] + step), 1)) / step
# yy_gg[1] <- yy_gg[1] + 0.5 / step
#
# ggplot() +
#   geom_col(aes(x = xx_gg, y = yy_gg), alpha = 0.5, fill = "red",
#            width = 0.1, just = 0) +
#   geom_histogram(aes(x = LR_stat_corrected, y = after_stat(density)), alpha = 0.5,
#                  breaks = xx_gg) +
#   geom_vline(xintercept = critical, linetype = "dashed") +
#   coord_cartesian(xlim = c(0, 5), ylim = c(0, 7.5)) +
#   labs(title = expression("H0:"~zeta~"= 0 (INGARCH only) vs. H1:"~zeta~"> 0"),
#        subtitle = bquote(Scenario == .(s)~~~Length == .(lgt)~~~"H0 holds"), x = "LR test statistic", y = "density") +
#   theme_bw()
#
# # "Zoom in" at the beginning
# ggplot() +
#   geom_col(aes(x = xx_gg, y = yy_gg), alpha = 0.5, fill = "red",
#            width = 0.1, just = 0) +
#   geom_histogram(aes(x = LR_stat_corrected, y = after_stat(density)), alpha = 0.5,
#                  breaks = xx_gg) +
#   geom_vline(xintercept = critical, linetype = "dashed") +
#   coord_cartesian(xlim = c(0, 5), ylim = c(0, 1.5)) +
#   labs(title = expression("H0:"~zeta~"= 0 (INGARCH only) vs. H1:"~zeta~"> 0"),
#        subtitle = bquote(Scenario == .(s)~~~Length == .(lgt)~~~"H0 holds"), x = "LR test statistic", y = "density") +
#   theme_bw()
#
# # "Zoom in" at the end
# ggplot() +
#   geom_col(aes(x = xx_gg, y = yy_gg), alpha = 0.5, fill = "red",
#            width = 0.1, just = 0) +
#   geom_histogram(aes(x = LR_stat_corrected, y = after_stat(density)), alpha = 0.5,
#                  breaks = xx_gg) +
#   geom_vline(xintercept = critical, linetype = "dashed") +
#   coord_cartesian(xlim = c(2, 8), ylim = c(0, 0.1)) +
#   labs(title = expression("H0:"~zeta~"= 0 (INGARCH only) vs. H1:"~zeta~"> 0"),
#        subtitle = bquote(Scenario == .(s)~~~Length == .(lgt)~~~"H0 holds"), x = "LR test statistic", y = "density") +
#   theme_bw()

# Calculate the p-values
p_val <- 1 - 0.5 * (ifelse(LR_stat_corrected > 0, 1, 0) + pchisq(LR_stat_corrected, 1))

# # We have too many p-values between 0.49 and 0.5
# hist(p_val, breaks = seq(0, 1, by = 0.01), ylim = c(0, 510))
# abline(h = 500, col = "red")

# Calculate the error of the first kind
alpha <- sum(p_val < 0.05) / 1000

print(paste("Type I error:", alpha * 100, "%"))

plot(x = res$zeta[LR_stat <  critical], y = LR_stat[LR_stat <  critical],
     xlim = c(0, 1), ylim = c(0, 12), xlab = expression(zeta),
     ylab = "LR Statistic")
abline(h = critical, col = "red")
points(x = res$zeta[LR_stat >=  critical], y = LR_stat[LR_stat >=  critical],
       col = "red")

########################################################
## H0: beta = 0 (INARCH) vs. H1: beta > 0 (INGARCH)
## H0 holds
########################################################

s <- 1  # Scenario number
lgt <- 500  # Length

res <- read.csv(paste0("inst/LR_test/Results/H0_INARCH_holds_results_s", s, "_", lgt, ".csv"))
res_null <- read.csv(paste0("inst/LR_test/Results/H0_INARCH_holds_results_null_s", s, "_", lgt, ".csv"))
# res_new <- read.csv(paste0("Results/H0_INGARCH_holds_results_s", s, "_", lgt, "_refit.csv"))

# Check the point estimates
mean(res$kappa) - vals_kappa[s]
mean(res_null$kappa) - vals_kappa[s]
mean(res$tau) - vals_tau[s]
mean(res_null$tau) - vals_tau[s]

# Calculate the LR statistic
LR_stat <- -2 * (res_null$max_lik - res$max_lik)
LR_stat_corrected <- ifelse(LR_stat < 0, 0, LR_stat)

# The 95% critical value of the mixture distribution
critical <- qchisq(0.9, 1)

# Display with ggplot
step <- 0.1
xx_gg <- seq(0, max(LR_stat_corrected), by = step)  # breaks for the histogram
yy_gg <- 0.5 * diff(pchisq(c(xx_gg, xx_gg[length(xx_gg)] + step), 1)) / step
yy_gg[1] <- yy_gg[1] + 0.5 / step

ggplot() +
  geom_col(aes(x = xx_gg, y = yy_gg), alpha = 0.5, fill = "red",
           width = 0.1, just = 0) +
  geom_histogram(aes(x = LR_stat_corrected, y = after_stat(density)), alpha = 0.5,
                 breaks = xx_gg) +
  geom_vline(xintercept = critical, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 7.5)) +
  labs(title = expression("H0:"~zeta~"= 0 (INGARCH only) vs. H1:"~zeta~"> 0"),
       subtitle = bquote(Scenario == .(s)~~~Length == .(lgt)~~~"H0 holds"), x = "LR test statistic", y = "density") +
  theme_bw()

# "Zoom in" at the beginning
ggplot() +
  geom_col(aes(x = xx_gg, y = yy_gg), alpha = 0.5, fill = "red",
           width = 0.1, just = 0) +
  geom_histogram(aes(x = LR_stat_corrected, y = after_stat(density)), alpha = 0.5,
                 breaks = xx_gg) +
  geom_vline(xintercept = critical, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 1.5)) +
  labs(title = expression("H0:"~zeta~"= 0 (INGARCH only) vs. H1:"~zeta~"> 0"),
       subtitle = bquote(Scenario == .(s)~~~Length == .(lgt)~~~"H0 holds"), x = "LR test statistic", y = "density") +
  theme_bw()

# "Zoom in" at the end
ggplot() +
  geom_col(aes(x = xx_gg, y = yy_gg), alpha = 0.5, fill = "red",
           width = 0.1, just = 0) +
  geom_histogram(aes(x = LR_stat_corrected, y = after_stat(density)), alpha = 0.5,
                 breaks = xx_gg) +
  geom_vline(xintercept = critical, linetype = "dashed") +
  coord_cartesian(xlim = c(2, 8), ylim = c(0, 0.1)) +
  labs(title = expression("H0:"~zeta~"= 0 (INGARCH only) vs. H1:"~zeta~"> 0"),
       subtitle = bquote(Scenario == .(s)~~~Length == .(lgt)~~~"H0 holds"), x = "LR test statistic", y = "density") +
  theme_bw()

# Calculate the p-values
p_val <- 1 - 0.5 * (ifelse(LR_stat_corrected > 0, 1, 0) + pchisq(LR_stat_corrected, 1))

# Display the histogram of p-values - THERE IS A SPIKE BEFORE 0.5
hist(p_val, breaks = seq(0, 1, by = 0.05), ylim = c(0, 550))
abline(h = 500, col = "red")

# We have still quite a lot p-values between 0.49 and 0.5
hist(p_val, breaks = seq(0, 1, by = 0.01), ylim = c(0, 510))
abline(h = 500, col = "red")

# Small values of the LR-statistic
summary(LR_stat_corrected[p_val >= 0.49 & p_val < 0.5])

# A lot of close to zero estimates of the additional parameter beta
length(res$beta[p_val >= 0.49 & p_val < 0.5])
summary(res$beta[p_val >= 0.49 & p_val < 0.5])

# Calculate the error of the first kind
alpha <- sum(p_val < 0.05) / length(LR_stat_corrected2)
# Alternative calculation (same value)
# alpha <- sum(LR_stat_corrected2 > critical) / length(LR_stat_corrected2)
print(paste("Type I error:", alpha * 100, "%"))



print(paste("Type I error:", alpha * 100, "%"))

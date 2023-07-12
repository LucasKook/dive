# Simple normal IV example with PIT residuals
# LK 2023

# DEPs --------------------------------------------------------------------

set.seed(1111)

library("tram")
library("coin")
library("dHSIC")

# GEN ---------------------------------------------------------------------

simple_dgp <- function(n = 1e5, doD = FALSE) {
  H <- rnorm(n)
  Z <- sample(0:1, n, TRUE)
  D <- as.numeric((1 - doD) * H + Z > rnorm(n))
  Y <- D + H + rnorm(n)
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

dobs <- simple_dgp()
dint <- simple_dgp(n = 1e5 - 1, doD = TRUE)

# FIT ---------------------------------------------------------------------

mobs <- Lm(Y ~ D, data = dobs, prob = c(0, 1))
mint <- Lm(Y ~ D, data = dint, prob = c(0, 1))

# RUN ---------------------------------------------------------------------

### Residuals observational data observational distribution
pit_obs <- predict(mobs, type = "distribution")

### Residuals observational data interventional distribution
pit_int <- predict(mint, newdata = dobs, type = "distribution")

### Plot residuals against instrument and check independence (obs)
boxplot(pit_obs ~ Z, data = dobs)
spearman_test(pit_obs ~ Z, data = dobs)
# dhsic.test(pit_obs, dobs$Z, method = "gamma")$p.value

# In-sample PIT are marginally uniform
hist(pit_obs, probability = TRUE)
abline(h = 1, lwd = 2, col = 2)

# but not conditional on Z
hist(pit_obs[dobs$Z == 0], col = rgb(.5, .1, .1, .1), main = "PIT given Z", xlab = "PIT", probability = TRUE)
hist(pit_obs[dobs$Z == 1], col = rgb(.1, .1, .5, .1), add = TRUE, probability = TRUE)
abline(h = 1, lwd = 2, col = 2)

### Plot residuals against instrument and check independence (int)
boxplot(pit_int ~ Z, data = dobs)
spearman_test(pit_int ~ Z, data = dobs)
# dhsic.test(pit_int, dobs$Z, method = "gamma")$p.value

# PIT from the interventional distribution are marginally uniform
hist(pit_int, probability = TRUE)
abline(h = 1, lwd = 2, col = 2)

# and remain so conditional on Z
hist(pit_int[dobs$Z == 0], col = rgb(.5, .1, .1, .1), main = "PIT given Z", xlab = "PIT", probability = TRUE)
hist(pit_int[dobs$Z == 1], col = rgb(.1, .1, .5, .1), add = TRUE, probability = TRUE)
abline(h = 1, lwd = 2, col = 2)

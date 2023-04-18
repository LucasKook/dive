# Distributional random forest with control function
# LK March 2023

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("ranger")
library("tram")

# FUNs --------------------------------------------------------------------

gen_dat <- function(n = 1e3, doD = FALSE, nfine = 1e6) {
  ### Instrument
  Z <- sample(c(-1, 1), n, TRUE) # rt(n, df = 5)
  ### Hidden
  fH <- rt(nfine, df = 5)
  H <- rt(n, df = 5)
  ### Treatment
  UD <- runif(n)
  D <- as.numeric(plogis(Z + (1 - doD) * H) >= UD)
  ### Response
  UY <- runif(n)
  Y <- qnorm(UY, mean = D + H, sd = 1 + abs(D + H))
  ### Return
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

# Oracle ------------------------------------------------------------------

oracle_distr <- Vectorize(\(y, d = 0) {
  integrate(\(x) pnorm(y, mean = d + x, sd = 1 + abs(d + x)) * dt(x, df = 5),
            -20, 20)$value
}, "y")

ys <- seq(-13, 13, length.out = 1e3)
F1 <- oracle_distr(ys, d = 1)
F0 <- oracle_distr(ys, d = 0)

plot(ys, F0, type = "l", lty = 2)
lines(ys, F1, type = "l", col = 2, lty = 2)

# Sanity checks -----------------------------------------------------------

# Data under intervention on D (d0) and observational (d)
d0 <- gen_dat(1e3, doD = TRUE)
d1 <- gen_dat(1e3, doD = FALSE)

# plot(ecdf(d0$Y[d0$D == 0]), add = TRUE, lty = 2, cex = 0.1)
# plot(ecdf(d0$Y[d0$D == 1]), add = TRUE, col = 2, lty = 2, cex = 0.1)

# Nonparametric control function ------------------------------------------

### Fit RF for control function
cf <- ranger(D ~ Z, data = d1, probability = TRUE)
preds <- predict(cf, data = d1)$predictions
d1$ps <- d1$D - preds[, 1]

### Fit RF with control function prediction and compute RF weights for prediction
rf <- ranger(Y ~ D + ps, data = d1, quantreg = TRUE)

# TRAM --------------------------------------------------------------------

# m0 <- BoxCox(Y ~ ps, data = d1, subset = d1$D == 0, order = 20)
# m1 <- BoxCox(Y ~ ps, data = d1, subset = d1$D == 1, order = 20)
#
# lines(ys, rowMeans(predict(m0, type = "distribution", newdata = d1[,-1], q = ys)), lty = 3)
# lines(ys, rowMeans(predict(m1, type = "distribution", newdata = d1[,-1], q = ys)), lty = 3, col = 2)

# Results -----------------------------------------------------------------

nd0 <- nd1 <- d1
nd0$D <- 0
nd1$D <- 1
p0 <- predict(rf, data = nd0, quantiles = qs <- seq(0, 1, length.out = 1e3), type = "quantiles")
p1 <- predict(rf, data = nd1, quantiles = qs, type = "quantiles")

p0 <- colMeans(p0$predictions)
p1 <- colMeans(p1$predictions)

lines(p0, qs)
lines(p1, qs, col = 2)

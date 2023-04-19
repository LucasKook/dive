# Distributional random forest with control function
# LK March 2023

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("ranger")
library("tram")
devtools::load_all()

# Data --------------------------------------------------------------------

n <- 1e5

# Data under intervention on D (d0) and observational (d)
d0 <- dgp_ex1_cont(n, doD = TRUE)
d1 <- dgp_ex1_cont(n, doD = FALSE)

# plot(ecdf(d0$Y[d0$D == 0]), add = TRUE, lty = 2, cex = 0.1)
# plot(ecdf(d0$Y[d0$D == 1]), add = TRUE, col = 2, lty = 2, cex = 0.1)

oracle_distr <- attr(d1, "odist")

ys <- seq(-13, 13, length.out = 1e3)
F1 <- oracle_distr(ys, d = 1)
F0 <- oracle_distr(ys, d = 0)

plot(ys, F0, type = "l", lty = 2)
lines(ys, F1, type = "l", col = 2, lty = 2)

# Nonparametric control function ------------------------------------------

### Fit RF for control function
cf <- ranger(factor(D) ~ Z, data = d1, probability = TRUE)
preds <- predict(cf, data = d1)$predictions
d1$ps <- d1$D - preds[, 2]

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
p0 <- predict(rf, data = nd0, quantiles = qs <- seq(0, 1, length.out = 1e3), type = "quantiles")$pred
p1 <- predict(rf, data = nd1, quantiles = qs, type = "quantiles")$pred

lines(colMeans(p0), qs)
lines(colMeans(p1), qs, col = 2)

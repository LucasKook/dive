# Distributional random forest with control function
# LK March 2023

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("ranger")
library("tram")
devtools::load_all()

# Data --------------------------------------------------------------------

n <- 1e3

# Data under intervention on D (d0) and observational (d)
# d0 <- dgp_ex1_cont(n, doD = TRUE)
# d1 <- dgp_ex1_cont(n, doD = FALSE)

dgp <- function(n = 1e3, doD = FALSE, cf = rnorm(5)) {
  ### Instrument
  Z <- rt(n, df = 5) # sample(0:1, n, TRUE) # rt(n, df = 5)
  ### Hidden
  H <- rt(n, df = 5)
  ### Treatment
  UD <- runif(n)
  D <- as.numeric(plogis(cf[1] + cf[2] * Z + cf[3] * (1 - doD) * H) >= UD)
  ### Response
  UY <- runif(n)
  Y <- qnorm(UY, mean = cf[4] * D + cf[5] * H, sd = 1 + abs(D + H))
  ### Return
  ret <- data.frame(Y = Y, D = D, Z = Z, H = H)

  ### Compute oracle
  oracle_distr <- Vectorize(\(y, d = 0) {
    integrate(\(x) pnorm(y, mean = d + x, sd = 1 + abs(d + x)) * dt(x, df = 5),
              -20, 20)$value
  }, "y")

  structure(ret, odist = oracle_distr, cf = cf)
}

# d0 <- dgp(n, doD = TRUE)
d1 <- dgp(n, doD = FALSE)
d1t <- dgp(n, doD = FALSE)

# plot(ecdf(d0$Y[d0$D == 0]), add = TRUE, lty = 2, cex = 0.1)
# plot(ecdf(d0$Y[d0$D == 1]), add = TRUE, col = 2, lty = 2, cex = 0.1)

oracle_distr <- attr(d1, "odist")

ys <- seq(min(d1$Y), max(d1$Y), length.out = 1e3)
F1 <- oracle_distr(ys, d = 1)
F0 <- oracle_distr(ys, d = 0)

plot(ys, F0, type = "l", lty = 2)
lines(ys, F1, type = "l", col = 2, lty = 2)

# Nonparametric control function ------------------------------------------

### Fit RF for control function
cf <- ranger(factor(D) ~ Z, data = d1, probability = TRUE)
preds <- predict(cf, data = d1t)$predictions
d1t$V <- preds[, 1] * (d1t$D == 0) + preds[, 2] * (d1t$D == 1)
# d1$V <- d1$D - preds[, 2]
# d1$V <- randomized_pit(preds[, 1], d1$D, trafo = qnorm)

### Fit RF with control function prediction and compute RF weights for prediction
rf <- ranger(Y ~ D + V, data = d1t, quantreg = TRUE)

# TRAM --------------------------------------------------------------------

# m0 <- BoxCox(Y ~ V, data = d1t, subset = d1$D == 0, order = 20)
# m1 <- BoxCox(Y ~ V, data = d1t, subset = d1$D == 1, order = 20)
#
# lines(ys, rowMeans(predict(m0, type = "distribution", newdata = d1t[,-1], q = ys)), lty = 1)
# lines(ys, rowMeans(predict(m1, type = "distribution", newdata = d1t[,-1], q = ys)), lty = 1, col = 2)

# Results -----------------------------------------------------------------

nd0 <- nd1 <- d1t
nd0$D <- 0
nd1$D <- 1
p0 <- predict(rf, data = nd0, quantiles = qs <- seq(0, 1, length.out = 1e3), type = "quantiles")$pred
p1 <- predict(rf, data = nd1, quantiles = qs, type = "quantiles")$pred

lines(colMeans(p0), qs)
lines(colMeans(p1), qs, col = 2)

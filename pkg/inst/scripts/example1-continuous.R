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
  # gD <- (1 - doD) * fH + rlogis(nfine)
  # UD <- dive:::.clip(ecdf(gD)((1 - doD) * H + rlogis(n)))
  UD <- runif(n)
  D <- as.numeric(plogis(Z + (1 - doD) * H) >= UD)
  ### Response
  # gY <- fH + rnorm(nfine)
  # UY <- dive:::.clip(ecdf(gY)(H + rnorm(n)))
  UY <- runif(n)
  Y <- qnorm(UY, mean = D + H, sd = 1 + abs(D + H))
  ### Return
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

# rfk <- function(rf, data) {
#   preds <- predict(rf, data = d, type = "terminalNodes")
#   inbag <- simplify2array(rf$inbag.counts)
#   rfw <- matrix(0, nrow = nrow(d), ncol = nrow(d))
#   for (i in 1:nrow(d)) {
#     for (j in 1:i) {
#       tree_idx <- inbag[i, ] == 0 & inbag[j, ] == 0
#       rfw[i, j] <- sum(preds$predictions[i, tree_idx] ==
#                          preds$predictions[j, tree_idx]) / sum(tree_idx)
#     }
#   }
#   rfw <- rfw + t(rfw - diag(nrow = nrow(d)))
#   rfw / matrix(pmax(colSums(rfw), .Machine$double.eps),
#                ncol = nrow(d), nrow = nrow(d), byrow = TRUE)
# }

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
d0 <- gen_dat(1e5, doD = TRUE)
d1 <- gen_dat(1e5, doD = FALSE)

# plot(ecdf(d0$Y[d0$D == 0]), add = TRUE, lty = 2)
# plot(ecdf(d0$Y[d0$D == 1]), add = TRUE, col = 2, lty = 2)

# Nonparametric control function ------------------------------------------

### Fit RF for control function
cf <- ranger(D ~ Z, data = d1, probability = TRUE)
preds <- predict(cf, data = d1)$predictions
d1$ps <- d1$D - preds[, 1]
# cf <- glm(D ~ Z, data = dZ, family = "binomial")
# preds <- predict(cf, newdata = d, type = "response")
# d$ps <- d$D - preds

### Fit RF with control function prediction and compute RF weights for prediction
rf <- ranger(Y ~ D + ps, data = d1, quantreg = TRUE)
# rfw <- rfk(rf, d)

# TRAM --------------------------------------------------------------------

m0 <- BoxCox(Y | ps ~ 1, data = d1, subset = d1$D == 0, order = 10)
m1 <- BoxCox(Y | ps ~ 1, data = d1, subset = d1$D == 1, order = 10)

lines(ys, rowMeans(predict(m0, type = "distribution", newdata = d[,-1], q = ys)), lty = 3)
lines(ys, rowMeans(predict(m1, type = "distribution", newdata = d[,-1], q = ys)), lty = 3, col = 2)

# Results -----------------------------------------------------------------

idx0 <- which(d$D == 0)
idx1 <- which(d$D == 1)

preds <- predict(rf, data = d1, quantiles = qs <- seq(0, 1, length.out = 1e3),
                 type = "quantiles")
p0 <- colMeans(preds$predictions[idx0,])
p1 <- colMeans(preds$predictions[idx1,])

lines(p0, qs)
lines(p1, qs, col = 2)

# Nonparametric DTE estimation
# LK March 2023

# set.seed(0)

# dependencies ------------------------------------------------------------

library("deeptrafo")
library("tram")

# DGP ---------------------------------------------------------------------

tp <- 2
tnz <- 1
bX <- rep(c(1, 0), c(tnz, tp - tnz))
bY <- rep(c(1, 0), c(tnz, tp - tnz))

gen_dat <- function(n = 5e4, p = tp, nz = tnz, betaX = bX, betaY = bY) {
  X <- matrix(runif(n * p), ncol = p, nrow = n)
  props <- plogis(X %*% betaX)
  D <- as.numeric(props >= runif(n))
  Y <- X %*% betaY + D + sqrt(exp(X %*% betaY)) * rnorm(n)
  structure(data.frame(Y = Y, D = D, X = X), betaX = betaX, betaY = betaY,
            props = props * D + (1 - props) * (1 - D))
}

# Generate data -----------------------------------------------------------

d <- gen_dat()
dp <- gen_dat()
dt <- gen_dat()

# Fit models --------------------------------------------------------------

### Model for untreated
m0 <- BoxCoxNN(Y | X.1 + X.2 ~ 1, data = d[d$D == 0, ])
fit(m0, epochs = 1e3, early_stopping = TRUE)
# m0 <- Lm(Y ~ X.1 + X.2 | X.1 + X.2, data = d[d$D == 0, ])

### Model for treated
m1 <- BoxCoxNN(Y | X.1 + X.2 ~ 1, data = d[d$D == 1, ])
fit(m1, epochs = 1e3, early_stopping = TRUE)
# m1 <- Lm(Y ~ X.1 + X.2 | X.1 + X.2, data = d[d$D == 1, ])

### Predictions
qs <- seq(max(min(d$Y[d$D == 0]), min(d$Y[d$D == 1])),
          min(max(d$Y[d$D == 0]), max(d$Y[d$D == 1])), length.out = 1e2)
p0 <- do.call("cbind", predict(m0, type = "cdf", newdata = dt[dt$D == 0, -1], q = qs))
p1 <- do.call("cbind", predict(m1, type = "cdf", newdata = dt[dt$D == 1, -1], q = qs))
# p0 <- t(predict(m0, newdata = dt[dt$D == 0, -1], q = qs, type = "distribution"))
# p1 <- t(predict(m1, newdata = dt[dt$D == 1, -1], q = qs, type = "distribution"))
#
# pp0 <- t(predict(m0, newdata = dp[dp$D == 0, -1], q = qs, type = "distribution"))
# pp1 <- t(predict(m1, newdata = dp[dp$D == 1, -1], q = qs, type = "distribution"))

### Propensity scores
# mp <- glm(D ~ X.1 + X.2, data = dp, family = "binomial")
# preds <- predict(mp, newdata = d[, -1], type = "response")
# ps <- d$D * preds + (1 - d$D) * (1 - preds)
ps <- attr(dt, "props")

### Influence function
IF <- function(y, preds, Y, prop) {
  inds <- as.numeric(Y <= y)
  (1 / prop) * (inds - preds)
}

IF0 <- sapply(seq_along(qs), \(idx) IF(qs[idx], p0[, idx], dt$Y[dt$D == 0], ps[dt$D == 0]))
IF1 <- sapply(seq_along(qs), \(idx) IF(qs[idx], p1[, idx], dt$Y[dt$D == 1], ps[dt$D == 1]))

# Plot --------------------------------------------------------------------

dd <- \(y, dat) pnorm(y, mean = dat$D + dat$X.1 * attr(dat, "betaY")[1], sd = sqrt(exp(dat$X.1 * attr(dat, "betaY")[1])))
plot(qs, sapply(qs, \(q) mean(dd(q, dt[dt$D == 1,])) - mean(dd(q, dt[dt$D == 0, ]))), col = 2, type = "l")
lines(qs, colMeans(p1) - colMeans(p0), type = "l", col = 3)
lines(qs, colMeans(p1) + colMeans(IF1) - colMeans(p0) - colMeans(IF0), type = "l")
legend("bottomright", c("One-step", "Oracle", "Plug-in"), col = c(1, 2, 3), lwd = 1)

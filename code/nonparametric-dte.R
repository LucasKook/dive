# Nonparametric DTE estimation
# LK March 2023

set.seed(0)

# dependencies ------------------------------------------------------------

library("deeptrafo")

# DGP ---------------------------------------------------------------------

gen_dat <- function(n = 1e4, p = 2, nz = 1, betaX = c(1, 0), betaY = c(-1, 0)) {
  X <- matrix(rnorm(n * p), ncol = p, nrow = n)
  props <- plogis(X %*% betaX)
  D <- as.numeric(props <= runif(n))
  Y <- X %*% betaY + D + rnorm(n)
  structure(data.frame(Y = Y, D = D, X = X), betaX = betaX, betaY = betaY,
            props = props * D + (1 - props) * (1 - D))
}

# Generate data -----------------------------------------------------------

d <- gen_dat()
dp <- gen_dat()
dt <- gen_dat()

# Fit models --------------------------------------------------------------

### Model for untreated
m0 <- BoxCoxNN(Y | te(X.1, X.2) ~ 1, data = d[d$D == 0, ])
fit(m0, epochs = 1e3, early_stopping = TRUE)

### Model for treated
m1 <- BoxCoxNN(Y | te(X.1, X.2) ~ 1, data = d[d$D == 1, ])
fit(m1, epochs = 1e3, early_stopping = TRUE)

### Predictions
qs <- seq(max(min(d$Y[d$D == 0]), min(d$Y[d$D == 1])),
          min(max(d$Y[d$D == 0]), max(d$Y[d$D == 1])), length.out = 1e2)
p0 <- do.call("cbind", predict(m0, type = "cdf", newdata = d[d$D == 0, -1], q = qs))
p1 <- do.call("cbind", predict(m1, type = "cdf", newdata = d[d$D == 1, -1], q = qs))

### Propensity scores
mp <- glm(D ~ X.1 + X.2, data = dp, family = "binomial")
preds <- predict(mp, newdata = d[, -1], type = "response")
ps <- d$D * preds + (1 - d$D) * (1 - preds)

### Influence function
IF <- function(y, preds, Y, prop) {
  inds <- as.numeric(Y <= y)
  (1 / prop) * (inds - preds)
}

IF0 <- sapply(seq_along(qs), \(idx) IF(qs[idx], p0[, idx], d$Y[d$D == 0], ps[d$D == 0]))
IF1 <- sapply(seq_along(qs), \(idx) IF(qs[idx], p1[, idx], d$Y[d$D == 1], ps[d$D == 1]))

# Plot --------------------------------------------------------------------

plot(qs, colMeans(p1) - colMeans(p0) + colMeans(IF1) - colMeans(IF0), type = "s")
dd <- \(y, dat) pnorm(y, mean = dat$D + dat$X.1 * attr(dat, "betaY")[1])
lines(qs, sapply(qs, \(q) mean(dd(q, d[d$D == 1,])) - mean(dd(q, d[d$D == 0, ]))), col = 2)
lines(qs, colMeans(p1) - colMeans(p0), type = "s", col = 3)

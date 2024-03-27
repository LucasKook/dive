### Example discrete breaks independence
### LK 2024

set.seed(12)

# Data --------------------------------------------------------------------

dgp <- function(n = 1e4) {
  z <- as.numeric(runif(n) > 0.5)
  x <- as.numeric(z + rlogis(n) > 0)
  y <- as.numeric(-1 + x + rlogis(n) > 0)
  data.frame(y = y, x = x, z = z)
}

# correctly specified by logistic regression, no hiddens, but
# correlation test rejects for PIT vs Z
d <- dgp()
m <- glm(y ~ x, data = d, family = "binomial")
pit <- (d$y == 0) * (1 - predict(m, type = "response")) + (d$y == 1)
rpit <- randomized_pit(1 - predict(m, type = "response"), d$y)
cor.test(pit, d$z)
cor.test(rpit, d$z)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot(pit ~ z, data = d)
plot(rpit ~ z, data = d)
par(opar)

coef(m)
indep_iv(y ~ x, ~ 0 + z, data = d, method = "DIVE")

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")
library("coin")

# Run ---------------------------------------------------------------------

dive <- BoxCoxDA(y | x ~ 1, data = d, anchor = ~ z, loss = "indep",
                 optimizer = optimizer_adam(0.01))
fit(dive, epochs = 1e4)

ipit <- predict(dive, type = "cdf")
plot(ecdf(ipit))
abline(0, 1)
cor.test(ipit, d$z)

### EX WITH NP
### Example discrete breaks independence
### LK 2024

set.seed(12)

# Data --------------------------------------------------------------------

dgp <- function(n = 1e2) {
  z <- as.numeric(runif(n) > 0.5)
  x <- as.numeric(z + rlogis(n) > 0)
  y <- as.numeric(-1 + x + rlogis(n) > 0)
  data.frame(y = y, x = x, z = z)
}

obj <- \(b, Y, X, E, lambda = 1) {
  Xb <- X %*% b
  p0 <- 1 - (Xb >= 0) * (Xb <= 1) * Xb
  R <- randomized_pit(p0, Y)
  cmv <- mean((R - ecdf(R)(R))^2)
  cor <- mean(fitted(lm(R ~ E))^2)
  max(lambda * cor, cmv)
}

d <- dgp()
ps <- seq(0, 1, length.out = 1e2)
bs <- data.frame(expand.grid(p0 = ps, p1 = ps))

loss <- apply(bs, 1, \(tb) obj(b = tb, Y = d$y, X = cbind(1, d$x), E = d$z, lambda = 0))
obj(c(plogis(-1), 0.5 - plogis(-1)), Y = d$y, X = cbind(1, d$x), E = d$z, lambda = 0)
min(loss)

c(plogis(-1), plogis(0))
bs[which.min(loss), ]

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

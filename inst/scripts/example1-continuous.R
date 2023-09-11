# Distributional random forest with control function
# LK March 2023

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("tidyverse")
library("ranger")
library("tram")
devtools::load_all()

# Data --------------------------------------------------------------------

n <- 1e3

# Data under intervention on D (d0) and observational (d)

# tcf <- c(-1.43743764887867, -0.798806347083827, -1.58585991850206)
# dgp <- function(n = 1e3, doD = FALSE, cf = rnorm(3)) {
#   ### Instrument
#   Z <- rt(n, df = 5)
#   # Z <- sample(0:1, n, TRUE)
#   if (doD) cop <- copula::indepCopula(2) else cop <- copula::claytonCopula(-0.5, 2)
#   U <- copula::rCopula(cop, n = n)
#   UD <- U[, 1]
#   UY <- U[, 2]
#   ### Treatment
#   D <- as.numeric(plogis(cf[1] + cf[2] * Z) >= UD)
#   ### Response
#   Y <- qnorm(UY, mean = cf[3] * D, sd = 1 + abs(D))
#   ### Return
#   ret <- data.frame(Y = Y, D = D, Z = Z, U = U)
#   structure(ret, cf = cf)
# }

tcf <- c(-1.43743764887867, -0.798806347083827, -0.197635299911542,
         -1.58585991850206, 1.81201327858338)
dgp <- function(n = 1e3, doD = FALSE, cf = rnorm(5)) {
  ### Instrument
  Z <- rt(n, df = 5)
  # Z <- sample(0:1, n, TRUE)
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
    integrate(\(x) pnorm(y, mean = cf[4] * d + cf[5] * x, sd = 1 + abs(d + x)) *
                dt(x, df = 5), -20, 20)$value
  }, "y")

  structure(ret, odist = oracle_distr, cf = cf)
}

d0 <- dgp(10 * n, doD = TRUE, cf = tcf)

# plot(ecdf(d0$Y[d0$D == 0]), lty = 2, cex = 0.1)
# plot(ecdf(d0$Y[d0$D == 1]), add = TRUE, col = 2, lty = 2, cex = 0.1)
# oracle_distr <- attr(d1, "odist")
# ys <- seq(min(d1$Y), max(d1$Y), length.out = 1e3)
# F1 <- oracle_distr(ys, d = 1)
# F0 <- oracle_distr(ys, d = 0)
# lines(ys, F0, type = "l", lty = 2)
# lines(ys, F1, type = "l", col = 2, lty = 2)

nsim <- 50
pb <- txtProgressBar(0, nsim, style = 3)
res <- lapply(1:nsim, \(iter) {
  setTxtProgressBar(pb, iter)

  ### Generate data
  d1 <- d1t <- dgp(n, doD = FALSE, cf = tcf)
  d1t <- dgp(n, doD = FALSE, cf = tcf)

  ### Fit RF for control function
  cf <- ranger(factor(D) ~ Z, data = d1, probability = TRUE)
  preds <- predict(cf, data = d1t)$predictions
  d1t$V <- preds[, 1] * (d1t$D == 0) + preds[, 2] * (d1t$D == 1)

  ### Fit RF with control function prediction and compute RF weights for prediction
  rf <- ranger(Y ~ D + V, data = d1t, quantreg = TRUE)

  ### Compute CDFs
  nd0 <- nd1 <- d1t
  nd0$D <- 0
  nd1$D <- 1
  p0 <- predict(rf, data = nd0, quantiles = qs <- seq(0, 1, length.out = 3e2),
                type = "quantiles")$pred
  p1 <- predict(rf, data = nd1, quantiles = qs, type = "quantiles")$pred

  data.frame(p0 = colMeans(p0), p1 = colMeans(p1), q = qs)
})

pdat <- res %>%
  bind_rows(.id = "iter") %>%
  pivot_longer(names_to = "group", values_to = "y", p0:p1)

mdat <- pdat %>% group_by(q, group) %>% summarise(y = mean(y)) %>% ungroup()

ggplot(pdat, aes(x = y, y = q, color = group, group = interaction(iter, group))) +
  geom_line(alpha = 0.1) +
  geom_line(aes(group = group), data = mdat, lwd = 0.9) +
  stat_ecdf(inherit.aes = FALSE, aes(x = Y, color = "p0"), data = d0[d0$D == 0, ],
            lty = 2) +
  stat_ecdf(inherit.aes = FALSE, aes(x = Y, color = "p1"), data = d0[d0$D == 1, ],
            lty = 2) +
  theme_bw() +
  scale_color_manual(values = c("p0" = "darkblue", p1 = "darkred"),
                     labels = c("p0" = "D = 0", "p1" = "D = 1")) +
  labs(color = element_blank(), title = "Conditional DGP")

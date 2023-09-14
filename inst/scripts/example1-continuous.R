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
dgp <- function(n = 1e3, doD = FALSE, cf = rnorm(5)) {
  ### Instrument
  Z <- rt(n, df = 5)
  # Z <- sample(0:1, n, TRUE)
  ### Hidden
  H <- rt(n, df = 5)
  ### Treatment
  UD <- runif(n)
  D <- as.numeric(plogis(cf[1] + cf[2] * Z + cf[3] * (1 - doD) * H) >= UD)
  ### Covariate
  X <- rnorm(n)
  ### Response
  NY <- rnorm(n)
  Y <- qchisq(pnorm(0.5 * X + cf[4] * D + cf[5] * H +
                      (1 + abs(0.5 * D + 0.3 * H + 0.3 * X)) * NY), df = 10)
  ### Return
  ret <- data.frame(Y = Y, D = D, X = X, Z = Z, H = H)
  structure(ret, cf = cf)
}

### Generate large interventional data set
tcf <- c(-1.43, -0.79, -1.19, -1.58, 0.81)
d0 <- dgp(10 * n, doD = TRUE, cf = tcf)

# Simulation --------------------------------------------------------------

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

  ### Fit RF with ctrl fn prediction and compute RF weights for prediction
  rf <- ranger(Y ~ D + X + V, data = d1t, quantreg = TRUE)

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
  labs(color = element_blank())

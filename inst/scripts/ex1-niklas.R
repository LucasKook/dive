# Distributional random forest with control function
# LK March 2023

set.seed(2410)

# Dependencies ------------------------------------------------------------

library("tidyverse")
library("ranger")
library("tram")
devtools::load_all()

# Data --------------------------------------------------------------------

n <- 3e3

bpath <- file.path("inst", "figures", Sys.Date())
if (!dir.exists(bpath))
  dir.create(bpath, recursive = TRUE, showWarnings = FALSE)
fname <- paste0("ex1_n", n)

pz <- \(x, location = 0) pnorm(x, mean = location) # plogis

dgp <- function(n = 1e3) {

  # Exogenous variables
  H <- rnorm(n)
  Z <- runif(n, 0, 0.5)
  eps <- rnorm(n)*0.05

  # Endogenous variables
  D0 <- sapply(Z, function(p) rbinom(1, 1, p))
  D1 <- sapply(Z, function(p) rbinom(1, 1, 2*p))
  D <- D0
  D[H > 0] <- D1[H > 0]

  # Two possible responses
  Y <- H + D + eps
  # Y <- H*(D + 1) + eps

  p1 <- Z + D * Z
  # ctrl <- D * p1 + (1 - D) * (1 - p1)
  ctrl <- D * (p1 > 0.5) - (1 - D) * (p1 <= 0.5)
  # ctrl <- (1 - p1)^(1 - D)

  ### Return
  data.frame(Y = Y, D = D, Z = Z, H = H, ctrl = ctrl)
}

### Oracle conditional distribution
# orc <- Vectorize(\(y, d) integrate(\(h) pnorm(y, mean = h * (1 + d), sd = 0.05) * dnorm(h), -Inf, Inf)$value, "y")
orc <- Vectorize(\(y, d) integrate(\(h) pnorm(y, mean = h + d, sd = 0.05) * dnorm(h), -Inf, Inf)$value, "y")

# ys <- seq(-7, 7, length.out = 1e3)
# plot(ys, orc(ys, d = 0), type = "l", col = "darkred")
# lines(ys, orc(ys, d = 1), type = "l", col = "darkblue")
# dd <- dgp(n = n)
# m <- Lm(Y ~ H * D, data = dd)
# dd0 <- dd1 <- dd
# dd0$D <- 0
# dd1$D <- 1
# p0 <- predict(m, type = "distribution", newdata = dd0[, -1], q = ys)
# p1 <- predict(m, type = "distribution", newdata = dd1[, -1], q = ys)
# lines(ys, rowMeans(p0), lty = 2, col = "darkred")
# lines(ys, rowMeans(p1), lty = 2, col = "darkblue")
#
# rf0 <- ranger(Y ~ H, data = dd[dd$D == 0, ], quantreg = TRUE, num.trees = 2e3)
# rf1 <- ranger(Y ~ H, data = dd[dd$D == 1, ], quantreg = TRUE, num.trees = 2e3)
# qs <- seq(0, 1, length.out = 3e2)
# nd0 <- nd1 <- dd; nd0$D <- 0; nd1$D <- 1
# p0r <- predict(rf0, type = "quantiles", data = nd0, quantiles = qs)$pred
# p1r <- predict(rf1, type = "quantiles", data = nd1, quantiles = qs)$pred
# cdf0 <- Vectorize(\(y) mean(apply(p0r, 1, \(x) qs[which.min(abs(x - y))])))
# cdf1 <- Vectorize(\(y) mean(apply(p1r, 1, \(x) qs[which.min(abs(x - y))])))
# lines(ys, cdf0(ys), lty = 3, col = "darkred", lwd = 1.5)
# lines(ys, cdf1(ys), lty = 3, col = "darkblue", lwd = 1.5)
#
# legend("topleft", c("ranger", "tram", "oracle"), lty = 3:1, bty = "n", lwd = 1.5)
# legend("left", c("D = 0", "D = 1"), col = c("darkred", "darkblue"), bty = "n", lwd = 1.5)

# Simulation --------------------------------------------------------------

nsim <- 10
pb <- txtProgressBar(0, nsim, style = 3)
res <- lapply(1:nsim, \(iter) {
  setTxtProgressBar(pb, iter)

  ### Generate data
  d1 <- dgp(n)
  d1t <- dgp(n)

  ### In-sample CTRL
  # cf <- ranger(factor(D) ~ Z, data = d1t, probability = TRUE)
  # preds <- predict(cf, data = d1t)$predictions
  # d1t$iV <- preds[, 1]^(1 - d1t$D)

  ### Out-of-sample CTRL
  # cf <- ranger(factor(D) ~ Z, data = d1, probability = TRUE)
  # preds <- predict(cf, data = d1t)$predictions
  # d1t$V <- preds[, 1]^(1 - d1t$D)

  d1t$iH <- as.numeric(d1t$H <= 0)

  d1t0 <- d1t[d1t$D == 0, ]
  d1t1 <- d1t[d1t$D == 1, ]

  # ### RF + CTRL out-of-sample
  # CTRL0 <- ranger(Y ~ V, data = d1t0, quantreg = TRUE)
  # CTRL1 <- ranger(Y ~ V, data = d1t1, quantreg = TRUE)
  #
  # ### RF + CTRL out-of-sample
  # INSA0 <- ranger(Y ~ iV, data = d1t0, quantreg = TRUE)
  # INSA1 <- ranger(Y ~ iV, data = d1t1, quantreg = TRUE)

  ### RF + oracle CTRL
  ORAC0 <- ranger(Y ~ ctrl, data = d1t0, quantreg = TRUE)
  ORAC1 <- ranger(Y ~ ctrl, data = d1t1, quantreg = TRUE)

  ### RF + confounder
  CONF0 <- ranger(Y ~ iH, data = d1t0, quantreg = TRUE)
  CONF1 <- ranger(Y ~ iH, data = d1t1, quantreg = TRUE)

  ### Compute CDFs forests
  all <- list(
    # "CTRL0" = CTRL0,
    # "CTRL1" = CTRL1,
    "ORAC0" = ORAC0,
    "ORAC1" = ORAC1,
    "CONF0" = CONF0,
    "CONF1" = CONF1#,
    # "INSA0" = INSA0,
    # "INSA1" = INSA1
  )
  qs <- seq(0, 1, length.out = 3e2)
  ys <- quantile(d1t$Y, probs = qs)
  lapply(seq_along(all), \(idx) {
    rf <- all[[idx]]
    ps <- predict(rf, data = d1t, quantiles = qs, type = "quantiles")$pred
    cdf <- Vectorize(\(y) mean(apply(ps, 1, \(x) qs[which.min(abs(x - y))])))
    data.frame(cdf = cdf(ys), y = ys, method = names(all)[idx], q = qs,
               group = ifelse(grepl("0", names(all)[idx]), "p0", "p1"))
  }) |> bind_rows()
})

pdat <- res %>%
  bind_rows(.id = "iter") |>
  mutate(method = str_remove(method, "[0-9]"))

mdat <- pdat %>% group_by(q, group, method) %>% summarise(cdf = mean(cdf), y = mean(y)) %>% ungroup()

ggplot(pdat, aes(x = y, y = cdf, color = group, group = interaction(iter, group, method))) +
  geom_line(alpha = 0.3) +
  facet_grid(~ method) +
  geom_line(aes(group = group), data = mdat, lwd = 0.9) +
  stat_function(fun = ~ orc(.x, d = 0), aes(color = "p0"), inherit.aes = FALSE, linetype = 2, linewidth = 1) +
  stat_function(fun = ~ orc(.x, d = 1), aes(color = "p1"), inherit.aes = FALSE, linetype = 2, linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = c("p0" = "darkblue", p1 = "darkred"),
                     labels = c("p0" = "D = 0", "p1" = "D = 1")) +
  labs(color = element_blank(), subtitle = fname, y = "CDF",
       caption = "Dashed lines: Oracle interventional marginal distribution")

ggsave(file.path(bpath, paste0(fname, ".pdf")))

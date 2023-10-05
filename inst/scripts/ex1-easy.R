# Distributional random forest with control function
# LK March 2023

set.seed(2410)

# Dependencies ------------------------------------------------------------

library("tidyverse")
library("ranger")
library("tram")
devtools::load_all()

# Data --------------------------------------------------------------------

n <- 1e3

bpath <- file.path("inst", "figures", Sys.Date())
if (!dir.exists(bpath))
  dir.create(bpath, recursive = TRUE, showWarnings = FALSE)
fname <- paste0("ex1_n", n)

dgp <- function(n = 1e3) {
  ### Instrument
  Z <- rt(n, df = 5) / 0.8
  ### Hidden
  H <- sample(c(-1, 1), n, TRUE)
  ### Treatment
  UD <- runif(n)
  D <- as.numeric(plogis(Z + H) <= UD)
  ctrl <- \(z) (plogis(z - 1) + plogis(z + 1))
  ### Response
  NY <- rlogis(n)
  Y <- 2 * (1 + D) * H + NY
  ### Return
  data.frame(Y = Y, D = D, Z = Z, H = H, ctrl = ctrl(Z)^(1 - D))
}

### Oracle conditional distribution
orc <- Vectorize(\(y, d) (plogis(y, location = - 2 * (1 + d)) +
                            plogis(y, location = 2 * (1 + d))) / 2, "y")

# ys <- seq(-7, 7, length.out = 1e3)
# plot(ys, orc(ys, d = 0), type = "l", col = "darkred")
# lines(ys, orc(ys, d = 1), type = "l", col = "darkblue")
# dd <- dgp(n = 3e2)
# m <- Survreg(Y ~ H * D, data = dd, dist = "logistic")
# dd0 <- dd1 <- dd
# dd0$D <- 0
# dd1$D <- 1
# p0 <- predict(m, type = "distribution", newdata = dd0[, -1], q = ys)
# p1 <- predict(m, type = "distribution", newdata = dd1[, -1], q = ys)
# lines(ys, rowMeans(p0), lty = 2, col = "darkred")
# lines(ys, rowMeans(p1), lty = 2, col = "darkblue")
#
# rf <- ranger(Y ~ D + H, data = dd, quantreg = TRUE)
# qs <- seq(0.001, 0.999, length.out = 3e2)
# nd0 <- nd1 <- dd; nd0$D <- 0; nd1$D <- 1
# p0r <- predict(rf, type = "quantiles", data = nd0, quantiles = qs)$pred
# p1r <- predict(rf, type = "quantiles", data = nd1, quantiles = qs)$pred
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
  cf <- ranger(factor(D) ~ Z, data = d1t, probability = TRUE)
  preds <- predict(cf, data = d1t)$predictions
  d1t$iV <- preds[, 1]^(1 - d1t$D)

  ### Out-of-sample CTRL
  cf <- ranger(factor(D) ~ Z, data = d1, probability = TRUE)
  preds <- predict(cf, data = d1t)$predictions
  d1t$V <- preds[, 1]^(1 - d1t$D)

  ### RF + CTRL out-of-sample
  CTRL <- ranger(Y ~ D + V, data = d1t, quantreg = TRUE)

  ### RF + CTRL out-of-sample
  INSA <- ranger(Y ~ D + iV, data = d1t, quantreg = TRUE)

  ### RF + oracle CTRL
  ORAC <- ranger(Y ~ D + ctrl, data = d1t, quantreg = TRUE)

  ### RF + confounder
  CONF <- ranger(Y ~ D + H, data = d1t, quantreg = TRUE)

  ### Compute CDFs forests
  all <- list("CTRL" = CTRL, "ORAC" = ORAC, "CONF" = CONF, "INSA" = INSA)
  nd0 <- nd1 <- d1t
  nd0$D <- 0
  nd1$D <- 1
  qs <- seq(0.001, 0.999, length.out = 3e2)
  ys <- quantile(d1t$Y, probs = qs)
  lapply(seq_along(all), \(idx) {
    rf <- all[[idx]]
    p0 <- predict(rf, data = nd0, quantiles = qs, type = "quantiles")$pred
    p1 <- predict(rf, data = nd1, quantiles = qs, type = "quantiles")$pred
    cdf0 <- Vectorize(\(y) mean(apply(p0, 1, \(x) qs[which.min(abs(x - y))])))
    cdf1 <- Vectorize(\(y) mean(apply(p1, 1, \(x) qs[which.min(abs(x - y))])))
    data.frame(p0 = cdf0(ys), p1 = cdf1(ys), y = ys, method = names(all)[idx], q = qs)
  }) |> bind_rows()
})

pdat <- res %>%
  bind_rows(.id = "iter") %>%
  pivot_longer(names_to = "group", values_to = "cdf", p0:p1)

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

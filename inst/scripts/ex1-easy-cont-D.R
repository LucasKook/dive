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
  Z <- rnorm(n)
  ### Hidden
  H <- rnorm(n)
  ### Treatment
  D <- Z + H + rnorm(n)
  ctrl <- Vectorize(\(d, z) integrate(\(h) pnorm(d, mean = z + h) * dnorm(h),
                                      -Inf, Inf)$value)
  ### Response
  Y <- D/2 + (2 + D) * H + rnorm(n)
  ### Return
  data.frame(Y = Y, D = D, Z = Z, H = H, ctrl = ctrl(D, Z))
}

# res <- replicate(1e2, {
#   dd <- dgp(n = 1e3)
#   c(oracle = unname(coef(lm(Y ~ D + H, data = dd))["D"]),
#     control = unname(coef(lm(Y ~ D + residuals(lm(D ~ Z, data = dd)), data = dd))["D"]),
#     newey = unname(coef(lm(Y ~ D + ctrl, data = dd))["D"]))
# })
# boxplot(t(res))

### Oracle conditional distribution
orc <- Vectorize(\(y, d) {
  integrate(\(h) pnorm(y, mean = d/2 + (2 + d) * h) * dnorm(h), -Inf, Inf)$value
})

# ys <- seq(-15, 15, length.out = 1e3)
# plot(ys, orc(ys, d = -3), type = "l", col = "darkred")
# lines(ys, orc(ys, d = 3), type = "l", col = "darkblue")
# dd <- dgp(n = 1e3)
# m <- Lm(Y ~ H * D, data = dd)
# dd0 <- dd1 <- dd
# dd0$D <- -3
# dd1$D <- 3
# p0 <- predict(m, type = "distribution", newdata = dd0[, -1], q = ys)
# p1 <- predict(m, type = "distribution", newdata = dd1[, -1], q = ys)
# lines(ys, rowMeans(p0), lty = 2, col = "darkred")
# lines(ys, rowMeans(p1), lty = 2, col = "darkblue")
#
# rf <- ranger(Y ~ D + H, data = dd, quantreg = TRUE, num.trees = 2e3, mtry = 2)
# qs <- seq(0.001, 0.999, length.out = 1e3)
# nd0 <- nd1 <- dd; nd0$D <- -3; nd1$D <- 3
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

gcdf <- function(preds, q = seq(0.001, 0.999, length.out = 3e2)) {
  Vectorize(\(y) mean(apply(preds, 1, \(x) q[which.min(abs(x - y))])))
}

nsim <- 1
pb <- txtProgressBar(0, nsim, style = 3)
res <- lapply(1:nsim, \(iter) {
  setTxtProgressBar(pb, iter)
  qs <- seq(0.001, 0.999, length.out = 3e2)

  ### Generate data
  d1 <- dgp(n)
  d1t <- dgp(n)

  ### In-sample CTRL
  # cf <- ranger(D ~ Z, data = d1t, quantreg = TRUE, num.trees = 2e3)
  # preds <- predict(cf, data = d1t, type = "quantiles", quantiles = qs)$predictions
  # d1t$iV <- gcdf(preds)(d1t$D)
  ism <- Lm(D ~ Z, data = d1t)
  d1t$iV <- predict(ism, type = "distribution")

  ### Out-of-sample CTRL
  # cf <- ranger(D ~ Z, data = d1, quantreg = TRUE, num.trees = 2e3)
  # preds <- predict(cf, data = d1t, type = "quantiles", quantiles = qs)$predictions
  # d1t$V <- gcdf(preds)(d1t$D)
  ism <- Lm(D ~ Z, data = d1)
  d1t$V <- predict(ism, newdata = d1t, type = "distribution")

  ### RF + CTRL out-of-sample
  CTRL <- ranger(Y ~ D + V, data = d1t, quantreg = TRUE, mtry = 2, num.trees = 2e3)

  ### RF + CTRL out-of-sample
  INSA <- ranger(Y ~ D + iV, data = d1t, quantreg = TRUE, mtry = 2, num.trees = 2e3)

  ### RF + oracle CTRL
  ORAC <- ranger(Y ~ D + ctrl, data = d1t, quantreg = TRUE, mtry = 2, num.trees = 2e3)

  ### RF + confounder
  CONF <- ranger(Y ~ D + H, data = d1t, quantreg = TRUE, mtry = 2, num.trees = 2e3)

  ### Compute CDFs forests
  all <- list("CTRL" = CTRL, "ORAC" = ORAC, "CONF" = CONF, "INSA" = INSA)
  ys <- quantile(d1t$Y, probs = qs)
  lapply(seq_along(all), \(idx) {
    rf <- all[[idx]]
    dd1 <- dd0 <- d1t; dd1$D <- 3; dd0$D <- -3
    ps0 <- predict(rf, data = dd0, quantiles = qs, type = "quantiles")$pred
    ps1 <- predict(rf, data = dd1, quantiles = qs, type = "quantiles")$pred
    data.frame(cdf0 = gcdf(ps0)(ys), cdf1 = gcdf(ps1)(ys), y = ys,
               method = names(all)[idx], q = qs)
  }) |> bind_rows()
})

pdat <- res %>%
  bind_rows(.id = "iter") |>
  pivot_longer(cdf0:cdf1, names_to = "group", values_to = "cdf")

mdat <- pdat %>% group_by(q, group, method) %>% summarise(cdf = mean(cdf), y = mean(y)) %>% ungroup()

ggplot(pdat, aes(x = y, y = cdf, color = group, group = interaction(iter, group, method))) +
  geom_line(alpha = 0.3) +
  facet_grid(~ method) +
  geom_line(aes(group = group), data = mdat, lwd = 0.9) +
  stat_function(fun = ~ orc(.x, d = -3), aes(color = "cdf0"), inherit.aes = FALSE, linetype = 2, linewidth = 1) +
  stat_function(fun = ~ orc(.x, d = 3), aes(color = "cdf1"), inherit.aes = FALSE, linetype = 2, linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = c("cdf0" = "darkblue", "cdf1" = "darkred"),
                     labels = c("cdf0" = "D = 0", "cdf1" = "D = 1")) +
  labs(color = element_blank(), subtitle = fname, y = "CDF",
       caption = "Dashed lines: Oracle interventional marginal distribution")

ggsave(file.path(bpath, paste0(fname, ".pdf")))

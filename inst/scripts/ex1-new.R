# Distributional random forest with control function
# LK March 2023

set.seed(2410)
tcf <- runif(5, min = 1, max = 2) * sample(c(-1, 1), 5, TRUE)

# Dependencies ------------------------------------------------------------

library("tidyverse")
library("ranger")
library("tram")
devtools::load_all()

# Data --------------------------------------------------------------------

n <- 1e5

bpath <- file.path("inst", "figures", Sys.Date())
if (!dir.exists(bpath))
  dir.create(bpath, recursive = TRUE, showWarnings = FALSE)
fname <- paste0("ex1_n", n)

.clip <- function(x, lower = 0.0001, upper = 0.999) {
  x[x < lower] <- lower
  x[x > upper] <- upper
  x
}

dgp <- function(n = 1e3, cf = rnorm(5)) {
  ### Instrument
  Z <- rt(n, df = 5) / 0.8
  # Z <- sample(0:1, n, TRUE)
  ### Hidden
  H <- rt(n, df = 5) / 0.8
  ### Treatment
  UD <- runif(n)
  D <- as.numeric(plogis(cf[1] + cf[2] * Z + cf[3] * H) <= UD)
  ctrl <- Vectorize(\(z) integrate(\(u) plogis(cf[1] + cf[2] * z + cf[3] * u) *
                                     dt(u, df = 5), lower = -Inf, upper = Inf)$value)
  ### Response
  NY <- rlogis(n)
  tshift <- cf[4] * D + cf[5] * H - 2
  Y <- qchisq(.clip(plogis(tshift + NY)), df = 10)
  ### Return
  ret <- data.frame(Y = Y, D = D, Z = Z, H = H, ctrl = ctrl(Z)^(1 - D))
  structure(ret, cf = cf)
}

### Oracle conditional distribution
orc <- Vectorize(\(u, d = 1) integrate(\(h) plogis(qlogis(pchisq(
  u, df = 10)) - d * tcf[4] - tcf[5] * h + 2) * dt(h, df = 5), -Inf, u)$value,
  vectorize.args = "u")

# Simulation --------------------------------------------------------------

nsim <- 1
pb <- txtProgressBar(0, nsim, style = 3)
res <- lapply(1:nsim, \(iter) {
  setTxtProgressBar(pb, iter)

  ### Generate data
  d1 <- dgp(n, cf = tcf)
  d1t <- dgp(n, cf = tcf)

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
  lapply(seq_along(all), \(idx) {
    rf <- all[[idx]]
    p0 <- predict(rf, data = nd0, quantiles = qs, type = "quantiles")$pred
    p1 <- predict(rf, data = nd1, quantiles = qs, type = "quantiles")$pred
    data.frame(p0 = colMeans(p0), p1 = colMeans(p1), q = qs,
               method = names(all)[idx])
  }) |> bind_rows()
})

pdat <- res %>%
  bind_rows(.id = "iter") %>%
  pivot_longer(names_to = "group", values_to = "y", p0:p1)

mdat <- pdat %>% group_by(q, group, method) %>% summarise(y = mean(y)) %>% ungroup()

ggplot(pdat, aes(x = y, y = q, color = group, group = interaction(iter, group, method))) +
  geom_line(alpha = 0.1) +
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

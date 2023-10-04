# Distributional random forest with control function
# LK March 2023

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, character(0)))
  args <- c(0, 0, 0)
args <- as.numeric(args)

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("tidyverse")
library("ranger")
library("tram")
devtools::load_all()

# Data --------------------------------------------------------------------

n <- 3e3
use_oracle_ctrl <- as.logical(args[1])
scale_effect <- as.logical(args[2])
split_sample <- as.logical(args[3])

bpath <- file.path("inst", "figures", Sys.Date())
if (!dir.exists(bpath))
  dir.create(bpath, recursive = TRUE, showWarnings = FALSE)
fname <- paste0("ex1-cont_use-oracle-ctrl-", use_oracle_ctrl, "_scale-effect-",
                scale_effect, "_split-sample-", split_sample, "_n-", n)

# Data under intervention on D (d0) and observational (d)
dgp <- function(n = 1e3, doD = FALSE, cf = rnorm(5), scale = scale_effect) {
  ### Instrument
  Z <- rt(n, df = 5)
  # Z <- sample(0:1, n, TRUE)
  ### Hidden
  H <- rt(n, df = 5)
  ### Treatment
  UD <- runif(n)
  lp <- plogis(cf[1] + cf[2] * Z + cf[3] * (1 - doD) * H)
  D <- as.numeric(lp <= UD)
  ctrl <- lp^(1-D) - D * lp
  ### Covariate
  X <- rnorm(n)
  ### Response
  NY <- rlogis(n)
  tshift <- 0.5 * X + cf[4] * D + cf[5] * H
  tscale <- (1 + abs(0.5 * D + 0.3 * H + 0.3 * X))^scale
  Y <- qchisq(plogis(tshift + tscale * NY), df = 10)
  ### Return
  ret <- data.frame(Y = Y, D = D, X = X, Z = Z, H = H, ctrl = ctrl)
  structure(ret, cf = cf)
}

### Generate large interventional data set
tcf <- c(-1.43, -0.79, -1.19, -1.58, 0.81)
d0 <- dgp(10 * n, doD = TRUE, cf = tcf)

### Interventional data fitted with Colr model gives correct coef/cdf
# m <- Colr(Y ~ D + X + H, data = d0, prob = c(0.001, 0.999), order = 10)
# m0 <- Colr(Y | D ~ 1, data = d0, prob = c(0.001, 0.999), order = 10)
# plot(m0, which = "distribution")
# lines(ecdf(d0$Y[d0$D == 0]), col = 2)
# lines(ecdf(d0$Y[d0$D == 1]), col = 2)

# Simulation --------------------------------------------------------------

nsim <- 50
pb <- txtProgressBar(0, nsim, style = 3)
res <- lapply(1:nsim, \(iter) {
  setTxtProgressBar(pb, iter)

  ### Generate data
  d1 <- dgp(n, doD = FALSE, cf = tcf)
  d1t <- if (split_sample) dgp(n, doD = FALSE, cf = tcf) else d1

  ### Fit RF for control function
  if (use_oracle_ctrl) {
    d1t$V <- d1t$ctrl
  } else {
    cf <- ranger(factor(D) ~ Z, data = d1, probability = TRUE)
    preds <- predict(cf, data = d1t)$predictions
    d1t$V <- preds[, 1]^(1 - d1t$D) - preds[, 1]^d1t$D
  }

  # plot(V ~ ctrl, data = d1t)
  # abline(lm(V ~ ctrl, data = d1t))

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
  labs(color = element_blank(), subtitle = fname,
       caption = "Dashed lines: Oracle interventional marginal distribution")

# Save --------------------------------------------------------------------

out <- pdat |>
  mutate(use_oracle_ctrl = use_oracle_ctrl,
         scale_effect = scale_effect,
         split_sample = split_sample)

ggsave(file.path(bpath, paste0(fname, ".pdf")))
write_csv(out, file.path(bpath, paste0(fname, ".csv")))

### Vis all
if (FALSE) {
  pdat <- tibble(file = list.files(bpath, "*.csv", full.names = TRUE)) |>
    mutate(dat = map(file, ~ read_csv(.x, show_col_types = FALSE))) |>
    unnest(dat)
  mdat <- pdat %>%
    group_by(q, group, use_oracle_ctrl, scale_effect, split_sample) %>%
    summarise(y = mean(y)) %>% ungroup()

  ggplot(pdat, aes(x = y, y = q, color = group, group = interaction(iter, group))) +
    facet_grid(use_oracle_ctrl ~ split_sample + scale_effect, labeller = label_both) +
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

  ggsave(file.path(bpath, "ex1-all.pdf"))
}

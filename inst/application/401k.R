### 401(k) DIVE
### LK 2024

set.seed(12)

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Data --------------------------------------------------------------------

d401k <- read_csv("inst/data/401k.csv")

dat <- data.frame(
  y = d401k$net_tfa / 1e3,
  d = factor(d401k$p401),
  z = d401k$e401
)[sample.int(nrow(d401k), 1e3), ]

idx <- which(dat$y == 0)
dat$y[idx] <- runif(length(idx), max(dat$y[dat$y < 0]), min(dat$y[dat$y > 0]))

# FUNs --------------------------------------------------------------------

fit_adaptive <- function(
    args, epochs, max_iter = 5, stepsize = 2, alpha = 0.1, ws = NULL, ...
) {
  for (iter in seq_len(max_iter)) {
    mod <- do.call("ColrDA", c(args, list(tf_seed = iter)))
    if (!is.null(ws))
      set_weights(mod$model, ws)
    browser()
    # plot(args$data$y, predict(mod, type = "cdf"), col = args$data$d)
    fit(mod, epochs = epochs, ...)
    iPIT <- predict(mod, type = "cdf")
    unif <- ks.test(iPIT, "punif")$p.value
    indep <- dHSIC::dhsic.test(iPIT, args$data$z, method = "gamma")$p.value
    mod$xi <- args$xi
    mod$p.unif <- unif
    mod$p.indep <- indep
    if (min(unif, indep) > alpha)
      return(mod)
    else
      args$xi <- ifelse(indep < unif, args$xi * (1 + stepsize),
                        args$xi / (1 + stepsize))
  }
  message("No solution for which uniformity and independence is not
          rejected at level alpha.")
  return(do.call("ColrDA", args))
}

# Run ---------------------------------------------------------------------

tord <- 10

### Fit
m0 <- BoxCox(y | 0 + d ~ 1, data = dat, support = range(dat$y), order = tord)

### Warmstart with conditional distribution
mtmp <- ColrNN(y | d ~ 1, data = dat, order = tord,
               optimizer = optimizer_adam(0.1))
fit(mtmp, epochs = 3e3, validation_split = 0, callbacks = list(
  callback_reduce_lr_on_plateau("loss", factor = 0.9, patience = 20),
  callback_early_stopping("loss", patience = 40)))
tmp <- get_weights(mtmp$model)
args <- list(formula = y | d ~ 1, data = dat, anchor = ~ z, loss = "indep",
             optimizer = optimizer_adam(0.01), order = tord, xi = 0.001)
cb <- list(callback_reduce_lr_on_plateau("loss", patience = 2e2, factor = 0.9),
           callback_early_stopping("loss", patience = 4e2))
m <- fit_adaptive(args, 1e4, callbacks = cb, ws = tmp)

### Predict
dat$TRAM <- predict(m0, which = "distribution", type = "distribution")
dat$DIVE <- predict(m, type = "cdf")

pdat <- dat |> pivot_longer(TRAM:DIVE, names_to = "model", values_to = "rank")

p1 <- ggplot(pdat, aes(x = rank, color = factor(z))) +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "gray40") +
  facet_wrap(~ model) +
  stat_ecdf() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "PIT rank", y = "ECDF") +
  theme_bw() +
  theme(text = element_text(size = 13.5))

nd <- expand_grid(d = sort(unique(dat$d)), y = seq(
  min(dat$y), max(dat$y), length.out = 1e3))
nd$TRAM <- predict(m0, type = "distribution", newdata = nd)
nd$DIVE <- predict(m, type = "cdf", newdata = nd)

p2 <- ggplot(nd |> pivot_longer(TRAM:DIVE, names_to = "model", values_to = "cdf"),
       aes(x = y, y = cdf, color = factor(d))) +
  geom_rug(aes(x = y), data = dat, inherit.aes = FALSE, color = "gray80",
           alpha = 0.1) +
  facet_wrap(~ model) +
  geom_line() +
  labs(x = "log(wage)", y = "CDF") +
  theme_bw() +
  theme(text = element_text(size = 13.5))

ggpubr::ggarrange(p2, p1, ncol = 1, align = "hv")

ggsave("inst/figures/401k.pdf", height = 6, width = 7)

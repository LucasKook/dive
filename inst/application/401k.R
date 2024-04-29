### 401(k) DIVE
### LK 2024

set.seed(12)
save <- TRUE

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
  z = factor(d401k$e401)
)[sample.int(nrow(d401k), 1e3), ]

idx <- which(dat$y == 0)
dat$y[idx] <- runif(length(idx), max(dat$y[dat$y < 0]), min(dat$y[dat$y > 0]))

dat$oy <- ordered(dat$y)

# FUNs --------------------------------------------------------------------

fit_adaptive <- function(
    args, epochs, max_iter = 5, stepsize = 2, alpha = 0.1, ws = NULL, ...
) {
  for (iter in seq_len(max_iter)) {
    mod <- do.call("PolrDA", c(args, list(tf_seed = iter)))
    if (!is.null(ws))
      set_weights(mod$model, ws)
    plot(args$data$y, predict(mod, type = "cdf"), col = args$data$d)
    fit(mod, epochs = epochs, ...)
    plot(args$data$y, predict(mod, type = "cdf"), col = args$data$d)
    iPIT <- predict(mod, type = "cdf")
    unif <- ks.test(iPIT, "punif")$p.value
    indep <- dHSIC::dhsic.test(iPIT, as.numeric(args$data$z),
                               method = "gamma")$p.value
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

tord <- 6

### Fit
m0 <- BoxCox(y | 0 + d ~ 1, data = dat, support = range(dat$y), order = tord)

### Initialization
mtmp <- PolrNN(oy | d ~ 1, data = dat)
tmp <- get_weights(mtmp$model)
tmp[[1]][] <- -4.5
tmp[[2]][] <- -4.5
tmp[[3]][] <- -1
args <- list(formula = oy | d ~ 1, data = dat, anchor = ~ z, loss = "indep",
             optimizer = optimizer_adam(0.1), xi = 1/3, tf_seed = 1)
cb <- list(callback_reduce_lr_on_plateau("loss", patience = 2e1, factor = 0.9),
           callback_early_stopping("loss", patience = 6e1))
m <- fit_adaptive(args, 1e4, callbacks = cb, ws = tmp)

### Predict
dat$Nonparametric <- predict(m0, which = "distribution", type = "distribution")
dat$DIVE <- predict(m, type = "cdf")

pdat <- dat |> pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "rank")

# Vis ---------------------------------------------------------------------

p1 <- ggplot(pdat, aes(x = rank, color = factor(z))) +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "gray40") +
  facet_wrap(~ model) +
  stat_ecdf() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Estimated iPIT", y = "ECDF", color = "401(k) eligibility") +
  theme_bw() +
  theme(text = element_text(size = 13.5))

nd <- expand_grid(d = sort(unique(dat$d)), y = unique(dat$y))
nd$oy <- ordered(nd$y)
nd$Nonparametric <- predict(m0, type = "distribution", newdata = nd)
nd$DIVE <- predict(m, type = "cdf", newdata = nd)

p2 <- ggplot(nd |> pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "cdf"),
       aes(x = y, y = cdf, color = factor(d))) +
  geom_rug(aes(x = y), data = dat, inherit.aes = FALSE, color = "gray80",
           alpha = 0.1) +
  facet_wrap(~ model) +
  geom_step() +
  labs(x = "Net total financial assets", y = "Estimated CDF", color = "401(k) participation") +
  theme_bw() +
  theme(text = element_text(size = 13.5)) +
  scale_color_manual(values = colorspace::diverge_hcl(2))

ggpubr::ggarrange(p2, p1, ncol = 1, align = "hv")

# Save --------------------------------------------------------------------

if (save)
  ggsave("inst/figures/401k.pdf", height = 6, width = 8)

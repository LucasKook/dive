### Return to schooling application
### LK 2023

set.seed(12)
save <- TRUE

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Data --------------------------------------------------------------------

data("SchoolingReturns", package = "ivreg")
SchoolingReturns$wage <- log(SchoolingReturns$wage)

run <- \(iter) {
  set.seed(iter)
  dat <- SchoolingReturns[sample.int(nrow(SchoolingReturns), 1e3), ]

  # Run ---------------------------------------------------------------------

  ### Nonparametric
  m0 <- BoxCox(wage | smsa ~ 1, data = dat, support = range(dat$wage))

  ### DIVE
  args <- list(formula = wage | smsa ~ 1, data = dat, anchor = ~ nearcollege,
                loss = "indep", optimizer = optimizer_adam(0.05),
                order = 10, xi = 1, tf_seed = iter)
  cb <- list(callback_reduce_lr_on_plateau("loss", factor = 0.9, patience = 20),
             callback_early_stopping("loss", patience = 40))
  m <- fit_adaptive(args, epochs = 1e4, max_iter = 5, stepsize = 2, alpha = 0.1,
                    ws = NULL, modFUN = "BoxCoxDA", callbacks = cb)

  dat$Nonparametric <- predict(m0, which = "distribution", type = "distribution")
  dat$DIVE <- predict(m, type = "cdf")

  pd <- dat |>
    pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "rank") |>
    mutate(iter = iter)

  nd <- expand_grid(smsa = sort(unique(dat$smsa)), wage = seq(
    min(dat$wage), max(dat$wage), length.out = 1e3))
  nd$Nonparametric <- predict(m0, type = "distribution", newdata = nd)
  nd$DIVE <- predict(m, type = "cdf", newdata = nd)
  nd$iter <- iter

  list(pd = pd, nd = nd)
}

nsim <- 50
ret <- lapply(seq_len(nsim), run)
pdat <- do.call("rbind", lapply(ret, \(x) x[["pd"]]))
nd <- do.call("rbind", lapply(ret, \(x) x[["nd"]]))

# Vis ---------------------------------------------------------------------

p1 <- ggplot(pdat, aes(x = rank, color = nearcollege, linetype = factor(iter))) +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "gray40") +
  facet_wrap(~ model) +
  stat_ecdf() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Estimated iPIT", y = "ECDF", color = "Near college") +
  theme_bw() +
  theme(text = element_text(size = 13.5)) +
  guides(linetype = "none")

p2 <- ggplot(
  nd |> pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "cdf"),
  aes(x = wage, y = cdf, color = smsa, linetype = factor(iter))) +
  # geom_rug(aes(x = wage), data = dat, inherit.aes = FALSE, color = "gray80",
  #          alpha = 0.1) +
  facet_wrap(~ model) +
  geom_line() +
  labs(x = "log(wage)", y = "Estimated CDF", color = "Metropolitan area") +
  theme_bw() +
  theme(text = element_text(size = 13.5)) +
  scale_color_manual(values = colorspace::diverge_hcl(2)) +
  guides(linetype = "none")

ggpubr::ggarrange(p2, p1, ncol = 1, align = "hv")

# Save --------------------------------------------------------------------

if (save)
  ggsave("inst/figures/schooling.pdf", height = 6, width = 7)

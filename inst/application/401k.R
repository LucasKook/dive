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

raw <- read_csv("inst/data/401k.csv")

d401k <- data.frame(
  y = raw$net_tfa / 1e3,
  d = factor(raw$p401),
  z = factor(raw$e401)
) |> filter(y > 0) |> mutate(y = log(y))

run <- \(iter) {
  cat("\nIteration:", iter)
  set.seed(iter)
  dat <- d401k[sample.int(nrow(d401k), 1e3), ]
  dat$oy <- dat$y

  # Run ---------------------------------------------------------------------

  ### Fit
  m0 <- BoxCox(y | 0 + d ~ 1, data = dat, support = range(dat$y), order = 20)

  ### Warmstart model
  mtmp <- BoxCoxNN(y | d ~ 1, data = dat, optimizer = optimizer_adam(1e-2),
                   order = 20)
  fit(mtmp, epochs = 3e3, validation_split = 0, verbose = FALSE, callbacks = list(
    callback_reduce_lr_on_plateau("loss", factor = 0.9, patience = 20, min_delta = 1e-3),
    callback_early_stopping("loss", patience = 60, min_delta = 1e-3)))
  tmp <- get_weights(mtmp$model)

  args <- list(formula = oy | d ~ 1, data = dat, anchor = ~ z, loss = "indep",
               xi = 1/3, tf_seed = iter, order = 20)
  cb <- \() list(callback_reduce_lr_on_plateau(
    "loss", patience = 2e1, factor = 0.9, min_delta = 1e-4),
    callback_early_stopping("loss", patience = 60, min_delta = 1e-4))
  m <- fit_adaptive(args, epochs = 1e4, max_iter = 10, stepsize = 5, alpha = 0.1,
                    cb = cb, ws = tmp, modFUN = "BoxCoxDA", lr = 0.1,
                    indep_over_unif = 1e2, start_xi = TRUE, verbose = FALSE)

  ### Predict
  dat$Nonparametric <- c(predict(m0, which = "distribution", type = "distribution"))
  dat$DIVE <- c(predict(m, type = "cdf"))

  pd <- dat |>
    pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "rank")
  pd$iter <- iter

  nd <- expand_grid(d = sort(unique(dat$d)), y = unique(dat$y))
  # nd$oy <- ordered(nd$y)
  nd$oy <- nd$y
  nd$Nonparametric <- c(predict(m0, type = "distribution", newdata = nd))
  nd$DIVE <- c(predict(m, type = "cdf", newdata = nd))
  nd$iter <- iter

  list(pd = pd, nd = nd)

}

nsim <- 50
ret <- lapply(seq_len(nsim), run)
pdat <- do.call("rbind", lapply(ret, \(x) x[["pd"]]))
nd <- do.call("rbind", lapply(ret, \(x) x[["nd"]]))

# Vis ---------------------------------------------------------------------

p1 <- ggplot(pdat, aes(x = rank, color = factor(z), linetype = factor(iter))) +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "gray40") +
  facet_wrap(~ model) +
  stat_ecdf(alpha = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Estimated iPIT", y = "ECDF", color = "401(k) eligibility") +
  theme_bw() +
  theme(text = element_text(size = 13.5)) +
  guides(linetype = "none")

p2 <- nd |>
  pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "cdf") |>
  ggplot(aes(x = y, y = cdf, color = factor(d), linetype = factor(iter))) +
  facet_wrap(~ model) +
  geom_line(alpha = 0.5) +
  labs(x = "Net total financial assets", y = "Estimated CDF", color = "401(k) participation") +
  theme_bw() +
  theme(text = element_text(size = 13.5)) +
  scale_color_manual(values = colorspace::diverge_hcl(2)) +
  guides(linetype = "none")

ggpubr::ggarrange(p2, p1, nrow = 1, align = "hv", legend = "top")

# Save --------------------------------------------------------------------

if (save) {
  write_csv(pdat, "inst/figures/401k-pdat.csv")
  write_csv(nd, "inst/figures/401k-nd.csv")
  ggsave("inst/figures/401k.pdf", height = 3.5, width = 12)
}

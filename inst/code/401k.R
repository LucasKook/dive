### 401(k) DIVE
### LK 2024

waggr <- as.numeric(commandArgs(TRUE)[1])
if (is.na(waggr)) waggr <- 1

set.seed(12)
save <- TRUE

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Data --------------------------------------------------------------------

raw <- read_csv("inst/data/401k.csv")
taggr <- c("k_sum", "k_max")[waggr]

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
               xi = 1/3, tf_seed = iter, order = 20, aggr = taggr)
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

# Save --------------------------------------------------------------------

if (save) {
  if (!dir.exists("inst/figures"))
    dir.create("inst/figures", recursive = TRUE)
  bp <- file.path("inst", "figures")
  pa1 <- file.path(bp, paste0("401k", ifelse(waggr == 2, "-max", ""), "-pdat.csv"))
  pa2 <- file.path(bp, paste0("401k", ifelse(waggr == 2, "-max", ""), "-nd.csv"))
  write_csv(pdat, pa1)
  write_csv(nd, pa2)
}

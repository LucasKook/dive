### Return to schooling application
### LK 2023

waggr <- as.numeric(commandArgs(TRUE)[1])
if (is.na(waggr)) waggr <- 1

set.seed(12)
save <- TRUE
warmstart <- TRUE

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Data --------------------------------------------------------------------

data("SchoolingReturns", package = "ivreg")
SchoolingReturns$wage <- log(SchoolingReturns$wage)
taggr <- c("k_sum", "k_max")[waggr]

run <- \(iter) {
  cat("\nIteration:", iter)
  set.seed(iter)
  dat <- SchoolingReturns[sample.int(nrow(SchoolingReturns), 1e3), ]

  # Run ---------------------------------------------------------------------

  ### Nonparametric
  m0 <- BoxCox(wage | smsa ~ 1, data = dat, support = range(dat$wage))

  ### Warmstart model
  tmp <- NULL
  if (warmstart) {
    mtmp <- BoxCoxNN(wage | smsa ~ 1, data = dat,
                     optimizer = optimizer_adam(1e-2))
    fit(mtmp, epochs = 3e3, validation_split = 0, verbose = FALSE, callbacks = list(
      callback_reduce_lr_on_plateau("loss", factor = 0.9, patience = 20, min_delta = 1e-3),
      callback_early_stopping("loss", patience = 60, min_delta = 1e-3)))
    tmp <- get_weights(mtmp$model)
  }

  ### DIVE
  args <- list(formula = wage | smsa ~ 1, data = dat, anchor = ~ nearcollege,
               loss = "indep",
               order = 10, xi = 1, tf_seed = iter, aggr = taggr)
  cb <- \() list(callback_reduce_lr_on_plateau(
    "loss", factor = 0.9, patience = 20, min_delta = 1e-4),
    callback_early_stopping("loss", patience = 60, min_delta = 1e-4))
  m <- fit_adaptive(args, epochs = 1e4, max_iter = 10, stepsize = 5,
                    alpha = 0.1, ws = tmp, modFUN = "BoxCoxDA", cb = cb,
                    start_xi = warmstart, lr = 0.05, indep_over_unif = 1e3,
                    verbose = FALSE)

  dat$Nonparametric <- c(predict(m0, which = "distribution", type = "distribution"))
  dat$DIVE <- c(predict(m, type = "cdf"))

  pd <- dat |>
    pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "rank") |>
    mutate(iter = iter)

  nd <- expand_grid(smsa = sort(unique(dat$smsa)), wage = seq(
    min(dat$wage), max(dat$wage), length.out = 1e3))
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
  bp <- file.path("inst", "results", "figures")
  if (!dir.exists(bp))
    dir.create(bp, recursive = TRUE)
  pa1 <- file.path(bp, paste0("schooling", ifelse(waggr == 2, "max", ""),
                              "-pdat.csv"))
  pa2 <- file.path(bp, paste0("schooling", ifelse(waggr == 2, "max", ""),
                              "-nd.csv"))
  write_csv(pdat, pa1)
  write_csv(nd, pa2)
}

### Coverage simulation DIVE
### LK 2025

set.seed(12)

### CLI args
args <- commandArgs(trailingOnly = TRUE)
scenario <- if (length(args) != 0) args[1] else 1
run <- if (length(args) != 0) as.numeric(args[2]) else 1
vb <- FALSE

### File names
save <- TRUE
odir <- file.path("inst/results/coverage")
fname <- paste0("sim-res-", scenario, "-run-", run)

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Settings ----------------------------------------------------------------

dgp <- switch(scenario,
  "1" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(4 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- -10 + 8 * D + 6 * H
    data.frame(Y = Y, H = H, D = D, Z = Z)
  },
  "2" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(4 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- 16 * D + 6 * H + H * rlogis(n)
    data.frame(Y = Y, H = H, D = D, Z = Z)
  },
  "3" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(4 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- log(1 + exp(18 + 8 * D + 6 * H))
    data.frame(Y = Y, H = H, D = D, Z = Z)
  },
  "4" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(4 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- -6 * D + (6 + 3 * D) * H
    data.frame(Y = Y, H = H, D = D, Z = Z)
  }
)

# ORACLE ------------------------------------------------------------------

dint <- dgp(1e6, do = TRUE)
F0 <- ecdf(dint$Y[dint$D == 0])
F1 <- ecdf(dint$Y[dint$D == 1])
oracle <- Vectorize(\(y, d) d * F1(y) + (1 - d) * F0(y))

# Params ------------------------------------------------------------------

nep <- 3e3
wep <- 3e3
rep <- 50
ns <- 100 * 2^(0:4)
tlr <- 0.1
tord <- 50
iou <- 1
alp <- 0.1
tss <- 10

# FUNs --------------------------------------------------------------------

check_unif_indep <- function(iPIT, Z) {
  plot(iPIT ~ Z)
  c(
    unif = ks.test(iPIT, "punif")$p.value,
    indep = dHSIC::dhsic.test(iPIT, Z, method = "gamma")$p.value
  )
}

# Run ---------------------------------------------------------------------

full_dat <- dgp(1e4)
grd <- seq(min(full_dat$Y), max(full_dat$Y), length.out = 3e2)
nd <- data.frame(expand.grid(Y = grd, D = 0:1))

res <- lapply(ns, \(tn) {
  cat("\nRunning with n =", tn, "\n")
  pb <- txtProgressBar(min = 0, max = rep, style = 3, width = 60)
  lapply(seq_len(rep), \(iter) {
    setTxtProgressBar(pb, iter)
    set.seed(1e4 + iter + 1e3 * run)

    ### Generate data
    idx <- sample.int(NROW(full_dat), tn)
    dat <- full_dat[idx, ]

    ### Compute oracle interventional CDF
    nd$ORACLE <- oracle(nd$Y, nd$D)

    ### Fit vanilla nonparametric CDF
    m0 <- BoxCox(Y | D ~ 1,
      data = dat, support = supp <- range(full_dat$Y),
      order = tord
    )

    ### Warmstart with conditional distribution
    mtmp <- BoxCoxNN(Y | D ~ 1,
      data = dat, order = tord,
      optimizer = optimizer_adam(1e-2),
      tf_seed = iter
    )
    fit(mtmp, epochs = wep, validation_split = 0, callbacks = list(
      callback_reduce_lr_on_plateau("loss", factor = 0.9, patience = 20, min_delta = 1e-3),
      callback_early_stopping("loss", patience = 60, min_delta = 1e-3)
    ), verbose = vb)
    tmp <- get_weights(mtmp$model)

    ### Fit DIVE
    args <- list(
      formula = Y | D ~ 1, data = dat, anchor = ~Z,
      loss = "indep", xi = 1, trafo_options = trafo_control(
        order_bsp = tord, support = supp
      )
    )
    cb <- \() list(
      callback_reduce_lr_on_plateau(
        "loss",
        patience = 20, factor = 0.9, min_delta = 1e-4
      ),
      callback_early_stopping("loss", patience = 60, min_delta = 1e-4)
    )
    m <- fit_adaptive(args, nep,
      max_iter = 10, ws = tmp,
      modFUN = "BoxCoxDA", verbose = vb, lr = tlr,
      cb = cb, start_xi = TRUE, stepsize = tss,
      indep_over_unif = iou, alpha = alp
    )

    ### Evaluate
    nd$TRAM <- c(predict(m0,
      which = "distribution",
      type = "distribution",
      newdata = nd
    ))
    nd$DIVE <- c(predict(m, type = "cdf", newdata = nd))

    nd |>
      pivot_longer(TRAM:DIVE,
        names_to = "method",
        values_to = "cdf"
      ) |>
      mutate(n = tn, order = tord, lr = tlr, iter = iter, scenario = scenario, run = run)
  }) |> bind_rows()
}) |> bind_rows()

out <- res |>
  group_by(Y, D, n, order, lr, method, scenario, run) |>
  summarize(
    lwr = quantile(cdf, 0.1),
    med = quantile(cdf, 0.5),
    upr = quantile(cdf, 0.9),
    ORACLE = median(ORACLE)
  )

# Save --------------------------------------------------------------------

if (save) {
  if (!dir.exists(odir)) {
    dir.create(odir, recursive = TRUE)
  }
  write_csv(out, file.path(odir, paste0(fname, ".csv")))
}

### Simulation DIVE
### LK 2024

set.seed(12)

### CLI args
args <- commandArgs(trailingOnly = TRUE)
scenario <- if (length(args) != 0) args[1] else 1
vb <- FALSE

### File names
save <- TRUE
odir <- file.path("inst/results/simulations")
fname <- paste0("sim-res-nonp", scenario)

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Settings ----------------------------------------------------------------

dgp <- switch(
  scenario,
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

# d <- dgp(3e3)
# plot(d$Y, oracle(d$Y, d$D))
# plot(ecdf(d$Y[d$D == 0]), add = TRUE)
# plot(ecdf(d$Y[d$D == 1]), add = TRUE)

# Params ------------------------------------------------------------------

nep <- 3e3
wep <- 3e3
rep <- 50
ns <- 100 * 2^(0:4)
lrs <- 0.01
iou <- 1
alp <- 0.1
tss <- 10

# FUNs --------------------------------------------------------------------

check_unif_indep <- function(iPIT, Z) {
  plot(iPIT ~ Z)
  c(unif = ks.test(iPIT, "punif")$p.value,
    indep = dHSIC::dhsic.test(iPIT, Z, method = "gamma")$p.value)
}

# Run ---------------------------------------------------------------------

res <- lapply(ns, \(tn) {
  lapply(lrs, \(tlr) {
    cat("\nRunning with n =", tn, ", lr =", tlr, "\n")
    pb <- txtProgressBar(min = 0, max = rep, style = 3)
    lapply(seq_len(rep), \(iter) {
      setTxtProgressBar(pb, iter)
      set.seed(1e4 + iter)
      ### Generate data
      dat <- dgp(tn)

      ### Compute oracle interventional CDF
      dat$ORACLE <- oracle(dat$Y, dat$D)
      dat$Y <- ordered(dat$Y)

      ### Checks
      # print(check_unif_indep(dat$ORACLE, dat$Z))
      # print(check_unif_indep(ecdf(dat$Y)(dat$Y), dat$Z))

      ### Fit vanilla nonparametric CDF
      m0 <- Polr(Y | D ~ 1, data = dat, method = "logistic")

      ### Warmstart with conditional distribution
      mtmp <- PolrNN(Y | D ~ 1, data = dat, optimizer = optimizer_adam(1e-2),
                     tf_seed = iter, latent_distr = "logistic")
      tmp <- get_weights(mtmp$model)
      tmp[[1]][] <- -1
      tmp[[2]][] <- -1
      tmp[[3]][] <- -1
      set_weights(mtmp$model, tmp)
      fit(mtmp, epochs = wep, validation_split = 0, callbacks = list(
        callback_reduce_lr_on_plateau("loss", factor = 0.9, patience = 20, min_delta = 1e-3),
        callback_early_stopping("loss", patience = 60, min_delta = 1e-3)), verbose = vb)
      tmp <- get_weights(mtmp$model)
      ### Fit DIVE
      args <- list(formula = Y | D ~ 1, data = dat, anchor = ~ Z,
                   loss = "indep", xi = 1, latent_distr = "logistic")
      cb <- \() list(callback_reduce_lr_on_plateau(
        "loss", patience = 20, factor = 0.9, min_delta = 1e-4),
        callback_early_stopping("loss", patience = 60, min_delta = 1e-4))
      m <- fit_adaptive(args, nep, max_iter = 10, ws = tmp,
                        modFUN = "PolrDA", verbose = vb, lr = tlr,
                        cb = cb, start_xi = TRUE, stepsize = tss,
                        indep_over_unif = iou, alpha = alp)

      ### Evaluate
      dat$TRAM <- c(predict(m0, which = "distribution",
                            type = "distribution"))
      dat$DIVE <- c(predict(m, type = "cdf"))

      dat |> pivot_longer(TRAM:DIVE, names_to = "method",
                          values_to = "cdf") |>
        group_by(method) |>
        summarize(MSE = mean((ORACLE - cdf)^2),
                  MAE = max(abs(ORACLE - cdf))) |>
        mutate(n = tn, lr = tlr, iter = iter)
    }) |> bind_rows()
  }) |> bind_rows()
}) |> bind_rows()

# Vis ---------------------------------------------------------------------

res |> pivot_longer(MSE:MAE, names_to = "metric", values_to = "value") |>
  ggplot(aes(x = ordered(n), y = value, fill = method)) +
  facet_wrap(~ metric) +
  geom_boxplot() +
  theme_bw()

# Save --------------------------------------------------------------------

if (save) {
  if (!dir.exists(odir))
    dir.create(odir, recursive = TRUE)
  write_csv(res, file.path(odir, paste0(fname, ".csv")))
  ggsave(file.path(odir, paste0(fname, ".pdf")))
}

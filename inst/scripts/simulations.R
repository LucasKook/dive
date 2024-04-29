### Simulation DIVE
### LK 2024

set.seed(12)

### CLI args
args <- commandArgs(trailingOnly = TRUE)
scenario <- if (length(args) != 0) args[1] else 1

### File names
save <- TRUE
odir <- file.path("inst/results/simulations", Sys.Date())
fname <- paste0("sim-res-", scenario)

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
    Y <- 16 * D + 6 * H + 2 * H * rlogis(n)
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
    Y <- - 4 * (1 + 3 * D) * H
    data.frame(Y = Y, H = H, D = D, Z = Z)
  }
)

# ORACLE ------------------------------------------------------------------

dint <- dgp(1e6, do = TRUE)
F0 <- ecdf(dint$Y[dint$D == 0])
F1 <- ecdf(dint$Y[dint$D == 1])
oracle <- Vectorize(\(y, d) d * F1(y) + (1 - d) * F0(y))

# Params ------------------------------------------------------------------

nep <- 1e4
rep <- 20
ords <- c(10, 30, 50)
lrs <- c(0.01, 0.05, 0.1)
ns <- c(1e2, 3e2, 7e2, 1e3)

# FUNs --------------------------------------------------------------------

check_unif_indep <- function(iPIT, Z) {
  plot(iPIT ~ Z)
  c(unif = ks.test(iPIT, "punif")$p.value,
    indep = dHSIC::dhsic.test(iPIT, Z, method = "gamma")$p.value)
}

# Run ---------------------------------------------------------------------

res <- lapply(ns, \(tn) {
  lapply(ords, \(tord) {
    lapply(lrs, \(tlr) {
      lapply(seq_len(rep), \(iter) {
        ### Generate data
        dat <- dgp(tn)

        ### Compute oracle interventional CDF
        dat$ORACLE <- oracle(dat$Y, dat$D)

        ### Checks
        # print(check_unif_indep(dat$ORACLE, dat$Z))
        # print(check_unif_indep(ecdf(dat$Y)(dat$Y), dat$Z))

        ### Fit vanilla nonparametric CDF
        m0 <- Colr(Y | D ~ 1, data = dat, support = supp <- range(dat$Y),
                   order = tord)

        ### Warmstart with conditional distribution
        mtmp <- ColrNN(Y | D ~ 1, data = dat, order = tord,
                       optimizer = optimizer_adam(1e-2))
        fit(mtmp, epochs = 3e3, validation_split = 0, callbacks = list(
          callback_reduce_lr_on_plateau("loss", factor = 0.9, patience = 20),
          callback_early_stopping("loss", patience = 200)))
        tmp <- get_weights(mtmp$model)
        ### Fit DIVE
        args <- list(formula = Y | D ~ 1, data = dat, anchor = ~ Z,
                     loss = "indep", optimizer = optimizer_adam(tlr),
                     xi = 1, trafo_options = trafo_control(
                       order_bsp = tord, support = supp))
        cb <- list(callback_reduce_lr_on_plateau("loss", patience = 2e2,
                                                 factor = 0.9),
                   callback_early_stopping("loss", patience = 4e2))
        m <- fit_adaptive(args, nep, callbacks = cb, ws = tmp, modFUN = ColrDA)

        ### Evaluate
        dat$TRAM <- c(predict(m0, which = "distribution",
                              type = "distribution"))
        dat$DIVE <- c(predict(m, type = "cdf"))

        dat |> pivot_longer(TRAM:DIVE, names_to = "method",
                            values_to = "cdf") |>
          group_by(method) |>
          summarize(CmV = mean((ORACLE - cdf)^2),
                    KS = max(abs(ORACLE - cdf))) |>
          mutate(n = tn, order = tord, lr = tlr, iter = iter)
      }) |> bind_rows()
    }) |> bind_rows()
  }) |> bind_rows()
}) |> bind_rows()

# Vis ---------------------------------------------------------------------

res |> pivot_longer(CmV:KS, names_to = "metric", values_to = "value") |>
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

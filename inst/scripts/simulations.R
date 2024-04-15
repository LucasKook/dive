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
    D <- as.numeric(3 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- 4 * D + 6 * H
    data.frame(Y = Y, H = H, D = D, Z = Z)
  },
  "2" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(3 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- 4 * D + 6 * H + 2 * H * rlogis(n)
    data.frame(Y = Y, H = H, D = D, Z = Z)
  },
  "3" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(3 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- log(1 + exp(18 + 8 * D + 6 * H))
    data.frame(Y = Y, H = H, D = D, Z = Z)
  },
  "4" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(3 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- - 4 * (1 + D) * H
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
lams <- c(0, 0.5, 1, 3, 10)

# Run ---------------------------------------------------------------------

res <- lapply(ns, \(tn) {
  lapply(ords, \(tord) {
    lapply(lrs, \(tlr) {
      lapply(lams, \(tlam) {
        lapply(seq_len(rep), \(iter) {
          dat <- dgp(tn)
          dat$ORACLE <- oracle(dat$Y, dat$D)
          m0 <- BoxCox(Y | D ~ 1, data = dat, support = range(dat$Y),
                       order = tord)
          m <- BoxCoxDA(Y | D ~ 1, data = dat, anchor = ~ Z, loss = "indep",
                        optimizer = optimizer_adam(tlr), order = tord,
                        xi = tlam)
          fit(m, epochs = nep)
          dat$TRAM <- c(predict(m0, which = "distribution",
                                type = "distribution"))
          dat$DIVE <- c(predict(m, type = "cdf"))

          dat |> pivot_longer(TRAM:DIVE, names_to = "method",
                              values_to = "cdf") |>
            group_by(method) |>
            summarize(CmV = mean((ORACLE - cdf)^2),
                      KS = max(abs(ORACLE - cdf))) |>
            mutate(n = tn, lam = tlam, order = tord, lr = tlr)
        }) |> bind_rows()
      }) |> bind_rows()
    }) |> bind_rows()
  }) |> bind_rows()
}) |> bind_rows()

# Vis ---------------------------------------------------------------------

ggplot(res |> pivot_longer(CmV:KS, names_to = "metric", values_to = "value"),
       aes(x = ordered(n), y = value, fill = method)) +
  facet_grid(metric ~ lam + lr, scales = "free", labeller = label_both) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw()

# Save --------------------------------------------------------------------

if (save) {
  if (!dir.exists(odir))
    dir.create(odir, recursive = TRUE)
  write_csv(res, file.path(odir, paste0(fname, ".csv")))
  ggsave(file.path(odir, paste0(fname, ".pdf")))
}

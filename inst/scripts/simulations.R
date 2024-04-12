### Simulation DIVE
### LK 2024

set.seed(12)
save <- TRUE
odir <- "inst/results/simulations"
fname <- "sim-res"

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Settings ----------------------------------------------------------------

dgp <- function(n = 1e2) {
  Z <- rlogis(n)
  H <- rnorm(n)
  D <- as.numeric(3 * Z + 0.5 * H > rlogis(n))
  Y <- 4 * D + 4 * H
  oracle <- Vectorize(\(y, d) pnorm(y, mean = 4 * d, sd = 2))
  data.frame(Y = Y, H = H, D = D, Z = Z, ORACLE = oracle(Y, D))
}

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
  facet_grid(lam ~ metric, scales = "free", labeller = label_both) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw()

# Save --------------------------------------------------------------------

if (save) {
  if (!dir.exists(odir))
    dir.create(odir, recursive = TRUE)
  write_csv(res, file.path(odir, paste0(fname, ".csv")))
  ggsave(file.path(odir, paste0(fname, ".pdf")))
}

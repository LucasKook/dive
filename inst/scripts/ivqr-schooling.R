
set.seed(12)

library("tidyverse")
library("ivreg")
library("IVQR")

data("SchoolingReturns", package = "ivreg")
SchoolingReturns$wage <- log(SchoolingReturns$wage)
SchoolingReturns$smsa <- as.numeric(SchoolingReturns$smsa) - 1
SchoolingReturns$nearcollege <- as.numeric(SchoolingReturns$nearcollege) - 1

pb <- txtProgressBar(0, 50, style = 3)
res <- do.call("rbind", lapply(1:50, \(iter) {
  setTxtProgressBar(pb, iter)
  dat <- SchoolingReturns[sample.int(nrow(SchoolingReturns), 1e3), ]
  lm <- lm(wage ~ smsa, data = dat)
  iv <- ivreg(wage ~ smsa | nearcollege, data = dat)
  grd <- seq(-0.5, 2, length.out = 1e2)
  fit <- ivqr(wage ~ smsa | nearcollege | 1,
              tt <- seq(0.05, 0.95, length.out = 10),
              grid = grd, data = dat, qrMethod = "br")
  tibble(ntau = tt, "OLS" = unname(coef(lm)[2]), "2SLS" = unname(coef(iv)[2]),
         "IVQR" = c(fit$coef$endg_var), iter = iter)
})) |> as.data.frame() |>
  pivot_longer(OLS:IVQR, names_to = "method", values_to = "qte")

nd <- read_csv("inst/results/figures/schooling-nd.csv")
pdat <- read_csv("inst/results/figures/schooling-pdat.csv") |>
  mutate(model = factor(model, levels = c("DIVE", "Nonparametric"),
                        labels = c("DIVE", "CCDF")))

f <- function(y, cdf, newx) {
  s <- spline(x = y, y = cdf, method = "hyman")
  approx(x = s$y, y = s$x, xout = newx, yleft = -Inf, yright = Inf)$y
}

pd <- nd |>
  mutate(smsa = ifelse(smsa == "no", 0, 1)) |>
  pivot_longer(Nonparametric:DIVE, names_to = "method", values_to = "tau") |>
  group_by(smsa, iter) |>
  mutate(ntau = seq(0.05, 0.95, length.out = length(smsa))) |>
  group_by(smsa, iter) |>
  mutate(ny = f(wage, tau, ntau)) |>
  select(smsa, ntau, ny) |>
  pivot_wider(names_from = "smsa", values_from = "ny") |>
  mutate(qte = `1` - `0`, method = "DIVE") |>
  select(ntau, iter, method, qte)

bind_rows(res, pd) |>
  ggplot(aes(x = ntau, y = qte, color = method, group = interaction(method, iter))) +
  geom_line(alpha = 0.4) +
  labs(x = bquote(tau), y = "Estimated QCE", color = "Estimator") +
  theme_bw() + theme(legend.position = "top")

ggsave("inst/figures/qce-comparison-schooling.pdf", height = 3, width = 5)


set.seed(12)

library("tidyverse")
library("ivreg")
library("IVQR")

raw <- read_csv("inst/data/401k.csv")

d401k <- data.frame(
  y = raw$net_tfa / 1e3,
  d = raw$p401,
  z = factor(raw$e401)
) |> filter(y > 0) |> mutate(y = log(y))

pb <- txtProgressBar(0, 50, style = 3)
res <- do.call("rbind", lapply(1:50, \(iter) {
  setTxtProgressBar(pb, iter)
  dat <- d401k[sample.int(nrow(d401k), 1e3), ]
  lm <- lm(y ~ d, data = dat)
  iv <- ivreg(y ~ d | z, data = dat)
  grd <- seq(-0.5, 4, length.out = 1e2)
  fit <- ivqr(y ~ d | z | 1, tt <- seq(0.05, 0.95, length.out = 10),
              grid = grd, data = dat)
  tibble(ntau = tt, "OLS" = unname(coef(lm)[2]), "2SLS" = unname(coef(iv)[2]),
         "IVQR" = c(fit$coef$endg_var), iter = iter)
})) |> as.data.frame() |>
  pivot_longer(OLS:IVQR, names_to = "method", values_to = "qte")

nd <- read_csv("inst/results/figures/401k-nd.csv")
pdat <- read_csv("inst/results/figures/401k-pdat.csv") |>
  mutate(model = factor(model, levels = c("DIVE", "Nonparametric"),
                        labels = c("DIVE", "CCDF")))

f <- function(y, cdf, newx) {
  s <- spline(x = y, y = cdf, method = "hyman")
  approx(x = s$y, y = s$x, xout = newx, yleft = -Inf, yright = Inf)$y
}

pd <- nd |>
  pivot_longer(Nonparametric:DIVE, names_to = "method", values_to = "tau") |>
  group_by(d, iter) |>
  mutate(ntau = seq(0.05, 0.95, length.out = length(d))) |>
  group_by(d, iter) |>
  mutate(ny = f(y, tau, ntau)) |>
  select(d, ntau, ny) |>
  pivot_wider(names_from = "d", values_from = "ny") |>
  mutate(qte = `1` - `0`, method = "DIVE") |>
  select(ntau, iter, method, qte)

bind_rows(res, pd) |>
  ggplot(aes(x = ntau, y = qte, color = method, group = interaction(method, iter))) +
  geom_line(alpha = 0.4) +
  labs(x = bquote(tau), y = "Estimated QCE", color = "Estimator") +
  theme_bw() + theme(legend.position = "top")

ggsave("inst/figures/qce-comparison-401k.pdf", height = 3, width = 5)

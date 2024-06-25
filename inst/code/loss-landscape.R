### Loss landscape for normal CDFs
### LK 2024

set.seed(12)

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")
library("dHSIC")

# FUNs --------------------------------------------------------------------

tg0 <- \(h, ny) h
tg1 <- \(h, ny) 1 + h
tqH <- \(h, z) plogis(h + z)

dgp_np <- function(n, g1 = tg1, g0 = tg0, fqH = tqH, intervene = FALSE) {
  Z <- rnorm(n)
  H <- rnorm(n)
  NY <- rnorm(n)
  qH <- if (intervene) fqH(0, Z) else fqH(H, Z)
  D <- rbinom(n, 1, prob = qH)
  Y0 <- g0(H, NY)
  Y1 <- g1(H, NY)
  Y <- D * Y1 + (1 - D) * Y0
  data.frame(Y0 = Y0, Y1 = Y1, H = H, D = D, Y = Y, qH = qH, Z = Z)
}

indep_unif_loss <- function(m0, m1, Y, D, Z, lam = 1) {
  dive <- Vectorize(\(y, d) pnorm(y, mean = d * m1 + (1 - d) * m0))
  pit <- dive(Y, D)
  unif <- goftest::cvm.test(pit)$statistic / NROW(Y) # sum((pit - ecdf(pit)(pit))^2)
  indep <- NROW(Y) * dhsic(list(qnorm(pit), Z), kernel = c("gaussian", "gaussian"))$dHSIC
  tibble(unif = unif, indep = indep, loss = sum(unif, lam * indep))
}

# Plot --------------------------------------------------------------------

d <- dgp_np(1e2)
ms <- seq(-3, 3, length.out = 1e2)
grd <- data.frame(expand.grid(m0 = ms, m1 = ms, lam = c(0.01, 0.1, 1)))
loss <- pmap(grd, indep_unif_loss, Y = d$Y, D = d$D, Z = d$Z) |>
  bind_rows()
pd <- bind_cols(grd, loss)

est <- pd |> group_by(lam) |> slice(which.min(loss)) |> ungroup()

ggplot(pd, aes(x = m0, y = m1, z = loss)) +
  facet_wrap(~ lam, nrow = 1, labeller = label_bquote(lambda==.(lam))) +
  labs(x = bquote(mu[0]), y = bquote(mu[1])) +
  geom_contour_filled(bins = 30, show.legend = FALSE) +
  annotate("point", x = 0, y = 1, color = "white", size = 3, pch = 4) +
  annotate("text", x = 0.8, y = 0.8, color = "white", size = 3, label = "ground truth") +
  geom_point(data = est, color = "white", size = 3, pch = 3) +
  geom_text(data = est, color = "white", size = 3, label = "estimate", nudge_y = 0.3) +
  theme_bw() + theme(text = element_text(size = 13.5))

ggsave("inst/figures/loss-landscape-gaussian.pdf", height = 3, width = 10)

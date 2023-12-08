### Visualize rank assumptions in examples
### LK 2023

set.seed(2410)
library("tidyverse")
library("patchwork")
theme_set(theme_bw() + theme(text = element_text(size = 13.5)))

# DGP ---------------------------------------------------------------------

setting <- c("np", "sqh", "cmrs-holds")[2]
if (setting == "np") {
  tg0 <- \(h, ny) h
  tg1 <- \(h, ny) -h
  tqH <- \(h) plogis(h)
} else if (setting == "sqh") {
  tg0 <- \(h, ny) h + ny / 3
  tg1 <- \(h, ny) h + 1.5 * ny / 3 - 1.5
  # tg0 <- \(h, ny) qlogis(pnorm(h + ny / 3))
  # tg1 <- \(h, ny) qlogis(pnorm(h + 1.5 * ny / 3 - 1.5))
  tqH <- \(h) 0.2 + 0.6 * as.numeric(h > 0)
} else if (setting == "cmrs-holds") {
  tg0 <- \(h, ny) h^2
  tg1 <- \(h, ny) -h^2
  tqH <- \(h) 0.2 + 0.6 * as.numeric(h > 0)
}

dgp_np <- function(n, g1 = tg1, g0 = tg0, fqH = tqH, intervene = FALSE) {
  ### Observational data
  H <- rnorm(n)
  NY <- rnorm(n)
  qH <- if (intervene) 0.5 else fqH(H)
  D <- rbinom(n, 1, prob = qH)
  Y0 <- g0(H, NY)
  Y1 <- g1(H, NY)
  Y <- D * Y1 + (1 - D) * Y0
  data.frame(Y0 = Y0, Y1 = Y1, H = H, D = D, Y = Y, qH = qH)
}

# Compute interventional distribution -------------------------------------

dint <- dgp_np(1e5, intervene = TRUE)
F0 <- ecdf(dint$Y0)
F1 <- ecdf(dint$Y1)

# Plot --------------------------------------------------------------------

### Generate data and compute ranks
pd <- dgp_np(1e2)
pd$R0 <- F0(pd$Y0)
pd$R1 <- F1(pd$Y1)
pd$W <- pd$D * pd$R1 + (1 - pd$D) * pd$R0

### Plot rank invariance
pri <- ggplot(pd, aes(x = Y0, y = Y1)) +
  geom_point() +
  labs(x = "Y(0)", y = "Y(1)")

### Plot rank similarity
prs <- ggplot(pd |> pivot_longer(R0:R1), aes(x = H, y = value, color = name)) +
  geom_point(show.legend = FALSE) +
  labs(x = "H", y = "Rank")

### Plot conditional mean rank similarity
pcmrs <- ggplot(pd |> pivot_longer(R0:R1), aes(x = qH, y = value, color = name)) +
  {if (length(unique(pd$qH)) > 2) geom_point() else ggbeeswarm::geom_quasirandom()} +
  labs(x = "q(H)", color = element_blank(), y = "Rank") +
  scale_color_discrete(labels = c("Y1" = "Y(1)", "Y0" = "Y(0)"))

### Plot uniformity condition
punif <- ggplot(pd, aes(x = W)) +
  stat_ecdf() +
  geom_abline(intercept = 0, slope = 1) +
  labs(y = "ECDF", x = "Interventional PIT")

(pri + labs(tag = "A", subtitle = "Rank invariance")) +
  (prs + labs(tag = "B", subtitle = "Rank similarity")) +
  (pcmrs + labs(tag = "C", subtitle = "Conditional mean\nrank similarity")) +
  (punif + labs(tag = "D", subtitle = "Uniformity condition")) +
  plot_layout(nrow = 1, guides = "collect")

ggsave(paste0("inst/figures/rank-assumptions", setting, ".pdf"),
       width = 12, height = 3.5)

### Visualize rank assumptions in examples
### LK 2023

set.seed(2410)
library("tidyverse")
library("patchwork")
theme_set(theme_bw() + theme(text = element_text(size = 13.5)))

# DGP ---------------------------------------------------------------------

setting <- c("rinv-violated", "rsim-violated", "cmrs-violated", "intro")[4]
if (setting == "cmrs-violated") {
  tg0 <- \(h, ny) h
  tg1 <- \(h, ny) -h
  tqH <- \(h) plogis(h)
} else if (setting == "rinv-violated") {
  tg0 <- \(h, ny) h + ny / 3
  tg1 <- \(h, ny) h + 1.5 * ny / 3 - 1.5
  # tg0 <- \(h, ny) qlogis(pnorm(h + ny / 3))
  # tg1 <- \(h, ny) qlogis(pnorm(h + 1.5 * ny / 3 - 1.5))
  tqH <- \(h) 0.2 + 0.6 * as.numeric(h > 0)
} else if (setting == "rsim-violated") {
  tg0 <- \(h, ny) h^2
  tg1 <- \(h, ny) -h^2
  tqH <- \(h) 0.2 + 0.6 * as.numeric(h > 0)
} else if (setting == "intro") {
  tg0 <- \(h, ny) 3 * cos(-1 + h + ny)
  tg1 <- \(h, ny) 3 * sin(-1 + h + ny)
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
  labs(x = parse(text = "g[0](H,N[Y])"), y = parse(text = "g[1](H,N[Y])"))

### Plot rank similarity
prs <- ggplot(pd |> pivot_longer(R0:R1), aes(x = H, y = value, color = name)) +
  geom_point(show.legend = FALSE) +
  labs(x = "H", y = parse(text = "Rank:~{{F[d]*'*'}}(g[d](H,N[Y]))")) +
  scale_color_discrete(labels = c("R1" = "1", "R0" = "0"))

### Plot conditional mean rank similarity
pcmrs <- ggplot(pd |> pivot_longer(R0:R1), aes(x = qH, y = value, color = name)) +
  {if (length(unique(pd$qH)) > 2) geom_point() else ggbeeswarm::geom_quasirandom()} +
  labs(x = parse(text = "pi(H)"), color = "", y = parse(text = "Rank:~{{F[d]*'*'}}(g[d](H,N[Y]))")) +
  scale_color_discrete(labels = c("R1" = parse(text = "d==1"), "R0" = parse(text = "d==0")))

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


# Intro plot --------------------------------------------------------------

do <- dgp_np(3e4, intervene = FALSE)
di <- dgp_np(3e4, intervene = TRUE)
pdat <- bind_rows(do |> mutate(setting = "Observational"),
                  di |> mutate(setting = "Interventional"))

p1 <- ggplot(pdat, aes(x = Y, color = factor(D))) +
  stat_ecdf() +
  facet_grid(~ setting) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Response", y = "CDF", color = "Treatment", tag = "A") +
  theme(legend.position = "top")

# qte <- \(p, d) quantile(d$Y1, probs = p) - quantile(d$Y0, probs = p)
# rte <- \(y, d) log(ecdf(d$Y1)(y)) - log(ecdf(d$Y0)(y))
# rpd <- \(p, d) ecdf(d$Y1)(quantile(d$Y0, probs = p)) - p
dte <- \(y, d) ecdf(d$Y[d$D == 1])(y) - ecdf(d$Y[d$D == 0])(y)
tte <- \(y, qZ = \(p) log(p) - log(1 - p), d) {
  ret <- qZ(ecdf(d$Y[d$D == 1])(y)) - qZ(ecdf(d$Y[d$D == 0])(y))
  ret[is.infinite(ret)] <- NA
  ret
}
dok <- \(y, d) quantile(d$Y[d$D == 0], probs = ecdf(d$Y[d$D == 1])(y)) - y
ys <- seq(-3 + sqrt(.Machine$double.eps), 3 - sqrt(.Machine$double.eps),
          length.out = 1e3)
FUNs <- c("Distributional" = "dte", "Logit-transformed" = "tte", "Doksum-type" = "dok")
make_plot_data <- \(td) {
  out <- bind_cols(lapply(FUNs, \(x) {
    do.call(x, list(y = ys, d = td))
  }))
  colnames(out) <- names(FUNs)
  out$y <-ys
  out |> pivot_longer(-y)
}

pdat2 <- bind_rows(make_plot_data(do) |> mutate(setting = "Observational"),
                   make_plot_data(di) |> mutate(setting = "Interventional"))

p2 <- ggplot(pdat2, aes(x = y, y = value, color = setting)) +
  geom_line() +
  facet_wrap(~ name, scales = "free") +
  scale_color_manual(values = colorspace::diverge_hcl(2)) +
  labs(x = "Response", y = "Treatment effect", color = "Data", tag = "B") +
  theme(legend.position = "top")

ggpubr::ggarrange(p1, p2, widths = c(2, 3))

ggsave("inst/figures/intro.pdf", height = 3, width = 9)

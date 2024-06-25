
waggr <- as.numeric(commandArgs(TRUE)[1])
if (is.na(waggr)) waggr <- 1

set.seed(12)

library("tidyverse")
library("ivreg")

raw <- read_csv("inst/data/401k.csv")

d401k <- data.frame(
  y = raw$net_tfa / 1e3,
  d = factor(raw$p401),
  z = factor(raw$e401)
) |> filter(y > 0) |> mutate(y = log(y))

res <- do.call("rbind", lapply(1:50, \(iter) {
  dat <- d401k[sample.int(nrow(d401k), 1e3), ]
  lm <- lm(y ~ d, data = dat)
  iv <- ivreg(y ~ d | z, data = dat)
  c("OLS" = unname(coef(lm)[2]), "2SLS" = unname(coef(iv)[2]))
})) |> as.data.frame()

itr <- ifelse(waggr == 2, "-max", "")
pa1 <- paste0("inst/results/figures/401k", itr, "-nd.csv")
pa2 <- paste0("inst/results/figures/401k", itr, "-pdat.csv")

nd <- read_csv(pa1)
pdat <- read_csv(pa2) |>
  mutate(model = factor(model, levels = c("DIVE", "Nonparametric"),
                        labels = c("DIVE", "CCDF")))

rmi <- pdat |> group_by(z, iter, model) |>
  summarize(p = sum(rank > 0.99) / length(rank)) |>
  filter(p > 0.5) |> pull(iter) |> unique()

nd <- nd |> filter(!iter %in% rmi)
pdat <- pdat |> filter(!iter %in% rmi)

c1 <- nd |>
  pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "cdf") |>
  group_by(d, iter, model) |>
  arrange(cdf) |>
  summarise(mean = sum(diff(cdf) * y[-1])) |>
  pivot_wider(names_from = "d", values_from = "mean") |>
  mutate(est = `1` - `0`) |>
  group_by(model) |>
  summarize(mest = mean(est), sd = sd(est) / sqrt(length(est)))

c2 <- res |> pivot_longer(everything(), names_to = "model", values_to = "est") |>
  group_by(model) |> summarise(mest = mean(est), sd = sd(est) / sqrt(length(est)))

pd <- rbind(c1, c2) |>
  mutate(model = factor(model, levels = c("OLS", "Nonparametric", "2SLS", "DIVE"),
                        labels = c("OLS", "CCDF", "2SLS", "DIVE")))

p0 <- ggplot(pd, aes(x = mest, y = model, xmin = mest - sd, xmax = mest + sd)) +
  geom_pointrange() +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  labs(y = element_blank(), x = "Estimated average causal effect") +
  theme(text = element_text(size = 13.5))

p1 <- ggplot(pdat, aes(x = rank, color = factor(z), group = interaction(z, iter))) +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "gray40") +
  facet_wrap(~ model) +
  stat_ecdf(alpha = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Estimated iPIT residual", y = "ECDF", color = "401(k) eligibility") +
  theme_bw() +
  theme(text = element_text(size = 13.5))

p2 <- nd |>
  pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "cdf") |>
  mutate(model = factor(model, levels = c("DIVE", "Nonparametric"),
                        labels = c("DIVE", "CCDF"))) |>
  ggplot(aes(x = y, y = cdf, color = factor(d), group = interaction(d, iter))) +
  facet_wrap(~ model) +
  geom_line(alpha = 0.5) +
  labs(x = "Log-transformed net total financial assets",
       y = "Estimated CDF", color = "401(k) participation") +
  theme_bw() +
  theme(text = element_text(size = 13.5)) +
  scale_color_manual(values = colorspace::diverge_hcl(2))

p3 <- nd |>
  pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "cdf") |>
  mutate(model = factor(model, levels = c("DIVE", "Nonparametric"),
                        labels = c("DIVE", "CCDF"))) |>
  pivot_wider(names_from = "d", values_from = "cdf") |>
  mutate(dce = `1` - `0`) |>
  ggplot(aes(x = y, y = dce, group = interaction(model, iter))) +
  facet_wrap(~ model) +
  geom_line(alpha = 0.3) +
  labs(x = "Log-transformed net total financial assets",
       y = "Estimated DCE", color = "Model") +
  theme_bw() +
  theme(text = element_text(size = 13.5)) +
  scale_color_manual(values = colorspace::diverge_hcl(2))

ggpubr::ggarrange(p2 + labs(tag = "A"), p1 + labs(tag = "B"),
                  p3 + labs(tag = "C"), p0 + labs(tag = "D"),
                  nrow = 2, ncol = 2, legend = "top", heights = c(0.53, 0.47))
pa3 <- paste0("inst/figures/401k-comparison", itr, ".pdf")
ggsave(pa3, height = 7, width = 12)

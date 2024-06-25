
waggr <- as.numeric(commandArgs(TRUE)[1])
if (is.na(waggr)) waggr <- 1

set.seed(12)

library("tidyverse")
library("ivreg")

data("SchoolingReturns", package = "ivreg")
SchoolingReturns$wage <- log(SchoolingReturns$wage)

res <- do.call("rbind", lapply(1:50, \(iter) {
  dat <- SchoolingReturns[sample.int(nrow(SchoolingReturns), 1e3), ]
  lm <- lm(wage ~ smsa, data = dat)
  iv <- ivreg(wage ~ smsa | nearcollege, data = dat)
  c("OLS" = unname(coef(lm)[2]), "2SLS" = unname(coef(iv)[2]))
})) |> as.data.frame()

itr <- ifelse(waggr == 2, "-max", "")
pa1 <- paste0("inst/results/figures/schooling", itr, "-nd.csv")
pa2 <- paste0("inst/results/figures/schooling", itr, "-pdat.csv")
nd <- read_csv(pa1)
pdat <- read_csv(pa2) |>
  mutate(model = factor(model, levels = c("DIVE", "Nonparametric"),
                        labels = c("DIVE", "CCDF")))

rmi <- pdat |> group_by(nearcollege, iter, model) |>
  summarize(p = sum(rank > 0.99) / length(rank)) |>
  filter(p > 0.5) |> pull(iter) |> unique()

nd <- nd |> filter(!iter %in% rmi)
pdat <- pdat |> filter(!iter %in% rmi)

c1 <- nd |>
  pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "cdf") |>
  group_by(smsa, iter, model) |>
  arrange(cdf) |>
  summarise(mean = sum(diff(cdf) * wage[-1])) |>
  pivot_wider(names_from = "smsa", values_from = "mean") |>
  mutate(est = yes - no) |>
  group_by(model) |>
  summarize(mest = mean(est), sd = sd(est)/ sqrt(length(est)))

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

p1 <- ggplot(pdat, aes(x = rank, color = nearcollege, group = interaction(nearcollege, iter))) +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "gray40") +
  facet_wrap(~ model) +
  stat_ecdf(alpha = 0.2) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Estimated iPIT residual", y = "ECDF", color = "Near college") +
  theme_bw() +
  theme(text = element_text(size = 13.5)) +
  guides(linetype = "none")

p2 <- ggplot(
  nd |> pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "cdf") |>
  mutate(model = factor(model, levels = c("DIVE", "Nonparametric"),
                        labels = c("DIVE", "CCDF"))),
  aes(x = wage, y = cdf, color = smsa, group = interaction(smsa, iter))) +
  facet_wrap(~ model) +
  geom_line(alpha = 0.2) +
  labs(x = "Log-transformed wage ", y = "Estimated CDF", color = "Metropolitan area") +
  theme_bw() +
  theme(text = element_text(size = 13.5)) +
  scale_color_manual(values = colorspace::diverge_hcl(2)) +
  guides(linetype = "none")

p3 <- nd |>
  pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "cdf") |>
  mutate(model = factor(model, levels = c("DIVE", "Nonparametric"),
                        labels = c("DIVE", "CCDF"))) |>
  pivot_wider(names_from = "smsa", values_from = "cdf") |>
  mutate(dce = yes - no) |>
  ggplot(aes(x = wage, y = dce, group = interaction(model, iter))) +
  facet_wrap(~ model) +
  geom_line(alpha = 0.3) +
  labs(x = "Log-transformed wage", y = "Estimated DCE", color = "Model") +
  theme_bw() +
  theme(text = element_text(size = 13.5)) +
  scale_color_manual(values = colorspace::diverge_hcl(2))

ggpubr::ggarrange(p2 + labs(tag = "A"), p1 + labs(tag = "B"),
                  p3 + labs(tag = "C"), p0 + labs(tag = "D"),
                  nrow = 2, ncol = 2, legend = "top", heights = c(0.53, 0.47))
pa3 <- paste0("inst/figures/schooling-comparison", itr, ".pdf")
ggsave(pa3, height = 7, width = 12)

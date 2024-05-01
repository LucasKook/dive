
set.seed(12)

library("tidyverse")
library("ivreg")

raw <- read_csv("inst/data/401k.csv")

d401k <- data.frame(
  y = raw$net_tfa / 1e3,
  d = factor(raw$p401),
  z = factor(raw$e401)
) |> filter(y > 0)

res <- do.call("rbind", lapply(1:50, \(iter) {
  dat <- d401k[sample.int(nrow(d401k), 1e3), ]
  lm <- lm(y ~ d, data = dat)
  iv <- ivreg(y ~ d | z, data = dat)
  c("OLS" = unname(coef(lm)[2]), "2SLS" = unname(coef(iv)[2]))
})) |> as.data.frame()

nd <- read_csv("inst/results/figures/401k-nd.csv")

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

ggplot(pd, aes(x = mest, y = model, xmin = mest - sd, xmax = mest + sd)) +
  geom_pointrange() +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  labs(y = element_blank(), x = "Estimated average causal effect") +
  theme(text = element_text(size = 13.5))

ggsave("inst/figures/401k-comparison.pdf", height = 3.5, width = 4.5)

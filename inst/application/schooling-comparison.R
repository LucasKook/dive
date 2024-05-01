
library("tidyverse")
library("ivreg")

data("SchoolingReturns", package = "ivreg")
SchoolingReturns$wage <- log(SchoolingReturns$wage)

res <- do.call("rbind", lapply(1:50, \(iter) {
  dat <- SchoolingReturns[sample.int(nrow(SchoolingReturns), 1e3), ]
  lm <- lm(wage ~ smsa, data = dat)
  iv <- ivreg(wage ~ smsa | nearcollege, data = dat)
  c("LM" = unname(coef(lm)[2]), "2SLS" = unname(coef(iv)[2]))
})) |> as.data.frame()

nd <- read_csv("inst/results/figures/schooling-nd.csv")

c1 <- nd |>
  pivot_longer(Nonparametric:DIVE, names_to = "model", values_to = "cdf") |>
  group_by(smsa, iter, model) |>
  summarise(mean = sum(diff(cdf) * wage[-1])) |>
  pivot_wider(names_from = "smsa", values_from = "mean") |>
  mutate(est = yes - no) |>
  group_by(model) |>
  summarize(mest = mean(est), sd = sd(est)) # / sqrt(length(est)))

c2 <- res |> pivot_longer(everything(), names_to = "model", values_to = "est") |>
  group_by(model) |> summarise(mest = mean(est), sd = sd(est)) # / sqrt(length(est)))

pd <- rbind(c1, c2) |>
  mutate(model = factor(model, levels = c("LM", "Nonparametric", "2SLS", "DIVE")))

ggplot(pd, aes(x = mest, y = model, xmin = mest - sd, xmax = mest + sd)) +
  geom_pointrange() +
  theme_bw() +
  xlim(c(0, 0.7)) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  labs(y = element_blank(), x = "Estimated average causal effect") +
  theme(text = element_text(size = 13.5))

ggsave("schooling-comparison.pdf", height = 3.5, width = 4.5)

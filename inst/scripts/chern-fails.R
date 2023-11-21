set.seed(42)

### Observational data
n <- 3e3
H <- rnorm(n)
D <- as.numeric(H > rlogis(n))
Y <- 2 * D * H - H

### Interventional data
doD <- as.numeric(rnorm(n) > rlogis(n))
doY <- 2 * doD * H - H

### Chernozhukov version
QQ <- Vectorize(\(t, dd) dd * quantile(doY[doD == 1], probs = t) +
                  (1-dd) * quantile(doY[doD == 0], probs = t))
outs <- sapply(probs <- (0:100)/100, \(x) {
  mean(Y <= QQ(x, c(0, 1))[D + 1])
})
plot(probs, outs)
abline(0, 1)

### Same results with theoretical quantiles
QQtheory <- \(t) qnorm(t)
outs <- sapply(probs, \(x) {
  mean(Y <= QQtheory(x))
})
plot(probs, outs)
abline(0, 1)

### Potential outcomes
Y0 <- -H
Y1 <- H
plot(rank(Y0)/n, rank(Y1)/n)

### Plots
library("tidyverse")
library("patchwork")
theme_set(theme_bw() + theme(text = element_text(size = 13.5)))

pd <- data.frame(Y0 = Y0, Y1 = Y1, H = H, D = D, Y = Y)[sample.int(n, 100), ]

pri <- ggplot(pd, aes(x = Y0, y = Y1)) +
  geom_point() +
  labs(x = "Y(0)", y = "Y(1)")

prs <- ggplot(pd |> pivot_longer(Y0:Y1), aes(x = name, y = pnorm(value))) +
  geom_boxplot() +
  labs(x = element_blank(), y = "Rank") +
  scale_x_discrete(labels = c("Y1" = "Y(1)", "Y0" = "Y(0)"))

pcmrs <- ggplot(pd |> pivot_longer(Y0:Y1), aes(x = plogis(H), y = pnorm(value), color = name)) +
  geom_point() +
  labs(x = "q(H)", color = element_blank(), y = "Rank") +
  scale_color_discrete(labels = c("Y1" = "Y(1)", "Y0" = "Y(0)"))

punif <- ggplot(pd, aes(x = D * pnorm(Y) + (1 - D) * pnorm(Y))) +
  stat_ecdf() +
  geom_abline(intercept = 0, slope = 1) +
  labs(y = "ECDF", x = "Interventional PIT")

(pri + labs(tag = "A", subtitle = "Rank invariance")) +
  (prs + labs(tag = "B", subtitle = "Rank similarity")) +
  (pcmrs + labs(tag = "C", subtitle = "Conditional mean rank similarity")) +
  (punif + labs(tag = "D", subtitle = "Uniformity condition")) +
  plot_layout(nrow = 1)

ggsave("inst/figures/rank-assumptions.pdf", width = 12, height = 3.5)

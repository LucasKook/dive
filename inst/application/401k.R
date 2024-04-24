### 401(k) DIVE
### LK 2024

set.seed(12)

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Data --------------------------------------------------------------------

d401k <- read_csv("inst/data/401k.csv")

dat <- data.frame(
  y = d401k$net_tfa / 1e3,
  d = d401k$p401,
  z = d401k$e401
)[sample.int(nrow(d401k), 1e3), ]

# Run ---------------------------------------------------------------------

### Fit
m0 <- BoxCox(y | d ~ 1, data = dat, support = range(dat$y), order = 10)
m <- BoxCoxDA(y | d ~ 1, data = dat, anchor = ~ z, loss = "indep",
              optimizer = optimizer_adam(0.1), order = 10, xi = 0.1)
fit(m, epochs = 1e4)

### Predict
dat$TRAM <- predict(m0, which = "distribution", type = "distribution")
dat$DIVE <- predict(m, type = "cdf")

pdat <- dat |> pivot_longer(TRAM:DIVE, names_to = "model", values_to = "rank")

p1 <- ggplot(pdat, aes(x = rank, color = factor(z))) +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "gray40") +
  facet_wrap(~ model) +
  stat_ecdf() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "PIT rank", y = "ECDF") +
  theme_bw() +
  theme(text = element_text(size = 13.5))

nd <- expand_grid(d = sort(unique(dat$d)), y = seq(
  min(dat$y), max(dat$y), length.out = 1e3))
nd$TRAM <- predict(m0, type = "distribution", newdata = nd)
nd$DIVE <- predict(m, type = "cdf", newdata = nd)

p2 <- ggplot(nd |> pivot_longer(TRAM:DIVE, names_to = "model", values_to = "cdf"),
       aes(x = y, y = cdf, color = factor(d))) +
  geom_rug(aes(x = y), data = dat, inherit.aes = FALSE, color = "gray80",
           alpha = 0.1) +
  facet_wrap(~ model) +
  geom_line() +
  labs(x = "log(wage)", y = "CDF") +
  theme_bw() +
  theme(text = element_text(size = 13.5))

ggpubr::ggarrange(p2, p1, ncol = 1, align = "hv")

ggsave("inst/figures/401k.pdf", height = 6, width = 7)

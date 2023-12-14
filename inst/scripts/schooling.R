### Return to schooling application
### LK 2023

set.seed(12)

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Data --------------------------------------------------------------------

data("SchoolingReturns", package = "ivreg")
SchoolingReturns$wage <- log(SchoolingReturns$wage)
dat <- SchoolingReturns[sample.int(nrow(SchoolingReturns), 1e3), ]

# Run ---------------------------------------------------------------------

m0 <- BoxCox(wage | smsa ~ 1, data = dat, support = range(dat$wage))
plot(m0, which = "distribution", newdata = data.frame(
  smsa = sort(unique(SchoolingReturns$smsa))), col = 1:2)

m <- BoxCoxDA(wage | smsa ~ 1, data = dat, anchor = ~ nearcollege,
              loss = "indep", optimizer = optimizer_adam(0.05),
              order = 30)
fit(m, epochs = 1e4)

plot(m, newdata = data.frame(smsa = sort(unique(SchoolingReturns$smsa))),
     type = "cdf")

dat$TRAM <- predict(m0, which = "distribution", type = "distribution")
dat$DIVE <- predict(m, type = "cdf")

pdat <- dat |> pivot_longer(TRAM:DIVE, names_to = "model", values_to = "rank")

p1 <- ggplot(pdat, aes(x = rank, color = nearcollege)) +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "gray40") +
  facet_wrap(~ model) +
  stat_ecdf() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "PIT rank", y = "ECDF") +
  theme_bw() +
  theme(text = element_text(size = 13.5))

nd <- expand_grid(smsa = sort(unique(dat$smsa)), wage = seq(
  min(dat$wage), max(dat$wage), length.out = 1e3))
nd$TRAM <- predict(m0, type = "distribution", newdata = nd)
nd$DIVE <- predict(m, type = "cdf", newdata = nd)

p2 <- ggplot(nd |> pivot_longer(TRAM:DIVE, names_to = "model", values_to = "cdf"),
       aes(x = wage, y = cdf, color = smsa)) +
  geom_rug(aes(x = wage), data = dat, inherit.aes = FALSE, color = "gray80",
           alpha = 0.1) +
  facet_wrap(~ model) +
  geom_line() +
  labs(x = "log(wage)", y = "CDF") +
  theme_bw() +
  theme(text = element_text(size = 13.5))

ggpubr::ggarrange(p2, p1, ncol = 1, align = "hv")

ggsave("inst/figures/case-study.pdf", height = 6, width = 7)

# ggplot(pdat, aes(x = model, y = rank, color = nearcollege)) +
#   geom_boxplot(position = position_dodge(0.8)) +
#   ggbeeswarm::geom_quasirandom(dodge.width = 0.8) +
#   theme_bw()
#
# ggplot(data.frame(pr = preds), aes(x = pr)) +
#   stat_ecdf() +
#   geom_abline(intercept = 0, slope = 1, linetype = 2) +
#   theme_bw() +
#   labs(subtitle = paste0("F(Y, D) indep Z: p-value = ", round(
#     dHSIC::dhsic.test(
#       preds, as.numeric(dat$nearcollege), method = "gamma")$p.value, 4)))
#
# preds <- predict(m0, which = "distribution", type = "distribution")
# ggplot(data.frame(pr = preds), aes(x = pr)) +
#   stat_ecdf() +
#   geom_abline(intercept = 0, slope = 1, linetype = 2) +
#   theme_bw() +
#   labs(subtitle = paste0("F(Y, D) indep Z: p-value = ", round(
#     dHSIC::dhsic.test(
#       preds, as.numeric(dat$nearcollege), method = "gamma")$p.value, 4)))

# Distributional causal effects
# LK 2023

# Distributional causal effect:
#   f(y, x) = \partial_x F_{Y | X = x}(y)

library("tram")
library("tidyverse")
theme_set(theme_bw() + theme(legend.position = "top", text = element_text(size = 13.5)))

# Normal linear regression ------------------------------------------------

nlrm <- function(y, x, beta = 1, sigma = 1) {
  - dnorm(y - x * beta, sd = sigma) * beta / sigma
}

curve(nlrm(x, 1, 0.5), from = -4, to = 4)

# nonparametric -----------------------------------------------------------

set.seed(1)

# g1 <- rchisq(1e3, df = 5)
# g2 <- qchisq(pnorm(rnorm(1e3) + 1), df = 5)

g1 <- rt(1e3, df = 3)
g2 <- rt(1e3, df = 30)

plot(ecdf(g1))
plot(ecdf(g2), add = TRUE, col = 2)

dte <- \(y) ecdf(g2)(y) - ecdf(g1)(y)
qte <- \(p) quantile(g2, probs = p) - quantile(g1, probs = p)
tte <- \(y, qZ = \(p) log(-log(1-p))) qZ(ecdf(g2)(y)) - qZ(ecdf(g1)(y))
rte <- \(y) log(ecdf(g2)(y)) - log(ecdf(g1)(y))
dok <- \(y) quantile(g1, probs = ecdf(g2)(y)) - y

trange <- c(max(range(g1)[1], range(g2)[1]), min(range(g1)[2], range(g2)[2]))
curve(dte(x), trange[1], trange[2])
curve(qte(x), 0, 1)
curve(tte(x), trange[1], trange[2])
curve(rte(x), trange[1], trange[2])
curve(dok(x), trange[1], trange[2])

# Survival example --------------------------------------------------------

data("GBSG2", package = "TH.data")
g1 <- GBSG2$time[GBSG2$horTh == "no"]
g2 <- GBSG2$time[GBSG2$horTh == "yes"]

plot(ecdf(g1))
plot(ecdf(g2), add = TRUE, col = 2)

trange <- c(max(range(g1)[1], range(g2)[1]), min(range(g1)[2], range(g2)[2]))
curve(dte(x), trange[1], trange[2])
curve(qte(x), 0, 1)
curve(tte(x), trange[1], trange[2])
curve(rte(x), trange[1], trange[2])
curve(dok(x), trange[1], trange[2])

# With tram ---------------------------------------------------------------

# Ignoring censoring! TODO: Update with KM estimator
GBSG2$surv <- with(GBSG2, survival::Surv(time, rep(1, nrow(GBSG2)))) # cens))
m <- Coxph(surv | 0 + horTh ~ 1, data = GBSG2, prob = c(0.001, 0.999), order = 10)
nd0 <- data.frame(
  surv = seq(0, max(GBSG2$time), length.out = 1e3),
  horTh = unique(GBSG2$horTh)[1]
)
nd1 <- data.frame(
  surv = seq(0, max(GBSG2$time), length.out = 1e3),
  horTh = unique(GBSG2$horTh)[2]
)

# Hazard
plot(nd0$surv, tte(nd0$surv), type = "l")
sTTE <- predict(m, type = "trafo", newdata = nd1) - predict(m, type = "trafo", newdata = nd0)
lines(nd0$surv, sTTE, type = "l")

# CDF
plot(nd0$surv, dte(nd0$surv), type = "l")
sDTE <- predict(m, type = "distribution", newdata = nd1) - predict(m, type = "distribution", newdata = nd0)
lines(nd0$surv, sDTE, type = "l")

# Quantile
ps <- seq(1e-6, 1 - 1e-6, length.out = 1e3)
plot(ps, qte(ps), type = "l")
sQTE <- as.double(predict(m, type = "quantile", newdata = nd1[1, ], prob = ps)) -
  as.double(predict(m, type = "quantile", newdata = nd0[1, ], prob = ps))
lines(ps, sQTE, type = "l")

# Risk ratio
plot(nd0$surv, rte(nd0$surv), type = "l")
sRTE <- log(predict(m, type = "distribution", newdata = nd1) /
              predict(m, type = "distribution", newdata = nd0))
lines(nd0$surv, sRTE, type = "l")

# Doksum
plot(nd0$surv, dok(nd0$surv), type = "l")
sDOK <- as.double(predict(m, type = "quantile", newdata = nd0[1, ],
  prob = predict(m, type = "distribution", newdata = nd1[-1, ]))) - nd0$surv[-1]
lines(nd0$surv[-1], sDOK, type = "l")

# ggplot ------------------------------------------------------------------

ggplot(GBSG2, aes(x = time, color = horTh)) +
  stat_ecdf() + labs(y = "ECDF", color = "Treatment")

ggsave("gbsg2.pdf", height = 4, width = 5)

pd <- nd0 %>%
  mutate(
    DTE = dte(surv),
    QTE = -qte(ecdf(surv)(surv)),
    RTE = rte(surv),
    TTE = tte(surv),
    DOK = dok(surv),
    method = "Nonparametric plug-in"
  ) %>% full_join(nd0 %>% mutate(
    DTE = sDTE, QTE = -sQTE,
    RTE = sRTE, TTE = sTTE, DOK = c(0, sDOK),
    method = "Transformation model"
  )) %>% pivot_longer(DTE:DOK, names_to = "type", values_to = "est") %>%
  mutate(type = factor(type, levels = c("DTE", "QTE", "TTE", "RTE", "DOK")))

ggplot(pd, aes(x = surv, y = est, color = method)) +
  geom_line() +
  facet_wrap(~ type, scales = "free", nrow = 2) +
  labs(x = "time", y = "Estimate") +
  scale_color_brewer(palette = "Dark2")

ggsave("dte.pdf", height = 6, width = 7)

# Distributional causal effects
# LK 2023

set.seed(0)

odir <- file.path("inst/figures")
if (!dir.exists(odir))
  dir.create(odir)

# Dependencies ------------------------------------------------------------

library("tram")
library("survival")
library("tidyverse")
theme_set(theme_bw() +
            theme(legend.position = "top", text = element_text(size = 13.5)))

# Normal linear regression ------------------------------------------------
# Distributional causal effect:
#   f(y, x) = \partial_x F_{Y | X = x}(y)

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
rpd <- \(p) ecdf(g2)(quantile(g1, probs = p)) - p

trange <- c(max(range(g1)[1], range(g2)[1]), min(range(g1)[2], range(g2)[2]))
curve(dte(x), trange[1], trange[2])
curve(qte(x), 0, 1)
curve(tte(x), trange[1], trange[2])
curve(rte(x), trange[1], trange[2])
curve(dok(x), trange[1], trange[2])
curve(rpd(x), 0, 1)

# Survival example --------------------------------------------------------

data("GBSG2", package = "TH.data")
g1 <- GBSG2$time[GBSG2$horTh == "no"]
g2 <- GBSG2$time[GBSG2$horTh == "yes"]

plot(ecdf(g1))
plot(ecdf(g2), add = TRUE, col = 2)

trange <- c(max(range(g1)[1], range(g2)[1]), min(range(g1)[2], range(g2)[2]))
ys <- seq(trange[1], trange[2], length.out = 1e3)
res <- lapply(1:2999, \(iter) {
  idx1 <- sample.int(length(g1), length(g1), replace = TRUE)
  idx2 <- sample.int(length(g2), length(g2), replace = TRUE)
  ecdf(g2[idx2])(ys) - ecdf(g1[idx1])(ys)
})
out <- do.call("rbind", res)

plot(ys, dte(ys), type = "s", ylim = c(-0.3, 0.1))
lines(ys, apply(out, 2, quantile, probs = 0.025), type = "s")
lines(ys, apply(out, 2, quantile, probs = 0.975), type = "s")
abline(h = 0, col = "gray80")

curve(qte(x), 0, 1, type = "s")

plot(ys, tte(ys), type = "s")
ys <- seq(trange[1], trange[2], length.out = 1e3)
res <- lapply(1:2999, \(iter) {
  idx1 <- sample.int(length(g1), length(g1), replace = TRUE)
  idx2 <- sample.int(length(g2), length(g2), replace = TRUE)
  log(-log(1-ecdf(g2[idx2])(ys))) - log(-log(1-ecdf(g1[idx1])(ys)))
})
out <- do.call("rbind", res)
lines(ys, apply(out, 2, quantile, probs = 0.025, na.rm = TRUE), type = "s")
lines(ys, apply(out, 2, quantile, probs = 0.975, na.rm = TRUE), type = "s")
abline(h = 0, col = "gray80")

curve(rte(x), trange[1], trange[2], type = "s")
curve(dok(x), trange[1], trange[2], type = "s")

# With tram ---------------------------------------------------------------

# Ignoring censoring! TODO: Update with KM estimator
GBSG2$surv <- with(GBSG2, Surv(time, rep(1, nrow(GBSG2)))) # cens))
m <- Coxph(surv | 0 + horTh ~ 1, data = GBSG2, prob = c(1e-7, 1 - 1e-7), order = 15)
nd0 <- data.frame(
  surv = seq(0, max(GBSG2$time), length.out = 1e3),
  horTh = unique(GBSG2$horTh)[1]
)
nd1 <- data.frame(
  surv = seq(0, max(GBSG2$time), length.out = 1e3),
  horTh = unique(GBSG2$horTh)[2]
)

# Hazard
plot(nd0$surv, tte(nd0$surv), type = "s")
sTTE <- predict(m, type = "trafo", newdata = nd1) - predict(m, type = "trafo", newdata = nd0)
lines(nd0$surv, sTTE, type = "l")

res <- lapply(1:999, \(iter) {
  cat("-", iter, "-")
  idx <- sample.int(nrow(GBSG2), nrow(GBSG2), replace = TRUE)
  mm <- Coxph(surv | 0 + horTh ~ 1, data = GBSG2[idx, ], prob = c(1e-6, 1 - 1e-6))
  predict(mm, type = "trafo", newdata = nd1) -
    predict(mm, type = "trafo", newdata = nd0)
})
out <- do.call("rbind", res)
lines(ys, apply(out, 2, quantile, probs = 0.025, na.rm = TRUE), type = "s")
lines(ys, apply(out, 2, quantile, probs = 0.975, na.rm = TRUE), type = "s")
abline(h = 0, col = "gray80")

# CDF
plot(nd0$surv, dte(nd0$surv), type = "s")
sDTE <- predict(m, type = "distribution", newdata = nd1) - predict(m, type = "distribution", newdata = nd0)
lines(nd0$surv, sDTE, type = "l")

# Quantile
ps <- seq(1e-6, 1 - 1e-6, length.out = 1e3)
plot(ps, qte(ps), type = "s")
sQTE <- as.double(predict(m, type = "quantile", newdata = nd1[1, ], prob = ps)) -
  as.double(predict(m, type = "quantile", newdata = nd0[1, ], prob = ps))
lines(ps, sQTE, type = "l")

# Risk ratio
plot(nd0$surv, rte(nd0$surv), type = "s")
sRTE <- log(predict(m, type = "distribution", newdata = nd1) /
              predict(m, type = "distribution", newdata = nd0))
lines(nd0$surv, sRTE, type = "l")

# Doksum
plot(nd0$surv, dok(nd0$surv), type = "s")
sDOK <- as.double(predict(m, type = "quantile", newdata = nd0[1, ],
  prob = predict(m, type = "distribution", newdata = nd1[-1, ]))) - nd0$surv[-1]
lines(nd0$surv[-1], sDOK, type = "l")

# Reparameterized DTE
nnd0 <- nd1
nnd0$p <- seq(1e-6, 1 - 1e-6, length.out = nrow(nd0))
plot(nnd0$p, rpd(nnd0$p), type = "l")
nnd0$surv <- as.double(predict(m, type = "quantile", newdata = nd0[1, ], prob = nnd0$p))
sRPD <- as.double(predict(m, type = "distribution", newdata = nnd0)) - nnd0$p
lines(nnd0$p, sRPD, type = "l")

# DTE with censoring ------------------------------------------------------

cecdf <- function(y) {
  f1 <- survfit(y ~ 1)
  Vectorize(\(y) 1 - f1$surv[rev(which(f1$time <= y))[1]])
}

GBSG2$tsurv <- with(GBSG2, Surv(time, cens))
g1 <- with(GBSG2, cecdf(tsurv[horTh == "no", ]))
g2 <- with(GBSG2, cecdf(tsurv[horTh == "yes", ]))
plot(nd0$surv, g1(nd0$surv), type = "s", ylim = c(0, 1))
lines(nd0$surv, g2(nd0$surv), type = "s", col = 2)
mm <- Coxph(tsurv | horTh ~ 1, data = GBSG2, prob = c(0.001, 0.999), order = 10)

lines(nd0$surv, p1 <- predict(mm, newdata = nd0 %>% rename(tsurv = surv), type = "distribution"))
lines(nd0$surv, p2 <- predict(mm, newdata = nd1 %>% rename(tsurv = surv), type = "distribution"), col = 2)

plot(nd0$surv, g2(nd0$surv) - g1(nd0$surv), type = "s")
lines(nd0$surv, p2 - p1, type = "s")

qMEV <- \(p) log(-log(1 - p))
plot(nd0$surv, qMEV(g2(nd0$surv)) - qMEV(g1(nd0$surv)), type = "s")
lines(nd0$surv, qMEV(p2) - qMEV(p1), type = "s")

# ggplot ------------------------------------------------------------------

ggplot(GBSG2, aes(x = time, color = horTh)) +
  stat_ecdf() + labs(y = "ECDF", color = "Treatment")

ggsave(file.path(odir, "gbsg2.pdf"), height = 4, width = 5)

g1 <- GBSG2$time[GBSG2$horTh == "no"]
g2 <- GBSG2$time[GBSG2$horTh == "yes"]

pd <- nd0 %>%
  mutate(
    DTE = dte(surv),
    QTE = -qte(ecdf(surv)(surv)),
    RTE = rte(surv),
    TTE = tte(surv),
    RPD = rpd(ecdf(surv)(surv)),
    DOK = dok(surv),
    method = "Nonparametric plug-in"
  ) %>% full_join(nd0 %>% mutate(
    DTE = sDTE, QTE = -sQTE,
    RTE = sRTE, TTE = sTTE, RPD = sRPD, DOK = c(0, sDOK),
    method = "Transformation model"
  )) %>% pivot_longer(DTE:DOK, names_to = "type", values_to = "est") %>%
  mutate(type = factor(type, levels = c("DTE", "QTE", "TTE", "RTE", "DOK", "RPD")))

ggplot(pd %>% filter(type %in% c("DTE", "QTE", "TTE", "DOK")), aes(x = surv, y = est, color = method)) +
  geom_step() +
  facet_wrap(~ type, scales = "free", nrow = 2, labeller = as_labeller(c(
    "DTE" = "Distributional", "QTE" = "Quantile", "TTE" = "Log-hazard",
    "DOK" = "Doksum-type"
  ))) +
  labs(x = "time", y = "Estimate") +
  scale_color_brewer(palette = "Dark2") +
  theme(text = element_text(size = 13.5))

ggsave(file.path(odir, "dte.pdf"), height = 6, width = 7)

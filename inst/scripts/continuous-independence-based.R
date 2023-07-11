# Continuous outcome, binary treatment, binary treatment, binary instrument
# LK 2023

set.seed(1234)

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tram")
library("deeptrafo")
library("coin")
library("ranger")
library("dHSIC")

# GEN ---------------------------------------------------------------------

df <- marginal_dgp_ex1_cont(n = 1e4)
dint <- marginal_dgp_ex1_cont(n = 1e4, doD = TRUE)

### Oracle from intervention data
m0i <- BoxCox(Y ~ 1, data = dint, prob = c(0.001, 0.999), subset = dint$D == 0, bounds = c(0, 1), order = 10)
m1i <- BoxCox(Y ~ 1, data = dint, prob = c(0.001, 0.999), subset = dint$D == 1, bounds = c(0, 1), order = 10)

### Evaluate residuals of counfounded data in interventional model
df$Ri <- NA
df$Ri[df$D == 0] <- predict(m0i, newdata = df[df$D == 0,], type = "trafo")
df$Ri[df$D == 1] <- predict(m1i, newdata = df[df$D == 1,], type = "trafo")

# boxplot(Ri ~ Z, data = df)
spearman_test(Ri ~ Z, data = df)
# dhsic.test(df$Ri, df$Z, method = "gamma")$p.value

odist <- attr(df, "odist")

# RUN ---------------------------------------------------------------------

### Fit separate models
m0 <- BoxCox(Y ~ 1, data = df, prob = c(0.001, 0.999), subset = df$D == 0, order = 15)
m1 <- BoxCox(Y ~ 1, data = df, prob = c(0.001, 0.999), subset = df$D == 1, order = 15)

### Inspect residuals
df$R <- NA
df$R[df$D == 0] <- residuals(m0)
df$R[df$D == 1] <- residuals(m1)

# boxplot(R ~ Z, data = df)
spearman_test(R ~ Z, data = df)
# dhsic.test(df$R, df$Z, method = "gamma")$p.value

### Plot models
plot(m0, which = "distribution", type = "distribution", lwd = 2, K = 3e2)
plot(m1, which = "distribution", type = "distribution", lwd = 2, col = 2, add = TRUE, K = 3e2)

### Add oracle
plot(m0i, which = "distribution", type = "distribution", lwd = 2, lty = 2, add = TRUE, K = 3e2)
plot(m1i, which = "distribution", type = "distribution", lwd = 2, col = 2, add = TRUE, lty = 2, K = 3e2)

# Control function --------------------------------------------------------

### Compute ctrl fun
cm <- glm(D ~ Z, data = df, family = "binomial")
df$ctrl <- df$D - predict(cm, type = "response")

### CTRL WITH RF
# tau <- seq(0, 1, length.out = 3e2)
# rf0 <- ranger(Y ~ ctrl, data = df[df$D == 0, ], quantreg = TRUE)
# rf1 <- ranger(Y ~ ctrl, data = df[df$D == 1, ], quantreg = TRUE)
# pr0 <- t(predict(rf0, data = df, quantiles = tau, type = "quantiles")$predictions)
# pr1 <- t(predict(rf1, data = df, quantiles = tau, type = "quantiles")$predictions)
# plot(rowMeans(pr0), tau, type = "l", lwd = 2, col = 1)
# lines(rowMeans(pr1), tau, type = "l", lwd = 2, col = 2)

ys <- seq(min(df$Y), max(df$Y), length.out = 3e2)

### CTRL WITH TRAM
# library("tidyverse")
# c0 <- deeptrafo(Y | s(ctrl) ~ 1, data = df[df$D == 0, ], trafo_options = trafo_control(support = c(0, 1)))
# fit(c0, epochs = 1e3, validation_split = 0)
# pr0 <- do.call("cbind", predict(c0, type = "cdf", q = ys, newdata = df[1:1e3, -1]))
#
# c1 <- deeptrafo(Y | s(ctrl) ~ 1, data = df[df$D == 1, ], trafo_options = trafo_control(support = c(0, 1)))
# fit(c1, epochs = 1e3, validation_split = 0)
# pr1 <- do.call("cbind", predict(c1, type = "cdf", q = ys, newdata = df[1:1e3, -1]))

c0 <- BoxCox(Y | ctrl ~ 1, data = df[df$D == 0, ], prob = c(0.001, 0.999), bounds = c(0, 1), order = 15)
c1 <- BoxCox(Y | ctrl ~ 1, data = df[df$D == 1, ], prob = c(0.001, 0.999), bounds = c(0, 1), order = 15)

df$Rc <- NA
df$Rc[df$D == 0] <- residuals(c0)
df$Rc[df$D == 1] <- residuals(c1)

# boxplot(Rc ~ Z, data = df)
spearman_test(Rc ~ Z, data = df)
# dhsic.test(df$Rc, df$Z, method = "gamma")$p.value

### Plot estimates # Does not work b/c need to extrapolate in ctrl
ys <- seq(min(df$Y), max(df$Y), length.out = 1e3)
pr0 <- predict(c0, which = "distribution", type = "distribution", newdata = droplevels(df), q = ys)
pr1 <- predict(c1, which = "distribution", type = "distribution", newdata = droplevels(df), q = ys)
plot(ys, rowMeans(pr0), col = 1, lwd = 2, type = "l")
lines(ys, rowMeans(pr1), col = 2, lwd = 2)
legend("topleft", c("D = 0", "D = 1"), col = c(1, 2), lwd = 2, bty = "n")

### Add interventional fit (dashed)
plot(m0i, which = "distribution", type = "distribution", lwd = 2, lty = 2, add = TRUE, K = 3e2)
plot(m1i, which = "distribution", type = "distribution", lwd = 2, col = 2, add = TRUE, lty = 2, K = 3e2)

### Add oracle (dotted)
lines(ys, odist(ys, d = 0), lty = 3, lwd = 2, col = 1)
lines(ys, odist(ys, d = 1), lty = 3, lwd = 2, col = 2)

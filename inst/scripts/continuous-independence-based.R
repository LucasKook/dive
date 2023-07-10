# Continuous outcome, binary treatment, binary treatment, binary instrument
# LK 2023

set.seed(1234)

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tram")
library("coin")

# GEN ---------------------------------------------------------------------

df <- marginal_dgp_ex1_cont(n = 1e3)
dint <- marginal_dgp_ex1_cont(n = 1e3, doD = TRUE)

# RUN ---------------------------------------------------------------------

### Fit separate models
m0 <- BoxCox(Y ~ 1, data = df, prob = c(0.001, 0.999), subset = df$D == 0)
m1 <- BoxCox(Y ~ 1, data = df, prob = c(0.001, 0.999), subset = df$D == 1)

### Inspect residuals
df$R <- NA
df$R[df$D == 0] <- residuals(m0)
df$R[df$D == 1] <- residuals(m1)

boxplot(R ~ Z, data = df)
independence_test(R ~ Z, data = df)

### Fit mdoels on interventional data
m0i <- BoxCox(Y ~ 1, data = dint, prob = c(0.001, 0.999), subset = dint$D == 0)
m1i <- BoxCox(Y ~ 1, data = dint, prob = c(0.001, 0.999), subset = dint$D == 1)

### Evaluate residuals of counfounded data in interventional model
df$Ri <- NA
df$Ri[df$D == 0] <- predict(m0i, newdata = df[df$D == 0,], type = "trafo")
df$Ri[df$D == 1] <- predict(m1i, newdata = df[df$D == 1,], type = "trafo")

boxplot(Ri ~ Z, data = df)
independence_test(Ri ~ Z, data = df)

### Plot estimates
plot(m0i, which = "distribution", lwd = 2)
plot(m1i, which = "distribution", add = TRUE, col = 2, lwd = 2)
legend("topleft", c("D = 0", "D = 1"), col = c(1, 2), lwd = 2, bty = "n")

### Add true distributions
ys <- seq(-4, 4, length.out = 1e3)
odist <- attr(df, "odist")
lines(ys, odist(ys), lty = 2)
lines(ys, attr(df, "odist")(ys, d = 1), lty = 2, col = 2)

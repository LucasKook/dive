
library("tram")

dgp <- function(n = 1e3, doD = FALSE) {
  Z <- rlogis(n)
  H <- rlogis(n)
  D <- as.numeric((Z + (1 - doD) * H) / 2 > rlogis(n))
  Y <- D + H + rnorm(n, sd = abs(H))
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

d <- dgp(3e3)

hist(rr <- residuals(Colr(Y | D ~ 1, data = d, order = 30)))
coin::independence_test(rr ~ d$Z, xtrafo = rank, ytrafo = rank)

PIT <- Vectorize(\(y, d, sd = 1) integrate(
  \(h)pnorm(y, mean = d + h, sd = abs(h) * sd) * dlogis(h), -Inf, Inf)$value,
  c("y", "d")
)

pit_pen <- \(r) max(ecdf(r)(sort(r)) - sort(r))
ind_pen <- \(r, e) coin::independence_test(r ~ e, xtrafo = rank, ytrafo = rank)

R <- PIT(d$Y, d$D)
R1 <- PIT(d$Y, d$D, sd = 1.5)

pit_pen(R)
pit_pen(R1)

tmp <- ind_pen(R, d$Z)
ind_pen(R1, d$Z)

hist(R)
hist(R1)

library("dare")
m <- ColrDA(Y | D ~ 1, anchor = ~ Z, data = d, order = 30, xi = 1e2,
            optimizer = optimizer_adam(0.1))
fit(m, epochs = 1e4)
plot(ecdf(residuals(m)))
abline(0.5, 0.5)

# library("tram")
# dd <- dgp(1e4, TRUE)
# coef(mm <- Colr(Y | D ~ 1, data = dd, order = 30), with_baseline = TRUE)
# hist(predict(mm, type = "distribution"))

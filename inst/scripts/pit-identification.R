set.seed(42)

library("tram")

dgp <- function(n = 1e3, doD = FALSE) {
  Z <- rlogis(n)
  H <- rlogis(n)
  D <- as.numeric((3 * Z + (1 - doD) * H) / 2 > rlogis(n))
  Y <- D + H + rnorm(n, sd = abs(H))
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

do <- dgp(1e5, doD = TRUE)
F0 <- ecdf(do$Y[do$D == 0])
F1 <- ecdf(do$Y[do$D == 1])
FF <- Vectorize(\(y, d) d * F1(y) + (1 - d) * F0(y))

d <- dgp(1e5)
F0o <- ecdf(d$Y[d$D == 0])
F1o <- ecdf(d$Y[d$D == 1])
FFo <- Vectorize(\(y, d) d * F1o(y) + (1 - d) * F0o(y))

plot(ecdf(FF(d$Y, d$D)))
abline(0, 1)

plot(ecdf(FFo(d$Y, d$D)))
abline(0, 1)

coin::independence_test(FF(d$Y, d$D) ~ d$Z, xtrafo = rank, ytrafo = rank)
coin::independence_test(FFo(d$Y, d$D) ~ d$Z, xtrafo = rank, ytrafo = rank)

# plot(ecdf(mean(d$D == 0) * ecdf(do$Y[do$D == 0])(d$Y) +
#             mean(d$D == 1) * ecdf(do$Y[do$D == 1])(d$Y)),
#      cex = 0.1)
# abline(0, 1)
#
# plot(ecdf(ecdf(do$Y[do$D == 0])(d$Y[d$D == 0])), add = TRUE)
# plot(ecdf(ecdf(do$Y[do$D == 1])(d$Y[d$D == 1])), add = TRUE)
#
# plot(ecdf(ecdf(do$Y)(d$Y[d$D == 0])), add = TRUE, cex = 0.1)
# plot(ecdf(ecdf(do$Y)(d$Y[d$D == 1])), add = TRUE, cex = 0.1)

plot(ecdf((rr <- predict(mm <- Colr(Y | D ~ 1, data = d, order = 30), type = "distribution"))))
abline(0, 1)
coin::independence_test(rr ~ d$Z, xtrafo = rank, ytrafo = rank)

PIT <- Vectorize(\(y, d, sd = 1) integrate(
  \(h) pnorm(y, mean = d + h, sd = abs(h) * sd) * dlogis(h), -Inf, Inf)$value,
  c("y", "d")
)

pit_pen <- \(r) max(ecdf(r)(sort(r)) - sort(r))
ind_pen <- \(r, e) coin::independence_test(r ~ e, xtrafo = rank, ytrafo = rank)

R <- PIT(d$Y, d$D)
R1 <- PIT(d$Y, d$D, sd = 1.5)

pit_pen(R)
pit_pen(R1)

ind_pen(R, d$Z)
ind_pen(R1, d$Z)

hist(R)
hist(R1)

library("dare")
m <- ColrDA(Y | D ~ 1, anchor = ~ Z, data = d, order = 30, xi = 1e2,
            optimizer = optimizer_adam(0.1), loss = "indep")
fit(m, epochs = 1e4)
plot(ecdf(predict(m, type = "cdf")))
abline(0, 1)

pit_pen(predict(m, type = "cdf"))
ind_pen(predict(m, type = "cdf"), d$Z)

dint <- dgp(n = 1e4, TRUE)
plot(d$Y, predict(m, type = "cdf"), cex = 0.1)
plot(ecdf(dint$Y[dint$D == 0]), col = "darkblue", add = TRUE, lwd = 2, cex = 0.1)
plot(ecdf(dint$Y[dint$D == 1]), col = "darkred", add = TRUE, lwd = 2)
plot(mm, type = "distribution", which = "distribution", K = 300, lty = 1, add = TRUE, col = "gray80")

# library("tram")
# dd <- dgp(1e4, TRUE)
# coef(mm <- Colr(Y | D ~ 1, data = dd, order = 30), with_baseline = TRUE)
# hist(predict(mm, type = "distribution"))

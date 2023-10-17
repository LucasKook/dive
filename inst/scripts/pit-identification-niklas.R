set.seed(42)

library("tram")

dgp <- function(n = 1e3, doD = FALSE) {

  # Exogenous variables
  H <- rnorm(n)
  Z <- runif(n, 0, 0.25)
  eps <- rnorm(n)*0.05

  # Endogenous variables
  D0 <- sapply(Z, function(p) rbinom(1, 1, p))
  D1 <- sapply(Z, function(p) rbinom(1, 1, (4^(1 - doD))*p))
  D <- D0
  D[H > 0] <- D1[H > 0]

  # Two possible responses
  # Y <- H + D + eps
  Y <- H*(D + 1) + eps

  ### Return
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

d <- dgp(3e3)

hist(rr <- predict(BoxCox(Y | D ~ 1, data = d, order = 30), type = "distribution"))
coin::independence_test(rr ~ d$Z, xtrafo = rank, ytrafo = rank)

PIT <- Vectorize(\(y, d, sd = 0.05) integrate(
  \(h) pnorm(y, mean = h * (1 + d), sd = sd) * dnorm(h), -Inf, Inf)$value,
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
m <- BoxCoxDA(Y | D ~ 1, anchor = ~ Z, data = d, order = 30, xi = 1e4,
              optimizer = optimizer_adam(0.1), loss = "indep")
fit(m, epochs = 1e4)
plot(ecdf(predict(m, type = "cdf")))
abline(0, 1)

pit_pen(predict(m, type = "cdf"))
ind_pen(predict(m, type = "cdf"), d$Z)

dint <- dgp(n = 1e4, TRUE)
plot(m, type = "cdf", cex = 0.1)
plot(ecdf(dint$Y[dint$D == 0]), col = "darkblue", add = TRUE, lwd = 2)
plot(ecdf(dint$Y[dint$D == 1]), col = "darkred", add = TRUE, lwd = 2)

# library("tram")
# dd <- dgp(1e4, TRUE)
# coef(mm <- Colr(Y | D ~ 1, data = dd, order = 30), with_baseline = TRUE)
# hist(predict(mm, type = "distribution"))

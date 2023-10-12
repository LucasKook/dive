
dgp <- function(n = 1e3, doD = FALSE) {
  Z <- rlogis(n)
  H <- rlogis(n)
  D <- as.numeric((Z + (1 - doD) * H) / 2 > rlogis(n))
  Y <- D + H + rnorm(n, sd = abs(H))
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

d <- dgp(1e3)

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
d$Y <- as.ordered(d$Y)
m <- PolrDA(Y | D ~ 1, data = d, xi = 1e2, anchor = ~ Z, loss = "indep",
          optimizer = optimizer_adam(0.00001))
tmp <- get_weights(m$model)
tmp[[1]][] <- c(qlogis(ecdf(d$Y)(sort(d$Y)))[-length(d$Y)], 1)
tmp[[2]][] <- c(qlogis(ecdf(d$Y)(sort(d$Y)))[-length(d$Y)], 1)
tmp[[3]][] <- 0
set_weights(m$model, tmp)
fit(m, epochs = 1e4)
c(unlist(coef(m, "int")), unlist(coef(m)))
hist(predict(m, type = "cdf"))

# library("tram")
# dd <- dgp(1e4, TRUE)
# coef(mm <- Colr(Y | D ~ 1, data = dd, order = 30), with_baseline = TRUE)
# hist(predict(mm, type = "distribution"))

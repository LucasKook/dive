
set.seed(42)

library("tram")

dgp <- function(n = 1e3, doD = FALSE) {
  Z <- rlogis(n)
  H <- rnorm(n)
  D <- as.numeric((2 * Z + (1 - doD) * 3 * H) / 2 > rlogis(n))
  Y <- D + 3 * H + rnorm(n, sd = 0.1)
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

ll <- \(parm, Y, X, Z) {
  tsd <- log(1 + exp(parm[1]))
  tbe <- parm[-1]
  nll <- mean(-dnorm(Y, mean = X %*% tbe, sd = tsd, log = TRUE))
  resids <- (Y - X %*% tbe) / tsd
  pit <- pnorm(resids)
  ppen <- mean((ecdf(pit)(pit) - pit)^2)
  ipen <- dHSIC::dhsic(resids, Z)$dHSIC
  weighted.mean(c(nll, ppen, ipen), w = c(0, 1, 1))
}

d <- dgp(1e3)
Y <- model.matrix(~ 0 + Y, data = d)
X <- model.matrix(~ D, data = d)
Z <- model.matrix(~ 0 + Z, data = d)

opt <- optim(c(0, 0, 0), ll, method = "L-BFGS-B", Y = Y, X = X, Z = Z)
opt$par

dd <- dgp(1e5, TRUE)
mm <- lm(Y ~ D, data = dd)
c(exp(1 - log(mean(residuals(mm)^2))), coef(mm))
ll(oparm <- c(exp(1 - log(sqrt(mean(residuals(mm)^2)))), coef(mm)), Y, X, Z)
opt$value

plot(ecdf(d$Y[d$D == 0]), col = "blue", main = "")
plot(ecdf(d$Y[d$D == 1]), add = TRUE, col = "red")
plot(ecdf(dd$Y[dd$D == 0]), col = "darkblue", add = TRUE)
plot(ecdf(dd$Y[dd$D == 1]), add = TRUE, col = "darkred")
ys <- seq(-30, 30, length.out = 1e3)
p0 <- \(y) pnorm(y, mean = opt$par[2], sd = log(1 + exp(opt$par[1])))
p1 <- \(y) pnorm(y, mean = opt$par[2] + opt$par[3], sd = log(1 + exp(opt$par[1])))
lines(ys, p0(ys), col = "cornflowerblue")
lines(ys, p1(ys), col = "pink")
legend("topleft", c("Observational", "Interventional", "Estimate"),
       col = c("blue", "darkblue", "cornflowerblue"), lwd = 1, title = "D = 0",
       bty = "n")
legend("left", c("Observational", "Interventional", "Estimate"),
       col = c("red", "darkred", "pink"), lwd = 1, title = "D = 1",
       bty = "n")
#
# plot(ecdf(PIT <- (1 - d$D) * p0(d$Y) + d$D * p1(d$Y)))
# abline(0, 1)
# coin::independence_test(PIT ~ Z, data = d)
#

oparm <- c(sqrt(mean(residuals(mm)^2)), coef(mm))
p0 <- Vectorize(\(y, d) pnorm(y, mean = oparm[2] + d * oparm[3], sd = oparm[1]))
# plot(ecdf(PIT <- p0(d$Y, d$D)))
# abline(0, 1)
coin::independence_test(PIT ~ Z, data = d)


# Test independence/uniform loss with NN
# LK 2023

set.seed(42)

# FUNs --------------------------------------------------------------------

dgp <- function(n = 1e3, doD = FALSE) {
  Z <- sample(0:1, n, TRUE)
  H <- rnorm(n)
  D <- as.numeric((3 * Z + 2 * (1 - doD) * H) / 2 > rlogis(n))
  # Y <- 4 * D + 4 * H + rnorm(n, sd = 1 + abs(H)/10)
  # Y <- log(1 + exp(10 + 8 * D + H * ( 4 + rlogis(n) / 10)))
  Y <- 2 * D * H - H
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

d <- dgp(1e5)

# <FIXME> These need to be conditional densities </FIXME>
D0Z0 <- ecdf(d$Y[d$D == 0 & d$Z == 0])
D1Z0 <- ecdf(d$Y[d$D == 1 & d$Z == 0])
D0Z1 <- ecdf(d$Y[d$D == 0 & d$Z == 1])
D1Z1 <- ecdf(d$Y[d$D == 1 & d$Z == 1])

p1Z0 <- mean(d$D[d$Z == 0])
p1Z1 <- mean(d$D[d$Z == 0])

ys <- seq(min(d$Y), max(d$Y))
plot(ys, (D1Z1(ys) * p1Z1)/(D0Z1(ys) * (1 - p1Z1)) -
       (D1Z0(ys) * p1Z0)/(D0Z0(ys) * (1 - p1Z0)), type = "s")
abline(h = 0)



tfun <- function(n) {
  qnorm(ecdf(smpl <- rnorm(n) + rnorm(n) + rt(n, df = 5) + rchisq(n, df = 1))(smpl))
}

plot(ecdf(tfun(1e3)))
ry <- seq(-4, 4, length.out = 1e3)
lines(ry, pnorm(ry))

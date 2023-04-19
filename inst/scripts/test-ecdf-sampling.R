
mclip <- function(x) {
  stopifnot(all(x <= 1) && all(x >= 0))
  x[x == 0] <- max(x[x != 0])
  x[x == 1] <- min(x[x != 1])
  x
}

tfun <- function(n) {
  smpl <- \(n) rnorm(n) + rnorm(n) + rt(n, df = 5) + rchisq(n, df = 1)
  tsmpl <- ecdf(smpl(n))(smpl(n))
  qnorm(mclip(tsmpl))
}

plot(ecdf(tfun(1e4)))
ry <- seq(-4, 4, length.out = 1e3)
lines(ry, pnorm(ry))

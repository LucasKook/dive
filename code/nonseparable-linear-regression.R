# Demo: Non-separable linear regression (sanity check)
# LK 03/2023

set.seed(1)

library("coin")

# FUNs --------------------------------------------------------------------

gen_dat <- function(n = 1e3, parmD = 0, parmY = c(-0.5, 1)) {
  H <- rt(n, df = 5)
  E <- sample(c(-1, 1), n, TRUE)
  ND <- rlogis(n)
  NY <- rlogis(n)
  gD <- H + ND
  gY <- H + NY
  g2u <- \(g) ecdf(g)
  UD <- g2u(gD)(gD)
  UY <- g2u(gY)(gY)
  D <- as.numeric(plogis(parmD + E) >= UD)
  UY[UY == 1] <- 1 - 1e-6
  UY[UY == 0] <- 1e-6
  Y <- qnorm(UY, mean = parmY[1] + D * parmY[2])
  data.frame(Y = Y, D = D, E = E, H = H)
}

ind_obj <- \(b, Y, X, E, tstat = "quadratic", trafo = identity) {
  R <- (Y - X %*% b)
  E <- factor(E)
  trafo(statistic(independence_test(R ~ E, teststat = tstat)))
}

ret <- replicate(1e2, {
  d <- gen_dat()
  lm(Y ~ D, data = d)
  c("LM" = optim(c(-0.5, 1), ind_obj, Y = d$Y, X = cbind(1, d$D), E = d$E)$par[2])
})

boxplot(ret)
abline(h = 1)

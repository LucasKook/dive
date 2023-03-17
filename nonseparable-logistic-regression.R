# Demo: Non-separable logistic regression
# LK 03/2023

set.seed(2)

library("coin")

# FUNs --------------------------------------------------------------------

gen_dat <- function(n = 1e3, parmD = 0, parmY = c(-0.5, 1), discreteD = TRUE,
                    discreteE = TRUE) {
  H <- rt(n, df = 5)
  E <- if (discreteE) sample(c(-1, 1), n, TRUE) else rnorm(n)
  ND <- rlogis(n)
  NY <- rlogis(n)
  gD <- H + ND
  gY <- H + NY
  g2u <- \(g) ecdf(g)
  UD <- g2u(gD)(gD)
  UY <- g2u(gY)(gY)
  D <- if (discreteD) {
    as.numeric(plogis(parmD + E) <= UD)
  } else {
    pnorm(UD, mean = parmD + E)
  }
  Y <- as.numeric(plogis(parmY[1] + D * parmY[2]) >= UY)
  data.frame(Y = Y, D = D, E = E, H = H)
}

cor_obj <- \(b, Y, X, E) {
  R <- (Y - plogis(X %*% b))
  # c(t(R) %*% prm %*% R)
  sum(fitted(lm(R ~ E))^2)
}

ind_obj <- \(b, Y, X, E, tstat = "quadratic", trafo = identity) {
  R <- (Y - plogis(X %*% b))
  trafo(statistic(independence_test(R ~ E, teststat = tstat)))
}

# d <- gen_dat(n = 1e3)
# ts <- unlist(lapply(bs <- seq(-2, 2, length.out = 1e2), \(bb) ind_obj(
#   b = c(-0.5, bb), Y = d$Y, X = cbind(1, d$D), E = d$E, tstat = "max", trafo = abs)))
# plot(bs, ts, type = "l")
hsic_obj <- \(b, Y, X, E) {
  R <- (Y - plogis(X %*% b))
  # R1 <- R[id1 <- (R > median(R))]
  # R2 <- R[id2 <- (R <= median(R))]
  # m1 <- mean(c(dHSIC:::median_bandwidth_rcpp(as.matrix(R1), length(R1), 1),
  #              dHSIC:::median_bandwidth_rcpp(as.matrix(R2), length(R2), 1)))
  # if (m1 == 0) m1 <- 0.001
  m1 <- 0.1
  m2 <- dHSIC:::median_bandwidth_rcpp(as.matrix(E), length(E), 1)
  dHSIC::dhsic(R, E, kernel = c("gaussian.fixed", "gaussian.fixed"), bandwidth = c(m1, m2))$dHSIC
}

res <- replicate(1e2, {
  d <- gen_dat()
  # d$R <- residuals(m0 <- lm(D ~ E, data = d))
  # d$R <- residuals(m0 <- glm(D ~ E, data = d, family = "binomial"), type = "response")
  # d$PR <- fitted(m0)
  # E <- cbind(1, d$E)
  # prm <- E %*% solve(t(E) %*% E) %*% t(E)
  c(
    YD = unname(coef(glm(Y ~ D, data = d, family = "binomial"))["D"]),
    # YDH = unname(coef(glm(Y ~ D + H, data = d, family = "binomial"))["D"]),
    # SRI = unname(coef(glm(Y ~ D + R, data = d, family = "binomial"))["D"]),
    # PR = unname(coef(glm(Y ~ PR, data = d, family = "binomial"))["PR"]),
    # COR = optim(c(-0.5, 0.5), cor_obj, Y = d$Y, X = cbind(1, d$D), E = E)$par[2],
    IND = optim(c(-0.5, 1), ind_obj, Y = d$Y, X = cbind(1, d$D), E = d$E)$par[2]
    # HSIC = optim(c(-0.5, 1), hsic_obj, Y = d$Y, X = cbind(1, d$D), E = d$E)$par[2]
  )
})

boxplot(t(res))
abline(h = 1, col = 2, lty = 3)

if (FALSE) {
  set.seed(12)
  d <- gen_dat()
  m <- glm(Y ~ D, data = d, family = "binomial")
  d$R <- residuals(m0 <- glm(D ~ E, data = d, family = "binomial"))
  d$RR <- residuals(m, type = "response")
  cf <- coef(glm(Y ~ I(R - D), data = d, family = "binomial"))
  d$mR <- with(d, Y - plogis(cf[1] + D * cf[2]))
  with(d, t(RR) %*% E %*% solve(t(E) %*% E) %*% t(E) %*% RR)
  with(d, t(mR) %*% E %*% solve(t(E) %*% E) %*% t(E) %*% mR)

  coin::independence_test(R ~ E, data = d)
  coin::independence_test(RR ~ E, data = d)
  coin::independence_test(mR ~ E, data = d)

  # TARGET ATE
  plogis(0.5) - plogis(-0.5)
  mean(d$Y[d$R > 0]) - mean(d$Y[d$R < 0])

  hist(replicate(1e2, {
    d <- gen_dat()
    d$R <- residuals(m0 <- glm(D ~ E, data = d, family = "binomial"))
    mean(d$Y[d$R > 0]) - mean(d$Y[d$R < 0])
  }), main = "", xlab = "")
}

# Demo: Non-separable logistic regression
# LK 03/2023

set.seed(2)

library("coin")

# FUNs --------------------------------------------------------------------

gen_dat <- function(n = 1e4, parmD = 0, parmY = c(-0.5, 1), discreteD = TRUE,
                    discreteE = TRUE, normalHDE = FALSE, conditional = FALSE,
                    discreteH = TRUE) {
  H <- if (discreteH) sample(c(-1, 1), n, TRUE) else if (normalHDE) rnorm(n)
  else rt(n, df = 5) / 2
  E <- if (discreteE) sample(c(-1, 1), n, TRUE) else rnorm(n)
  ND <- if (normalHDE) rnorm(n) else rlogis(n)
  NY <- if (normalHDE) rnorm(n) else rlogis(n)
  gD <- H + ND
  gY <- H + NY
  g2u <- if (normalHDE) \(...) {\(g) pnorm(g, sd = sqrt(2))} else \(g) ecdf(g)
  UD <- g2u(gD)(gD)
  UY <- g2u(gY)(gY)
  D <- if (discreteD) {
    if (conditional) as.numeric(plogis(parmD + E + H) <= ND)
    else as.numeric(plogis(parmD + E) <= UD)
  } else {
    pnorm(UD, mean = parmD + E)
  }
  Y <- if (conditional) as.numeric(parmY[1] + D * parmY[2] + H >= NY)
  else as.numeric(plogis(parmY[1] + D * parmY[2]) >= UY)
  data.frame(Y = Y, D = D, E = E, H = H)
}

cor_obj <- \(b, Y, X, E) {
  R <- (Y - plogis(X %*% b))
  sum(fitted(lm(R ~ E))^2)
}

ind_obj <- \(b, Y, X, E, tstat = "quadratic", trafo = identity) {
  R <- (Y - plogis(X %*% b))
  trafo(statistic(independence_test(R ~ E, teststat = tstat)))
}

hsic_obj <- \(b, Y, X, E) {
  R <- (Y - plogis(X %*% b))
  m1 <- 0.001
  m2 <- dHSIC:::median_bandwidth_rcpp(as.matrix(E), length(E), 1)
  dHSIC::dhsic(R, E, kernel = c("gaussian.fixed", "gaussian.fixed"), bandwidth = c(m1, m2))$dHSIC
}

res <- replicate(1e2, {
  d <- gen_dat(discreteD = TRUE, discreteE = TRUE, conditional = TRUE)
  d$R <- residuals(m0 <- glm(D ~ E, data = d, family = "binomial"), type = "response")
  d$PR <- fitted(m0)
  c(
    YD = unname(coef(glm(Y ~ D, data = d, family = "binomial"))["D"]),
    YDH = unname(coef(glm(Y ~ D + H, data = d, family = "binomial"))["D"]),
    SRI = unname(coef(glm(Y ~ D + R, data = d, family = "binomial"))["D"]),
    PR = unname(coef(glm(Y ~ PR, data = d, family = "binomial"))["PR"]),
    COR = optim(c(-0.5, 0.5), cor_obj, Y = d$Y, X = cbind(1, d$D), E = d$E)$par[2],
    IND = optim(c(-0.5, 1), ind_obj, Y = d$Y, X = cbind(1, d$D), E = d$E)$par[2]
  )
})

boxplot(t(res))
abline(h = 1, col = 2, lty = 3)

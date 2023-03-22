# Demo: Non-separable logistic regression
# LK 03/2023

set.seed(2)

library("coin")

# FUNs --------------------------------------------------------------------

gen_dat <- function(n = 5e3, parmD = 0, parmY = c(-0.5, 1), discreteD = TRUE,
                    discreteE = TRUE, normalHDE = FALSE, conditional = TRUE,
                    discreteH = FALSE) {
  E <- if (discreteE) sample(c(-1, 1), n, TRUE) else rnorm(n)
  H <- if (discreteH) sample(c(-1, 1), n, TRUE) else if (normalHDE) rnorm(n)
  else rt(n, df = 5) / 2
  ND <- if (normalHDE) rnorm(n) else rlogis(n)
  NY <- if (normalHDE) rnorm(n) else rlogis(n)
  gD <- parmD * H + ND
  gY <- H + NY

  g2u <- if (normalHDE) \(...) {\(g) pnorm(g, sd = sqrt(2))} else \(g) ecdf(g)
  UD <- g2u(gD)(gD)
  UY <- g2u(gY)(gY)
  D <- if (discreteD) {
    if (conditional) as.numeric(plogis(E + parmD * H) <= ND)
    else as.numeric(plogis(E) <= UD)
  } else {
    pnorm(UD, mean = parmD + E)
  }
  Y <- if (conditional) as.numeric(parmY[1] + D * parmY[2] + H >= NY)
  else as.numeric(plogis(parmY[1] + D * parmY[2]) >= UY)
  data.frame(Y = Y, D = D, E = E, H = H)
}

ind_obj <- \(b, Y, X, E, tstat = "quadratic", trafo = identity) {
  R <- (Y - plogis(X %*% b))
  trafo(statistic(independence_test(R ~ E, teststat = tstat)))
}

cor_obj <- \(b, Y, X, E) {
  R <- (Y - plogis(X %*% b))
  sum(fitted(lm(R ~ E))^2)
}

res <- replicate(5e1, {
  d0 <- gen_dat(parmD = 0)
  d1 <- gen_dat(parmD = 1)

  c(
    GLM1 = mean(coef(glm(Y ~ D, data = d1, family = "binomial")) - c(-0.5, 1)),
    GLM0 = mean(coef(glm(Y ~ D, data = d0, family = "binomial")) - c(-0.5, 1)),
    # ORC1 = mean(coef(glm(Y ~ D + H, data = d1, family = "binomial")) - c(-0.5, 1, 1)),
    # ORC0 = mean(coef(glm(Y ~ D + H, data = d0, family = "binomial")) - c(-0.5, 1, 1)),
    IND1 = mean(optim(c(-0.5, 1), ind_obj, Y = d1$Y, X = cbind(1, d1$D), E = d1$E)$par - c(-0.5, 1)),
    IND0 = mean(optim(c(-0.5, 1), ind_obj, Y = d0$Y, X = cbind(1, d0$D), E = d0$E)$par - c(-0.5, 1)),
    COR1 = mean(optim(c(-0.5, 1), cor_obj, Y = d1$Y, X = cbind(1, d1$D), E = d1$E)$par - c(-0.5, 1)),
    COR0 = mean(optim(c(-0.5, 1), cor_obj, Y = d0$Y, X = cbind(1, d0$D), E = d0$E)$par - c(-0.5, 1))
  )
})

boxplot(t(res), ylab = "bias")
abline(h = 0, lty = 2, col = 2)

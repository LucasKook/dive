# Demo: Non-separable logistic regression
# LK 03/2023

set.seed(2)

library("coin")

tn <- 1e3
nsim <- 1e2

### Defaults
oparm <- c("1" = -0.5, "D" = 1, "H" = 1)
tparmD <- 1
tcond <- TRUE
tdiscrD <- TRUE
tdiscrE <- TRUE
tnormalH <- FALSE

# Should remain TRUE for computing the oracle (conditional case)
tdiscrH <- TRUE

### Oracle
if (tcond) {
  oracle <- plogis(sum(oparm * c(1, 1, 1))) * 0.5 +
    plogis(sum(oparm * c(1, 1, 0))) * 0.5 -
    plogis(sum(oparm * c(1, 0, 1))) * 0.5 -
    plogis(sum(oparm * c(1, 0, 0))) * 0.5
} else {
  oracle <- plogis(sum(oparm * c(1, 1, 0))) -
    plogis(sum(oparm * c(1, 0, 0)))
}

# FUNs --------------------------------------------------------------------

gen_dat <- function(n = tn, parmD = tparmD, parmY = oparm[1:2],
                    discreteD = tdiscrD, discreteE = tdiscrE,
                    normalHDE = tnormalH, conditional = tcond,
                    discreteH = tdiscrH, parmH = oparm[3]) {
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
    if (conditional) as.numeric(plogis(parmD * E + parmH * H) <= ND)
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

res <- replicate(nsim, {
  ### Generate data and fit models for two-stage procedures
  d <- gen_dat(discreteD = tdiscrD, discreteE = tdiscrE, conditional = tcond,
               discreteH = tdiscrH)
  d$R <- residuals(m0 <- glm(D ~ E, data = d, family = "binomial"),
                   type = "response")
  d$PR <- fitted(m0)

  ### Marginal model
  nd <- data.frame(D = c(0, 1))
  mD <- glm(Y ~ D, data = d, family = "binomial")
  YD <- unname(diff(predict(mD, newdata = nd, type = "response")))

  ### Oracle conditional model
  ndH <- data.frame(expand.grid(D = c(0, 1), H = sort(unique(d$H))))
  ndH$pH <- rep(table(d$H) / nrow(d), each = 2)
  mDH <- glm(Y ~ D * H, data = d, family = "binomial")
  YDH <- ndH$pH * predict(mDH, newdata = ndH, type = "response")
  YDH <- diff(c(sum(YDH[ndH$D == 0]), sum(YDH[ndH$D == 1])))

  ### Residual inclusion
  ndR <- data.frame(expand.grid(D = c(0, 1), R = sort(unique(d$R))))
  ndR$pR <- rep(table(d$R) / nrow(d), each = 2)
  mSRI <- glm(Y ~ D + R, data = d, family = "binomial")
  SRI <- ndR$pR * predict(mSRI, newdata = ndR, type = "response")
  SRI <- diff(c(sum(SRI[ndR$D == 0]), sum(SRI[ndR$D == 1])))

  ### Two-stage GLM
  # ndP <- data.frame(expand.grid(PR = sort(unique(d$PR))))
  # mPR <- glm(Y ~ PR, data = d, family = "binomial")
  # PR <- diff(unname(predict(mPR, newdata = ndP, type = "response")))
  PR <- NA

  ### Foster correlation
  pCOR <- optim(c(-0.5, 0.5), cor_obj, Y = d$Y, X = cbind(1, d$D), E = d$E)$par
  nX <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
  COR <- diff(c(plogis(nX %*% pCOR)))

  ### Independence
  pIND <- optim(c(-0.5, 1), ind_obj, Y = d$Y, X = cbind(1, d$D), E = d$E)$par
  IND <- diff(c(plogis(nX %*% pIND)))

  ### Return
  c(YD = YD, YDH = YDH, SRI = SRI, PR = PR, COR = COR, IND = IND)
})

boxplot(t(res))
abline(h = oracle, col = 2, lty = 3)

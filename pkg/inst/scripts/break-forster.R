# Distributional random forest with control function
# LK March 2023

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("ranger")
library("randomForest")
library("coin")

# FUNs --------------------------------------------------------------------

gen_dat <- function(n = 1e3, doD = FALSE, nfine = 1e6) {
  ### Instrument
  Z <- sample(c(-1, 1), n, TRUE) # rt(n, df = 5)
  ### Hidden
  fH <- rt(nfine, df = 5)
  H <- rt(n, df = 5)
  ### Another covariate
  X <- rt(n, df = 10)
  ### Treatment
  UD <- runif(n)
  D <- as.numeric(plogis(Z + (1 - doD) * H + X) >= UD)
  ### Response
  UY <- runif(n)
  Y <- as.numeric(plogis(D + H + X) >= UY)
  ### Return
  data.frame(Y = Y, D = D, Z = Z, X = X, H = H)
}

cor_obj <- \(b, Y, X, E) {
  R <- (Y - plogis(X %*% b))
  sum(fitted(lm(R ~ E))^2)
}

ind_obj <- \(b, Y, X, E, tstat = "quadratic", trafo = identity) {
  R <- (Y - plogis(X %*% b))
  trafo(statistic(independence_test(R ~ E, teststat = tstat)))
}

# Sanity checks -----------------------------------------------------------

# Data under intervention on D (d0) and observational (d)
d0 <- gen_dat(1e5, doD = TRUE)
d1 <- gen_dat(1e5, doD = FALSE)

### Oracle
mDXH <- glm(Y ~ D + X + H, data = d0, family = "binomial")
d0$resp <- predict(mDXH, newdata = d0, type = "response")
op0 <- mean(d0$resp[d0$D == 0])
op1 <- mean(d0$resp[d0$D == 1])
(oOR <- log((op1 * (1 - op0)) / ((1 - op1) * op0)))

# not causal conditional OR (b/c non-collapsible)
(coef(mDX <- glm(Y ~ D + X, data = d0, family = "binomial"))["D"])

# not causal marginal OR (b/c non-collapsible)
(coef(mD <- glm(Y ~ D, data = d0, family = "binomial"))["D"])
# conditional causal OR
# exp(coef(mDH <- glm(Y ~ D + X + H, data = d0, family = "binomial"))["D"])

### Under doD should return the causal OR
(COR0 <- optim(c(0, 0), cor_obj, Y = d0$Y, X = cbind(1, d0$D), E = d0$Z)$par[2])
(IND0 <- optim(c(0, 0), ind_obj, Y = d0$Y, X = cbind(1, d0$D), E = d0$Z)$par[2])

### Even under obs should return the causal OR
(COR <- optim(c(0, 0), cor_obj, Y = d1$Y, X = cbind(1, d1$D), E = d1$Z)$par[2])
(IND <- optim(c(0, 0), ind_obj, Y = d1$Y, X = cbind(1, d1$D), E = d1$Z)$par[2])

### Naive control function (parametric, breaks down for more complex examples)
S1 <- glm(D ~ Z, data = d1, family = "binomial")
d1$R <- d1$D - predict(S1, type = "response")
S2 <- glm(Y ~ D + R, data = d1, family = "binomial")
(NCTL <- unname(coef(S2)[2]))

# Nonparametric control function ------------------------------------------

### Fit RF for control function
d <- d1
cf <- ranger(D ~ Z, data = d, probability = TRUE)
preds <- predict(cf, data = d)$predictions
d$ps <- d$D - preds[, 1]
rf <- ranger(Y ~ D + ps, data = d, probability = TRUE)
rfp <- predict(rf, data = d)$pred[, 1]
idx0 <- which(d$D == 0)
idx1 <- which(d$D == 1)
p0 <- mean(rfp[idx0])
p1 <- mean(rfp[idx1])

# Output ------------------------------------------------------------------

c(
  ORACLE = oOR, mD = unname(coef(mD)["D"]), mDX = unname(coef(mDX)["D"]),
  COR = COR, IND = IND, NCTL = NCTL, RFCF = log((p1 * (1 - p0)) / ((1 - p1) * p0))
)

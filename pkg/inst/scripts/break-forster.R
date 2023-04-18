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

OR <- \(p1, p2, cf = identity) {
  cf((p1 * (1 - p2)) / ((1 - p1) * p2))
}

# Sanity checks -----------------------------------------------------------

# Data under intervention on D (d0) and observational (d)
d0 <- gen_dat(1e5, doD = TRUE)
d1 <- gen_dat(1e5, doD = FALSE)

### Oracle
mDXH0 <- glm(Y ~ D + X + H, data = d0, family = "binomial")
nd0 <- nd1 <- d1
nd0$D <- 0
nd1$D <- 1
(OR(mean(predict(mDXH0, newdata = nd1, type = "response")),
           mean(predict(mDXH0, newdata = nd0, type = "response")), log))

# Include hidden confounder and marginalize
mDXH <- glm(Y ~ D + X + H, data = d1, family = "binomial")
nd0 <- nd1 <- d1
nd0$D <- 0
nd1$D <- 1
(oOR <- OR(mean(predict(mDXH, newdata = nd1, type = "response")),
           mean(predict(mDXH, newdata = nd0, type = "response")), log))

# not causal conditional OR (b/c non-collapsible)
(coef(mDX <- glm(Y ~ D + X, data = d0, family = "binomial"))["D"])

# not causal marginal OR (b/c non-collapsible)
(coef(mD <- glm(Y ~ D, data = d0, family = "binomial"))["D"])
# conditional causal OR
# exp(coef(mDH <- glm(Y ~ D + X + H, data = d0, family = "binomial"))["D"])

### Under obs should return the causal conditional OR
(COR <- optim(c(0, 0), cor_obj, Y = d1$Y, X = cbind(1, d1$D), E = d1$Z)$par[2])
(IND <- optim(c(0, 0), ind_obj, Y = d1$Y, X = cbind(1, d1$D), E = d1$Z)$par[2])

### Naive control function (parametric, breaks down for more complex examples)
S1 <- glm(D ~ Z, data = d1, family = "binomial")
d1$R <- d1$D - predict(S1, type = "response")
S2 <- glm(Y ~ D + X + R, data = d1, family = "binomial")
nd0 <- nd1 <- d1
nd0$D <- 0
nd1$D <- 1
S2p0 <- predict(S2, newdata = nd0, type = "response")
S2p1 <- predict(S2, newdata = nd1, type = "response")
(NCTL <- OR(mean(S2p1), mean(S2p0), log))
# (NCTL <- unname(coef(S2)[2])) # non-collapsible!

# Nonparametric control function ------------------------------------------

cf <- ranger(D ~ Z, data = d1, probability = TRUE)
preds <- predict(cf, data = d1)$predictions
d1$ps <- d1$D - preds[, 1]
rf <- ranger(Y ~ D + X + ps, data = d1, probability = TRUE)
nd0 <- nd1 <- d1
nd0$D <- 0
nd1$D <- 1
rfp0 <- predict(rf, data = nd0)$pred[, 1]
rfp1 <- predict(rf, data = nd1)$pred[, 1]
p0 <- mean(rfp0)
p1 <- mean(rfp1)

# Output ------------------------------------------------------------------

c(
  ORACLE = oOR, mD = unname(coef(mD)["D"]), mDX = unname(coef(mDX)["D"]),
  COR = COR, IND = IND, NCTL = NCTL, RFCF = OR(p1, p0, log)
)

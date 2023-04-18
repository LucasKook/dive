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
  ### Treatment
  # gD <- (1 - doD) * fH + rlogis(nfine)
  # UD <- dive:::.clip(ecdf(gD)((1 - doD) * H + rlogis(n)))
  UD <- runif(n)
  D <- as.numeric(plogis(Z + (1 - doD) * H) >= UD)
  ### Response
  # gY <- fH + rnorm(nfine)
  # UY <- dive:::.clip(ecdf(gY)(H + rnorm(n)))
  UY <- runif(n)
  Y <- as.numeric(plogis(D + H) >= UY)
  ### Return
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

cor_obj <- \(b, Y, X, E) {
  R <- (Y - plogis(X %*% b))
  sum(fitted(lm(R ~ E))^2)
}

ind_obj <- \(b, Y, X, E, tstat = "quadratic", trafo = identity) {
  R <- (Y - plogis(X %*% b))
  trafo(statistic(independence_test(R ~ E, teststat = tstat)))
}

# Oracle ------------------------------------------------------------------

op1 <- integrate(\(x) plogis(1 + x) * dt(x, df = 5), -20, 20)$value
op0 <- integrate(\(x) plogis(x) * dt(x, df = 5), -20, 20)$value

### Oracle ATE (collapsible)
(oATE <- op1 - op0)

### Oracle OR (non-collapsible)
(oOR <- (op1 * (1 - op0)) / ((1 - op1) * op0))

# Sanity checks -----------------------------------------------------------

# Data under intervention on D (d0) and observational (d)
d0 <- gen_dat(1e5, doD = TRUE)
d1 <- gen_dat(1e5, doD = FALSE)

# causal OR b/c nonparametric (everything binary)
exp(coef(mD <- glm(Y ~ D, data = d0, family = "binomial"))["D"])
diff(predict(mD, newdata = data.frame(D = c(0, 1)), type = "response"))
# conditional causal OR
# exp(coef(mDH <- glm(Y ~ D + H, data = d0, family = "binomial"))["D"])
# diff(predict(mDH, newdata = data.frame(D = c(1, 0), H = 0), type = "response"))

### Under doD should return the causal OR
exp(COR0 <- optim(c(0, 0), cor_obj, Y = d0$Y, X = cbind(1, d0$D), E = d0$Z)$par)[2]
exp(IND0 <- optim(c(0, 0), ind_obj, Y = d0$Y, X = cbind(1, d0$D), E = d0$Z)$par)[2]

### Even under obs should return the causal OR
exp(COR <- optim(c(0, 0), cor_obj, Y = d1$Y, X = cbind(1, d1$D), E = d1$Z)$par)[2]
exp(IND <- optim(c(0, 0), ind_obj, Y = d1$Y, X = cbind(1, d1$D), E = d1$Z)$par)[2]

### Naive control function (parametric, breaks down for more complex examples)
S1 <- glm(D ~ Z, data = d1, family = "binomial")
d1$R <- d1$D - predict(S1, type = "response")
S2 <- glm(Y ~ D + R, data = d1, family = "binomial")
NCTL <- unname(exp(coef(S2)[2]))

# Nonparametric control function ------------------------------------------

### Fit RF for control function
d <- d1
cf <- ranger(D ~ Z, data = d, probability = TRUE)
preds <- predict(cf, data = d)$predictions
d$ps <- d$D - preds[, 1]
# Virtually identical to RF-based PS
# cf <- glm(D ~ Z, data = dZ, family = "binomial")
# preds <- predict(cf, newdata = d, type = "response")
# d$ps <- d$D - preds

### Fit RF with control function prediction and compute RF weights for prediction
rf <- ranger(Y ~ D + ps, data = d)
nd0 <- nd1 <-  d
nd0$D <- 0
nd1$D <- 1
p0 <- mean(predict(rf, data = nd0)$pred)
p1 <- mean(predict(rf, data = nd1)$pred)

# Results -----------------------------------------------------------------

### ATE
DM <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
c(ORACLE = oATE, RFCF = p1 - p0, COR = diff(plogis(DM %*% COR)),
  IND = diff(plogis(DM %*% IND)))

### OR
c(ORACLE = oOR, RFCF = ((p1 * (1 - p0)) / ((1 - p1) * p0)),
  COR = exp(COR[2]), IND = exp(IND[2]), NCTL = NCTL)


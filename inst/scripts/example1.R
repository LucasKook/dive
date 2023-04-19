# Distributional random forest with control function
# LK March 2023

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("ranger")
library("randomForest")
library("coin")
devtools::load_all()

# Data --------------------------------------------------------------------

n <- 1e4

### Data under intervention on D (d0) and observational (d1)
d0 <- dgp_ex1_binary(n, doD = TRUE)
d1 <- dgp_ex1_binary(n, doD = FALSE)

# Oracle ------------------------------------------------------------------

op0 <- attr(d1, "p0")
op1 <- attr(d1, "p1")

### Oracle ATE (collapsible)
(oATE <- op1 - op0)

### Oracle OR (non-collapsible)
(oOR <- OR(op1, op0))

### causal OR b/c nonparametric (everything binary)
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

# New data for prediction
nd0 <- nd1 <-  d1
nd0$D <- 0
nd1$D <- 1

S2 <- glm(Y ~ D + R, data = d1, family = "binomial")
NCTL <- OR(mean(predict(S2, newdata = nd1, type = "response")),
           mean(predict(S2, newdata = nd0, type = "response")))

# Nonparametric control function ------------------------------------------

### Fit RF for control function
cf <- ranger(D ~ Z, data = d1, probability = TRUE)
preds <- predict(cf, data = d1)$predictions
d1$ps <- d1$D - preds[, 1]

# Virtually identical to RF-based PS
# cf <- glm(D ~ Z, data = dZ, family = "binomial")
# preds <- predict(cf, newdata = d, type = "response")
# d$ps <- d$D - preds

# New data for prediction
nd0 <- nd1 <-  d1
nd0$D <- 0
nd1$D <- 1

### Fit RF with control function prediction and compute RF weights for prediction
rf <- ranger(Y ~ D + ps, data = d1)
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

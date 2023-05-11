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
dtune <- dgp_ex1_binary(n, doD = FALSE)

# Oracle ------------------------------------------------------------------

op0 <- attr(d1, "p0")
op1 <- attr(d1, "p1")

### Oracle ATE (collapsible)
(oATE <- ATE(op1, op0))
# (oATE <- mean(d0$Y[d0$D == 1]) - mean(d0$Y[d0$D == 0]))

### Oracle OR (non-collapsible)
(oOR <- OR(op1, op0, log))

### causal OR b/c nonparametric (everything binary)
coef(mD <- glm(Y ~ D, data = d0, family = "binomial"))["D"]
diff(predict(mD, newdata = data.frame(D = c(0, 1)), type = "response"))
# conditional causal OR
# exp(coef(mDH <- glm(Y ~ D + H, data = d0, family = "binomial"))["D"])
# diff(predict(mDH, newdata = data.frame(D = c(1, 0), H = 0), type = "response"))

### Under doD should return the causal OR
parCOR0 <- indep_iv(Y ~ D, ~ Z, d0, "COR")
pCOR0 <- indep_marginal_predictions(parCOR0, d0)
COR0 <- OR(pCOR0[, "p1"], pCOR0[, "p0"], log)
parIND0 <- indep_iv(Y ~ D, ~ 0 + Z, d0, "IND", ytrafo = rank)
pIND0 <- indep_marginal_predictions(parIND0, d0)
IND0 <- OR(pIND0[, "p1"], pIND0[, "p0"], log)

### Even under obs should return the causal OR
parCOR <- indep_iv(Y ~ D, ~ Z, d1, "COR")
pCOR <- indep_marginal_predictions(parCOR, d1)
COR <- OR(pCOR[, "p1"], pCOR[, "p0"], log)
parIND <- indep_iv(Y ~ D, ~ 0 + Z, d1, "IND", ytrafo = rank)
pIND <- indep_marginal_predictions(parIND, d1)
IND <- OR(pIND[, "p1"], pIND[, "p0"], log)

### Naive control function (parametric, breaks down for more complex examples)
S1 <- glm(D ~ Z, data = d1, family = "binomial")
d1$R <- d1$D - predict(S1, type = "response")
S2 <- glm_marginal_predictions(Y ~ D + R, data = d1)
NCTL <- OR(S2[, "p1"], S2[, "p0"], log)

# Nonparametric control function ------------------------------------------

### Fit RF for control function
cf <- ranger(factor(D) ~ Z, data = dtune, probability = TRUE)
preds <- predict(cf, data = d1)$predictions
d1$ps <- d1$D - preds[, 2]
# Virtually identical to RF-based PS
# cf <- glm(D ~ Z, data = dZ, family = "binomial")
# preds <- predict(cf, newdata = d, type = "response")
# d$ps <- d$D - preds

pRF <- ranger_marginal_predictions(factor(Y) ~ D + ps, d1)
RF <- OR(pRF[, "p1"], pRF[, "p0"], log)

# Results -----------------------------------------------------------------

### ATE
c(ORACLE = oATE, RFCF = diff(-c(pRF)), COR = diff(-c(pCOR)[1:2]),
  IND = diff(-c(pIND)[1:2]), NCTL = diff(-c(S2)[1:2]))

### OR
c(ORACLE = oOR, RFCF = RF, COR = COR, IND = IND, NCTL = NCTL)

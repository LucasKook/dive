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

rfk <- function(rf, data) {
  preds <- predict(rf, data = d, type = "terminalNodes")
  inbag <- simplify2array(rf$inbag.counts)
  rfw <- matrix(nrow = nrow(d), ncol = nrow(d))
  for (i in 1:nrow(d)) {
    for (j in (1:nrow(d))) {
      tree_idx <- inbag[i, ] == 0 & inbag[j, ] == 0
      rfw[i, j] <- sum(preds$predictions[i, tree_idx] ==
                          preds$predictions[j, tree_idx]) / sum(tree_idx)
    }
  }
  rfw / matrix(pmax(colSums(rfw), .Machine$double.eps),
               ncol = nrow(d), nrow = nrow(d), byrow = TRUE)
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

### Generate data
# learn <- seq_len(floor(nrow(d) / 2) + 1)
# test <- seq_len(nrow(d))[-learn]

### Fit RF for control function
d <- d1[1:1e3,]
dZ <- d1[-1:-1e3,]
cf <- ranger(D ~ Z, data = dZ, probability = TRUE)
preds <- predict(cf, data = d)$predictions
d$ps <- d$D - preds[, 1]
# cf <- glm(D ~ Z, data = dZ, family = "binomial")
# preds <- predict(cf, newdata = d, type = "response")
# d$ps <- d$D - preds

### Fit RF with control function prediction and compute RF weights for prediction
rf <- ranger(Y ~ D + ps, data = d, keep.inbag = TRUE)
rfw <- rfk(rf, d)

# Results -----------------------------------------------------------------

idx0 <- which(d$D == 0)
idx1 <- which(d$D == 1)

pcdf <- Vectorize(\(y, idx) c(t(rfw[, idx]) %*% as.numeric(d$Y <= y)), "y")

p0 <- colMeans(pcdf(0, idx0))
p1 <- colMeans(pcdf(0, idx1))

### ATE
DM <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
c(ORACLE = oATE, CF = p0 - p1, COR = diff(plogis(DM %*% COR)),
  IND = diff(plogis(DM %*% IND)))

### OR
c(ORACLE = oOR, CF = 1 / ((p1 * (1 - p0)) / ((1 - p1) * p0)),
  COR = exp(COR[2]), IND = exp(IND[2]), NCTL = NCTL)

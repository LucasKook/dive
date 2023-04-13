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
  Z <- sample(0:1, n, TRUE) # rt(n, df = 5)
  ### Hidden
  fH <- rt(nfine, df = 5)
  H <- rt(n, df = 5)
  ### Treatment
  # gD <- (1 - doD) * fH + rlogis(nfine)
  # UD <- dive:::.clip(ecdf(gD)((1 - doD) * H + rlogis(n)))
  UD <- runif(n)
  D <- as.numeric(plogis(Z + doD * H) >= UD)
  ### Response
  # gY <- fH + rnorm(nfine)
  # UY <- dive:::.clip(ecdf(gY)(H + rnorm(n)))
  UY <- runif(n)
  Y <- as.numeric(plogis(D + H) >= UY)
  ### Return
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

rfk <- function(rf, data, idx) {
  preds <- predict(rf, data = d, type = "terminalNodes")
  prox <- matrix(nrow = nrow(d), ncol = nrow(d))
  for (i in idx) {
    for (j in (1:nrow(d))[-idx]) {
      prox[i, j] <- mean(preds$predictions[i, ] == preds$predictions[j, ])
    }
  }
  rfw <- prox[idx,][,-idx]
  rfw / matrix(pmax(colSums(rfw), .Machine$double.eps),
               ncol = nrow(d) - length(idx),
               nrow = length(idx), byrow = TRUE)
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
oATE <- op1 - op0
oOR <- (op1 * (1 - op0)) / ((1 - op1) * op0)

### Non-collapsibility
# coef(glm(Y ~ D, data = gen_dat(1e5, doD = TRUE), family = "binomial"))["D"]

# Run ---------------------------------------------------------------------

### Generate data
d <- gen_dat()
learn <- seq_len(floor(nrow(d) / 2) + 1)
test <- seq_len(nrow(d))[-learn]

### COR and IND
optim(c(0, 0), cor_obj, Y = d$Y, X = cbind(1, d$D), E = d$Z)$par

### Fit RF for control function
cf <- ranger(D ~ Z, data = d[learn,], probability = TRUE)
preds <- dive:::.clip(predict(cf, data = d)$predictions)
d$ps <- d$D - 1 - preds[, 2]

### Fit RF with control function prediction and compute RF weights for prediction
rf <- ranger(Y ~ D + ps, data = d[test,], probability = TRUE)
rfw <- rfk(rf, d, test)

# Predict and plot --------------------------------------------------------

idx0 <- which(d$D[learn] == 0)
idx1 <- which(d$D[learn] == 1)

pcdf <- Vectorize(\(y, idx) c(t(rfw[, idx]) %*% as.numeric(d$Y[test] <= y)), "y")

p0 <- colMeans(pcdf(0, idx0))
p1 <- colMeans(pcdf(0, idx1))

### ATE
c(ORACLE = oATE, CF = p0 - p1, COR = NA, IND = NA)

### OR
c(ORACLE = oOR, CF = 1 / ((p1 * (1 - p0)) / ((1 - p1) * p0)),
  COR = NA, IND = NA)

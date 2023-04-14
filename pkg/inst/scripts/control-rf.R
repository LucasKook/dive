# Distributional random forest with control function
# LK March 2023

set.seed(1)

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
  D <- as.numeric(plogis(Z + H) >= UD)
  ### Response
  # gY <- fH + rnorm(nfine)
  # UY <- dive:::.clip(ecdf(gY)(H + rnorm(n)))
  UY <- runif(n)
  Y <- qnorm(UY, mean = D + H) # , sd = 1 + D)
  ### Return
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

# rfk <- function(rf, data, idx) {
#   preds <- predict(rf, data = d, type = "terminalNodes")
#   prox <- matrix(nrow = nrow(d), ncol = nrow(d))
#   for (i in idx) {
#     for (j in (1:nrow(d))[-idx]) {
#       prox[i, j] <- mean(preds$predictions[i, ] == preds$predictions[j, ])
#     }
#   }
#   rfw <- prox[idx,][,-idx]
#   rfw / matrix(pmax(colSums(rfw), .Machine$double.eps),
#                ncol = nrow(d) - length(idx),
#                nrow = length(idx), byrow = TRUE)
# }
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

# Run ---------------------------------------------------------------------

### Generate data
d <- gen_dat()
learn <- seq_len(floor(nrow(d) / 2) + 1)
test <- seq_len(nrow(d))[-learn]

### Fit RF for control function
# cf <- glm(D ~ Z, data = d[learn,], family = "binomial")
# preds <- predict(cf, type = "response", newdata = d)
# d$ps <- d$D - 1 - preds
cf <- ranger(D ~ Z, data = d, probability = TRUE)
preds <- dive:::.clip(predict(cf, data = d)$predictions)
d$ps <- d$D - 1 - preds[, 2]

# coef(lm(Y ~ D + ps, data = d[-learn,]))["D"]

### Fit RF with control function prediction and compute RF weights for prediction
rf <- ranger(Y ~ D + ps, data = d, keep.inbag = TRUE)
rfw <- rfk(rf, d)

# Predict and plot --------------------------------------------------------

idx0 <- which(d$D == 0)
idx1 <- which(d$D == 1)

pcdf <- Vectorize(\(y, idx) c(t(rfw[, idx]) %*% as.numeric(d$Y <= y)), "y")
ts <- seq(min(d$Y), max(d$Y), length.out = 1e2)

gt <- Vectorize(\(y, idx, m = 0) pnorm(y, mean = m + d$H[idx]), "y")

plot(ts, colMeans(pcdf(ts, idx0)), type = "l")
lines(ts, colMeans(pcdf(ts, idx1)), type = "l", col = 2)
lines(ts, colMeans(gt(ts, idx0)), lty = 2)
lines(ts, colMeans(gt(ts, idx1, m = 1)), #, sd = 2),
      lty = 2, col = 2)

legend("topleft", c("Truth", "Estimated"), lty = c(2, 1), bty = "n")

# R <- Vectorize(\(y) mean(t(rfw) %*% as.numeric(d$Y[learn] <= y)))(d$Y[-learn])
# independence_test(R ~ Z, data = d[-learn,])

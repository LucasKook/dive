# Distributional random forest with control function
# LK March 2023

set.seed(1)

# Dependencies ------------------------------------------------------------

library("ranger")
library("randomForest")
library("coin")

# FUNs --------------------------------------------------------------------

clipU <- function(x, eps = 1e-6) {
  x[x == 0] <- eps
  x[x == 1] <- 1 - eps
  x
}

gen_dat <- function(n = 2e3, doD = FALSE) {
  ### Instrument
  Z <- rt(n, df = 5)
  ### Hidden
  H <- rt(n, df = 5)
  ### Treatment
  ND <- rlogis(n)
  gHD <- (1 - doD) * H + ND
  UD <- clipU(ecdf(gHD)(gHD))
  D <- as.numeric(plogis(Z) >= UD)
  ### Response
  NY <- rnorm(n)
  gHY <- H + NY
  UY <- clipU(ecdf(gHY)(gHY))
  UY[UY == 0] <- 1e-6
  UY[UY == 1] <- 1 - 1e-6
  Y <- qnorm(UY, mean = D, sd = 1 + D)
  ### Return
  data.frame(Y = Y, D = D, Z = Z, H = H)
}

# Run ---------------------------------------------------------------------

### Generate data
d <- gen_dat()
learn <- seq_len(floor(nrow(d) / 2) + 1)

### Fit RF for control function
cf <- ranger(D ~ Z, data = d[learn,], probability = TRUE)
preds <- clipU(predict(cf, data = d)$predictions)
d$ps <- d$D - preds[, 2]

### Fit RF with control function prediction and compute RF weights for prediction
rf <- randomForest(Y ~ D + ps, data = d[-learn,])
rfw <- predict(rf, newdata = d, proximity = TRUE)$proximity[learn,][,-learn]
rfw <- rfw / matrix(pmax(colSums(rfw), .Machine$double.eps),
                    ncol = nrow(d) - length(learn),
                    nrow = length(learn), byrow = TRUE)

# Predict and plot --------------------------------------------------------

idx0 <- which(d$D[-learn] == 0)
idx1 <- which(d$D[-learn] == 1)

pcdf <- Vectorize(\(y, idx) c(t(rfw[, idx]) %*% as.numeric(d$Y[learn] <= y)), "y")
ts <- seq(min(d$Y), max(d$Y), length.out = 1e2)

plot(ts, colMeans(pcdf(ts, idx0)), type = "l")
lines(ts, colMeans(pcdf(ts, idx1)), type = "l", col = 2)
lines(ts, pnorm(ts), lty = 2)
lines(ts, pnorm(ts, mean = 1, sd = 2), lty = 2, col = 2)

# R <- Vectorize(\(y) mean(t(rfw) %*% as.numeric(d$Y[learn] <= y)))(d$Y[-learn])
# independence_test(R ~ Z, data = d[-learn,])

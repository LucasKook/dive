
library("ranger")

rfk <- function(rf, data) {
  preds <- predict(rf, data = d, type = "terminalNodes")
  inbag <- simplify2array(rf$inbag.counts)
  rfw <- matrix(0, nrow = nrow(d), ncol = nrow(d))
  for (i in 1:nrow(d)) {
    for (j in 1:i) {
      tree_idx <- inbag[i, ] == 0 & inbag[j, ] == 0
      rfw[i, j] <- sum(preds$predictions[i, tree_idx] ==
                         preds$predictions[j, tree_idx]) / sum(tree_idx)
    }
  }
  rfw <- rfw + t(rfw - diag(nrow = nrow(d)))
  rfw / matrix(pmax(colSums(rfw), .Machine$double.eps),
               ncol = nrow(d), nrow = nrow(d), byrow = TRUE)
}

x <- runif(1e3)
y <- sin(pi * x) + rnorm(1e3, sd = 0.1)

rf <- ranger(y ~ x, data = d <- data.frame(y = y, x = x), keep.inbag = TRUE)
out <- rfk(rf, d)

plot(d$y, preds <- c(t(out) %*% d$y))
abline(0, 1)

mean((d$y - sin(pi * d$x))^2)
mean((d$y - preds)^2)

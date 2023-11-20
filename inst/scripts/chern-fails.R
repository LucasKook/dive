set.seed(42)

### Observational data
n <- 3e3
H <- rnorm(n)
D <- as.numeric(H > rlogis(n))
Y <- 2 * D * H - H

### Interventional data
doD <- as.numeric(rnorm(n) > rlogis(n))
doY <- 2 * doD * H - H

### Chernozhukov version
QQ <- Vectorize(\(t, dd) dd * quantile(doY[doD == 1], probs = t) +
                  (1-dd) * quantile(doY[doD == 0], probs = t))
outs <- sapply(probs <- (0:100)/100, \(x) {
  mean(Y <= QQ(x, c(0, 1))[D + 1])
})
plot(probs, outs)
abline(0, 1)

### Same results with theoretical quantiles
QQtheory <- \(t) qnorm(t)
outs <- sapply(probs, \(x) {
  mean(Y <= QQtheory(x))
})
plot(probs, outs)
abline(0, 1)


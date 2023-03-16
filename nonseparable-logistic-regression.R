# Demo: Non-separable logistic regression
# LK 03/2023

set.seed(-42)

# FUNs --------------------------------------------------------------------

gen_dat <- function(n = 3e3, parmD = 0, parmY = c(-0.5, 1)) {
  H <- rt(n, df = 5)
  E <- sample(c(-1, 1), n, TRUE)
  ND <- rlogis(n)
  NY <- rlogis(n)
  gD <- H + ND
  gY <- H + NY
  g2u <- \(g) ecdf(g)
  UD <- g2u(gD)(gD)
  UY <- g2u(gY)(gY)
  D <- as.numeric(plogis(parmD + E) >= UD)
  Y <- as.numeric(plogis(parmY[1] + D * parmY[2]) >= UY)
  data.frame(Y = Y, D = D, E = E, H = H)
}

res <- replicate(1e2, {
  d <- gen_dat()
  d$R <- residuals(m0 <- glm(D ~ E, data = d, family = "binomial"), type = "response")
  d$PR <- fitted(m0)
  c(
    YD = unname(coef(glm(Y ~ D, data = d, family = "binomial"))["D"]),
    YDH = unname(coef(glm(Y ~ D + H, data = d, family = "binomial"))["D"]),
    SRI = unname(coef(glm(Y ~ D + R, data = d, family = "binomial"))["D"]),
    PR = unname(coef(glm(Y ~ PR, data = d, family = "binomial"))["PR"])
  )
})

boxplot(t(res))
abline(h = 1, col = 2, lty = 3)

if (FALSE) {
  set.seed(12)
  d <- gen_dat()
  m <- glm(Y ~ D, data = d, family = "binomial")
  d$R <- residuals(m0 <- glm(D ~ E, data = d, family = "binomial"))
  d$RR <- residuals(m, type = "response")
  cf <- coef(glm(Y ~ I(R - D), data = d, family = "binomial"))
  d$mR <- with(d, Y - plogis(cf[1] + D * cf[2]))
  with(d, t(RR) %*% E %*% solve(t(E) %*% E) %*% t(E) %*% RR)
  with(d, t(mR) %*% E %*% solve(t(E) %*% E) %*% t(E) %*% mR)

  coin::independence_test(R ~ E, data = d)
  coin::independence_test(RR ~ E, data = d)
  coin::independence_test(mR ~ E, data = d)

  # TARGET ATE
  plogis(0.5) - plogis(-0.5)
  mean(d$Y[d$R > 0]) - mean(d$Y[d$R < 0])

  hist(replicate(1e2, {
    d <- gen_dat()
    d$R <- residuals(m0 <- glm(D ~ E, data = d, family = "binomial"))
    mean(d$Y[d$R > 0]) - mean(d$Y[d$R < 0])
  }), main = "", xlab = "")
}

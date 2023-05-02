
set.seed(3)
devtools::load_all()

dat <- marginal_dgp_ex1_binary(n = 1e4, doD = TRUE)
# dat <- dgp_ex1_binary(n = 1e4)

dat$Y <- factor(dat$Y)
dat$D <- factor(dat$D)
dat$Z <- factor(dat$Z)

m <- glm(D ~ Z, data = dat, family = "binomial")
nD <- simulate(m, nsim = 1e2)

cfx0 <- coef(glm(Y ~ D, data = dat, family = "binomial"))

res <- lapply(nD, \(nd) {
  tmp <- dat
  tmp$D <- nd
  cfx <- coef(glm(Y ~ D, data = tmp, family = "binomial"))
  plogis(sum(cfx)) - plogis(cfx[1])
})

rres <- do.call("rbind", res)
boxplot(rres)

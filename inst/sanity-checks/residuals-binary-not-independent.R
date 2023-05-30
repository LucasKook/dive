
devtools::load_all()

dgp <- function(n = 1e3, doD = FALSE) {
  ### Instrument
  # Z <- sample(c(-1, 1), n, TRUE) # rt(n, df = 5)
  Z <- rt(n, df = 5)
  ### Hidden
  if (doD) cop <- copula::indepCopula(2) else cop <- copula::claytonCopula(-0.5, 2)
  U <- copula::rCopula(cop, n = n)
  ### Treatment
  D <- as.numeric(plogis(Z) >= U[, 1])
  # D <- Z + rnorm(n, sd = 1 + abs(Z))
  ### Response
  UY <- runif(n)
  Y <- as.numeric(plogis(D) >= U[, 2])
  ### Return
  data.frame(Y = Y, D = D, Z = Z, UD = U[, 1], UY = U[, 2])
}

d <- dgp(1e4, TRUE)
m <- glm(Y ~ D, data = d, family = "binomial")
preds <- predict(m, newdata = d, type = "response")

# Classical residuals -----------------------------------------------------

d$cR <- d$Y - preds

### Against treatment
plot(factor(cR) ~ factor(D), data = d)
cor.test(d$cR, d$D)
coin::independence_test(cR ~ D, data = d, ytrafo = rank)
dHSIC::dhsic.test(d$cR, d$D, method = "gamma", kernel = c("discrete", "discrete"))

### Against instrument
cor.test(d$cR, d$Z)
coin::independence_test(cR ~ Z, data = d, ytrafo = rank)
dHSIC::dhsic.test(d$cR, d$Z, method = "gamma", kernel = c("discrete", "discrete"))

# rPIT residuals ----------------------------------------------------------

d$R <- randomized_pit(1 - preds, d$Y)

### Against treatment
plot(R ~ factor(D), data = d)
cor.test(d$R, d$D)
coin::independence_test(R ~ D, data = d)
coin::independence_test(R ~ D, data = d, ytrafo = rank)
dHSIC::dhsic.test(d$R, d$D, method = "gamma", kernel = c("discrete", "discrete"))

### Against instrument
plot(R ~ Z, data = d)
cor.test(d$R, d$Z)
coin::independence_test(R ~ Z, data = d, ytrafo = rank)
dHSIC::dhsic.test(d$R, d$Z, method = "gamma", kernel = c("discrete", "discrete"))
summary(multcomp::glht(glm(Y ~ D*Z, data = d, family = "binomial"),
                       linfct = c("Z == 0", "D:Z == 0")))

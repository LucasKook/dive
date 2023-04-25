#
# joint <- function(a = rnorm(1), b0 = rnorm(1), b = rnorm(1), c0 = rnorm(1),
#                   c = rnorm(1)) {
#   Z <- c("Z0" = plogis(a), "Z1" = 1 - plogis(a))
#   DZ <- matrix(c(plogis(b0), 1 - plogis(b0),
#                  plogis(b0 + b), 1 - plogis(b0 + b)),
#                nrow = 2, byrow = FALSE)
#   colnames(DZ) <- c("Z0", "Z1")
#   rownames(DZ) <- c("D0", "D1")
#   D <- (DZ %*% Z)[, 1]
#   YD <- matrix(c(plogis(c0), 1 - plogis(c0),
#                  plogis(c0 + c), 1 - plogis(c0 + c)),
#                nrow = 2, byrow = FALSE)
#   colnames(YD) <- c("D0", "D1")
#   rownames(YD) <- c("Y0", "Y1")
#   YD * rbind(D, D)
#   YD * rbind(DZ[, 1], DZ[, 1]) * Z[1] + YD * rbind(DZ[, 2], DZ[, 2]) * Z[2]
#   list("Z0" = YD * rbind(DZ[, 1], DZ[, 1]) * Z[1], "Z1" = YD * rbind(DZ[, 2], DZ[, 2]) * Z[2])
# }
#
# joint()
# sum(unlist(joint()))

devtools::load_all()
dgp <- function(n = 1e3, doD = FALSE) {
  ### Instrument
  # Z <- sample(c(-1, 1), n, TRUE) # rt(n, df = 5)
  Z <- rt(n, df = 5)
  ### Hidden
  if (doD) cop <- copula::indepCopula(2) else cop <- copula::claytonCopula(-0.5, 2)
  U <- copula::rCopula(cop, n = n)
  ### Treatment
  # D <- as.numeric(plogis(Z) >= U[, 1])
  D <- Z + rnorm(n, sd = 1 + abs(Z))
  ### Response
  UY <- runif(n)
  Y <- as.numeric(plogis(D) >= U[, 2])
  ### Return
  data.frame(Y = Y, D = D, Z = Z, UD = U[, 1], UY = U[, 2])
}
d <- dgp(1e4, TRUE)
m <- glm(Y ~ D, data = d, family = "binomial")
d$R <- d$Y - predict(m, newdata = d, type = "response")
# plot(factor(R) ~ factor(Z), data = d)
cor.test(d$R, d$Z)
coin::independence_test(R ~ Z, data = d)
coin::independence_test(R ~ Z, data = d, ytrafo = rank)
plot(R ~ Z, data = d)
dHSIC::dhsic.test(d$R, d$Z, method = "gamma", kernel = c("discrete", "discrete"))
summary(multcomp::glht(glm(Y ~ D*Z, data = d, family = "binomial"),
                       linfct = c("Z == 0", "D:Z == 0")))


dgp_ex1_binary <- function(n = 1e3, doD = FALSE) {
  ### Instrument
  Z <- sample(c(-1, 1), n, TRUE) # rt(n, df = 5)
  ### Hidden
  H <- rt(n, df = 5)
  ### Treatment
  UD <- runif(n)
  D <- as.numeric(plogis(Z + (1 - doD) * H) >= UD)
  ### Response
  UY <- runif(n)
  Y <- as.numeric(plogis(D + H) >= UY)
  ### Return
  ret <- data.frame(Y = Y, D = D, Z = Z, H = H)

  ### Compute oracle
  op1 <- integrate(\(x) plogis(1 + x) * dt(x, df = 5), -20, 20)$value
  op0 <- integrate(\(x) plogis(x) * dt(x, df = 5), -20, 20)$value

  structure(ret, p1 = op1, p0 = op0)
}

dgp_ex1_cont <- function(n = 1e3, doD = FALSE) {
  ### Instrument
  Z <- sample(c(-1, 1), n, TRUE) # rt(n, df = 5)
  ### Hidden
  H <- rt(n, df = 5)
  ### Treatment
  UD <- runif(n)
  D <- as.numeric(plogis(Z + (1 - doD) * H) >= UD)
  ### Response
  UY <- runif(n)
  Y <- qnorm(UY, mean = D + H, sd = 1 + abs(D + H))
  ### Return
  ret <- data.frame(Y = Y, D = D, Z = Z, H = H)

  ### Compute oracle
  oracle_distr <- Vectorize(\(y, d = 0) {
    integrate(\(x) pnorm(y, mean = d + x, sd = 1 + abs(d + x)) * dt(x, df = 5),
              -20, 20)$value
  }, "y")

  structure(ret, odist = oracle_distr)
}

dgp_foster <- function(n = 1e3, doD = FALSE, prs = rnorm(6)) {
  ### Instrument
  Z <- sample(c(-1, 1), n, TRUE) # rt(n, df = 5)
  ### Hidden
  H <- rt(n, df = 5)
  ### Another covariate
  X <- rt(n, df = 10)
  ### Treatment
  UD <- runif(n)
  D <- as.numeric(plogis(prs[1] * Z + prs[2] * (1 - doD) * H + prs[3] * X) >= UD)
  ### Response
  UY <- runif(n)
  Y <- as.numeric(plogis(prs[4] * D + prs[5] * H + prs[6] * X) >= UY)
  ### Return
  data.frame(Y = Y, D = D, Z = Z, X = X, H = H)

  ### Compute oracle
  # op <- glm_marginal_predictions(Y ~ D + X + H, data = ret, trt = "D")
  # structure(ret, p1 = op["p1"], p0 = op["p0"], prs = prs)
}

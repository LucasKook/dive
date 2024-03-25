
# Clip uniform away from 0, 1
.clip <- function(x) {
  stopifnot(all(x <= 1) && all(x >= 0))
  x[x == 0] <- max(x[x != 0])
  x[x == 1] <- min(x[x != 1])
  x
}

randomized_pit <- function(p0, y, trafo = identity) {
  lwr <- c(y * p0)
  upr <- c(p0^(1 - y))
  trafo(runif(NROW(y), lwr, upr))
}

# Takes function to generate random sample and turns it into a uniform RV
# using ECDF transform
# hist(.sample_to_uniform(\(n) rlogis(n) + rnorm(n), 1e3))
.sample_to_uniform <- function(sfun, n, nfine = 1e5) {
  S2U <- ecdf(sfun(nfine))
  .clip(S2U(sfun(n)))
}

cor_obj <- \(b, Y, X, E, ytrafo = NULL) {
  R <- (Y - plogis(X %*% b))
  sum(fitted(lm(R ~ E))^2)
}

ind_obj <- \(b, Y, X, E, tstat = "quadratic", trafo = identity, ytrafo = trafo) {
  p0 <- (1 - plogis(X %*% b))
  R <- randomized_pit(p0, Y)
  trafo(statistic(independence_test(R ~ E, teststat = tstat, ytrafo = ytrafo)))
}

ind_unif_obj <- \(b, Y, X, E, tstat = "quadratic", trafo = identity,
                  ytrafo = trafo, lambda = 1) {
  p0 <- (1 - plogis(X %*% b))
  R <- randomized_pit(p0, Y)
  cmv <- mean((R - ecdf(R)(R))^2)
  cor <- mean(fitted(lm(R ~ E))^2)
  max(cor, lambda * cmv)
}

OR <- \(p1, p2, cf = identity) {
  unname(cf((p1 * (1 - p2)) / ((1 - p1) * p2)))
}

ATE <- \(p1, p2, cf = identity) {
  unname(cf(p1 - p2))
}

glm_marginal_predictions <- function(fml, data, trt = "D", trafo = mean) {
  m0 <- glm(fml, data[data[[trt]] == 0, ], family = "binomial")
  m1 <- glm(fml, data[data[[trt]] == 1, ], family = "binomial")
  cbind(p1 = trafo(predict(m1, newdata = data, type = "response")),
        p0 = trafo(predict(m0, newdata = data, type = "response")))
}

indep_marginal_predictions <- function(par, data, trt = "D", trafo = mean) {
  X <- model.matrix(attr(par, "fml"), data)
  X0 <- X1 <- X
  X0[, trt] <- 0
  X1[, trt] <- 1
  cbind(p1 = trafo(plogis(X1 %*% par)), p0 = trafo(plogis(X0 %*% par)),
        cfx = unname(par[trt]))
}

ranger_marginal_predictions <- function(fml, data, trt = "D", trafo = mean) {
  rf0 <- ranger(fml, data = data[data[[trt]] == 0, ], probability = TRUE)
  rf1 <- ranger(fml, data = data[data[[trt]] == 1, ], probability = TRUE)
  rfp0 <- predict(rf0, data = data)$pred[, 2]
  rfp1 <- predict(rf1, data = data)$pred[, 2]
  cbind(p1 = trafo(rfp1), p0 = trafo(rfp0))
}

indep_iv <- function(formula, instrument, data, method = c("COR", "IND", "DIVE"), ytrafo = trafo) {
  method <- match.arg(method)
  obj <- switch(method, "COR" = cor_obj, "IND" = ind_obj, "DIVE" = ind_unif_obj)
  ### Set up model matrices
  Y <- model.response(model.frame(formula, data))
  X <- model.matrix(formula, data)
  Z <- model.matrix(instrument, data)
  ### Initial value
  b0 <- rep(1, ncol(X))
  ### Optim
  ret <- optim(b0, obj, Y = Y, X = X, E = Z, ytrafo = ytrafo)$par
  structure(ret, fml = formula, names = colnames(X))
}

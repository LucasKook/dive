
# Clip uniform away from 0, 1
.clip <- function(x) {
  stopifnot(all(x <= 1) && all(x >= 0))
  x[x == 0] <- max(x[x != 0])
  x[x == 1] <- min(x[x != 1])
  x
}

# Takes function to generate random sample and turns it into a uniform RV
# using ECDF transform
# hist(.sample_to_uniform(\(n) rlogis(n) + rnorm(n), 1e3))
.sample_to_uniform <- function(sfun, n, nfine = 1e5) {
  S2U <- ecdf(sfun(nfine))
  .clip(S2U(sfun(n)))
}

cor_obj <- \(b, Y, X, E) {
  R <- (Y - plogis(X %*% b))
  sum(fitted(lm(R ~ E))^2)
}

ind_obj <- \(b, Y, X, E, tstat = "quadratic", trafo = identity) {
  R <- (Y - plogis(X %*% b))
  trafo(statistic(independence_test(R ~ E, teststat = tstat)))
}

OR <- \(p1, p2, cf = identity) {
  unname(cf((p1 * (1 - p2)) / ((1 - p1) * p2)))
}

ATE <- \(p1, p2, cf = identity) {
  unname(cf(p1 - p2))
}

glm_marginal_predictions <- function(fml, data, trt = "D", trafo = mean) {
  m <- glm(fml, data, family = "binomial")
  nd0 <- nd1 <- data
  nd0[trt] <- 0
  nd1[trt] <- 1
  cbind(p1 = trafo(predict(m, newdata = nd1, type = "response")),
        p0 = trafo(predict(m, newdata = nd0, type = "response")),
        cfx = unname(coef(m)[trt]))
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
  rf <- ranger(fml, data = data, probability = TRUE)
  nd0 <- nd1 <- data
  nd0[trt] <- 0
  nd1[trt] <- 1
  rfp0 <- predict(rf, data = nd0)$pred[, 1]
  rfp1 <- predict(rf, data = nd1)$pred[, 1]
  cbind(p1 = trafo(rfp1), p0 = trafo(rfp0))
}

indep_iv <- function(formula, instrument, data, method = c("COR", "IND")) {
  method <- match.arg(method)
  obj <- switch(method, "COR" = cor_obj, "IND" = ind_obj)
  ### Set up model matrices
  Y <- model.response(model.frame(formula, data))
  X <- model.matrix(formula, data)
  Z <- model.matrix(instrument, data)
  ### Initial value
  b0 <- rep(0, ncol(X))
  ### Optim
  ret <- optim(b0, obj, Y = Y, X = X, E = Z)$par
  structure(ret, fml = formula, names = colnames(X))
}

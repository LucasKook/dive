# Distributional causal effects
# LK 2023

# Distributional causal effect:
#   f(y, x) = \partial_x F_{Y | X = x}(y)

# Normal linear regression ------------------------------------------------

nlrm <- function(y, x, beta = 1, sigma = 1) {
  - dnorm(y - x * beta, sd = sigma) * beta / sigma
}

curve(nlrm(x, 1, 0.5), from = -4, to = 4)

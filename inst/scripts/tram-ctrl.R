# TRAM control function
# LK 2023

library("tram")
library("tidyverse")

dgp <- function(n = 1e3) {
  Z <- runif(n, -2, 2)
  H <- runif(n, -1, 1)
  X <- as.numeric(plogis(Z + H) >= runif(n))
  Y <- qchisq(plogis(X + H + rlogis(n)), df = 10)
  data.frame(Y = Y, X = X, Z = Z, H = H)
}

nsim <- 3e2
pb <- txtProgressBar(min = 0, max = nsim, style = 3)
out <- lapply(1:nsim, \(iter) {
  setTxtProgressBar(pb, iter)
  dd <- dgp(3e2)
  m <- Colr(Y ~ X + H, data = dd, prob = c(0.001, 0.999))
  ctrl <- glm(X ~ Z, data = dd)
  dd$ctrl <- residuals(ctrl, type = "response")
  mc <- Colr(Y ~ X + ctrl, data = dd, prob = c(0.001, 0.999))
  dd$xhat <- predict(ctrl, type = "response")
  m2 <- Colr(Y ~ xhat, data = dd, prob = c(0.001, 0.999))
  c(oracle = unname(coef(m)["X"]), control = unname(coef(mc)["X"]),
    tsls = unname(coef(m2)["xhat"]))
}) |> bind_rows()

out |>
  pivot_longer(everything(), names_to = "method", values_to = "estimate") |>
ggplot(aes(x = method, y = estimate)) +
  geom_boxplot() +
  stat_summary() +
  geom_hline(yintercept = -1) +
  theme_bw()

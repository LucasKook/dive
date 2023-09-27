# Distributional random forest with control function
# LK March 2023

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("tidyverse")
library("dare")

# Data --------------------------------------------------------------------

n <- 5e2

# Data under intervention on D (d0) and observational (d)
dgp <- function(n = 1e3, doD = FALSE, cf = rnorm(5)) {
  ### Instrument
  Z <- rt(n, df = 5)
  # Z <- sample(0:1, n, TRUE)
  ### Hidden
  H <- rt(n, df = 5)
  ### Treatment
  UD <- runif(n)
  D <- as.numeric(plogis(cf[1] + cf[2] * Z + cf[3] * (1 - doD) * H) >= UD)
  ### Covariate
  X <- rnorm(n)
  ### Response
  NY <- rlogis(n)
  Y <- qchisq(plogis(0.5 * X + cf[4] * D + cf[5] * H + NY), df = 10)
  ### Return
  ret <- data.frame(Y = Y, D = D, X = X, Z = Z, H = H)
  structure(ret, cf = cf)
}

### Generate large interventional data set
tcf <- c(-1.43, -0.79, -1.19, -1.58, 0.81)
d0 <- dgp(100 * n, doD = TRUE, cf = tcf)

# Simulation --------------------------------------------------------------

nsim <- 10
pb <- txtProgressBar(0, nsim, style = 3)
res <- lapply(1:nsim, \(iter) {
  setTxtProgressBar(pb, iter)

  ### Generate data
  d1 <- dgp(n, doD = FALSE, cf = tcf)

  ### Fit RF for control function
  m <- ColrDA(Y ~ D + X, data = d1, anchor = ~ Z, xi = 1e3,
              optimizer = optimizer_adam(0.1), loss = "anchor",
              trafo_options = trafo_control(response_type = "continuous",
                                            support = \(y) c(min(d0$Y), max(d0$Y))))
  fit(m, epochs = 1e4, verbose = TRUE)

  ### Compute CDFs
  ys <- seq(min(d0$Y), max(d0$Y), length.out = 3e2)
  nd0 <- nd1 <- d1
  nd0 <- nd0 |> mutate(Y = list(ys), D = 0) |> unnest(Y)
  nd1 <- nd1 |> mutate(Y = list(ys), D = 1) |> unnest(Y)
  nd0$cdf <- predict(m, newdata = nd0, type = "cdf")
  nd1$cdf <- predict(m, newdata = nd1, type = "cdf")

  data.frame(q = ys, p0 = nd0 |> group_by(Y) |> summarize(cdf = mean(cdf)) |> pull(cdf),
             p1 = nd1 |> group_by(Y) |> summarize(cdf = mean(cdf)) |> pull(cdf))
})

pdat <- res %>%
  bind_rows(.id = "iter") |>
  pivot_longer(p0:p1, names_to = "group", values_to = "cdf")

mdat <- pdat %>% group_by(q, group) %>% summarise(cdf = mean(cdf)) %>% ungroup()

ggplot(pdat, aes(x = q, y = cdf, color = group, group = interaction(iter, group))) +
  geom_line(alpha = 0.1) +
  geom_line(aes(group = group), data = mdat, lwd = 0.9) +
  stat_ecdf(inherit.aes = FALSE, aes(x = Y, color = "p0"), data = d0[d0$D == 0, ],
            lty = 2) +
  stat_ecdf(inherit.aes = FALSE, aes(x = Y, color = "p1"), data = d0[d0$D == 1, ],
            lty = 2) +
  theme_bw() +
  scale_color_manual(values = c("p0" = "darkblue", p1 = "darkred"),
                     labels = c("p0" = "D = 0", "p1" = "D = 1")) +
  labs(color = element_blank())

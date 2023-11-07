# Test independence/uniform loss with NN
# LK 2023

set.seed(42)

# DEPs --------------------------------------------------------------------

library("tram")
library("tidyverse")
library("dare")

# FUNs --------------------------------------------------------------------

dgp <- function(n = 1e3, doD = FALSE) {
  Z <- rlogis(n)
  H <- rlogis(n)
  D <- as.numeric((3 * Z + 2 * (1 - doD) * H) / 2 > rlogis(n))
  Y <- 4 * D + 4 * H + rnorm(n, sd = 1 + abs(H)/10)
  # Y <- log(1 + exp(10 + 8 * D + H * ( 4 + rlogis(n) / 10)))
  data.frame(Y = Y, D = factor(D), Z = Z, H = H)
}

.to_gamma <- function(thetas) {
  gammas <- c(thetas[1L], log(exp(diff(thetas)) - 1))
  if (any(!is.finite(gammas))) {
    gammas[!is.finite(gammas)] <- 1e-20
  }
  gammas
}

# RUN ---------------------------------------------------------------------

ws <- dgp(1e4)
ys <- seq(min(ws$Y), max(ws$Y), length.out = 3e2)
nd0 <- data.frame(Y = ys, D = factor(0, levels = c(0, 1)))
nd1 <- data.frame(Y = ys, D = factor(1, levels = c(0, 1)))

n <- 6e2
nsim <- 2
out <- list()
for(iter in 1:nsim) {

  set.seed(111 * iter)

  ### Generate data
  d <- dgp(n)

  topt <- trafo_control(order_bsp = 30L, support = range(ws$Y),
                        response_type = "continuous")
  # m0 <- ColrNN(Y | D ~ 1, data = d, optimizer = optimizer_adam(0.1),
  #              trafo_options = topt, tf_seed = iter)
  # fit(m0, validation_split = 0, batch_size = 1e4, epochs = 1e4, callbacks = list(
  #   callback_reduce_lr_on_plateau("loss", patience = 20, factor = 0.5),
  #   callback_early_stopping("loss", patience = 50)))

  ### Setup model with independence loss and warm start
  m <- ColrDA(Y | D ~ 1, anchor = ~ Z, data = d, loss = "indep",
              optimizer = optimizer_adam(0.1), trafo_options = topt,
              tf_seed = 11 * iter)
  # set_weights(m$model, get_weights(m0$model))

  ### Fit
  fit(m, epochs = 1e4)

  ### Diagnostics
  preds <- predict(m, type = "cdf")
  plot(ecdf(preds))
  abline(0, 1)
  legend("topleft", legend = round(dHSIC::dhsic.test(preds, d$Z, method = "gamma")$p.value, 3),
         title = "Independence p-value")

  ### Return
  out[[iter]] <- data.frame(
    iter = iter, ys = ys,
    p0 = predict(m, newdata = nd0, type = "cdf"),
    p1 = predict(m, newdata = nd1, type = "cdf"))
}

dint <- dgp(n = 2e4, TRUE)
dobs <- dgp(n = 2e4)
pdat <- bind_rows(out) |>
  pivot_longer(p0:p1, names_to = "trt", values_to = "cdf") |>
  mutate(trt = factor(trt, levels = c("p0", "p1"), labels = c(0, 1)))
mdat <- pdat |> group_by(ys, trt) |> summarize(cdf = mean(cdf))

ggplot(pdat, aes(x = ys, y = cdf, color = trt, group = interaction(trt, iter))) +
  geom_line(alpha = 0.3, aes(linetype = "Est."), show.legend = FALSE) +
  geom_line(aes(linetype = "Est.", group = trt), data = mdat, lwd = 0.8) +
  theme_bw() +
  stat_ecdf(aes(x = Y, color = D, linetype = "Interv."), data = dint,
            inherit.aes = FALSE, lwd = 0.8, geom = "line") +
  stat_ecdf(aes(x = Y, color = D, linetype = "Observ."), data = dobs,
            inherit.aes = FALSE, lwd = 0.8, geom = "line") +
  scale_color_manual(values = c("0" = "darkblue", "1" = "darkred")) +
  labs(x = "Y", y = "CDF", color = "Treatment", linetype = "Type") +
  theme(legend.position = "top", text = element_text(size = 12.5)) +
  coord_cartesian(xlim = c(-30, 30))

ggsave("inst/ignore/indep-unif-loss.pdf", height = 5, width = 7)

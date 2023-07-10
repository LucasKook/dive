### Small example binary uncorrelated but dependent
### LK 2023

probs <- c(0.42, 0.28,
           0.18, 0.12)

res <- replicate(1e3, {
  smpl <- t(rmultinom(1e3, 1, prob = probs))

  dat <- do.call("rbind", apply(smpl, 1, \(obs) {
    data.frame(
      X = 0 * (obs[1] + obs[2]) + 1 * (obs[3] + obs[4]),
      Y = 0 * (obs[1] + obs[3]) + 1 * (obs[2] + obs[4])
    )
  }))

  cor(dat$X, dat$Y)
})

hist(res)
abline(v = 0)

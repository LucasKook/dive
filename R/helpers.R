
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

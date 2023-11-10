set.seed(1)

# True ECDFs
n <- 100000
H <- rnorm(n)
Y0 <- - H
Y1 <- H
F0 <- ecdf(Y0)
F1 <- ecdf(Y1)

# Data
n <- 10000

# Exogenous variable
H <- rnorm(n)

# Endogenous variables
D0 <- rbinom(n, 1, 0.1)
D1 <- rbinom(n, 1, 0.6)
D <- D0
D[H > 0] <- D1[H > 0]

# Response
# Y <- H * D - H * (1-D) 
Y <- 2 * H * D - H 

# tmp=F_D(Y), where F_d is the CDF of Y under do(D=d)
tmp <- F0(Y)
tmp[D == 1] <- F1(Y[D == 1])
hist(tmp)

plot(ecdf(D * F1(Y) + (1 - D) * F0(Y)))
abline(0, 1)


## code to prepare `generate_data_BMW` dataset goes here

## Based on DGP1 with correlated endogenous and exogenous regressors (delta = 1):
## y = beta 0 + beta1 * x + gamma * P + u
## P = delta * X + e , where e ~ gamma(1,1), x ~ gamma(1,1)
## u = rho * eta + epsilon
## where eta = Phi^{-1} (F_e (e)), epsilon ~ N(0,1)

## True parameters are:
## beta0 = 1
## beta1 = -1
## gamma = 1
## delta =1 (correlation between x and P), rho = 0.5 (moderate endogeneity)

library(MASS)

set.seed(123)

n <- 1000
beta0 <- 1
beta1 <- -1
gamma <- 1
delta <- 1
rho <- 0.5


#Non-normal exogenous regressor x ~ gamma(1,1)
x <- rgamma(n, shape=1, rate =1)

# non-normal first stage error e~ gamma(1,1)
e <- rgamma(n, shape = 1, rate = 1)

#eta = phi^{-1} (F_e (e)), where F_e is Gamma(1,1) CDF
#the control function from eq. 2.2

eta <- qnorm(pgamma(e, shape = 1, rate = 1))

# endogenous regressor correlated with x P = delta*x + e
P <- delta * x + e

# U = rho * eta + epsilon
epsilon <- rnorm(n)
u <- rho * eta + epsilon

#outcome equation

y = beta0 + beta1 * x + gamma * P + u

dataCopBMW <- data.frame(y=y, x=x, P=P)
usethis::use_data(dataCopBMW, overwrite = TRUE)



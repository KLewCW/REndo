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


## Dataset for dataCopMultiEndo
## Each endogenous regressor has its own first-stage regression on X
## and its own correction term eta_j = Phit^{-1} (F_{e_j} (e_j))

set.seed(123)
n <- 1000

#True parameters
#beta0 : 1
#beta1 : -1
#gamma1 :  1
#gamma2: 1

beta0 <- 1
beta1 <- -1
gamma1 <- 1
gamma2 <- 1


#From remark 2.1: P_j = delta_j * x + e_j , j in {1,2}
# e_1, e_2 ~ gamma(1,1). They are independent of each other and of x (assumption A4)
#u = rho_1 * eta_1 + rho_2 *eta_2 + epsilon

delta1 <- 1
delta2 <- 1
rho1 <- 0.5
rho2 <- 0.5

x <- rgamma(n, shape = 1, rate =1)
e1 <- rgamma(n, shape =1, rate =1)
e2 <- rgamma(n, shape = 1, rate = 1)

P1 <- delta1*x + e1
P2 <- delta2*x + e2

#Equation 2.3 control functions using true gamma(1,1) CDF
eta1 <- qnorm(pgamma(e1, shape = 1, rate = 1))
eta2 <- qnorm(pgamma(e2, shape = 1, rate =1))

epsilon <- rnorm(n)

u <- rho1 *eta1 + rho2 *eta2 + epsilon

y <- beta0 + beta1 * x + gamma1 * P1 + gamma2 * P2 + u

dataCopBMWMultiEndo <- data.frame(y =y , x =x, P1 = P1, P2 = P2)

usethis::use_data(dataCopBMWMultiEndo, overwrite = TRUE)

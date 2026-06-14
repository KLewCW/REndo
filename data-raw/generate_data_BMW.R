## code to prepare `generate_data_BMW` dataset goes here

## Based on DGP1 with correlated endogenous and exogenous regressors (delta = 1):
## y = mu + beta * x + alpha * P + epsilon
## P = delta * X + e , where e ~ gamma(1,1), x ~ gamma(1,1)
## epsilon = rho * eta + xi
## where eta = Phi^{-1} (F_e (e)), xi ~ N(0,1)

## True parameters are:
## mu = 1 (intercept)
## beta = -1 (x)
## alpha = 1 (P)
## delta = 1 (correlation between x and P), rho = 0.5 (moderate endogeneity)

library(MASS)

set.seed(123)

n <- 1000
mu <- 1
beta <- -1
alpha <- 1
delta <- 1
rho <- 0.5


#Non-normal exogenous regressor x ~ gamma(1,1)
x <- rgamma(n, shape=1, rate =1)

# non-normal first stage error e~ gamma(1,1)
e <- rgamma(n, shape = 1, rate = 1) #must be nonnormal as per assumption A5

#eta = phi^{-1} (F_e (e)), where F_e is Gamma(1,1) CDF
#the control function from eq. 2.2
eta <- qnorm(pgamma(e, shape = 1, rate = 1))

# endogenous regressor correlated with x P = delta*x + e
P <- delta * x + e

#structural error: epsilon = rho * eta + xi
xi <- rnorm(n)
epsilon <- rho * eta + xi

#outcome equation

y = mu + beta * x + alpha * P + epsilon

dataCopBMW <- data.frame(y=y, x=x, P=P)
usethis::use_data(dataCopBMW, overwrite = TRUE)


## Dataset for dataCopMultiEndo
## Each endogenous regressor has its own first-stage regression on X
## and its own correction term eta_j = Phi^{-1} (F_{e_j} (e_j))

## DGP:
## y = mu + beta * x + alpha1 * P1 + alpha2 * P2 + epsilon
## P_k = delta_k * x + e_k,  k in {1, 2}
## e1, e2 ~ Gamma(1,1) independent of each other and of x (assumption A4 of BMW 2024)
## epsilon = rho1 * eta1 + rho2 * eta2 + xi
## eta_k = Phi^{-1}(F_{e_k}(e_k))  (BMW 2024, eq. 2.3)
## xi ~ N(0,1)

set.seed(456)
n <- 5000 #larger sample (to get a more reliable estimate)

#True parameters
#mu : 1 (intercept)
#beta : -1 (x)
#alpha1 :  1 (P1)
#alpha2: 1 (P2)

mu <- 1
beta <- -1
alpha1 <- 1
alpha2 <- 1


#From remark 2.1: P_j = delta_j * x + e_j , j in {1,2}
# e_1, e_2 ~ gamma(1,1). They are independent of each other and of x (assumption A4)
#u = rho_1 * eta_1 + rho_2 *eta_2 + epsilon

delta1 <- 1
delta2 <- 1
rho1 <- 0.5
rho2 <- 0.5

x <- rgamma(n, shape = 1, rate =1) #nonnormal exo reg
e1 <- rgamma(n, shape =1, rate =1) #nonnormal first-stage errors independent of each other and of x
e2 <- rgamma(n, shape = 1, rate = 1)

P1 <- delta1*x + e1
P2 <- delta2*x + e2

#Equation 2.3 control functions using true gamma(1,1) CDF
eta1 <- qnorm(pgamma(e1, shape = 1, rate = 1))
eta2 <- qnorm(pgamma(e2, shape = 1, rate =1))

xi <- rnorm(n)

epsilon <- rho1 *eta1 + rho2 *eta2 + xi

y <- mu + beta * x + alpha1 * P1 + alpha2 * P2 + epsilon

dataCopBMWMultiEndo <- data.frame(y =y , x =x, P1 = P1, P2 = P2)

usethis::use_data(dataCopBMWMultiEndo, overwrite = TRUE)

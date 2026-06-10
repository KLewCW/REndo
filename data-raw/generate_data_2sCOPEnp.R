## code to prepare `generate_data_2sCOPEnp` dataset goes here

## Dataset 1: dataCopula2sCOPEnpCont - based on case 3 of Hu et al. 2025, continuous P

library(MASS)
library(truncnorm)

#Hu et al. 2025 section 4.3 neither joint Gaussian copula nor mean-dependence model
# This is a situation where (P, X) does not follow a joint Gaussian copula
#and where higher moments of P|X depend on X
#this is 'unique' case compared to other 2sCOPE , where bias is being eliminated better

set.seed(123)
n <- 1000

#true values: mu = 1, alpha =1 , beta = 2
mu <- 1
alpha <- 1
beta <- 2

# X ~ student-t(3) from eq. 34

X <- rt(n, df = 3)

# bivariate normal for (P*, Xi*), rho = 0.5 from eq. 37

Sigma <- matrix(c(1 , 0.5,
                  0.5, 1),
                nrow = 2, ncol = 2)

latent <- MASS:: mvrnorm(n = n , mu = c(0,0), Sigma = Sigma)

Pstar <- latent[,1]
Xistar  <- latent[,2]

# P|X ~ truncated normal with X-dependent bounds from eq. 35
# a = min(0, -2X + 2) , b = max(2, -2X + 2)

a <- pmin(0, -2 *X + 2)
b <- pmax(2, -2 * X + 2)

#P*_i = phi^{-1} (H_{TN(a,b)}(P|X)) from eq. 36

P <- truncnorm::qtruncnorm(pnorm(Pstar), a = a, b =b , mean = 0, sd = 1)

#Xi = Xi* , error is normal eq. 38
Xi <- Xistar

#outcome equation from eq. 39 : Y_i = mu + alpha * P_i + beta * X_i + xi_i
# Y_i = 1 + 1*P_i +2 * X_i + xi_i

y <- mu + alpha * P + beta * X + Xi

dataCopula2sCOPEnpCont <- data.frame(y = y, P = P, X = X)
usethis::use_data(dataCopula2sCOPEnpCont, overwrite = TRUE)

# Dataset 2: Binary endogenous regressor (from case 5)
# P ~ Bernoulli (binary treatment), X ~ t(3)
# rho_{px} = 0.5, rho_{pxi} = 0.5

#this dataset shows a discrete P. P being binary (0 or 1)

set.seed(123)
n <- 2000 #from paper

#true values: mu = 0, alpha = 1, beta = 2
mu <- 0
alpha <- 1
beta <- 2

#(P*, X* , Xi*) normally distributed eq. 44

Sigma <- matrix(c(1,   0.5, 0.5,
                  0.5, 1,   0,
                  0.5, 0,   1),
                nrow = 3, ncol = 3)

latent <- MASS::mvrnorm(n = n, mu = c(0, 0, 0), Sigma = Sigma)

Pstar <- latent[, 1]
Xstar <- latent[, 2]
Xistar <- latent[,3]

#Binary P: P = I{Phi(P*) > p1} = I{P* > Phi^{-1}(p1)} (eq. 46)
# p1 = 0.5
p1 <- 0.5
P <- as.integer(pnorm(Pstar) > p1) # 0 or 1

# X ~ t(3) : L(X) is student-t CDF eq.47
# X_i = L^{-1} (U_{X,i}) = L^{-1} (psi (X*_i))

X <- qt(pnorm(Xstar), df = 3)

# Xi = Xi* eq 4.5 error is normal
Xi <- Xistar

#outcome equation (Eq. 48): Y_i =  mu + alpha * P_i + beta * X_i + Xi
# Y_i = 0 + 1 * P_i + 2 * X_i + Xi_i

y <-  mu + alpha * P + beta * X + Xi

dataCopula2sCOPEnpBi <- data.frame(y = y, P = P, X = X)
usethis::use_data(dataCopula2sCOPEnpBi, overwrite = TRUE)

# Dataset 3: dataCopula2sCOPEnpMulti: extension with 2 endogenous and 2 exogenous regressors
# Based on case 3 and is extended to show that method works with multiple endogenous
#regressors.

set.seed(123)

n <- 1000

#True values: mu = 1, alpha1 = 1 (P1), alpha2 = 1 (P2), beta1 = 2 (X1), beta2 = -1 (X2)

mu <- 1
alpha1 <- 1
alpha2 <- 1
beta1 <- 2
beta2 <- -1

# Latent Gaussian dependence structure:
# rho(P1*, P2*) = 0.4 representing the correlation between the endogenous regressors
# rho(P1*, X1*) = 0.5 P1 correlated with X1
# rho(P2*, X2*) = 0.5 P2 correlated with X2
# rho(P1*, xi*) = 0.5 enogeneity of P1
# rho(P2*, xi*) = 0.5 endogeneity of P2

Sigma <- matrix(c(1, 0.4, 0.5, 0, 0.5,
                  0.4, 1, 0, 0.5, 0.5,
                  0.5, 0, 1, 0, 0,
                  0, 0.5, 0, 1, 0,
                  0.5, 0.5, 0, 0, 1),
                nrow = 5, ncol= 5, byrow = TRUE)

latent <- MASS::mvnorm(n = n, mu = rep(0,5), Sigma = Sigma)
P1star <- latent[,1]
P2star <- latent[,2]
X1star <- latent[,3]
X2star <- latent[,4]
Xistar <- latent[,5]

#X1 ~ t(3) nonnormal exo regressor (same as in case 3)
X1 <- qt(pnorm(X1star), df = 3)

#X2 ~ N(0,1) normal exo regressor
X2 <- X2star

#P1|X1 and X2 truncated normal with X1 dependent bounds
a1 <- pmin(0, -2*X1 + 2)
b1 <- pmax(2, -2 * X1 +2)
P1 <- truncnorm::qtruncnorm(pnorm(P1star), a = a1, b = b1, mean = 0, sd = 1)

#P2|X1 and X2 truncated normal with X2 depedent boungs

a2 <- pmin(0, -2 * X2 + 2)
b2 <- pmax(2, -2 * X2 + 2)
P2 <- truncnorm::qtruncnorm(pnorm(P2star), a = a2, b = b2, mean = 0, sd = 1)

#Xi = Xi star and is normally distributed
Xi <- Xistar

#outcome equation: Y_i = mu + alpha1*P1_i + alpha2*P2_i + beta1*X1_i + beta2*X2_i + Xi_i
y <- mu + alpha1 * P1 + alpha2 * P2 + beta1 * X1 + beta2 * X2 + Xi

dataCopula2sCOPEnpMulti <- data.frame(y = y, P1 = P1, P2 = P2, X1 = X1, X2 = X2)

usethis::use_data(dataCopula2sCOPEnpMulti, overwrite = TRUE)




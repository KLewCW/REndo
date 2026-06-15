## code to prepare `generate_data_JAMS` dataset goes here
library(MASS)
##Generating dataCOPJAMSSingle
#this dataset is to show how adj.ecdf is more advantageous
#compared to standard exdf in models with an intercept
# From Liengaard et al. 2024, eq. 10-11 : simulation of study 1 DGP
# y = mu + alpha * P + epsilon, where (epsilon, P*) ~ N(0, Sigma)
# sigma = [[1, 0.5], [0.5,1]]
#P = Phi(P*) ~ U(0,1)

set.seed(123)

#true parameters: mu = 3 (intercept), alpha = -1 (P) , rho = 0.5 (moderate endogeneity)
n <- 1000
mu <- 3
alpha <- -1


#joint distribution of (epsilon, P*) with endogeneity rho = 0.5
Sigma <- matrix(c(1, 0.5,
                0.5, 1),
                nrow = 2, ncol =2)

latent  <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
epsilon <- latent[, 1]
P.star  <- latent[, 2]

# P = Phi(P*) ~ U(0,1): non-normal endogenous regressor
P <- pnorm(P.star)

#outcome equation (eq. 10): y = mu + alpha * P + epsilon
y <- mu + alpha * P + epsilon

dataCopJAMSSingle <- data.frame(y = y, P = P)
usethis::use_data(dataCopJAMSSingle, overwrite = TRUE)

## dataset 2: dataCopJAMS
#based on the simulation study 2 DGP of Liengaard et al. 2024, eq. 23 & 24

#DGP : y = mu + alpha1*P1 + alpha2*P2 + beta1*W + beta2*Z +
#     delta1*P1*P2 + delta2*P1*W + delta3*P2*Z + epsilon

#The copula structure varies by Z (from eq. 24)

# Z=0: (epsilon, P1*, P2*, W*) ~ N(0, Sigma0) with rho = 0.4
# Z=1: (epsilon, P1*, P2*, W*) ~ N(0, Sigma1) with rho = 0.6
#
# Marginal distributions:
# P1, P2 ~ Gamma(2,1) are non-normal endogenous regressors
# W ~ Beta(2,2) is slightly non-normal continuous exogenous
# Z ~ Bernoulli(0.5) is binary exogenous
#
# This dataset has been generated to show the full generality of the JAMS method:
# 1) Multiple correlated endogenous regressors
# 2) Mixed discrete and continuous exogenous regressors
# 3) Nonlinear interaction terms
# 4) Copula structure that varies across Z categories

set.seed(456)

# True parameters:
# mu = 1 (intercept)
# alpha1 = -1 (P1)
# alpha2 = 1 (P2)
# beta1 = -1 (W)
# beta2 = -2 (Z)
# delta1 = -1 (P1:P2)
# delta2 = -1 (P1:W)
# delta3 = -2 (P2:Z)
n      <- 2000
mu     <-  1
alpha1 <- -1
alpha2 <-  1
beta1  <- -1
beta2  <- -2
delta1 <- -1   # P1:P2 interaction
delta2 <- -1   # P1:W interaction
delta3 <- -2   # P2:Z interaction

# Binary exogenous Z ~ Bernoulli(0.5)
Z <- rbinom(n, size = 1, prob = 0.5)

# Copula structure Z=0, rho = 0.4 (eq. 24)
Sigma1 <- matrix(c( 1,   0.4, 0.4, 0,
                    0.4, 1,   0.4, 0.4,
                    0.4, 0.4, 1,   0.4,
                    0,   0.4, 0.4, 1),
                 nrow = 4, ncol = 4)

# Copula structure Z=1, rho = 0.6 (eq. 24)
Sigma2 <- matrix(c( 1,   0.6, 0.6, 0,
                    0.6, 1,   0.6, 0.6,
                    0.6, 0.6, 1,   0.6,
                    0,   0.6, 0.6, 1),
                 nrow = 4, ncol = 4)

# Drawing latent normal variables separately per Z group
# (epsilon*, P1*, P2*, W*) ~ N(0, Sigma_z) for each z
latent <- matrix(NA, nrow = n, ncol = 4)

index1 <- Z == 0
index2 <- Z == 1

latent[index1, ] <- MASS::mvrnorm(n = sum(index1), mu = rep(0, 4), Sigma = Sigma1)
latent[index2, ] <- MASS::mvrnorm(n = sum(index2), mu = rep(0, 4), Sigma = Sigma2)

epsilon.star <- latent[, 1]
P1.star      <- latent[, 2]
P2.star      <- latent[, 3]
W.star       <- latent[, 4]

# Transforming to target marginal distributions using probability integral
# P1, P2 ~ Gamma(2,1) are non-normal endogenous regressors
P1      <- qgamma(pnorm(P1.star),  shape = 2, rate = 1)
P2      <- qgamma(pnorm(P2.star),  shape = 2, rate = 1)
# W ~ Beta(2,2) is slightly non-normal continuous exogenous regressor
W       <- qbeta(pnorm(W.star),    shape1 = 2, shape2 = 2)

epsilon <- epsilon.star  # epsilon ~ N(0,1) (Gaussian copula requirement)

#outcome equation (eq. 23):

y <- mu + alpha1 * P1 + alpha2 * P2 + beta1 * W + beta2 * Z + delta1 * P1 * P2 +
  delta2 * P1 * W + delta3 * P2 * Z + epsilon

dataCopJAMS <- data.frame( y  = y, P1 = P1, P2 = P2, W  = W, Z  = as.integer(Z))
usethis::use_data(dataCopJAMS, overwrite = TRUE)




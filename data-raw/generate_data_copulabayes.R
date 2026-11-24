## code to prepare `generate_data_CopulaBayes` dataset goes here

library(mvtnorm)
#This dataset is inspired from the section 4.1 of Haschka 2026
#based on equations 12 to 16
#one endogenous regressor (z_i) and one exogenous variable (x_i) being correlated with each other
#and an intercept

# y_i = alpga + x_i*beta  + z_i *delta + e_i; i = 1,.., N
# with alpha = 2 (intercept), beta = 6(x) and delta = -4 (z)#the true parameters
#sigam^2 = 5


#(xi_e, xi_z, xi_x) ~ N(0,Sigma)  #eq. 13
#x_i = phi^{-1} (Phi (xi_x)) ~ N (0,1) # eq. 14
#e_i = phi^{-1} (Phi(xi_e, sigma^2 = 5)) ~ N(0,5) # eq. 15
#z_i = F^{-1}_{lognormal(0,1)}(Phi(xi_z)) #eq. 16


#rho_xz = 0.3
#rho_ze = 0.7

#correlation matrix Sigma ordered as (xi_e, xi_z, xi_x):
#rho_ze = 0.7 endogeneity strength
#rho_xz =0.3 representing endo-exo correlation
#rho_xe = 0 (exogeneity restriction)

set.seed(123)
N <- 1000
alpha <- 2
beta <- 6
delta <- -4
sigma2 <- 5
rho_ze <- 0.7
rho_xz <- 0.3

Sigma <- matrix(
  c(1, rho_ze, 0,
     rho_ze, 1, rho_xz,
     0, rho_xz, 1),
  nrow = 3, ncol = 3, byrow = TRUE
)

eps <- mvtnorm::rmvnorm(N, mean = c(0,0,0), sigma = Sigma)
xi_e <- eps[,1]
xi_z <- eps[,2]
xi_x <- eps[,3]

#eq. 14:
x <- qnorm(pnorm(xi_x)) #normally distributed exo regressor

#eq. 15:
e <- qnorm(pnorm(xi_e), sd = sqrt(sigma2)) #normally distributed streuctural error

#q.16 z ~ lognormal(0,1)
z <- qlnorm(pnorm(xi_z)) #nonnormal for identification (page 520 as per chapter 2)

y <- alpha + beta * x + delta * z + e

dataCopulaBayes <- data.frame(y =y, z=z, x=x)
usethis::use_data(dataCopulaBayes, overwrite = TRUE)

## code to prepare `dataCopulaBayesMulti` dataset
## Extension of dataCopulaBayes to multiple endogenous and exogenous regressors
## K=2 endogenous, L=2 exogenous
##
## DGP is an extension of  Haschka (2026) Section 4.1
##
## y_i = alpha + beta1*x1_i + beta2*x2_i + delta1*z1_i + delta2*z2_i + e_i

library(mvtnorm)

set.seed(456)
N <- 2000

#true param
alpha  <-  1 #intercept
beta1  <-  0.5 #x1 exo reg
beta2  <- -1 #x2 exo reg
delta1 <- -2 #z1 endogenous reg
delta2 <-  1.5 #z2 endo reg
sigma2 <-  3 #error variance

rho_z1e  <- 0.6   # endogeneity of z1
rho_z2e  <- 0.5   # endogeneity of z2
rho_z1z2 <- 0.3   # correlation between z1 and z2
rho_x1z1 <- 0.2   # x1 correlated with z1
rho_x2z2 <- 0.2   # x2 correlated with z2

# Ordering: (xi_e, xi_z1, xi_z2, xi_x1, xi_x2)
Sigma.multi <- matrix(c(1,        rho_z1e,  rho_z2e,  0,        0,
                        rho_z1e,  1,        rho_z1z2, rho_x1z1, 0,
                        rho_z2e,  rho_z1z2, 1,        0,        rho_x2z2,
                        0,        rho_x1z1, 0,        1,        0,
                        0,        0,        rho_x2z2, 0,        1),
                      nrow = 5, ncol = 5, byrow = TRUE)

# check for positive definiteness
stopifnot(all(eigen(Sigma.multi)$values > 0))

eps.multi <- mvtnorm::rmvnorm(N, mean = rep(0, 5), sigma = Sigma.multi)
xi_e  <- eps.multi[, 1]
xi_z1 <- eps.multi[, 2]
xi_z2 <- eps.multi[, 3]
xi_x1 <- eps.multi[, 4]
xi_x2 <- eps.multi[, 5]

# x1, x2 ~ N(0,1): eq. 14 (normal exogenous regressors)
x1 <- qnorm(pnorm(xi_x1))
x2 <- qnorm(pnorm(xi_x2))

# e ~ N(0, sigma2)
e <- qnorm(pnorm(xi_e), sd = sqrt(sigma2))

# z1 ~ lognormal(0,1): eq. 16 (non-normal)
z1 <- qlnorm(pnorm(xi_z1))

# z2 ~ Gamma(2,1): eq. 16 (non-normal, different distribution from z1
# to show that the method can work with different marginal families)
z2 <- qgamma(pnorm(xi_z2), shape = 2, rate = 1)

y <- alpha + beta1*x1 + beta2*x2 + delta1*z1 + delta2*z2 + e


dataCopulaBayesMulti <- data.frame(y = y, z1 = z1, z2 = z2, x1 = x1, x2 = x2)
usethis::use_data(dataCopulaBayesMulti, overwrite = TRUE)



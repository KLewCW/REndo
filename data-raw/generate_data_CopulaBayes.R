## code to prepare `generate_data_CopulaBayes` dataset goes here

library(mvtnorm)
#This dataset is inspired from the section 4.1 of Haschka 2025
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
e <- qnorm(pnorm(xi_e, sd = sqrt(sigma2))) #normally distributed streuctural error

#q.16 z ~ lognormal(0,1)
z <- qlnorm(pnorm(xi_z)) #nonnormal for identification (page 520 as per chapter 2)

y <- alpha + beta * x + delta * z + e

dataCopulaBayes <- data.frame(y =y, z=z, x=x)
usethis::use_data(dataCopulaBayes, overwrite = TRUE)

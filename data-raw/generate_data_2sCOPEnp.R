## code to prepare `generate_data_2sCOPEnp` dataset goes here

## Dataset 1: data2sCOPEnoCont - based on case 3 of Hu et al. 2025, continuous P

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

#outcome equation from eq. 39 : Y_i = mu + alpha . P_i + beta . X_i + xi_i
# Y_i = 1 + 1.P_i +2 . X_i + xi_i

y <- mu + alpha * P + beta * X + Xi

data2sCOPEnpCont <- data.frame(y = y, P = P, X = X)
usethis::use_data(data2sCOPEnpCont, overwrite = TRUE)

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

#outcome equation (Eq. 48): Y_i =  mu + alpha . P_i + beta . X_i + Xi
# Y_i = 0 + 1 . P_i + 2 . X_i + Xi_i

y <-  mu + alpha * P + beta * X + Xi

data2sCOPEnpBi <- data.frame(y = y, P = P, X = X)
usethis::use_data(data2sCOPEnpBi, overwrite = TRUE)


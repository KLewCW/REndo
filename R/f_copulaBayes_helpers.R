#' @importFrom stats qnorm pnorm qgamma rgamma dnorm dgamma lgamma
#' @importFrom MCMCpack rdirichlet

#building margin structure for one variable

#First step: building per-unique-value margin structure for 1 regressor variable
#This is called once per regressor before the MCMC loop:
#This implements the categorial variable structure assumed in section 3 of Haschka 2025
# page 522, 'assume that vs exhibit categorical distribution, i.e., v_varpi ~ Cat (m_{varpi}, lambda_{varpi})
copulaBayesMargin <- function(v){
  uv <- sort(unique(v)) #sorted unique values (UV)
  grp <- match(v, uv) # integer index per observation into UV
  m <- length(uv) #number of UV
  cnt <- tabulate(grp, nbins = m) #count of observations per value
  list(uv = uv, grp = grp, m = m, cnt = cnt)

}


#converting a Dirichlet probability mass vector (lambda) to a normal score (xi) for
#one regressor. The nonparametric CDF step from section 3 is being implemented.
# The steps are:
# 1)  corresponding margins in the copula functions are obtained by taking cumulative sums of
# lambda_{varpi} of ascending ordered varpi, i.e., u_{varpi, j} = summation_{j=1}^{m_{varpi}} lambda_{varpi, j}
# 2) elements in copula are now xi_{varpi,i} = psi^{-1} (u_{varpi,i} (lambda_{varpi}))
# with u_{varpi, i}(lambda_varpi) assigning cumulative probability masses similar to a cdf
#this uses midpoint CDF rescaled to keep values in (0,1)
#The qnorm is then applied
copulaBayesConverter <- function(lambda, mg){ #lambda is the probability masses over
  #unique values (sum to 1) (from eq. 9)
  #mg is the margin structure from copulaBayesMargin

  Fmid <- (cumsum(lambda) - lambda/2) * (mg$m/ (mg$m + 0.01)) #rescaling by m/(m + 0.01 to keep values strictly inside (0,1))
  qnorm (Fmid[mg$grp])

}

# Gibbs update for Dirichlet masses (from Appendix B algorithm page 5) Eq. 9 (Haschka 2025)
#Drawing new probability masses lambda from the full conditional. The Gibbs step
#updates the nonparametric marginal distribution of one regressor


#The full conditional derivation is from W6, W7 and W9 (from appendix)
#lambda_varpi | omega_varpi = lambda_1 ,..., lambda_{m_varpi} ~ Dir(m_varpi, omega_varpi)

#Algorithm from appendix B (from Ng et al. 2011):
# 1) drawing K + L + 1 dimensional vector of multivariate normal variates v ~ N(0, sigma)
# generating dependent normal margins according to the copula representation
# 2) Applying univariate probability integral transform for the first K + L elements
# in v s.t. psi(v) element [0, 1]^{K + L}. To notw that ordering of model components
# for Sigma is z, then x, then e.

# 3) samples from W7 are then given by plugging elements of psi(v) into quantile functions
# of the univariate Dir(m_varpi, omega_{varpi} + v_{varpi}) distribution.
# converting u to Gamma(1,1) through quantile transform

copulaBayesDrawLambda <- function(u.channel, mg){
  #convert probability integral transform uniforms to Gamma(1,1) through quantile transform
  g.obs <- qgamma(u.channel, shape = 1, rate = 1)

  G.data <- as.numeric(tapply( g.obs, factor(mg$grp, levels = seq_len(mg$m)), sum))

  G.data[is.na(G.data)] <- 0

  #Dirichlet(1,...,1) prior: add one Gamma(1,1) pseudo count per cell
  #posterior : Gamma (1 + n_j, 1) per cell and normalising gives Dir(m; 1 + n_1,..., 1 +n_m)
  #posterior concentration = Dir(1,...,1) prior + data contribution
  as.vector(MCMCpack::rdirichlet(1, 1 + G.data))

}

#Log posterior: equation 4 + priors from equation 5 to 9

#Evaluating log (p (theta | y, z, x)) up to a normalising constant. This is the acceptance ratio
#numerator/denorminator in the MH step.


copulaBayeslogpost <- function(alpha, delta, beta, sigma2, Phi, xi.z, xi.x, sa, sb.delta, sb.beta, y, z, x){

  K <- ncol(xi.z)
  L <- ncol(xi.x)
  dim <- K + L + 1L

  #structural residuals
  e <- y - alpha - z %*% delta - x %*% beta
  xi.e <- pmin(pmax(e/sqrt(sigma2), -8), 8)

  #stacking all normal scores  (Appendix B ordering convention) (z_1,..,z_K, x_1,.., x_L, e)
  xi.mat <- pmin(pmax(cbind(xi.z, xi.x, xi.e), -8), 8)

  # Gaussian copula log density (eq.3)
  # equation: log c = -0.5 * log|Phi| - 0.5 * sum_i xi_i' (Phi^{-1} -I) xi_i
  A <- solve(Phi) - diag(dim)
  log.c <- -0.5 * log(det(Phi)) - 0.5 * sum(apply(xi.mat, 1, function(xi) as.numeric(t(xi) %*% A %*% xi)))


  #normal structural error log density
  log.e <- sum(dnorm(e, mean =0, sd = sqrt(sigma2), log =TRUE))

  #Inverse gamma (0.001, 0.001) prior on sigma2 from eq.7
  # log f_IG(x; a,b) = log f_G(1/x; a,b) - 2 log(x)
  log.s <- dgamma(1/sigma2, shape = 0.001, rate = 0.001, log = TRUE) - 2 * log(sigma2)

  # 'horseshoe' type shrinkage mentioned in eq. 5 and 6 prior on regression coefficients
  # gamma | phi^2 ~ N(0, phi^2)
  log.a <- dnorm(alpha, mean =0, sd = sqrt(sa), log = TRUE)
  log.d<- sum(dnorm(delta, mean = 0, sd = sqrt(sb.delta), log = TRUE))

  log.b <-if (length(beta) >0){
    sum(dnorm(beta, mean = 0, sd = sqrt(sb.beta), log = TRUE))
  } else{
    0
  }

  log.c + log.e + log.s + log.a + log.d + log.b
}

#Reparametrising the logposterior for the random walk MH step
#theta = c(alpha, delta_1,.., delta_K, beta_1,..., beta_L, log(sigma^2))
#log jacobian log(sigma^2) = log.sig2 converting from log-scale to natural scale
copulaBayeslogpostparam <- function(theta, phi, xi.z, xi.x, sa, sb.delta, sb.beta, y,z,x, K,L){
  alpha <- theta[1L]
  delta <- theta[seq(2L, 1L + K)]
  beta < - if (L > 0L){
    theta[seq(2L + K, 1L + K +L)]
  } else{
    numeric(0L)
  }

  log.s2 <- theta[1L + K + L + 1L]
  sigma2 <- exp(log.s2)

  #Jacobian = delta sigma^2/ delta log(sigma^2)
  copulaBayeslogpost(alpha,delta, beta, sigma2, Phi, xi.z, xi.x, sa, sb.delta, sb.beta, y, z, x) + log.s2

}

#Extracting upper triangle of correlation matrix for chain storage
copulaBayesMatrixtoVector <- function(Phi){
  Phi[upper.tri(Phi)]
}

#Reconstructing symmetic correlation matric from upper triangle vector
copulaBayesVectortoMatrix <- function(vec, d){
  Phi <- diag(d)
  Phi[upper.tri(Phi)] <- vec
  Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)]
  Phi
}

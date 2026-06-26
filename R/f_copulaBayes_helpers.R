#' @importFrom stats qnorm pnorm qgamma rgamma dnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom LaplacesDemon rdirichlet
#' @importFrom invgamma dinvgamma

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

  G.prior <- rgamma(mg$m, shape = 1, rate = 1)

  G <- G.prior + G.data #posterior : Gamma (1 + n_j, 1) per cell and normalising gives Dir(m; 1 + n_1,..., 1 +n_m)

  G/sum(G)
}

# Score vectors from Appendix C (W11)
# This should return N vector of per observation scores for the IWLS proposal
#eta_i = alpha + x'_i beta + z'_i delta is the linear predictor.


#As per W11, score for eta_i [ 0 ... 0 1 /sigma]' (Sigma^{-1} - I ) [ xi_{x,i}, xi_{z,i}, e_i/sigma] + e_i/sigma^2
#where [0 ... 0 1/sigma]' takes out the last row of (Sigma^{-1} -I) scaled by 1/sigma.
#this will be represented by A[dim, .]/sigma

#inverse of correlation matrix phi invPhi
# A = invPhi - I_{dim}
# xi.z = N x K matrix of normal scores for endogenous regressors
# xi.x = N x L matrix of normal scores for exogenous regressors
# e_i = y_i - eta_i, where eta_i = alpha + x'_i beta + z'_i delta (It is a N-vector)
#sigma^2 is the current error variance
#dimension(dim) = K+L+1

copulaBayesScoreEta <- function(A, xi.z, xi.x, e, sigma2, dim){

  xi.e <- e/sqrt(sigma2) #standardised error (last element of the full xi vector)

  # let t3_i be the dim-th row of A dotted with the full xi vector
  # t3_i = A[dim, ]' %*% xi_i
  #This represents the copula part of the score (the last element from W11)
  xi.all<- cbind(xi.z, xi.x, xi.e)
  t3 <- as.vector(xi.all %*% A[dim, ])

  #full score =  1/(sigma) * t3 + e /sigma^2 (element dim from W11)

  (1/sqrt(sigma2)) * t3 + e / sigma2

}

#score and approximate Hessian for log(sigma^2) (Appendix C W12 to W14)

#The function returns score (first derivative) and f2 (second derivative) of
#log L with respect to log(sigma^2) summed over observations. This is needed to build the Laplace proposal
#for log(sigma^2)

#From W12, delta (log L_i)/ delta (log sigma^2) has element [dim, dim] = 0.5 * (cross_i * xi.e_i + invPhi33 * xi.e_i ^2) - 0.5
#where cross_i = sum_{k} A [ dim, k] * xi.z[i,k] + sum_{l} A [dim, K + l] * xi.x[i,l] is the off-diagonal contribution

copulaBayesScorelogsigma2 <- function(A, invPhi33, xi.z, xi.x, e, sigma2){
  #invPhi33 represents (Xi^{-1})_{33}

  dim <- ncol(xi.z) + ncol(xi.x) + 1

  xi.e <- e/sqrt(sigma2)

  #off diagonal cross terms in row dim of A
  xi.regs <- cbind(xi.z, xi.x) #works for L=0 because xi.x has 0 columns
  cross <- as.vector(xi.regs %*% A[dim, seq_len(ncol(xi.regs))])


  #score = ( delta (log L ) ) / (delta (log sigma^2)) summed over observations (W12)
  #score derivative
  #from eq.7 in Haschka 2025, variance of Structural Error with common hyperparameters a = b = 0.001
  #where sigma^2~ inverse Gamma prior to the error variance (IG) (a,b)
  Score <- sum(0.5 * (cross * xi.e + invPhi33 * xi.e^2) - 0.5) - 0.001 + 0.001/sigma2 #sum over obs of W12 + IG prior derivative


  #approximate hessian (element [dim, dim] from W14) + IG prior curvature

  f2 <- sum (-0.25 * cross * xi.e - 0.5 * invPhi33 * xi.e^2) - 0.001/sigma2

  list(Score = Score, f2 = f2 )
}

#Log posterior: equation 4 + priors from equation 5 to 9

#Evaluating log (p (theta | y, z, x)) up to a normalising constant. This is the acceptance ratio
#numerator/denorminator in the MH step.


copulaBayeslogpost <- function( alpha, delta, beta, sigma2, Phi, xi.z, xi.x, sa, sb.delta, sb.beta, y, z, x){

  dim <- ncol(xi.z) + ncol(xi.x) + 1

  #structural residuals
  e <- y - alpha - z %*% delta - x %*% beta
  xi.e <- pmin(pmax(e/sqrt(sigma2), -8), 8)

  #stacking all normal scores  (Appendix B ordering convention)
  xi.mat <- pmin(pmax(cbind(xi.z, xi.x, xi.e), -8), 8)

  # Gaussian copula log density (eq.3)
  # equation: = -0.5 * log|Phi| - 0.5 * sum_i xi_i' (Phi^{-1} -I) xi_i
  A <- solve(Phi) - diag(dim)
  log.c <- -0.5 * log(det (Phi)) - 0.5 * sum(apply(xi.mat, 1, function(xi) as.numeric(t(xi) %*% A %*% xi)))


  #normal structural error log density
  log.e <- sum(dnorm(e, mean =0, sd = sqrt(sigma2), log =TRUE))

  #Inverse gamma (0.001, 0.001) prior on sigma2 from eq.7
  log.s <- invgamma::dinvgamma(sigma2, shape = 0.001, rate = 0.001, log = TRUE)

  # 'horseshoe' type shrinkage mentioned in eq. 5 and 6 prior on regression coefficients
  # gamma | phi^2 ~ N(0, phi^2)
  log.a <- dnorm(alpha, mean =0, sd = sqrt(sa), log = TRUE)
  log.d<- sum(dnorm(delta, mean = 0, sd = sqrt(sb.delta), log = TRUE))
  log.b <-sum(dnorm(beta, mean = 0, sd = sqrt(sb.beta), log = TRUE))

  log.c + log.e + log.s + log.a + log.d + log.b
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

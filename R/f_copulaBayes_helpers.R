#' @importFrom stats qnorm pnorm qgamma rgamma dnorm
#' @importFrom mvtnorm rmnorm
#' @importFrom LaplacesDemon rdirichlet
#' @importFrom invgamma dinvgamma

#building margin structure for one variable
#auxiliary functions
#This is called once per regressor before the MCMC loop:

copulaBayesMargin <- function(v){
  uv <- sort(unique(v)) #sorted unique values (UV)
  grp <- match(v, uv) #per observation index into UV
  m <- length(uv) #number of UV
  cnt <- tabulate(grp, nbins = m) #count per value
  list(uv = uv, grp = grp, m = m, cnt = cnt)

}

#now the Dirichlet mass vector lambda is converted to a normal score (xi)
# this uses midpoint CDF rescaled to keep values in (0,1)
#The qnorm is then applied
copulaBayesConverter <- function(lambda, mg){ #lambda is the probability masses over
  #unique values (sum to 1)
  #mg is the margin structure from copulaBayesMargin

  Fmid <- (cumsum(lambda) - lambda/2) * (mg$m/ (mg$m + 0.01))
  qnorm (Fmid[mg$grp])

}

# Gibbs update for Dirichlet masses (from Appendix B algorithm page 5)
copulaBayesDrawLambda <- function(u.channel, mg){
  #convert probability integral transform uniforms to Gamma(1,1) through quantile transform
  g.obs <- qgamma(u.channel, shape = 1, rate = 1)

  G.data <- as.numeric(tapply( g.obs, factor(mg$grp, levels = seq_len(mg$m)), sum))

  G.data[is.na(G.data)] <- 0

  #Dirichlet(1,...,1) prior: add one Gamma(1,1) pseudo count per cell

  G.prior <- rgamma(mg$m, shape = 1, rate = 1)

  G <- G.prior + G.data #posterior : Dir(m; 1 + n_1,..., 1 +n_m)

  G/sum(G)
}

# Score vectors from Appendix C (W11)
# This should return N vector of per observation scores for the IWLS proposal

#inverse of correlation matrix phi invPhi
# A = invPhi - I_{dim}
# xi.z = N x K matrix of normal scores for endogenous regressors
# xi.x = N x L matrix of normal scores for exogenous regressors
# e_i = y_i - eta_i, where eta_i = alpha + x'_i beta + z'_i delta (here we will
#write it as e a N vector )
#sigma^2 is the current error variance
#dimension(dim) = K+L+1

copulaBayesScoreEta <- function(invPhi, A, xi.z, xi.x, e, sigma2, dim){
  N <- length(e)
  K <- ncol(xi.z)
  L <- if(is.matrix(xi.x)){
    ncol(xi.x)
  } else{
    0
  }

  xi.e <- e/sqrt(sigma2)

  # let t3_i be the dim-th row of A dotted with the full xi vector

  t3 <- rep(0,N)

  for (k in seq_len(K))
    t3 <- t3 + A[dim, k] * xi.z[,k]

  if(L >0)
    for (l in seq_len(L))
      t3 <- t3 + A[dim, K + l] * xi.x[, l]
  t3 <- t3 + A[dim, dim] * xi.e

  #score is then (1/sigma) * t3 + e/sigma^2

  (1 / sqrt(sigma2)) * t3 + e / sigma2 #W11
}

#score and approximate Hessian for log(sigma^2) (Appendix C W12 to W14)

copulaBayesScorelogsigma2 <- function(A, invPhi33, xi.z, xi.z, e, sigma2){
  #invPhi33 represents (Xi^{-1})_{33}

  N <- length(e)
  K <- ncol(xi.z)
  L <- if(is.matrix(xi.x)){
    ncol(xi.x)
  } else{
    0
  }

  dim <- K + L + 1

  xi.e <- e/sqrt(sigma2)

  #off iagonal cross terms in row dim of A
  cross <- rep(0,N)
  for (k  in seq_len(K))
    cross <- cross + A[dim, k] * xi.z[,k]
  if(L > 0 )
    for (l in seq_len(L))
      cross <- cross + A[dim, K + l] * xi.x[, l]

  #score = ( delta (log L ) ) / (delta (log sigma^2)) summed over observations (W12)
  #score derivative
  Score <- sum(0.5 * (cross * xi.e + invPhi33 * xi.e^2) - 0.5) - 0.001 + 0.01/sigma^2


  #approximate hessian (W14). Second derivative (negative near the mode; P = -1/f2 below)

  f2 <- sum (-0.25 * cross * xi.e - 0.5 * invPhi33 * xi.e^2) - 0.001/sigma2

  list(Score = Score, f2 = f2 )
}

#Log posterior: equation 4 + prior
# sa, sb.delta and sb.beta are hyperprior variances

copulaBayeslogpost <- function( alpha, delta, beta, sigma2, Phi, xi.z, xi.x, sa, sb.delta, sb.beta, y, x, z){

  K <- ncol(z)
  L <- if(is.matrix(x) && ncol(x) >0){
    ncol(x)
  } else{
    0
  }

  dim <- K + L + 1

  e <- y - alpha - z %*% delta - if ( L >0 ) x %*% beta else 0
  xi.e <- e/sqrt(sigma2)
  xi.e <- pmin(pmax(xi.e, -8), 8)

  xi.mat <- cbind(xi.z, if (L>0) xi.x else NULL, xi.e)
  xi.mat <- pmin(pmax(xi.mat, -8), 8)

  # Gaussian copula log density (eq.3)

  A <- solve(Phi) - diag(dim)
  log.c <- -0.5 * log(det (Phi)) - 0.5 * sum(apply(xi.mat, 1, function(xi) as.numeric(t(xi) %*% A %*% xi)))


  #normal structural error density
  log.e <- sum(dnorm(e, mean =0, sd = sqrt(sigma2), log =TRUE))

  #Inverse gamma (0.001, 0.001) prior on sigma2 from eq.7

  log.s <- invgamma::dinvgamma(sigma2, shape = 0.001, rate = 0.001, log = TRUE)

  # 'horseshoe' type shrinkage mentioned in eq. 5 and 6 prior on regression coefficients
  log.a <- dnorm(alpha, mean =0, sd = sqrt(sa), log = TRUE)
  log.d<- sum(dnorm(delta, mean = 0, sd = sqrt(sb.delta), log = TRUE))
  log.b <- if (L >0) sum(dnorm(beta, mean = 0, sd = sqrt(sb.beta), log = TRUE)) else 0

  log.c + log.e + log.s + log.a + log.d + log.b
}



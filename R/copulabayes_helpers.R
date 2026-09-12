#' @importFrom stats qnorm pnorm qgamma rgamma dnorm dgamma
#' @importFrom MCMCpack riwish rdirichlet
copulabayes_margin <- function(v) {
  #building margin structure for one variable

  #First step: building per-unique-value margin structure for 1 regressor variable
  #This is called once per regressor before the MCMC loop:
  #This implements the categorial variable structure assumed in section 3 of Haschka 2025
  # page 522, 'assume that vs exhibit categorical distribution, i.e., v_varpi ~ Cat (m_{varpi}, lambda_{varpi})

  unique.vals <- sort(unique(v)) #sorted unique values (unique.vals)
  value.index <- match(v, unique.vals) # integer index per observation into unique.vals
  n.unique <- length(unique.vals) #number of unique.vals
  cnt <- tabulate(value.index, nbins = n.unique) #count of observations per value
  return(list(unique.vals = unique.vals, value.index = value.index, n.unique = n.unique, cnt = cnt))
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
copulabayes_converter <- function(lambda, mg) {
  #lambda is the probability masses over
  #unique values (sum to 1) (from eq. 9)
  #mg is the margin structure from copulabayes_margin

  cdf.midpoint <- (cumsum(lambda) - lambda / 2) * (mg$n.unique / (mg$n.unique + 0.01)) #rescaling by m/(m + 0.01 to keep values strictly inside (0,1))
  return(qnorm(cdf.midpoint[mg$value.index]))
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

copulabayes_drawlambda <- function(u.channel, mg) {
  #convert probability integral transform uniforms to Gamma(1,1) through quantile transform
  gamma.draw <- qgamma(u.channel, shape = 1, rate = 1)

  #per unique value sum of Gamma draws
  gamma.sum.per.val <- as.numeric(tapply(gamma.draw, factor(mg$value.index, levels = seq_len(mg$n.unique)), sum))

  gamma.sum.per.val[is.na(gamma.sum.per.val)] <- 0

  #Dirichlet(1,...,1) prior: add one Gamma(1,1) pseudo count per cell
  #posterior : Gamma (1 + n_j, 1) per cell and normalising gives Dir(m; 1 + n_1,..., 1 +n_m)
  #posterior concentration = Dir(1,...,1) prior + data contribution
  return(as.vector(mCmCpack::rdirichlet(1, 1 + gamma.sum.per.val)))
}

#Log posterior: equation 4 + priors from equation 5 to 9

#Evaluating log (p (theta | y, z, x)) up to a normalising constant. This is the acceptance ratio
#numerator/denorminator in the mH step.

copulabayes_logpost <- function(
  intercept,
  coef.endo,
  coef.exo,
  error.var,
  copula.cor,
  scores.endo,
  scores.exo,
  var.alpha,
  var.delta,
  var.beta,
  y,
  z,
  x
) {
  K <- ncol(scores.endo)
  L <- ncol(scores.exo)
  copula.dim <- K + L + 1L

  #structural residuals
  resid <- y - intercept - z %*% coef.endo - x %*% coef.exo
  scores.error <- pmin(pmax(qnorm(pnorm(as.vector(resid) / sqrt(error.var))), -8), 8)

  #stacking all normal scores  (Appendix B ordering convention) (z_1,..,z_K, x_1,.., x_L, e)
  scores.all <- pmin(pmax(cbind(scores.endo, scores.exo, scores.error), -8), 8)

  # Gaussian copula log density (eq.3)
  # equation: log c = -0.5 * log|copula.cor| - 0.5 * sum_i xi_i' (copula.cor^{-1} -I) xi_i
  A <- solve(copula.cor) - diag(copula.dim)
  log.copula <- -0.5 *
    log(det(copula.cor)) -
    0.5 * sum(apply(scores.all, 1, function(xi) as.numeric(t(xi) %*% A %*% xi)))

  #normal structural error log density
  log.error.density <- sum(dnorm(e, mean = 0, sd = sqrt(error.var), log = TRUE))

  #Inverse gamma (0.001, 0.001) prior on error.var from eq.7
  # log f_IG(x; a,b) = log f_G(1/x; a,b) - 2 log(x)
  log.prior.sigma2 <- dgamma(1 / error.var, shape = 0.001, rate = 0.001, log = TRUE) - 2 * log(error.var)

  # 'horseshoe' type shrinkage mentioned in eq. 5 and 6 prior on regression coefficients
  # gamma | copula.cor^2 ~ N(0, copula.cor^2)
  log.prior.alpha <- dnorm(intercept, mean = 0, sd = sqrt(var.alpha), log = TRUE)
  log.prior.delta <- sum(dnorm(coef.endo, mean = 0, sd = sqrt(var.delta), log = TRUE))

  if (length(coef.exo) > 0) {
    log.prior.beta <- sum(dnorm(coef.exo, mean = 0, sd = sqrt(var.beta), log = TRUE))
  } else {
    log.prior.beta <- 0
  }

  return(log.copula + log.error.density + log.prior.sigma2 + log.prior.alpha + log.prior.delta + log.prior.beta)
}

#Reparametrising the logposterior for the random walk mH step
#theta = c(intercept, coef.endo_1,.., coef.endo_K, coef.exo_1,..., coef.exo_L, log(sigma^2))
#log jacobian log(sigma^2) = log.prior.sigma2ig2 converting from log-scale to natural scale
copulabayes_logpostparam <- function(
  theta,
  copula.cor.cur,
  scores.endo,
  scores.exo,
  var.alpha,
  var.delta,
  var.beta,
  y,
  z,
  x,
  K,
  L
) {
  intercept <- theta[1L]
  coef.endo <- theta[seq(2L, 1L + K)]

  if (L > 0L) {
    coef.exo <- theta[seq(2L + K, 1L + K + L)]
  } else {
    coef.exo <- numeric(0L)
  }

  log.prior.sigma22 <- theta[1L + K + L + 1L]
  error.var <- exp(log.prior.sigma22)

  #Jacobian = coef.endo sigma^2/ coef.endo log(sigma^2)
  return(copulabayes_logpost(
    intercept,
    coef.endo,
    coef.exo,
    error.var,
    copula.cor.cur,
    scores.endo,
    scores.exo,
    var.alpha,
    var.delta,
    var.beta,
    y,
    z,
    x
  ) +
    log.prior.sigma22)
}

#Extracting upper triangle of correlation matrix for chain storage
copulabayes_matrix2vector <- function(copula.cor) {
  return(copula.cor[upper.tri(copula.cor)])
}

#Reconstructing symmetic correlation matric from upper triangle vector
copulaBayesVectortomatrix <- function(vec, copula.dim) {
# TODO: Unused??
  copula.cor <- diag(copula.dim)
  copula.cor[upper.tri(copula.cor)] <- vec
  copula.cor[lower.tri(copula.cor)] <- t(copula.cor)[lower.tri(copula.cor)]
  return(copula.cor)
}

copulabayes_logpost_coef <- function( #used in the step 1 for the IWLS for log-posterior of coeff
    coef.vec, #c(intercept, coef.endo, coef.exo)
    error.var,  #current sigma^2
    copula.cor, #current copula correlation matrix copula.cor
    scores.endo,  #current normal scores of endogenous regressors (N x K)
    scores.exo, #current normal scores of exogenous regressors (N x L)
    var.intercept,  #horseshoe variance for intercept
    var.coef.endo,  #horseshoe variances for endo coefficients (length K)
    var.coef.exo,   #horseshoe variances for exo coefficients (length L)
    y, z, x, K, L
) {X.design  <- cbind(1, z, x)
resid <- y - X.design %*% coef.vec
scores.error <- pmin(pmax(qnorm(pnorm(as.vector(resid) / sqrt(error.var))), -8), 8)
scores.all <- pmin(pmax(cbind(scores.endo, scores.exo, scores.error), -8), 8)

copula.dim <- K + L + 1L
A <- solve(copula.cor) - diag(copula.dim)
log.copulaopula <- -0.5 * log(det(copula.cor)) - 0.5 * sum(apply(scores.all, 1, function(xi) as.numeric(t(xi) %*% A %*% xi)))
log.error.densityrror <- sum(dnorm(resid, mean = 0, sd = sqrt(error.var), log = TRUE))
log.prior <- dnorm(coef.vec[1L], mean = 0, sd = sqrt(var.intercept), log = TRUE) +
  sum(dnorm(coef.vec[seq(2L, 1L + K)], mean = 0, sd = sqrt(var.coef.endo), log = TRUE))

if (L > 0L) {
  log.prior <- log.prior + sum(dnorm(coef.vec[seq(2L + K, 1L + K + L)], mean = 0, sd = sqrt(var.coef.exo), log = TRUE))
}
return(log.copulaopula + log.error.densityrror + log.prior)
}

# Log-posterior n for sigma^2 (in step 2 IWLS)
# Used for the Laplace proposal
copulabayes_logpost_sigma2 <- function(
    error.var,  #proposed sigma^2
    copula.cor,
    scores.endo,
    scores.exo,
    resid,  #current structural residuals
    y, z, x, K, L
) {
  if (error.var <= 0) return(-Inf)
  copula.dim <- K + L + 1L
  scores.error  <- pmin(pmax(qnorm(pnorm(as.vector(resid) / sqrt(error.var))), -8), 8)
  scores.all <- pmin(pmax(cbind(scores.endo, scores.exo, scores.error), -8), 8)
  A <- solve(copula.cor) - diag(copula.dim)
  log.copulaopula <- -0.5 * log(det(copula.cor)) -
    0.5 * sum(apply(scores.all, 1, function(xi) as.numeric(t(xi) %*% A %*% xi)))
  log.error.densityrror <- sum(dnorm(resid, mean = 0, sd = sqrt(error.var), log = TRUE))
  # IG(0.001, 0.001) prior on sigma^2
  log.prior <- dgamma(1 / error.var, shape = 0.001, rate = 0.001, log = TRUE) - 2 * log(error.var)
  return(log.copulaopula + log.error.densityrror + log.prior)
}


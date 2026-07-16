#' @importFrom stats lm residuals var runif rnorm rgamma
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack riwish rdirichlet
# MCMC sampler of the copula model from Appendix D

#sampler cycles through 4 blocks per iteration (eq.11):
# 1. Metropolis-Hastings algorithm (MH) step for alpha, delta, beta using IWLS proposal from appendix C
# 2. MH step for sigma^s using Laplace proposal (Appendix C)
# 3. Gibbs step for correlation matrix W (Appendix A W5, Wishart distribution)
# 4. Gibbs step for Dirichlet masses lambda (Appendix W7)

#Note: for steps one and 2 :  the MH (random walk) replaces the IWLS and Laplace proposals.
#From MCMCpack, random walk and Hastings ratio have been used.

#modification done from previous code because was a bit lost:
# X has now been kept as a N x L matrix (when L =0, it has 0 columns, i.e., X %*% beta = 0 )


copulaBayesMCMC <- function(y, z, x, num.iterations, verbose){

  N <- length(y)
  K <- ncol(z)
  L <- ncol(x)
  cop.dim <- K + L + 1L  #copula dimension = K endogenous + L exogenous + error


  #building margin structures per regressor (from Appendix D step 2)
  #This describe the unique values (UV) of each z_k and x_l and
  #this is fixed for the entire MCMC run.
  #note: observed data never change.
  mgz.list <- lapply(seq_len(K), function(k) copulaBayesMargin(z[,k]))
  mgx.list <- lapply(seq_len(L), function(l) copulaBayesMargin(x[,l]))

  # From appendix D: ols starting values
  # Research paper says starting value can be arbitrary. OLS is used here for conveniency
  #and close to the posterior mode when endogeneity is weak
  dat.ols <- data.frame(y =y, z =z, x=x)
  mod.ols <- lm(y ~ ., data= dat.ols)
  coef.ols <- coef(mod.ols)

  alpha.cur <- coef.ols["(Intercept)"]
  delta.cur <- coef.ols[1L + seq_len(K)]
  if( L > 0L){
    beta.cur <- coef.ols[1L + K + seq_len(L)]
  } else {
    beta.cur <- numeric(0L)
  }

  sigma2.cur <- var(residuals(mod.ols))

  # setting sigma^[0] = identiy matrix (I)
  # Here we start with no correlation assumed.
  Phi.cur <- diag(cop.dim)

  #initial lambda from Dir(1,...,1) are drawn
  lambdaz.list <- lapply(mgz.list, function(mg) as.vector(MCMCpack::rdirichlet(1, rep(1, mg$m))))
  lambdax.list <- lapply(mgx.list, function(mg) as.vector(MCMCpack::rdirichlet(1, rep(1, mg$m))))

  #initial normal scores from the starting lambda
  # xi_{varpi, i} = Phi^{-1} (u_{varpi, i} (lambda_varpi)) from section 3.1 of Haschka 2025, page 522
  xi.z <- matrix(NA_real_, N,K)
  for (k in seq_len(K))
    xi.z[,k] <- copulaBayesConverter(lambdaz.list[[k]], mgz.list[[k]]) #converting lambda to xi

  xi.x <- matrix(NA_real_, N, L)
  for (l in seq_len(L))
      xi.x[, l] <- copulaBayesConverter(lambdax.list[[l]], mgx.list[[l]])

  # Chain storage: one row per iteration.
  # Column layout (fixed throughout)
  n.rho <- cop.dim * (cop.dim - 1L)/2L
  n.hyper <- 1L + K + L
  n.coef <- 1L + K + L  #alpha + delta+ beta

  col.alpha <- 1L #intercept alpha
  col.delta <- 1L + seq_len(K) #delta_1,...,delta_K (endo coefficients)

  if (L > 0){
    col.beta <- 1L + K + seq_len(L)
  } else{
    col.beta <- integer(0L)  # exo coefficients beta_1,...,beta_L. Empty if L=0
  }

  col.sigma2 <- 1L + K + L + 1L #error variance sigma^2

  col.rho <- col.sigma2 + seq_len(n.rho) #upper triangle of Phi

  col.sa <- tail(col.rho, 1L) + 1L #'horseshoe' hyperprior variances for alpha

  col.sb.d <- col.sa + seq_len(K) #'horseshoe' hyperprior variances for delta_1,..., delta_K

  col.sb.b <- if (L >0L ){ #'horseshoe' hyperprior variances for  beta_1,..,beta_L
    col.sa + K + seq_len(L)
  } else{
    integer(0L)
  }

  n.cols <- tail(col.rho, 1L) + n.hyper

  chain <- matrix(NA_real_, nrow = num.iterations + 1L , ncol = n.cols)

  #Input first row with the starting values
  chain[1L, col.alpha] <- alpha.cur
  chain[1L, col.delta] <- delta.cur
  if ( L > 0L ) chain[1L, col.beta] <- beta.cur
  chain [1L, col.sigma2] <- sigma2.cur
  chain[1L, col.rho] <-copulaBayesMatrixtoVector(Phi.cur)

  chain[1L, col.sa] <- 1000
  chain[1L, col.sb.d] <- rep (1000, K)
  if (L > 0L ) chain[1L, col.sb.b] <- rep(1000, L)

  tune <- 1.0
  n.accepted <- 0L


  #MCMC loop

  for (i in seq_len(num.iterations)){

    if (verbose && i %% 500L == 0L){
      message( "Iteration: ", i, " | alpha = ", round(chain[i, col.alpha], 3L), " | delta[1] =",
                round(chain[i, col.delta[1L]], 3L), " | sigma2 = ", round(chain[i, col.sigma2], 3L),
               " | accept rate = ", round(n.accepted/i, 2L))
    }

    #Retrieving current hyperprior variances from chain
    sa.cur <- chain[i, col.sa]
    sb.d.cur <- chain[i, col.sb.d]

    if (L > 0L) {
      sb.b.cur <- chain[i, col.sb.b]
    } else{
        sb.b.cur <- numeric(0L)
      }

    #Step 1 and 2 of the iteration: random walk mH for alpha, delta, beta, log sigma^2
    #Using the MCMCpack metropolis principle: propose, then evalute ratio and accept.
    #Parameters are sampled jointly in the reparametrised space (log sigma^2) so that
    #the proposal is unconstrained. Replaces the custom IWLS and Laplace proposals from
    #the repository without any score computation.

    theta.cur <- c(alpha.cur, delta.cur, beta.cur, log(sigma2.cur))
    theta.prop <- theta.cur + rnorm(length(theta.cur), mean = 0, sd = tune)

    logpost.cur <- copulaBayeslogpostparam(theta.cur, Phi.cur, xi.z, xi.x, sa.cur, sb.d.cur, sb.b.cur, y, z, x, K, L)
    logpost.prop <- copulaBayeslogpostparam(theta.prop, Phi.cur,xi.z, xi.x, sa.cur, sb.d.cur, sb.b.cur, y, z, x, K, L )

    log.mh <- logpost.prop - logpost.cur
    if(!is.finite(log.mh)) log.mh <- -Inf

    if(log(runif(1L)) < log.mh){
      theta.cur <- theta.prop
      n.accepted <- n.accepted + 1L
    }

    alpha.cur  <- theta.cur[1L]
    delta.cur  <- theta.cur[seq(2L, 1L + K)]

    if(L > 0L){
      beta.cur <- theta.cur[seq(2L + K, 1L + K + L)]
    } else{
      beta.cur <- numeric(0L)
    }

    sigma2.cur <- exp(theta.cur[1L + K + L + 1L])

    chain[i + 1L, col.alpha]  <- alpha.cur
    chain[i + 1L, col.delta]  <- delta.cur
    if (L > 0L) chain[i + 1L, col.beta] <- beta.cur
    chain[i + 1L, col.sigma2] <- sigma2.cur

    #tuning every 50 iterations (Robbins-Monro)
    if (i %% 50L == 0L) {
      rate <- n.accepted / i
      if (rate < 0.20) tune <- tune * 0.9
      if (rate > 0.40) tune <- tune * 1.1
    }

    #Step 3: the Gibbs for copula correction matrix W (appendix A W5)
    e.cur <- y - alpha.cur - z %*% delta.cur - x %*% beta.cur
    xi.e <- pmin(pmax(qnorm(pnorm(as.vector(e.cur) /sqrt(sigma2.cur))), -8), 8)
    xi.all <- cbind(xi.z, xi.x, xi.e)

    W.new <- MCMCpack::riwish(N + cop.dim, diag(cop.dim) + crossprod(xi.all))
    sds <- sqrt(diag(W.new))
    D.inv <- diag(1/sds)
    Phi.cur <- D.inv %*% W.new %*% D.inv

    # Enforce exogeneity
    for (l in seq_len(L)) {
      Phi.cur[K + l, cop.dim] <- 0  #zero out (x_l, e) correlations
      Phi.cur[cop.dim, K + l] <- 0
    }

    chain[i + 1L , col.rho] <- copulaBayesMatrixtoVector(Phi.cur)

    #Step 4: Gibbs for Dirichlet masses lambda (Appendix B W7)
    #Drawing correlated normals from current phi, then copulaBayesDrawLambda updates
    #each variable's mass vector

    eps <- mvtnorm::rmvnorm(N,mean = rep(0,cop.dim), sigma = Phi.cur)

    for (k in seq_len(K)) {
      lambdaz.list[[k]] <- copulaBayesDrawLambda(pnorm(eps[, k]),mgz.list[[k]])
      xi.z[, k] <- copulaBayesConverter(lambdaz.list[[k]],mgz.list[[k]])
    }
    for (l in seq_len(L)) {
      lambdax.list[[l]] <- copulaBayesDrawLambda(pnorm(eps[, K + l]), mgx.list[[l]])
      xi.x[, l] <- copulaBayesConverter(lambdax.list[[l]], mgx.list[[l]])
    }

    #Horseshoe hyperprior updates (From equation 5 of Haschka 2025)

    chain[i + 1L, col.sa] <- 1 / rgamma(1L, shape = 0.501, rate = (alpha.cur^2 + 0.002) / 2)

    for (k in seq_len(K)) {
      chain[i + 1L, col.sb.d[k]] <- 1 / rgamma(1L, shape = 0.501, rate = (delta.cur[k]^2 + 0.002) / 2)
    }
    if (L > 0L) {
      for (l in seq_len(L)) {
        chain[i + 1L, col.sb.b[l]] <- 1 / rgamma(1L, shape = 0.501, rate = (beta.cur[l]^2 + 0.002) / 2)
      }
    }

  }

  attr(chain, "col.alpha")  <- col.alpha
  attr(chain, "col.delta")  <- col.delta
  attr(chain, "col.beta")   <- col.beta
  attr(chain, "col.sigma2") <- col.sigma2
  attr(chain, "col.rho")    <- col.rho
  attr(chain, "K")          <- K
  attr(chain, "L")          <- L
  attr(chain, "cop.dim")        <- cop.dim

  chain
}

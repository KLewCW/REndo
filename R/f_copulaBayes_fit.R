#' @importFrom stats lm residuals var runif rnorm
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom LaplacesDemon rinvwishart rdirichlet
#' @importFrom invgamma rinvgamma
#'

# MCMC sampler of the copula model from Appendix D

#sampler cycles through 4 blocks per iteration (eq.11):
# 1. Metropolis-Hastings algorithm (MH) step for alpha, delta, beta using IWLS proposal from appendix C
# 2. MH step for sigma^s using Laplace proposal (Appendix C)
# 3. Gibbs step for correlation matrix W (Appendix A W5, Wishart distribution)
# 4. Gibbs step for Dirichlet masses lambda (Appendix W7)

#modification done from previous code because was a bit lost:
# X has now been kept as a N x L matrix (when L =0, it has 0 columns, i.e., X %*% beta = 0 )


copulaBayesMCMC <- function(y, z, x, num.iterations, verbose){

  N <- length(y)
  K <- ncol(z)
  L <- ncol(x) # 0 when no exogenous regressors are present
  dim <- K + L + 1  #dimension of the copula correlation matrix


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
  delta.cur <- coef.ols[1 + seq_len(K)]
  beta.cur <- if( L >0){
    coef.ols[1 + K + seq_len(L)]
  } else {
    numeric(0)
  }

  sigma2.cur <- var(residuals(mod.ols))

  # setting sigma^[0] = identiy matrix (I)
  # Here we start with no correlation assumed.
  Phi.cur <- diag(dim)

  #initial lambda from Dir(1,...,1) are drawn
  lambdaz.list <- lapply(mgz.list, function(mg) as.vector(LaplacesDemon::rdirichlet(1, rep(1, mg$m))))
  lambdax.list <- lapply(mgx.list, function(mg) as.vector(LaplacesDemon::rdirichlet(1, rep(1, mg$m))))

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
  n.rho <- dim * (dim - 1)/2
  n.hyper <- 1 + K + L

  col.alpha <- 1 #intercept alpha
  col.delta <- 1 + seq_len(K) #delta_1,...,delta_K (endo coefficients)
  col.beta <- if (L > 0 ){
    1 + K + seq_len(L)
  } else{ # exo coefficients beta_1,...,beta_L. Empty if L=0
    integer (0)
  }

  col.sigma2 <- 1 + K + L + 1 #error variance sigma^2

  col.rho <- col.sigma2 + seq_len(n.rho) #upper triangle of Phi

  col.sa <- tail(col.rho, 1) + 1 #'horseshoe' hyperprior variances for alpha

  col.sb.d <- col.sa + seq_len(K) #'horseshoe' hyperprior variances for delta_1,..., delta_K

  col.sb.b <- if (L >0 ){ #'horseshoe' hyperprior variances for  beta_1,..,beta_L
    col.sa + K + seq_len(L)
  } else{
    integer(0)
  }

  n.cols <- tail(col.rho, 1) + n.hyper

  chain <- matrix(NA_real_, nrow = num.iterations + 1 , ncol = n.cols)

  #Input first row with the starting values
  chain[1, col.alpha] <- alpha.cur
  chain[1, col.delta] <- delta.cur
  if ( L > 0 ) chain[1, col.beta] <- beta.cur
  chain [1, col.sigma2] <- sigma2.cur
  chain[1, col.rho] <-copulaBayesMatrixtoVector(Phi.cur)
  chain[1, col.sa] <- 1000
  chain[1, col.sb.d] <- rep (1000, K)
  if (L > 0 ) chain[1, col.sb.b] <- rep(1000, L)

  #MCMC loop

  for (i in seq_len(num.iterations)){

    if (verbose && i %% 500 == 0){
      message( "Iteration: ", i, " | alpha = ", round(chain[i, col.alpha], 3), " | delta[1] =",
                round(chain[i, col.delta[1]], 3), " | sigma2 = ", round(chain[i, col.sigma2], 3))
    }

    #Retrieving current hyperprior variances from chain
    sa.cur <- chain[i, col.sa]
    sb.d.cur <- chain[i, col.sb.d]
    sb.b.cur <- if (L >0){
      chain[i, col.sb.b]
    } else {
      numeric(0)
    }

    #pre computing precision amtrix quantities used in both MH steps
    invPhi <- solve(Phi.cur)
    A <- invPhi - diag(dim) #sigma^{-1} - I from eq. 3
    invPhi33 <- invPhi[dim, dim] #(e,e) entry


    #MH for (alpha, delta, beta) - IWLS proposal
    #working weight M_I = 2/sigm^2 (W13 from Annex)
    #sigma_prop = (sigma^2 / 2) * solve (X'X)


    #MH step for alpha, delta and beta
    #From appendix D algorithm 2nd step


    e.cur <- y - alpha.cur - z %*% delta.cur - x %*% beta.cur

    #Score vecotr from W11
    nu.I <- copulaBayesScoreEta(A, xi.z, xi.x, e.cur, sigma2.cur, dim)

    #Design matrix X = (1, z, x)
    XX.aug <- cbind(1, z, x)
    sigma.prop <- (sigma2.cur/2) * solve(crossprod(XX.aug))
    sigma.prop <- (sigma.prop + t (sigma.prop))/2

    coef.cur <- c(alpha.cur, delta.cur, beta.cur)
    mu.cur <- as.vector(coef.cur + sigma.prop %*% ( t (XX.aug) %*% nu.I))
    coef.prop <- as.vector (mvtnorm::rmvnorm(1, mu.cur, sigma.prop))

    alpha.prop <- coef.prop[1]
    delta.prop <- coef.prop[1 + seq_len(K)]
    beta.prop <- if (L > 0) coef.prop[1 + K + seq_len(L)] else numeric(0)

    e.prop <- y - alpha.prop - z %*% delta.prop - x %*% beta.prop # reverse proposal mean (from the proposed values back to current)

    nu.prop <- copulaBayesScoreEta( A, xi.z, xi.x, e.prop, sigma2.cur, dim)
    mu.prop <- as.vector(coef.prop + sigma.prop %*% (t (XX.aug) %*% nu.prop))

    #Log Hastings ratio= log q(cur|prop) - log q(prop|cur)
    log.q.fwd <- mvtnorm :: dmvnorm(coef.prop, mu.cur, sigma.prop, log = TRUE)
    log.q.rev <- mvtnorm::dmvnorm(coef.cur, mu.prop, sigma.prop, log = TRUE)

    #log posterior at current anf proposed values
    logpost.cur <- copulaBayeslogpost(alpha.cur, delta.cur, beta.cur, sigma2.cur, Phi.cur, xi.z, xi.x,
                                      sa.cur, sb.d.cur, sb.b.cur, y, z, x)

    logpost.prop <- copulaBayeslogpost(alpha.prop, delta.prop, beta.prop, sigma2.cur,
                                       Phi.cur, xi.z, xi.x, sa.cur, sb.d.cur, sb.b.cur, y, z, x)

    log.mh <- logpost.prop - logpost.cur + log.q.rev - log.q.fwd #MH acceptance step

    if ( !is.finite(log.mh)) log.mh <- - Inf

    if (log(runif(1)) < log.mh){
      alpha.cur <- alpha.prop
      delta.cur <- delta.prop
      beta.cur <- beta.prop
    }

    chain[i + 1, col.alpha] <- alpha.cur
    chain [i + 1, col.delta ] <- delta.cur

    if ( L > 0) chain[i +1, col.beta] <- beta.cur

    #MH for log(sigma^2 ) - Laplace proposal

    e.cur <- y - alpha.cur - z %*% delta.cur - x %*% beta.cur
    sc.cur <- copulaBayesScorelogsigma2(A, invPhi33, xi.z, xi.x, e.cur, sigma2.cur)

    tau.cur <- log(sigma2.cur)
    P.cur <- if (sc.cur$f2 < 0 ) - 1/sc.cur$f2 else 1.0
    mu.tau <- P.cur * sc.cur$Score + tau.cur

    tau.new <- rnorm(1, mu.tau, sqrt(P.cur))
    sig2.new <- exp(tau.new)

    sc.new <- copulaBayesScorelogsigma2(A, invPhi33, xi.z, xi.x, e.cur, sig2.new) #reverse proposal from new back to current
    #evalutaed at sig2.new

    P.new <- if (sc.new$f2 < 0 ) -1/sc.new$f2 else 1.0
    mu.new.t <- P.new * sc.new$Score + tau.new

    q.ratio <- dnorm(tau.cur, mu.new.t, sqrt(P.new), log = TRUE) -
      dnorm(tau.new, mu.tau, sqrt(P.cur), log = TRUE)

    logpost.cur2 <- copulaBayeslogpost(alpha.cur, delta.cur, beta.cur,
                                    sigma2.cur, Phi.cur, xi.z, xi.x,
                                    sa.cur, sb.d.cur, sb.b.cur, y, z, x)
    logpost.new2 <- copulaBayeslogpost(alpha.cur, delta.cur, beta.cur,
                                    sig2.new, Phi.cur, xi.z, xi.x,
                                    sa.cur, sb.d.cur, sb.b.cur, y, z, x)

    log.mh.s <- logpost.new2 - logpost.cur2 + q.ratio
    if (!is.finite(log.mh.s)) log.mh.s <- -Inf

    if (log(runif(1)) < log.mh.s) sigma2.cur <- sig2.new
    chain[i+1, col.sigma2] <- sigma2.cur


    # Gibbs for W (Appendix A, W5)
    #W has na inverse Wishart prior from eq. 8 and the likelihood (W1 to W3 from appendix A) also
    #has a Wishart shape.
    # W | A ~ W^{-1}(sum(xi_i xi_i') + I, N + K+L+1)

    #Gibbs step draw directly from the exact full conditional. Afterwards, it is converted to a correlation
    #matrix Phi (from Appendix A) Cholesky factorisation.

    #exogeneiy restriction will also be enforced by setting Phi[x_l, e] = 0 for all l-x uncorrelated with e (from Haschka section 3, page 523)

    e.cur2  <- y - alpha.cur - z %*% delta.cur - x %*% beta.cur

    #normal score of standardised error. The last element of the xi vector
    xi.e    <- qnorm(pnorm(as.vector(e.cur2) / sqrt(sigma2.cur)))
    xi.e    <- pmin(pmax(xi.e, -8), 8)

    #stacking all normal scores and ordering
    xi.all  <- cbind(xi.z,xi.x , matrix(xi.e, N, 1))

    W.new    <- LaplacesDemon::rinvwishart( #inverse Wishart from W5 Appendix A
      nu = N + dim,
      S  = diag(dim) + crossprod(xi.all)
    )
    #coverting covatiance W to correlation Phi
    sds      <- sqrt(diag(W.new))
    D.inv    <- diag(1 / sds)
    Phi.cur  <- D.inv %*% W.new %*% D.inv

    # Enforce exogeneity
    for (l in seq_len(L)) {
      Phi.cur[K + l, dim] <- 0  #zero out (x_l, e) correlations
      Phi.cur[dim, K + l] <- 0
    }

    chain[i +1, col.rho] <- copulaBayesMatrixtoVector(Phi.cur)


    #  Gibbs for Dirichlet masses lambda (Appendix B)

    epsilon <- mvtnorm::rmvnorm(N, mean = rep(0, dim), sigma = Phi.cur)

    for (k in seq_len(K)) {
      lambdaz.list[[k]] <- copulaBayesDrawLambda(pnorm(epsilon[, k]),
                                                mgz.list[[k]])
      xi.z[, k]      <- copulaBayesConverter(lambdaz.list[[k]],
                                                   mgz.list[[k]])
    }

    for (l in seq_len(L)) {
        lambdax.list[[l]] <- copulaBayesDrawLambda(pnorm(epsilon[, K + l]),
                                                  mgx.list[[l]])
        xi.x[, l]      <- copulaBayesConverter(lambdax.list[[l]],
                                                     mgx.list[[l]])
    }


    # Horseshoe hyperprior updates (eq. 5)

    chain[i+1, col.sa] <- invgamma::rinvgamma(
      1, shape = 0.501, rate = (alpha.cur^2 + 0.002) / 2)

    for (k in seq_len(K))
      chain[i+1, col.sb.d[k]] <- invgamma::rinvgamma(
        1, shape = 0.501, rate = (delta.cur[k]^2 + 0.002) / 2)

    if (L > 0)
      for (l in seq_len(L))
        chain[i+1, col.sb.b[l]] <- invgamma::rinvgamma(
          1, shape = 0.501, rate = (beta.cur[l]^2 + 0.002) / 2)

  }

  attr(chain, "col.alpha")  <- col.alpha
  attr(chain, "col.delta")  <- col.delta
  attr(chain, "col.beta")   <- col.beta
  attr(chain, "col.sigma2") <- col.sigma2
  attr(chain, "col.rho")    <- col.rho
  attr(chain, "K")          <- K
  attr(chain, "L")          <- L
  attr(chain, "dim")        <- dim

  chain
}

#' @importFrom stats lm residuals var runif rnorm
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom LaplacesDemon rinvwishart
#' @importFrom invgamma rinvgamma
#' @importFrom copula P2p
#'

# MCMC sampler of the copula model from Appendix D

copulaBayesMCMC <- function(y, z, x, num.iterations, verbose){

  N <- length(y)
  K <- ncol(z)
  L <- if(!is.null(x) && ncol(x) > 0){
    ncol(x)
  } else {
    0
  }

  #building margin structures per regressor
  #this is fixed for the entire run

  mgz.list <- lapply(seq_len(K), function(k) copulaBayesMargin(z[,k]))

  mgx.list <- if(L>0){
    lapply(seq_len(L), function(l) copulaBayesMargin(x[, l]))
  } else {
    list()
  }

  # ols starting values

  dat.ols <- if (L >0){
    data.frame(y =y, z =z, x=x)
  } else {
    data.frame (y=y, z =z)
  }

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

  # setting sigma^[0] = I

  Phi.cur <- diag(dim)

  #initial lambda from Dir(1,...,1)
  lambdaz.list <- lapply(mgz.list, function(mg) as.vector(LaplacesDemon::rdirichlet(1, rep(1, mg$m))))
  lambdax.list <- if (L >0){
    lapply(mgx.list, function(mg) as.vector(LaplacesDemon::rdirichlet(1, rep(1, mg$m))))
  } else{
    list()
  }

  #initial normal scores from the starting lambda
  xi.z <- matrix(NA_real_, N,K)
  for (k in seq_len(K))
    xi.z[,k] <- copulaBayesConverter(lambdaz.list[[k]], mgz.list[[k]]) #converting lambda to xi

  xi.x <- if(L > 0 ){
    m <- matrix(NA_real_, N, L)
    for (l in seq_len(L))
      m[, l] <- copulaBayesConverter(lambdax.list[[l]], mgx.list[[l]])
    m
  } else{
    matrix(numeric(0), N, 0)
  }



}

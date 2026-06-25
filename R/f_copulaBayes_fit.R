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

  n.rho <- dim * (dim - 1)/2
  n.hyper <- 1 + K + L

  col.alpha <- 1
  col.delta <- 1 + seq_len(K)
  col.beta <- if (L > 0 ){
    1 + K + L + seq_len(L)
  } else{
    integer (0)
  }

  col.sigma2 <- 1 + K + L + 1

  col.rho <- col.sigma2 + seq_len(n.rho)

  col.sa <- tail(col.rho, 1) + 1

  col.sb.d <- col.sa + seq_len(K)

  col.sb.b <- if (L >0 ){
    col.sa + K + seq_len(L)
  } else{
    integer(0)
  }

  n.cols <- tail(col.rho, 1) + n.hyper

  chain <- matrix(NA_real_, nrow = num.iterations + 1 , ncol = n.cols)

  chain[1, col.alpha] <- alpha.cur
  chain[1, col.delta] <- delta.cur
  if ( L > 0 ) chain[1, col.beta] <- beta.cur
  chain [1, col.sigma2] <- sigma2.cur
  chain[1, col.rho] <- copula::P2p(Phi.cur) ##maybe change
  chain[1, col.sa] <- 1000
  chain[1, col.sb.d] <- rep (1000, K)
  if (L > 0 ) chain[1, col.sb.b] <- rep(1000, L)

  #MCMC loop

  for (i in seq_len(num.iterations)){

    if (verbose && i %% 500 == 0){
      message( "Iteration: ", i, "| alpha = ", round(chain[i, col.alpha], 3), "| delta[1] =",
               , round(chain[iter, col.delta[1]], 3), "| sigma2 = ", round(chain[i, col.sigma2], 3))
    }

    sa.cur <- chain[i, col.sa]
    sb.d.cur <- chain[i, col.sb.d]
    sb.b.cur <- if (L >0){
      chain[i, col.sb.b]
    } else {
      numeric(0)
    }

    invPhi <- solve(Phi.cur)
    A <- invPhi - diag(dim)
    invPhi33 <- invPhi[dim, dim]


    #MH for (alpha, delta, beta) - IWLS proposal
    #working weight M_I = 2/sigm^2 (W13 from Annex)
    #sigma_prop = (sigma^2 / 2) * solve (X'X)

    e.cur <- y - alpha.cur - z %*% delta.cur - if (L > 0 ) x %*% beta.cur else 0

    nu.I <- copulaBayesScoreEta(invPhi, A, xi.z, xi.z, e.cur, sigma2.cur, dim)

    XX.aug <- if ( L > 0){
      cbind(1, z, x)
    } else {
      cbind (1, z)
    }

    sigma.prop <- (sigma2.cur/2) * solve(crossprod(XX.aug))
    sigma.prop <- (sigma.prop + t (sigma.prop))/2

    coef.cur <- c(alpha.cur, delta.cur, if (L > 0) beta.cur else NULL)
    mu.cur <- as.vector(coef.cur + sigma.prop %*% ( t (XX.aug) %*% nu.I))
    coef.prop <- as.vector (mvtnorm::rmvnorm(1, mu.cur, sigma.propr))

    alpha.prop <- coef.prop[1]
    delta.prop <- coef.prop[1 + seq_len(K)]
    beta.prop <- if (L > 0) coef.prop[1 + K + seq_len(L)] else numeric(0)

    e.prop <- y - alpha.prop - z %*% delta.prop - if (L > 0) x %*% beta.prop else 0

    nu.prop <- copulaBayesScoreEta(invPhi, A, xi.z, xi.x, e.prop, sigma2.cur, dim)
    mu.prop <- as.vector(coef.prop + sigma.prop %*% (t (XX.aug) %*% nu.prop))

    log.q.fwd <- mvtnorm :: dmvnorm(coef.prop, mu.cur, sigma.prop, log = TRUE)
    log.q.rev <- mvtnorm::dmvnorm(coef.cur, mu.prop, sigma.prop, log = TRUE)

    logpost.cur <- copulaBayeslogpost(alpha.cur, delta.cur, beta.cur, sigma2.cur, Phi.cur, xi.z, xi.x,
                                      sa.cur, sb.d.cur, sb.b.cur, y, z, x)

    logpost.prop <- copulaBayeslogpost(alpha.prop, delta.prop, beta.prop, sigma2.cur,
                                       Phi.cur, xi.z, xi.x, sa.cur, sb.d.cur, sb.b.cur, y, z, x)

    log.mh <- lp.prop - logpost.cur + log.q.rev - log.q.fwd

    if ( !is.finite(log.mh)) log.mh <- -infert

    if (log(runif(1)) < log.mh){
      alpha.cur <- alpha.prop
      delta.cur <- delta.prop
      if (L > 0) beta.cur <- beta.prop
    }

    chain[i + 1, col.alpha] <- alpha.cur
    chain [i + 1, col.delta ] <- delta.cur

    if ( L > 0) chain[i +1, col.beta] <- beta.cur

    #MH for log(sigma^2 ) - Laplace proposal

    e.cur <- y - alpha.cur - z %*% delta.cur - if (L > 0 ) x %*% beta.cur else 0

    sc.cur <- copulaBayesScorelogsigma2(A, invPhi33, xi.z, xi.x, e.cur, sigma2.cur)

    tau.cur <- log(sigma2.cur)

    P.cur <- if (sc.cur$ f2 < 0 ) - 1/sc.cur$f2 else 1.0
    mu.tau <- P.cur * sc.cur$Score + tau.cur

    tau.new <- rnorm(1, mu.tau, sqrt(P.cur))
    sig2.new <- exp(tau.new)

    sc.new <- copulaBayesScorelogsigma2(A, invPhi33, xi.z, xi.x, e.cur, sig2.new)

    P.new <- if (sc.new$f2 < 0 ) -1/sc.new$f2 else 1.0
    mu.new.t <- P.new * sc.new$Score + tau.new

    q.ratio <- dnorm(tau.cur, mu.new.t, sqrt(P.new), log = TRUE) -
      dnorm(tau.new, mu.tau, sqrt(P.cur), log = TRUE)

    logpost.cur2 <- copulaBayes_log_post(alpha.cur, delta.cur, beta.cur,
                                    sigma2.cur, Phi.cur, xi.z, xi.x,
                                    sa.cur, sb.d.cur, sb.b.cur, y, z, x)
    logpost.new2 <- copulaBayes_log_post(alpha.cur, delta.cur, beta.cur,
                                    sig2.new, Phi.cur, xi.z, xi.x,
                                    sa.cur, sb.d.cur, sb.b.cur, y, z, x)

    log.mh.s <- logpost.new2 - logpost.cur2 + q.ratio
    if (!is.finite(log.mh.s)) log.mh.s <- -Inf

    if (log(runif(1)) < log.mh.s) sigma2.cur <- sig2.new
    chain[iter+1, col.sigma2] <- sigma2.cur


    # Gibbs for W (Appendix A, W5)
    # W | A ~ W^{-1}(sum(xi_i xi_i') + I, N + K+L+1)

    e.cur2  <- y - alpha.cur - z %*% delta.cur - if (L > 0) x %*% beta.cur else 0
    xi.e    <- qnorm(pnorm(as.vector(e.cur2) / sqrt(sigma2.cur)))
    xi.e    <- pmin(pmax(xi.e, -8), 8)
    xi.all  <- cbind(xi.z, if (L > 0) xi.x else NULL, matrix(xi.e, N, 1))

    W.new    <- LaplacesDemon::rinvwishart(
      nu = N + dim,
      S  = diag(dim) + crossprod(xi.all)
    )
    sds      <- sqrt(diag(W.new))
    D.inv    <- diag(1 / sds)
    Phi.cur  <- D.inv %*% W.new %*% D.inv

    # Enforce exogeneity
    for (l in seq_len(L)) {
      Phi.cur[K + l, dim] <- 0  #zero out (x_l, e) correlations
      Phi.cur[dim, K + l] <- 0
    }

    chain[i +1, col.rho] <- copula::P2p(Phi.cur)


    #  Gibbs for lambda (Appendix B)

    epsilon <- mvtnorm::rmvnorm(N, mean = rep(0, dim), sigma = Phi.cur)

    for (k in seq_len(K)) {
      lamz.list[[k]] <- copulaBayes_draw_lambda(pnorm(epsilon[, k]),
                                                mgZ.list[[k]])
      xi.z[, k]      <- copulaBayes_xi_from_lambda(lamZ.list[[k]],
                                                   mgZ.list[[k]])
    }

    if (L > 0)
      for (l in seq_len(L)) {
        lamX.list[[l]] <- copulaBayes_draw_lambda(pnorm(epsilon[, K + l]),
                                                  mgX.list[[l]])
        xi.x[, l]      <- copulaBayes_xi_from_lambda(lamX.list[[l]],
                                                     mgX.list[[l]])
      }

    # Horseshoe hyperprior updates (eq. 5)

    chain[i+1, col.sa] <- invgamma::rinvgamma(
      1, shape = 0.501, rate = (alpha.cur^2 + 0.002) / 2)

    for (k in seq_len(K))
      chain[i+1, col.sb.d[k]] <- invgamma::rinvgamma(
        1, shape = 0.501, rate = (delta.cur[k]^2 + 0.002) / 2)

    if (L > 0)
      for (l in seq_len(L))
        chain[iter+1, col.sb.b[l]] <- invgamma::rinvgamma(
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

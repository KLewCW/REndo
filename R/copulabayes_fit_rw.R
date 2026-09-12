
# MCMC sampler of the copula model from Appendix D
#
# sampler cycles through 4 blocks per iteration (eq.11):
# 1. Metropolis-Hastings algorithm (MH) step for alpha, delta, beta using IWLS proposal from appendix C
# 2. MH step for sigma^s using Laplace proposal (Appendix C)
# 3. Gibbs step for correlation matrix W (Appendix A W5, Wishart distribution)
# 4. Gibbs step for Dirichlet masses lambda (Appendix W7)
#
#Note: for stcopula.eps.draw one and 2 :  the MH (random walk) replaces the IWLS and Laplace proposals.
#From MCMCpack, random walk and Hastings ratio have been used.
#
#modification done from previous code because was a bit lost:
# X has now been kept as a N x L matrix (when L =0, it has 0 columns, i.e., X %*% beta = 0 )
#
#
#' @importFrom utils tail
#' @importFrom stats lm residuals var runif rnorm rgamma
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack riwish rdirichlet
copulabayes_mcmc_rw <- function(y, z, x, num.iterations, verbose) {
  N <- length(y)
  K <- ncol(z)
  L <- ncol(x)
  copula.dim <- K + L + 1L #copula dimension = K endogenous + L exogenous + error

  #building margin structures per regressor (from Appendix D step 2)
  #This describe the unique values (UV) of each z_k and x_l and
  #this is fixed for the entire MCMC run.
  #note: observed data never change.
  margin.endo.list <- lapply(seq_len(K), function(k) copulabayes_margin(z[, k]))
  margin.exo.list <- lapply(seq_len(L), function(l) copulabayes_margin(x[, l]))

  # From appendix D: ols starting values
  # Research paper says starting value can be arbitrary. OLS is used here for conveniency
  #and close to the posterior mode when endogeneity is weak
  dat.ols <- data.frame(y = y, z = z, x = x)
  mod.ols <- lm(y ~ ., data = dat.ols)
  coef.ols <- coef(mod.ols)

  intercept.cur <- coef.ols["(Intercept)"]
  coef.endo.cur <- coef.ols[1L + seq_len(K)]
  if (L > 0L) {
    coef.exo.cur <- coef.ols[1L + K + seq_len(L)]
  } else {
    coef.exo.cur <- numeric(0L)
  }

  error.var.cur <- var(residuals(mod.ols))

  # setting sigma^[0] = identiy matrix (I)
  # Here we start with no correlation assumed.
  copula.cor.cur <- diag(copula.dim)

  #initial lambda from Dir(1,...,1) are drawn
  masses.endo.list <- lapply(margin.endo.list, function(mg) {
    as.vector(MCMCpack::rdirichlet(1, rep(1, mg$m)))
  })
  masses.exo.list<- lapply(margin.exo.list, function(mg) {
    as.vector(MCMCpack::rdirichlet(1, rep(1, mg$m)))
  })

  #initial normal scores from the starting lambda
  # xi_{varpi, i} = Phi^{-1} (u_{varpi, i} (lambda_varpi)) from section 3.1 of Haschka 2025, page 522
  scores.endo<- matrix(NA_real_, N, K)
  for (k in seq_len(K)) {
    scores.endo[, k] <- copulabayes_converter(masses.endo.list[[k]], margin.endo.list[[k]])
  } #converting lambda to xi

  scores.exo<- matrix(NA_real_, N, L)
  for (l in seq_len(L)) {
    scores.exo[, l] <- copulabayes_converter(lambdax.list[[l]], margin.exo.list[[l]])
  }

  # Chain storage: one row per iteration.
  # Column layout (fixed throughout)
  n.copula.cor <- copula.dim * (copula.dim - 1L) / 2L
  n.hyper <- 1L + K + L
  n.coef <- 1L + K + L #alpha + delta+ beta

  col.intercept <- 1L #intercept alpha
  col.coef.endo <- 1L + seq_len(K) #delta_1,...,delta_K (endo coefficients)

  if (L > 0) {
    col.coef.exo <- 1L + K + seq_len(L)
  } else {
    col.coef.exo <- integer(0L) # exo coefficients beta_1,...,beta_L. Empty if L=0
  }

  col.error.var <- 1L + K + L + 1L #error variance sigma^2

  col.copula.cor <- col.error.var + seq_len(n.copula.cor) #upper triangle of Phi

  col.var.alpha <- tail(col.copula.cor, 1L) + 1L #'horseshoe' hyperprior variances for alpha

  col.var.delta <- col.var.alpha + seq_len(K) #'horseshoe' hyperprior variances for delta_1,..., delta_K

  col.var.beta <- if (L > 0L) {
    # 'horseshoe' hyperprior variances for  beta_1,..,beta_L
    col.var.alpha + K + seq_len(L)
  } else {
    integer(0L)
  }

  n.cols <- tail(col.copula.cor, 1L) + n.hyper

  chain <- matrix(NA_real_, nrow = num.iterations + 1L, ncol = n.cols)

  #Input first row with the starting values
  chain[1L, col.intercept] <- intercept.cur
  chain[1L, col.coef.endo] <- coef.endo.cur
  if (L > 0L) {
    chain[1L, col.coef.exo] <- coef.exo.cur
  }
  chain[1L, col.error.var] <- error.var.cur
  chain[1L, col.copula.cor] <- copulabayes_matrix2vector(copula.cor.cur)

  chain[1L, col.var.alpha] <- 1000
  chain[1L, col.var.delta] <- rep(1000, K)
  if (L > 0L) {
    chain[1L, col.var.beta] <- rep(1000, L)
  }

  mh.scale <- 1.0
  n.mh.accepted <- 0L

  #MCMC loop

  for (i in seq_len(num.iterations)) {
    if (verbose && i %% 500L == 0L) {
      message(
        "Iteration: ",
        i,
        " | intercept = ",
        round(chain[i, col.intercept], 3L),
        " | coef.endo[1] =",
        round(chain[i, col.coef.endo[1L]], 3L),
        " | error.var = ",
        round(chain[i, col.error.var], 3L),
        " | accept rate = ",
        round(n.mh.accepted / i, 2L)
      )
    }

    #Retrieving current hyperprior variances from chain
    var.alpha.cur <- chain[i, col.var.alpha]
    var.delta.cur <- chain[i, col.var.delta]

    if (L > 0L) {
      var.beta.cur <- chain[i, col.var.beta]
    } else {
      var.beta.cur <- numeric(0L)
    }

    #Step 1 and 2 of the iteration: random walk mH for alpha, delta, beta, log sigma^2
    #Using the MCMCpack metropolis principle: propose, then evalute ratio and accept.
    #Parameters are sampled jointly in the reparametrised space (log sigma^2) so that
    #the proposal is unconstrained. Replaces the custom IWLS and Laplace proposals from
    #the repository without any score computation.

    theta.cur <- c(intercept.cur, coef.endo.cur, coef.exo.cur, log(error.var.cur))
    theta.prop <- theta.cur + rnorm(length(theta.cur), mean = 0, sd = mh.scale)

    logpost.cur <- copulabayes_logpostparam(
      theta.cur,
      copula.cor.cur,
      scores.endo,
      scores.exo,
      var.alpha.cur,
      var.delta.cur,
      var.beta.cur,
      y,
      z,
      x,
      K,
      L
    )
    logpost.prop <- copulabayes_logpostparam(
      theta.prop,
      copula.cor.cur,
      scores.endo,
      scores.exo,
      var.alpha.cur,
      var.delta.cur,
      var.beta.cur,
      y,
      z,
      x,
      K,
      L
    )

    log.mh <- logpost.prop - logpost.cur
    if (!is.finite(log.mh)) {
      log.mh <- -Inf
    }

    if (log(runif(1L)) < log.mh) {
      theta.cur <- theta.prop
      n.mh.accepted <- n.mh.accepted + 1L
    }

    intercept.cur <- theta.cur[1L]
    coef.endo.cur <- theta.cur[seq(2L, 1L + K)]

    if (L > 0L) {
      coef.exo.cur <- theta.cur[seq(2L + K, 1L + K + L)]
    } else {
      coef.exo.cur <- numeric(0L)
    }

    error.var.cur <- exp(theta.cur[1L + K + L + 1L])

    chain[i + 1L, col.intercept] <- intercept.cur
    chain[i + 1L, col.coef.endo] <- coef.endo.cur
    if (L > 0L) {
      chain[i + 1L, col.coef.exo] <- coef.exo.cur
    }
    chain[i + 1L, col.error.var] <- error.var.cur

    #tuning every 50 iterations (Robbins-Monro)
    if (i %% 50L == 0L) {
      rate <- n.mh.accepted / i
      if (rate < 0.20) { # 20%-40% acceptance rate
        mh.scale <- mh.scale * 0.9
      }
      if (rate > 0.40) mh.scale <- mh.scale * 1.1
    }

    #Step 3: the Gibbs for copula correction matrix W (appendix A W5)
    resid.cur <- y - intercept.cur - z %*% coef.endo.cur - x %*% coef.exo.cur
    scores.error <- pmin(pmax(qnorm(pnorm(as.vector(resid.cur) / sqrt(error.var.cur))), -8), 8)
    scores.all <- cbind(scores.endo, scores.exo, scores.error)

    copula.cov.draw <- MCMCpack::riwish(N + copula.dim, diag(copula.dim) + crossprod(scores.all))
    cor.sd <- sqrt(diag(copula.cov.draw))
    cor.normaliser <- diag(1 / cor.sd)
    copula.cor.cur <- cor.normaliser %*% copula.cov.draw %*% cor.normaliser

    # Enforce exogeneity
    for (l in seq_len(L)) {
      copula.cor.cur[K + l, copula.dim] <- 0 #zero out (x_l, e) correlations (restriction mentioned in Haschka 2026 section 3.1)
      copula.cor.cur[copula.dim, K + l] <- 0
    }

    chain[i + 1L, col.copula.cor] <- copulabayes_matrix2vector(copula.cor.cur)

    #Step 4: Gibbs for Dirichlet masses lambda (Appendix B W7)
    #Drawing correlated normals from current phi, then copulabayes_drawlambda updates
    #each variable's mass vector

    copula.eps.draw <- mvtnorm::rmvnorm(N, mean = rep(0, copula.dim), sigma = copula.cor.cur)

    for (k in seq_len(K)) {
      masses.endo.list[[k]] <- copulabayes_drawlambda(pnorm(copula.eps.draw[, k]), margin.endo.list[[k]])
      scores.endo[, k] <- copulabayes_converter(masses.endo.list[[k]], margin.endo.list[[k]])
    }
    for (l in seq_len(L)) {
      lambdax.list[[l]] <- copulabayes_drawlambda(pnorm(copula.eps.draw[, K + l]), margin.exo.list[[l]])
      scores.exo[, l] <- copulabayes_converter(lambdax.list[[l]], margin.exo.list[[l]])
    }

    #Horseshoe hyperprior updates (From equation 5 of Haschka 2025)

    chain[i + 1L, col.var.alpha] <- 1 /
      rgamma(1L, shape = 0.501, rate = (intercept.cur^2 + 0.002) / 2)

    for (k in seq_len(K)) {
      chain[i + 1L, col.var.delta[k]] <- 1 /
        rgamma(1L, shape = 0.501, rate = (coef.endo.cur[k]^2 + 0.002) / 2)
    }
    if (L > 0L) {
      for (l in seq_len(L)) {
        chain[i + 1L, col.var.beta[l]] <- 1 /
          rgamma(1L, shape = 0.501, rate = (coef.exo.cur[l]^2 + 0.002) / 2)
      }
    }
  }

  attr(chain, "col.intercept") <- col.intercept
  attr(chain, "col.coef.endo") <- col.coef.endo
  attr(chain, "col.coef.exo") <- col.coef.exo
  attr(chain, "col.error.var") <- col.error.var
  attr(chain, "col.copula.cor") <- col.copula.cor
  attr(chain, "K") <- K
  attr(chain, "L") <- L
  attr(chain, "copula.dim") <- copula.dim

  return(chain)
}

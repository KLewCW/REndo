#' @importFrom stats lm residuals var runif rnorm rgamma qnorm pnorm
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom MCMCpack riwish rdirichlet

copulabayes_mcmc_iwls <- function(y,z,x,num.iterations,verbose){

  N <- length(y)
  K <- ncol(z)
  L <- ncol(x)

  cop.dim <- K + L + 1L

  ##this part is same as random walk
  margin.endo.list <- lapply(seq_len(K), function(k) copulabayes_margin(z[, k]))
  margin.exo.list <- lapply(seq_len(L), function(l) copulabayes_margin(x[,l]))

  #ols starting values same as RW

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

  copula.cor.cur <- diag(cop.dim)

  ## initial Dirichlet masses same as RW

  masses.endo.list <- lapply(margin.endo.list, function(mg){
    as.vector(MCMCpack::rdirichlet(1, rep(1, mg$m)))
  })

  masses.exo.list <- lapply(margin.exo.list, function(mg){
    as.vector(MCMCpack::rdirichlet(1, rep(1, mg$m)))
  })

  scores.endo <- matrix(NA_real_, N, K)
  for(k in seq_len(K)){
    scores.endo[, k] <- copulabayes_converter(masses.endo.list[[k]], margin.endo.list[[k]])
  }

  scores.exo <- matrix(NA_real_, N,L)
  for(l in seq_len(L)){
    scores.exo[, l] <- copulabayes_converter(masses.exo.list[[l]], margin.exo.list[[l]])
  }

  # chain storage (same as RW)

  n.copula.cor <- cop.dim * (cop.dim - 1L) / 2L
  n.hyper <- 1L + K + L
  col.intercept <- 1L
  col.coef.endo <- 1L + seq_len(K)

  if (L > 0L){
    col.coef.exo <- 1L + K + seq_len(L)
  } else {
    col.coef.exo <- integer(0L)
  }

  col.error.var <- 1L + K + L + 1L
  col.copula.cor <- col.error.var + seq_len(n.copula.cor)
  col.var.alpha <- tail(col.copula.cor, 1L) + 1L
  col.var.delta <- col.var.alpha + seq_len(K)


  if (L > 0L){
    col.var.beta <- col.var.alpha + K + seq_len(L)
  } else{
    col.var.beta <- integer(0L)
  }

  n.cols <- tail(col.copula.cor, 1L) + n.hyper

  chain <- matrix(NA_real_, nrow = num.iterations + 1L, ncol = n.cols)
  chain[1L, col.intercept] <- intercept.cur
  chain[1L, col.coef.endo] <- coef.endo.cur

  if (L > 0L) chain[1L, col.coef.exo]  <- coef.exo.cur
  chain[1L, col.error.var] <- error.var.cur
  chain[1L, col.copula.cor] <- copulabayes_matrix2vector(copula.cor.cur)
  chain[1L, col.var.alpha] <- 1000
  chain[1L, col.var.delta] <- rep(1000, K)
  if (L > 0L) chain[1L, col.var.beta] <- rep(1000, L)

  n.mh.accepted.coef <- 0L  # accepted proposals for coefficients
  n.mh.accepted.sigma <- 0L  # accepted proposals for sigma^2

  #IWLS

  #X.design = [ 1 z x]
  # (N X ( 1 + K + L)) full design matrix

  X.design <- cbind(1, z, x)

  #MCMC loop

  # ── MCMC loop ─────────────────────────────────────────────────────────────
  for (i in seq_len(num.iterations)) {

    if (verbose && i %% 500L == 0L) {
      message("Iteration: ", i,
              " | intercept = ", round(chain[i, col.intercept], 3L),
              " | coef.endo[1] = ", round(chain[i, col.coef.endo[1L]], 3L),
              " | error.var = ", round(chain[i, col.error.var], 3L),
              " | accept (coef) = ", round(n.mh.accepted.coef  / i, 2L),
              " | accept (sigma) = ", round(n.mh.accepted.sigma / i, 2L))
      }

    #current horseshoe variances
    var.alpha.cur <- chain[i, col.var.alpha]
    var.delta.cur <- chain[i, col.var.delta]

    if (L > 0L){
      var.beta.cur <- chain[i, col.var.beta]
    } else{
      var.beta.cur <- numeric(0L)
    }

    # Step 1&2 IWLS proposal for regression coefficients
    # Haschka (2025) Appendix C, W11 and W13

    # The IWLS step treats the copula log-density as a generalised linear
    # model with working weights M_i = 2/sigma^2 (negative of second derivative in W13)
    # and working responses nu_i (W11).
    #  covariance is Sigma_prop = (sigma^2 / 2) * (X'X)^{-1}
    # proposal mean is the IWLS update: mu_prop = coef.cur + Sigma_prop * X' * nu


    #current residuals and normal scores
    resid.cur  <- y - X.design %*% c(intercept.cur, coef.endo.cur, coef.exo.cur)
    scores.error.cur <- pmin(pmax( qnorm(pnorm(as.vector(resid.cur) / sqrt(error.var.cur))), -8), 8)
    scores.all.cur <- cbind(scores.endo, scores.exo, scores.error.cur)

    #copula inverse A = Phi^{-1} - I  (in score computation)
    A.cur <- solve(copula.cor.cur) - diag(cop.dim)

    #working weight (W13): M_i = 2 / sigma^2 (scalar)
    working.weight <- 2 / error.var.cur

    # Score vector for coefficients (need last element of W11 from Appendix C,
    #which is scalar):
    # last element is (1/sigma) * (A xi_i)[cop.dim] + e_i / sigma^2
    # where A = Sigma^{-1} - I and cop.dim = L + K + 1 (index of error dimension)
    A.d   <- A.cur[cop.dim, ]  # d-th row of A (error row). It is the last row of A = A[L+K+1, :]

    # the N-vector xi_i'
    nu.vec <- as.vector( scores.all.cur %*% A.d / sqrt(error.var.cur) +
        as.vector(resid.cur) / error.var.cur)

    #IWLS proposal covariance (W13)
    Sigma.prop.coef <- solve(t(X.design) %*% (working.weight * X.design))

    #IWLS proposal mean (W11)
    coef.cur.vec <- c(intercept.cur, coef.endo.cur, coef.exo.cur)
    mu.prop.coef <- coef.cur.vec + Sigma.prop.coef %*% t(X.design) %*% nu.vec

    # Draw proposal
    coef.prop.vec <- as.vector(mvtnorm::rmvnorm(1, mean = mu.prop.coef, sigma = Sigma.prop.coef))

    # Evaluating the log-posterior at current and proposed
    var.sb.cur <- c(var.alpha.cur, var.delta.cur, var.beta.cur)

    logpost.coef.cur  <- copulabayes_logpost_coef(
      coef.cur.vec,  error.var.cur, copula.cor.cur,
      scores.endo, scores.exo, var.alpha.cur, var.delta.cur, var.beta.cur,
      y, z, x, K, L)

    logpost.coef.prop <- copulabayes_logpost_coef(
      coef.prop.vec, error.var.cur, copula.cor.cur,
      scores.endo, scores.exo, var.alpha.cur, var.delta.cur, var.beta.cur,
      y, z, x, K, L)

    # Log proposal densities (needed for Hastings ratio with asymmetric proposal)
    log.q.prop.given.cur <- mvtnorm::dmvnorm(
      coef.prop.vec, mean = mu.prop.coef, sigma = Sigma.prop.coef, log = TRUE)

    # recomputing the proposal mean at the proposed point (for reverse proposal)
    resid.prop   <- y - X.design %*% coef.prop.vec

    scores.e.prop <- pmin(pmax(qnorm(pnorm(as.vector(resid.prop) / sqrt(error.var.cur))), -8), 8)

    scores.all.prop <- cbind(scores.endo, scores.exo, scores.e.prop)

    nu.vec.rev <- as.vector( scores.all.prop %*% A.d / sqrt(error.var.cur) +
        as.vector(resid.prop) / error.var.cur
    )

    mu.rev.coef <- coef.prop.vec + Sigma.prop.coef %*% t(X.design) %*% nu.vec.rev

    log.q.cur.given.prop <- mvtnorm::dmvnorm(
      coef.cur.vec, mean = mu.rev.coef, sigma = Sigma.prop.coef, log = TRUE)


    # Metropolis-Hastings acceptance ratio (asymmetric proposal)
    log.mh.coef <- (logpost.coef.prop + log.q.cur.given.prop) -
      (logpost.coef.cur  + log.q.prop.given.cur)

    if (!is.finite(log.mh.coef)) log.mh.coef <- -Inf

    if (log(runif(1L)) < log.mh.coef) {
      coef.cur.vec <- coef.prop.vec
      n.mh.accepted.coef <- n.mh.accepted.coef + 1L
    }

    intercept.cur <- coef.cur.vec[1L]
    coef.endo.cur <- coef.cur.vec[seq(2L, 1L + K)]

    if (L > 0L){
      coef.exo.cur <- coef.cur.vec[seq(2L + K, 1L + K + L)]
    } else{
      coef.exo.cur <-numeric(0L)
    }

    # Step 3: Laplace proposal for log(sigma^2)
    # Haschka (2025) Appendix C, W12 (first derivative) and W14 (second derivative)

    # The Laplace proposal uses a second-order Taylor expansion of the
    # log-posterior around log(sigma^2).

    #Note: W12 and W14 are matrix derivatives of dimension (L+K+1) x 1 and
    # (L+K+1) x (L+K+1) respectively.
    # We needthe scalar element:
    # For W12 the (L+K+1)-th element = score  for log(sigma^2) per obs
    # For W14 the (L+K+1, L+K+1) element = Hessian for log(sigma^2) per obs

    # Let d = cop.dim = K+L+1 (index of the error dimension in xi)
    # A = Sigma^{-1} - I
    # A[d,:] = last row of A (1 x cop.dim vector)
    # A[d,d] = last diagonal element of A (scalar)

    # W12 last element per observation i (derived from log L_i, from W12):
    #   s_i = (e_i / (2*sigma)) * (A[d,:] %*% xi_i) - 1/2 + e_i^2 / (2*sigma^2)

    # W14 last diagonal element per observation i (derived from W12, from W14):
    # h_i = -(e_i / (4*sigma)) * (A[d,:] %*% xi_i) - A[d,d] * e_i^2 / (4*sigma^2)- e_i^2 / (2*sigma^2)

    # Note: the paper gives d log L_i (likelihood only).
    # prior contributions are added separately (IG(0.001, 0.001) on sigma^2).
    # Prior score d log p / d log(sigma^2) = -1.001 + 0.001/sigma^2
    # Prior Hessian d^2 log p / d(log sigma^2)^2  = -0.001/sigma^2

    log.sigma2.cur <- log(error.var.cur)
    sigma.cur <- sqrt(error.var.cur)

    # Current residuals
    resid.cur <- y - X.design %*% c(intercept.cur, coef.endo.cur, coef.exo.cur)
    resid.vec <- as.vector(resid.cur)

    # Current normal score of the structural error
    scores.error.cur <- pmin(pmax(qnorm(pnorm(resid.vec / sigma.cur)), -8), 8)
    scores.all.cur <- cbind(scores.endo, scores.exo, scores.error.cur)

    # A = Sigma^{-1} - I; extract last row and last diagonal element
    A.cur <- solve(copula.cor.cur) - diag(cop.dim)
    A.d.row <- A.cur[cop.dim, ] # (1 x cop.dim): last row of A
    A.dd  <- A.cur[cop.dim, cop.dim]  # scalar: last diagonal element of A

    # N-vector: A[d,:] %*% xi_i for each observation i
    A.xi.cur <- as.vector(scores.all.cur %*% A.d.row)

    #score (W12: summed over obs + prior)
    score.obs.cur <- (resid.vec / (2 * sigma.cur)) * A.xi.cur - 0.5 +resid.vec^2 / (2 * error.var.cur)

    score.sigma <- sum(score.obs.cur) + (-1.001 + 0.001 / error.var.cur) # IG(0.001,0.001) prior

    #Hessian (W14: summed over obs + prior)
    hessian.obs.cur <- -(resid.vec / (4 * sigma.cur)) * A.xi.cur - A.dd * resid.vec^2 / (4 * error.var.cur) - resid.vec^2 / (2 * error.var.cur)
    hessian.sigma <- sum(hessian.obs.cur) + (-0.001 / error.var.cur) # IG(0.001,0.001) prior

    # Laplace proposal mean and variance
    # Proposal log(sigma^2)_proposed ~ N(mu.prop, prop.var)
    # mu.prop = log.sigma2.cur - score / hessian (Newton step)
    # prop.var = -1 / hessian (must be positive)

    if (hessian.sigma >= 0 || !is.finite(hessian.sigma)) {
      # Hessian not negative definite, Laplace approximation invalid here
      # Fall back to small random walk centered at current value
      sigma2.prop.sd  <- 0.1
      log.sigma2.prop.mean <- log.sigma2.cur
    } else {
      sigma2.prop.sd <- sqrt(-1 / hessian.sigma)
      log.sigma2.prop.mean <- log.sigma2.cur - score.sigma / hessian.sigma
    }

    log.sigma2.prop <- rnorm(1, mean = log.sigma2.prop.mean, sd = sigma2.prop.sd)
    sigma2.prop <- exp(log.sigma2.prop)
    sigma.prop <- sqrt(sigma2.prop)

    # Reverse proposal (requirements: Laplace proposal is asymmetric)
    # At the proposed value, re calculate W12 and W14 to find the reverse proposal mean
    scores.error.prop <- pmin(pmax(qnorm(pnorm(resid.vec / sigma.prop)), -8), 8)
    scores.all.prop <- cbind(scores.endo, scores.exo, scores.error.prop)
    A.xi.prop   <- as.vector(scores.all.prop %*% A.d.row)

    score.obs.rev   <- (resid.vec / (2 * sigma.prop)) * A.xi.prop - 0.5 + resid.vec^2 / (2 * sigma2.prop)
    score.sigma.rev <- sum(score.obs.rev) + (-1.001 + 0.001 / sigma2.prop)

    hessian.obs.rev <- -(resid.vec / (4 * sigma.prop)) * A.xi.prop - A.dd * resid.vec^2 / (4 * sigma2.prop) - resid.vec^2 / (2 * sigma2.prop)
    hessian.sigma.rev <- sum(hessian.obs.rev) + (-0.001 / sigma2.prop)

    if (hessian.sigma.rev >= 0 || !is.finite(hessian.sigma.rev)) {
      sigma2.rev.sd <- 0.1
      log.sigma2.rev.mean  <- log.sigma2.prop
    } else {
      sigma2.rev.sd <- sqrt(-1 / hessian.sigma.rev)
      log.sigma2.rev.mean <- log.sigma2.prop - score.sigma.rev / hessian.sigma.rev
    }

    # MH acceptance step
    # Log-posterior at current and proposed sigma^2
    # Jacobian: d sigma^2 / d log(sigma^2) = sigma^2
    # log p(log sigma^2) = log p(sigma^2) + log(sigma^2)
    lp.sigma.cur  <- copulabayes_logpost_sigma2( error.var.cur, copula.cor.cur, scores.endo, scores.exo,
                                                 resid.cur, y, z, x, K, L) + log.sigma2.cur

    lp.sigma.prop <- copulabayes_logpost_sigma2(sigma2.prop, copula.cor.cur, scores.endo, scores.exo,
                                                resid.cur, y, z, x, K, L) + log.sigma2.prop

    # Log proposal densities (both directions needed for asymmetric proposal)
    log.q.prop.given.cur <- dnorm(log.sigma2.prop, mean = log.sigma2.prop.mean,
                                  sd   = sigma2.prop.sd, log = TRUE)

    log.q.cur.given.prop <- dnorm(log.sigma2.cur, mean = log.sigma2.rev.mean,
                                  sd   = sigma2.rev.sd,  log = TRUE)

    # MH ratio: [p(prop) * q(cur|prop)] / [p(cur) * q(prop|cur)]
    log.mh.sigma <- (lp.sigma.prop + log.q.cur.given.prop) - (lp.sigma.cur  + log.q.prop.given.cur)

    if (!is.finite(log.mh.sigma)) log.mh.sigma <- -Inf

    if (log(runif(1L)) < log.mh.sigma) {
      error.var.cur <- sigma2.prop
      n.mh.accepted.sigma <- n.mh.accepted.sigma + 1L
    }

    chain[i + 1L, col.intercept] <- intercept.cur
    chain[i + 1L, col.coef.endo] <- coef.endo.cur
    if (L > 0L) chain[i + 1L, col.coef.exo] <- coef.exo.cur
    chain[i + 1L, col.error.var]  <- error.var.cur

    # Same as RW
    resid.cur <- y - X.design %*% c(intercept.cur, coef.endo.cur, coef.exo.cur)

    scores.error <- pmin(pmax(qnorm(pnorm(as.vector(resid.cur) / sqrt(error.var.cur))), -8), 8)
    scores.all <- cbind(scores.endo, scores.exo, scores.error)

    # Step 4: Gibbs for copula covariance W (Appendix A W5)
    W.draw <- MCMCpack::riwish(N + cop.dim, diag(cop.dim) + crossprod(scores.all))
    cor.normaliser <- diag(1 / sqrt(diag(W.draw)))
    copula.cor.cur <- cor.normaliser %*% W.draw %*% cor.normaliser
    for (l in seq_len(L)) {
      copula.cor.cur[K + l, cop.dim] <- 0
      copula.cor.cur[cop.dim, K + l] <- 0
    }
    chain[i + 1L, col.copula.cor] <- copulabayes_matrix2vector(copula.cor.cur)

    # Step 5: Gibbs for Dirichlet masses (Appendix B W7)
    eps.draw <- mvtnorm::rmvnorm(N, mean = rep(0, cop.dim), sigma = copula.cor.cur)
    for (k in seq_len(K)) {
      masses.endo.list[[k]] <- copulabayes_drawlambda(
        pnorm(eps.draw[, k]), margin.endo.list[[k]]
      )
      scores.endo[, k] <- copulabayes_converter(
        masses.endo.list[[k]], margin.endo.list[[k]]
      )
    }
    for (l in seq_len(L)) {
      masses.exo.list[[l]] <- copulabayes_drawlambda(
        pnorm(eps.draw[, K + l]), margin.exo.list[[l]]
      )
      scores.exo[, l] <- copulabayes_converter(
        masses.exo.list[[l]], margin.exo.list[[l]]
      )
    }

    # Horseshoe hyperprior updates (eq. 5-6 from main paper)
    chain[i + 1L, col.var.alpha] <- 1 /rgamma(1L, shape = 0.501, rate = (intercept.cur^2 + 0.002) / 2)

    for (k in seq_len(K)) {
      chain[i + 1L, col.var.delta[k]] <- 1 / rgamma(1L, shape = 0.501, rate = (coef.endo.cur[k]^2 + 0.002) / 2)
    }

    if (L > 0L) {
      for (l in seq_len(L)) {
        chain[i + 1L, col.var.beta[l]] <- 1 /rgamma(1L, shape = 0.501, rate = (coef.exo.cur[l]^2 + 0.002) / 2)
      }
    }
  }

  attr(chain, "col.intercept") <- col.intercept
  attr(chain, "col.coef.endo") <- col.coef.endo
  attr(chain, "col.coef.exo")  <- col.coef.exo
  attr(chain, "col.error.var") <- col.error.var
  attr(chain, "col.copula.cor") <- col.copula.cor
  attr(chain, "K") <- K
  attr(chain, "L") <- L
  attr(chain, "cop.dim") <- cop.dim
  return(chain)
}

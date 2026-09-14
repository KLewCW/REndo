#' @importFrom stats lm residuals var runif rnorm rgamma qnorm pnorm
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom MCMCpack riwish rdirichlet

copulabayes_mcmc_iwls <- function(y,z,x,num.iterations,verbose){

  N <- length(y)
  K <- ncol(z)
  L <- ncol(x)

  copula.dim <- K + L + 1L

  # Full design matrix including intercept
  X.design <- cbind(1, z, x)   # N x (1+K+L)
  p <- ncol(X.design)

  #column indices into copula matrix (except error dimension)
  col.indices <- seq_len(copula.dim - 1L)

  ##this part is same as random walk
  margin.endo.list <- lapply(seq_len(K), function(k) copulabayes_margin(z[, k]))
  margin.exo.list <- lapply(seq_len(L), function(l) copulabayes_margin(x[,l]))

  #ols starting values

  mod.ols <- lm(y ~ X.design - 1) # minus 1 because X.design already includes the intercept
  intercept.cur <- coef(mod.ols)[1L]
  coef.endo.cur <- coef(mod.ols)[1L + seq_len(K)]

  if(L > 0L){
    coef.exo.cur <-  coef(mod.ols)[1L + K + seq_len(L)]
  } else{
    numeric(0L)
  }
  coef.cur.vec   <- c(intercept.cur, coef.endo.cur, coef.exo.cur)
  error.var.cur  <- var(residuals(mod.ols))
  copula.cor.cur <- diag(copula.dim)

  ## initial Dirichlet masses and normal scores

  masses.endo.list <- lapply(margin.endo.list, function(mg) {g <- rgamma(mg$n.unique, 1); g / sum(g)})
  masses.exo.list <- lapply(margin.exo.list, function(mg) {g <- rgamma(mg$n.unique, 1); g / sum(g)})

  scores.endo <- matrix(NA_real_, N, K)
  for(k in seq_len(K)){
    scores.endo[, k] <- copulabayes_converter(masses.endo.list[[k]], margin.endo.list[[k]])
  }

  scores.exo <- matrix(NA_real_, N,L)
  for(l in seq_len(L)){
    scores.exo[, l] <- copulabayes_converter(masses.exo.list[[l]], margin.exo.list[[l]])
  }

  # chain storage (same as RW)

  n.copula.cor <- copula.dim * (copula.dim - 1L) / 2L
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


  #Omega represent prior for endo block of W (K X K identity)

  Omega.cur <- diag(K)

  n.mh.accepted.coef <- 0L  # accepted proposals for coefficients
  n.mh.accepted.sigma <- 0L  # accepted proposals for sigma^2

  #MCMC loop

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
    coef.cur.vec <- c(intercept.cur, coef.endo.cur, coef.exo.cur)
    resid.cur  <- y - X.design %*% coef.cur.vec
    scores.error.cur <- pmin(pmax( qnorm(pnorm(as.vector(resid.cur) / sqrt(error.var.cur))), -8), 8)
    scores.all.cur <- cbind(scores.endo, scores.exo, scores.error.cur)

    # A = Sigma^{-1} -I
    Sinv <- chol2inv(chol(copula.cor.cur))
    A <- Sinv - diag(copula.dim)

    #Copula score contribution (fixed during coefficient update)
    #q_i = A[1:d-1, d]' * xi [ 1: d-1,i] for each obs. i

    q <- as.vector(scores.all.cur[, col.indices, drop = FALSE] %*% A[col.indices, copula.dim])
    aee <- A[copula.dim, copula.dim] #scalar A[d,d] including the error dimension

    #First block: IWLS for regression coeff (Haschka 2026, Appendix C. W11 and W13)
    #Working weights = Sinv[d,d]/sigma^2
    # Proposal covariance is Sigma_prop = (Sinv[d,d]/s2 * X'X)^{-1}
    # Proposal meanis mu = cf + Sigma_prop * X' * nu

    e.cur  <- as.vector(resid.cur)
    sg.cur <- sqrt(error.var.cur)

    #Score per observation (last element of W11)
    nu.cur <- (q + aee * e.cur / sg.cur) / sg.cur + e.cur / error.var.cur

    # Proposal covariance through inverse Cholesky
    working.weight <- Sinv[copula.dim, copula.dim] / error.var.cur
    R.prop <- chol(chol2inv(chol(working.weight * crossprod(X.design)))) # upper Cholesky of Sigma_prop

    # Proposal mean
    mu.fwd <- coef.cur.vec + as.vector(crossprod(R.prop, R.prop %*% crossprod(X.design, nu.cur)))

    # Draw proposal using Cholesky
    coef.prop.vec <- as.vector(mu.fwd + crossprod(R.prop, rnorm(p)))

    # Reverse proposal mean (at the proposed coefficient values)
    e.prop  <- as.vector(y - X.design %*% coef.prop.vec)
    nu.prop <- (q + aee * e.prop / sg.cur) / sg.cur + e.prop / error.var.cur
    mu.rev  <- coef.prop.vec + as.vector(crossprod(R.prop, R.prop %*% crossprod(X.design, nu.prop)))

    logpost.coef.cur  <- copulabayes_logpost_coef(
      coef.cur.vec,  error.var.cur, copula.cor.cur,
      scores.endo, scores.exo, var.alpha.cur, var.delta.cur, var.beta.cur,
      y, z, x, K, L)

    logpost.coef.prop <- copulabayes_logpost_coef(
      coef.prop.vec, error.var.cur, copula.cor.cur,
      scores.endo, scores.exo, var.alpha.cur, var.delta.cur, var.beta.cur,
      y, z, x, K, L)

    # Log proposal densities using Cholesky log
    log.q.fwd <- copulabayes_dmvn_chol(coef.prop.vec, mu.fwd, R.prop)
    log.q.rev <- copulabayes_dmvn_chol(coef.cur.vec,  mu.rev, R.prop)

    # Metropolis-Hastings acceptance ratio
    log.mh.coef <- (logpost.coef.prop + log.q.rev) - (logpost.coef.cur + log.q.fwd)
    if (!is.finite(log.mh.coef)) log.mh.coef <- -Inf

    if (log(runif(1L)) < log.mh.coef) {
      coef.cur.vec <- coef.prop.vec
      n.mh.accepted.coef <- n.mh.accepted.coef + 1L
    }

    intercept.cur <- coef.cur.vec[1L]
    coef.endo.cur <- coef.cur.vec[seq(2L, 1L + K)]
    coef.exo.cur <- if (L > 0L) coef.cur.vec[seq(2L + K, 1L + K + L)] else numeric(0L)

    chain[i + 1L, col.intercept]  <- intercept.cur
    chain[i + 1L, col.coef.endo]  <- coef.endo.cur
    if (L > 0L) chain[i + 1L, col.coef.exo] <- coef.exo.cur

    # Block 2: Laplace proposal for log(sigma^2)
    # Haschka (2026) Appendix C, W12 (first derivative) and W14 (second derivative)

    # The Laplace proposal uses a second-order Taylor expansion of the
    # log-posterior around log(sigma^2).

    #Note: W12 and W14 are matrix derivatives of dimension (L+K+1) x 1 and
    # (L+K+1) x (L+K+1) respectively.
    # We needthe scalar element:
    # For W12 the (L+K+1)-th element = score  for log(sigma^2) per obs
    # For W14 the (L+K+1, L+K+1) element = Hessian for log(sigma^2) per obs

    # W12 last element per observation i (derived from log L_i, from W12):
    #   s_i = (e_i / (2*sigma)) * (A[d,:] %*% xi_i) - 1/2 + e_i^2 / (2*sigma^2)

    # W14 last diagonal element per observation i (derived from W12, from W14):
    # h_i = -(e_i / (4*sigma)) * (A[d,:] %*% xi_i) - A[d,d] * e_i^2 / (4*sigma^2)- e_i^2 / (2*sigma^2)

    # Note: the paper gives d log L_i (likelihood only).
    # prior contributions are added separately (IG(0.001, 0.001) on sigma^2).


    log.sigma2.cur <- log(error.var.cur)

    #recalculating the residuals with updated coefficients
    e.cur  <- as.vector(y - X.design %*% coef.cur.vec)
    xe.cur <- e.cur / sqrt(error.var.cur)   # standardised residuals = xi_e
    qxe <- sum(q * xe.cur)
    q2  <- sum(xe.cur^2)

    # last element of W12 + prior contribution
    # Prior: IG(0.001, 0.001) → contributes -a + b/sigma^2 after Jacobian
    score.sigma <- 0.5 * qxe + 0.5 * (aee + 1) * q2 - 0.5 * length(e.cur) - 0.001 + 0.001 / error.var.cur

    # (W14, last diagonal element) + prior contribution
    hessian.sigma <- -0.25 * qxe - 0.5 * (aee + 1) * q2 - 0.001 / error.var.cur

    # Laplace proposal variance: -1/Hessian (must be positive)
    if (hessian.sigma < -1e-12) {
      sigma.prop.var <- -1 / hessian.sigma
    } else {
      sigma.prop.var <- 1.0   # fallback if Hessian not negative definite
    }
    log.sigma2.prop.mean <- log.sigma2.cur + sigma.prop.var * score.sigma
    log.sigma2.prop <- rnorm(1, log.sigma2.prop.mean, sqrt(sigma.prop.var))
    sigma2.prop <- exp(log.sigma2.prop)

    # Reverse proposal (at proposed sigma^2)
    xe.prop <- e.cur / sqrt(sigma2.prop)
    qxe.rev <- sum(q * xe.prop)
    q2.rev  <- sum(xe.prop^2)
    score.sigma.rev   <- 0.5 * qxe.rev + 0.5 * (aee + 1) * q2.rev - 0.5 * length(e.cur) - 0.001 + 0.001 / sigma2.prop
    hessian.sigma.rev <- -0.25 * qxe.rev - 0.5 * (aee + 1) * q2.rev -0.001 / sigma2.prop

    sigma.rev.var <- if (hessian.sigma.rev < -1e-12) {
      -1 / hessian.sigma.rev
    } else {
      sigma.prop.var
    }
    log.sigma2.rev.mean <- log.sigma2.prop + sigma.rev.var * score.sigma.rev

    # MH ratio for log(sigma^2)
    # Log-posterior evaluated at sigma^2 directly
    lp.sigma.cur  <- copulabayes_logpost_sigma2(
      error.var.cur, copula.cor.cur, scores.endo, scores.exo,e.cur, y, z, x, K, L
    )

    lp.sigma.prop <- copulabayes_logpost_sigma2(
      sigma2.prop, copula.cor.cur, scores.endo, scores.exo,e.cur, y, z, x, K, L
    )

    # Jacobian: d(sigma^2)/d(log sigma^2) = sigma^2 (add log.sigma2.prop - log.sigma2.cur)
    log.mh.sigma <- (lp.sigma.prop - lp.sigma.cur) +
      dnorm(log.sigma2.cur, log.sigma2.rev.mean,sqrt(sigma.rev.var), log = TRUE)
    - dnorm(log.sigma2.prop, log.sigma2.prop.mean, sqrt(sigma.prop.var), log = TRUE)
    + (log.sigma2.prop - log.sigma2.cur)   #Jacobian

    if (!is.finite(log.mh.sigma)) log.mh.sigma <- -Inf

    if (log(runif(1L)) < log.mh.sigma) {
      error.var.cur <- sigma2.prop
      n.mh.accepted.sigma <- n.mh.accepted.sigma + 1L
    }
    chain[i + 1L, col.error.var] <- error.var.cur

    # Block 3: Gibbs for copula covariance W (Appendix A W5)
    # Recomputing xi_e with updated coefficients and sigma^2

    xi.e <- as.vector( y - X.design %*% coef.cur.vec) / sqrt(error.var.cur)

    dw.result <- copulabayes_draw_W(scores.endo, scores.exo, xi.e, Omega.cur, K, L)
    copula.cor.cur <- dw.result$Sigma
    Omega.cur <- dw.result$Omega
    chain[i + 1L, col.copula.cor] <- copulabayes_matrix2vector(copula.cor.cur)

    # Block 4: Gibbs for Dirichlet masses (Appendix B W7)
    #using Cholesky

    Rs <- chol(copula.cor.cur)
    U.mat <- pnorm(matrix(rnorm(N * copula.dim), N, copula.dim) %*% Rs)
    for (k in seq_len(K)) {
      masses.endo.list[[k]] <- copulabayes_drawlambda(U.mat[, k], margin.endo.list[[k]])
      scores.endo[, k]  <- copulabayes_converter(masses.endo.list[[k]], margin.endo.list[[k]])
    }
    for (l in seq_len(L)) {
      masses.exo.list[[l]] <- copulabayes_drawlambda(U.mat[, K + l], margin.exo.list[[l]])
      scores.exo[, l] <- copulabayes_converter(masses.exo.list[[l]], margin.exo.list[[l]])
    }

    # Horseshoe hyperprior updates
    chain[i + 1L, col.var.alpha] <- 1 / rgamma(1L, shape = 0.501, rate = (intercept.cur^2 + 0.002) / 2)
    for (k in seq_len(K))
      chain[i + 1L, col.var.delta[k]] <- 1 / rgamma(1L, shape = 0.501, rate = (coef.endo.cur[k]^2 + 0.002) / 2)
    if (L > 0L)
      for (l in seq_len(L))
        chain[i + 1L, col.var.beta[l]] <- 1 / rgamma(1L, shape = 0.501, rate = (coef.exo.cur[l]^2 + 0.002) / 2)
  }

  attr(chain, "col.intercept") <- col.intercept
  attr(chain, "col.coef.endo") <- col.coef.endo
  attr(chain, "col.coef.exo")  <- col.coef.exo
  attr(chain, "col.error.var") <- col.error.var
  attr(chain, "col.copula.cor") <- col.copula.cor
  attr(chain, "K") <- K
  attr(chain, "L") <- L
  attr(chain, "copula.dim") <- copula.dim
  return(chain)
}

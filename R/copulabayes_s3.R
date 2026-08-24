#' @importFrom coda as.mcmc effectiveSize geweke.diag
#' @importFrom stats sd quantile
new_rendo_copulabayes <- function(
  call,
  F.formula,
  chain, #post burnin thinned chain (all columns)
  chain.struct, #structural parameter draws: alpha, delta, beta, sigma2
  chain.rho, #endo-error copula correlation draws: rho_1,...,rho_K
  post.mean,
  post.sd,
  post.lo,
  post.hi,
  fitted.values,
  residuals,
  names.endo.regs,
  n.iterations,
  burnin,
  thin,
  n.draws
) {
  structure(
    list(
      call = call,
      F.formula = F.formula,
      chain = chain,
      chain.struct = chain.struct,
      chain.rho = chain.rho,
      post.mean = post.mean,
      post.sd = post.sd,
      post.lo = post.lo,
      post.hi = post.hi,
      fitted.values = fitted.values,
      residuals = residuals,
      names.endo.regs = names.endo.regs,
      n.iterations = n.iterations,
      burnin = burnin,
      thin = thin,
      n.draws = n.draws
    ),
    class = "rendo.copula.bayes"
  )
}


#' @export
print.rendo.copula.bayes <- function(x, ...) {
  cat("Bayesian Gaussian Copula Endogeneity Correction \n")
  cat(rep("-", 60), "\n", sep = "")

  cat("\nCall:\n")
  print(x$call)

  cat("\nPosterior means (structural parameters):\n")
  print(round(x$post.mean, 4))

  cat(
    "\n",
    x$n.draws,
    " draws retained from ",
    x$n.iterations,
    " iterations (burnin = ",
    x$burnin,
    ", thin = ",
    x$thin,
    ").\n",
    sep = ""
  )
  cat("Run summary() for full posterior summaries and convergence diagnostics.\n")
  invisible(x)
}


#' @export
summary.rendo.copula.bayes <- function(object, ...) {
  #summary method: full posterior summaries and convergence diagnostics

  #Structural parameters
  mcmc.struct <- coda::as.mcmc(object$chain.struct)
  ess.struct <- round(coda::effectiveSize(mcmc.struct)) #Effective sample size
  gw.struct <- round(coda::geweke.diag(mcmc.struct)$z, 2) #Geweke convergence diagnostic
  #Geweke: is a test for nonconvergence of a MCMC chain: diffenrece of means test that compares the mean
  #of early in the chain to the mean late in the chain.

  table.struct <- cbind(
    Mean = round(object$post.mean, 4),
    SD = round(object$post.sd, 4),
    `2.5%` = round(object$post.lo, 4),
    `97.5%` = round(object$post.hi, 4),
    ESS = ess.struct,
    `Geweke z` = gw.struct
  )

  #Endo-error copula correlations
  rho.mean <- colMeans(object$chain.rho)
  rho.sd <- apply(object$chain.rho, 2, sd)
  rho.low <- apply(object$chain.rho, 2, quantile, probs = 0.025)
  rho.high <- apply(object$chain.rho, 2, quantile, probs = 0.975)

  mcmc.rho <- coda::as.mcmc(object$chain.rho)
  ess.rho <- round(coda::effectiveSize(mcmc.rho))
  gw.rho <- round(coda::geweke.diag(mcmc.rho)$z, 2)

  table.rho <- cbind(
    Mean = round(rho.mean, 4),
    SD = round(rho.sd, 4),
    `2.5%` = round(rho.low, 4),
    `97.5%` = round(rho.high, 4),
    ESS = ess.rho,
    `Geweke z` = gw.rho
  )

  output <- list(
    call = object$call,
    table.struct = table.struct,
    table.rho = table.rho,
    n.draws = object$n.draws,
    n.iterations = object$n.iterations,
    burnin = object$burnin,
    thin = object$thin
  )
  class(output) <- "summary.rendo.copula.bayes"
  output
}

#' @export
print.summary.rendo.copula.bayes <- function(x, ...) {
  cat("Bayesian Gaussian Copula Endogeneity Correction\n")
  cat(rep("-", 65), "\n", sep = "")

  cat("\nCall:\n")
  print(x$call)

  cat(
    "\nMCMC settings:\n",
    "  Iterations: ",
    x$n.iterations,
    "  |  Burn-in: ",
    x$burnin,
    "  |  Thinning: ",
    x$thin,
    "  |  Draws retained: ",
    x$n.draws,
    "\n",
    sep = ""
  )

  # Structural parameters
  cat("\nStructural parameters:\n")
  print(x$table.struct, quote = FALSE)

  #Endogeneity strength
  cat("\nEndogeneity strength (copula correlation with structural error):\n")
  cat("rho > 0: positive endogeneity bias; rho < 0: negative bias.\n")
  cat("if the 95% credible interval excludes 0, endogeneity is significant.\n\n")
  print(x$table.rho, quote = FALSE)

  # Convergence guidance
  cat("\nConvergence guidance:\n")
  cat("ESS: effective sample size after autocorrelation correction.\n") #Recommended that ESS is around
  #400 (number of chains * 100)
  cat("Geweke: z-score comparing first 10% to last 50% of the chain.\n") #If the means are significantly different from each other,
  #then this is evidence that the chain has not converged.

  invisible(x)
}

# plot method
#' @export
#' @importFrom coda as.mcmc
plot.rendo.copula.bayes <- function(x, which = c("structural", "rho", "both"), ...) {
  which <- match.arg(
    which,
    choices = c("structural", "rho", "both"),
    several.ok = FALSE
  )

  if (which %in% c("structural", "both")) {
    cat("Plotting structural parameters (trace + density)...\n")
    plot(coda::as.mcmc(x$chain.struct), main = "Structural parameters", ...)
  }
  if (which %in% c("rho", "both")) {
    cat("Plotting endogeneity correlations (trace + density)...\n")
    plot(coda::as.mcmc(x$chain.rho), main = "Endogeneity (rho)", ...)
  }
  invisible(x)
}


#' @export
coef.rendo.copula.bayes <- function(object, ...) {
  object$post.mean
}

#' @export
residuals.rendo.copula.bayes <- function(object, ...) {
  object$residuals
}

#' @export
fitted.rendo.copula.bayes <- function(object, ...) {
  object$fitted.values
}

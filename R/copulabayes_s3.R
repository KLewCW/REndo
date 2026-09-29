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
  n.obs,
  n.iterations,
  burnin,
  thin,
  n.draws,
  method
) {
  return(structure(
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
      n.obs = n.obs,
      n.iterations = n.iterations,
      burnin = burnin,
      thin = thin,
      n.draws = n.draws,
      method = method
    ),
    class = "rendo.copula.bayes"
  ))
}


#' @export
print.rendo.copula.bayes <- function(x, ...) {
  cat("Bayesian Gaussian Copula Endogeneity Correction \n")
  cat(rep("-", 60), "\n", sep = "")

  cat("\nCall:\n")
  print(x$call)

  cat("\nSampler:", if (x$method == "RW")
    "Adaptive random walk MH"
    else "IWLS", "\n")

  cat("Endogenous regressor(s):", paste(x$names.endo.regs, collapse = ", "), "\n")
  cat("Observations:", x$n.obs, "\n")

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
  which <- match.arg(which)

  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par), add = TRUE)

  clean.param.names <- function(nms) {
    nms <- gsub("_endo$", " (endogenous)", nms)
    nms <- gsub("_exo$",  " (exogenous)",  nms)
    nms <- gsub("^sigma2$", "sigma^2 (error variance)", nms)
    return(nms)
  }

  clean.rho.names <- function(nms) {
    gsub("^rho_", "rho: endogeneity of ", nms)
  }

  plot_mcmc_label <- function(chain.mat, display.names, section.label){
    n.params <- ncol(chain.mat)

    info.line <- paste0("Sampler:", x$method, "| Draws retained:", x$n.draws,
                        " | Burnin:", x$burnin, "| Thin:", x$thin)

    par(mfrow = c(n.params, 2),
        mar = c(4, 4, 3, 1),
        oma = c(0, 0, 3, 0))

    mcmc.obj <- coda::as.mcmc(chain.mat)

    for (p in seq_len(n.params)){
      pname <- display.names[p]

      #left panel: the trace plot
      coda::traceplot(mcmc.obj[, p],main  = paste0("Trace: ", pname),
                      ylab  = pname, xlab  = "Retained draw", col   = "steelblue", ...)

      #right panel: density plot

      coda::densplot(mcmc.obj[, p], main  = paste0("Posterior: ", pname),
                     xlab  = pname,col   = "steelblue", ... )

      #title
      mtext( text  = paste0(section.label, "\n", info.line), outer = TRUE,
             cex   = 0.85, font  = 2, line  = 1)
    }
  }

    if (which %in% c("structural", "both")) {
      chain.plot  <- x$chain.struct
      param.names <- clean.param.names(colnames(chain.plot))
      plot_mcmc_label(chain.plot, param.names, "Structural parameters")
    }

    if (which %in% c("rho", "both")) {
      chain.rho.plot <- x$chain.rho
      rho.names <- clean.rho.names(colnames(chain.rho.plot))
      plot_mcmc_label( chain.rho.plot, rho.names,
        paste0( "Endogeneity strength (rho)", "\nrho > 0: OLS overestimates  |  rho < 0: OLS underestimates",
          "\nCredible interval excluding 0 confirms endogeneity"
        )
      )
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

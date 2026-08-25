#' Bayesian Gaussian Copula Endogeneity Correction
#'
#' @description
#' Fits the Bayesian in a one-step Gaussian copula endogeneity correction.
#' The method treats the marginal distributions of all regressors, the full copula
#' correlation matrix, and the regression coefficients as random variables and
#' samples them jointly via Markov chain Monte Carlo (MCMC). This avoids the use
#' of plug-in estimates for CDFs and correlations which allows precise finite-sample
#' inference without having to rely on asymptotic assumptions or any bootstrap procedure.
#'
#' @template template_param_formuladataverbose
#'
#' @param num.iterations a positive integer giving the total number of MCMC iterations.
#' From Haschka (2025), it is suggested to use 102000 as number of iterations, with a burnin of 2000
#' and a thinning by 100, in order to have 1000 posterior draws. For a 'quick' exploratory run, a smaller
#' value can be used for the number of iteration (e.g. 12000 with a thinning of 10). However, always check
#' convergence diagnostics first.
#'
#' @param burnin a non-negative integer giving the number of initial MCMC draws to discard as burn-in.
#' Default is 2000. Burn-in draws show the chain's dependence on the starting values and should not be used
#' for inference.
#'
#' @param thin a positive integer giving the thinning interval. Only keep the \code{thin} th draw after burn-in.
#' Default is 100 (for an iteration of 102000 and burning of 2000 in order to have 1000 posterior draws).
#' Thinning is used to reduce autocorrelation between consecutive draws.
#'
#'
#' @details
#' ## Model
#'
#' Consider the following linear regression model
#' \deqn{Y_i = \alpha + X_i' \beta + Z_i' \delta + \varepsilon_i}
#'
#' where
#' \eqn{i = 1, \dots, N} indexes observations,
#' \eqn{Z_i} is a \eqn{(K \times 1)} vector of continuous
#' endogenous regressors correlated with \eqn{\varepsilon_i},
#' \eqn{X_i} is an \eqn{(L \times 1)} vector of exogenous regressors
#' uncorrelated with \eqn{\varepsilon_i},
#' \eqn{\varepsilon_i \sim N(0, \sigma^2)}.
#'
#'
#' ## Methodology
#'
#' The method jointly samples all unknowns in one step via MCMC
#' (from Appendix D online, see algorithm 1). Each iteration would go through the
#' following steps:
#'
#' ##note: for the regression coeff and error var, since I tried to modify it and use
#' other package, it is slightly different from original paper Haschka (2015) where he
#' uses IWLS and Laplace. It is no more this but just a random walk MH, should this be changed ???
#'
#' \enumerate{
#'   \item {Regression coefficients} \eqn{(\alpha, \beta, \delta)} with
#'        'horseshoe'-type hierarchical shrinkage priors. They are updated by
#'         Metropolis-Hastings (MH) algorithm with an iteratively weighted least
#'         squares (IWLS) proposal (Appendix C, W11 and W13). The working weight
#'         \eqn{M_I = 2/\sigma^2} per observation (W13) gives proposal covariance
#'         \eqn{\Sigma_{\text{prop}} = (\sigma^2/2)(X'X)^{-1}}.
#'   \item {Error variance} \eqn{\sigma^2 \sim \text{IG}(0.001, 0.001)} (inverse Gamma, Haschka eq. 7).
#'         They are updated by Metropolis-Hastings with a Laplace (second-order
#'         Taylor) proposal for \eqn{\log \sigma^2} (Appendix C, W12/W14).
#'   \item {Copula correlation matrix} \eqn{W \sim W^{-1}(I, K+L+1)} ( Haschka 2025, eq. 8).
#'         This is updated by Gibbs sampling via the inverse Wishart full conditional (Appendix A, W5):
#'         \eqn{W | \cdot \sim W^{-1}(\sum_i \xi_i \xi_i' + I,\, N + K + L + 1)}.
#'         Converted to a correlation matrix \eqn{\Sigma} via Cholesky factorisation.
#'         Correlations between exogenous regressors and the structural error are
#'         set to zero, enforcing exogeneity of \eqn{X}.
#'   \item {Marginal distributions} of all regressors modelled nonparametrically via
#'         Dirichlet probability masses \eqn{\lambda_\varpi \sim \text{Dir}(1,\ldots,1)} (from Haschka 2025 eq. 9).
#'         This is updated by Gibbs sampling (Appendix B, W7):
#'         \eqn{\lambda_\varpi | v_\varpi \sim \text{Dir}(m_\varpi, \,1 + n_1,\ldots, 1 + n_{m_\varpi})}.
#'         The marginal CDF is never fixed but re-estimated at every iteration.
#' }
#'
#' ### Convergence diagnostics
#' The acceptance rate is printed every 500 iterations when \code{verbose = TRUE}.
#' A rate between 20\% and 40\% shows that the proposal scale is well-tuned (the sampler
#' adapts automatically every 50 iterations via Robbins-Monro scaling). This convergence
#' diagnostics should be verified before any conclusion.
#'
#' ### Identification requirements
#' Endogenous regressors \eqn{Z} needs to be non-normal, otherwise the chain will fail to converge.
#' Only continuous endogenous regressors are supported.
#'
#'
#' ## Formula interface
#' Formula follows a two-part notation:
#'
#' \preformatted{y ~ X + Z | continuous(Z) #one endogenous regressor}
#' \preformatted{y ~ X + Z1 + Z2 | continuous(Z1) + continuous(Z2) #two endogenous regressors}
#'
#'
#' @template template_references_parkgupta2012
#'
#'
#' @references
#' Haschka, R. E (2025) Bayesian Inference for Joint Estimation Models Using Copulas
#' to Handle Endogenous Regressors.
#' \emph{Oxford Bulletin of Economics and Statistics} 88(3), 519--534
#' \doi{10.1111/obes.70023}
#'
#'
#' @examples
#' # Example: Bayesian Gaussian copula endogeneity correction
#' # based on Section 4.1, eq. (12) to (16)
#' # N = 1000
#' # True parameters are: alpha=2, beta=6 (x), delta=-4 (z), sigma^2 =5
#' # z ~ lognormal(0,1): non-normal, required for identification as per section 2
#' # x ~ N(0,1) for exogenous regressor, correlated with z (rho_xz=0.3)
#' # Endogeneity strength is rho_ze = 0.7
#' # There is only one endogenous regressor(z_i) and one exogenous regressor (x_i) correlated
#' #with each other
#' #------------------------------------------------------------------------
#' data("dataCopulaBayes")
#' res_bayes <- copulaBayes(
#'   y ~ x + z | continuous(z),
#'   data           = dataCopulaBayes,
#'   num.iterations = 12000,
#'   burnin         = 2000,
#'   thin           = 10,
#'   verbose        = TRUE
#' )
#' summary(res_bayes)
#' plot(res_bayes, which = "both")
#'
#'
#' @md
#' @export
#'
#' @importFrom stats coef model.frame model.matrix model.response formula sd quantile
#' @importFrom Formula as.Formula
copulaBayes <- function(
  formula,
  data,
  num.iterations = 102000, #paper default but for quick testing 12000 iterations with 10 thin ?
  burnin = 2000,
  thin = 100,
  verbose = TRUE
) {
  cl <- match.call()

  check_err_msg(checkinput_copulashared_data_basics(data))
  check_err_msg(checkinput_copulashared_formula(formula))
  # check_err_msg(checkinput_copulabayes_formula_data(formula = formula, data = data))
  check_err_msg(checkinput_copulabayes_numiterations_burnin_thin(
    num.iterations = num.iterations,
    burnin = burnin,
    thin = thin
  ))
  check_err_msg(checkinput_copulashared_verbose(verbose))

  F.formula <- Formula::as.Formula(formula)

  names.endo.regs <- formula_readout_special(
    F.formula = F.formula,
    name.special = "continuous",
    from.rhs = 2,
    params.as.chars.only = TRUE
  )

  if (length(names.endo.regs) == 0) {
    stop(
      "No endogenous regressors found. Declare at least one using ",
      "continuous() in the second part of the formula, ",
      "e.g. y ~ X + Z | continuous(Z).",
      call. = FALSE
    )
  }

  f.main <- formula(F.formula, lhs = 1, rhs = 1)
  mf <- model.frame(f.main, data = data)
  y <- model.response(mf)

  X.main <- model.matrix(f.main, data = mf)
  X.main <- X.main[, colnames(X.main) != "(Intercept)", drop = FALSE]

  is.endo <- colnames(X.main) %in% names.endo.regs
  z <- X.main[, is.endo, drop = FALSE] # N x K endogenous
  x <- X.main[, !is.endo, drop = FALSE] # N x L exogenous (N x 0 if L=0)

  if (is.null(colnames(z))) {
    colnames(z) <- paste0("z", seq_len(ncol(z)))
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }

  if (ncol(z) < length(names.endo.regs)) {
    stop(
      "Could not match all declared endogenous regressors in the ",
      "design matrix. Check that continuous() arguments match ",
      "variable names exactly as they appear in the structural model.",
      call. = FALSE
    )
  }

  # MCMC
  if (verbose) {
    # fmt: skip
    message("Fitting Bayesian copula model for ",ncol(z)," endogenous and ",ncol(x)," exogenous regressor(s).")
    # fmt: skip
    message("Running a total of ",num.iterations," MCMC iterations (burnin = ",burnin,", thin = ",thin,").")
  }

  chain.full <- copulabayes_mcmc(
    y = y,
    z = z,
    x = x,
    num.iterations = num.iterations,
    verbose = verbose
  )

  # Burn-in and thinning
  idx.keep <- seq(burnin + 2L, num.iterations + 1L, by = thin)
  chain <- chain.full[idx.keep, , drop = FALSE]

  col.alpha <- attr(chain.full, "col.alpha")
  col.delta <- attr(chain.full, "col.delta")
  col.beta <- attr(chain.full, "col.beta")
  col.rho <- attr(chain.full, "col.rho")
  col.sigma2 <- attr(chain.full, "col.sigma2")
  K <- attr(chain.full, "K")
  L <- attr(chain.full, "L")
  cop.dim <- attr(chain.full, "cop.dim")

  coef.names <- c(
    "(Intercept)",
    paste0(colnames(z), "_endo"),
    if (L > 0) paste0(colnames(x), "_exo") else character(0),
    "sigma2"
  )
  structure.cols <- c(col.alpha, col.delta, col.beta, col.sigma2)

  chain.struct <- chain[, structure.cols, drop = FALSE]
  colnames(chain.struct) <- coef.names

  #endogenous and error copula correction draws

  idx.rho.endo <- (cop.dim - 1L) * (cop.dim - 2L) / 2L + seq_len(K)
  chain.rho <- chain[, col.rho[idx.rho.endo], drop = FALSE]
  colnames(chain.rho) <- paste0("rho_", colnames(z))

  # Posterior summaries
  post.mean <- colMeans(chain.struct)
  post.sd <- apply(chain.struct, 2, sd)
  post.lo <- apply(chain.struct, 2, quantile, probs = 0.025)
  post.hi <- apply(chain.struct, 2, quantile, probs = 0.975)

  # Structural fitted values and residuals
  alpha.pm <- post.mean["(Intercept)"]
  delta.pm <- post.mean[paste0(colnames(z), "_endo")]
  beta.pm <- if (L > 0) post.mean[paste0(colnames(x), "_exo")] else numeric(0)

  fitted.values <- as.vector(alpha.pm + z %*% delta.pm + x %*% beta.pm)
  residuals <- y - fitted.values

  return(new_rendo_copulabayes(
    call = cl,
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
    n.iterations = num.iterations,
    burnin = burnin,
    thin = thin,
    n.draws = nrow(chain)
  ))
}

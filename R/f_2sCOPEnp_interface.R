#' Two-Stage Nonparametric Copula Control Function Approach (2sCOPEnp)
#'
#' @description
#' Fitting the two-stage nonparametric copula control function (2sCOPEnp)
#' estimator of Hu et al. (2025) to address endogeneity without external
#' instruments.
#'
#' The 2sCOPEnp method generalises existing copula correction methods by
#' replacing the marginal CDF transformation with a nonparametric kernel
#' estimate of the \emph{conditional} CDF of each endogenous regressor given
#' the exogenous regressors, \eqn{\hat{F}(P_k | X)}. This eliminates the need
#' for a first-stage regression (as in \code{\link{copula2sCOPE}}) and relaxes
#' the Gaussian copula assumption on the regressor-error dependence structure.
#' 2sCOPEnp also supports \emph{discrete} endogenous regressors by leveraging
#' the exogenous control variables to smooth the discrete conditional CDF.
#'
#' @template template_param_formuladataverbose
#' @param data A data.frame containing the data of all parts specified in the formula
#' parameter. Variables must be correctly classed (\code{numeric}, \code{factor}, or
#'  \code{ordered}) otherwise it will silently yield wrong results!
#' @template template_param_numboots
#' @param npcdistbw.args A named list of arguments which are passed to
#' \code{\link[np]{npcdistbw}}. To tweak the bandwidth selection for the kernel
#' conditional CDF. Runs with the method's defaults if not specified.
#' See \code{\link[np]{npcdistbw}} for valid inputs.
#' @param bws A named list of estimated \code{condbandwidth} objects (outputs of
#' \code{np::npcdistbw()}) used in place of estimating the bandwidth from \code{data}.
#' If supplied, bandwidth estimation is skipped entirely. See Details.
#'
#'
#' @details
#' \subsection{Model}{
#'
#' Consider the structural regression model with \eqn{K} endogenous regressors:
#'
#' \deqn{Y_i = \mu + \sum_{k=1}^{K} P_{i,k} \alpha_k + X_i' \beta + \varepsilon_i}
#'
#' where \eqn{i=1, \ldots, n } indexes observations,
#' \eqn{Y_i} is the dependent variable,
#' \eqn{P_{i,k}} are endogenous regressors (continuous or discrete) that may be correlated
#' with the structural error \eqn{\varepsilon_i},
#' \eqn{X_i} is a vector of exogenous regressors uncorrelated with \eqn{\varepsilon_i}, and
#' \eqn{\mu, \alpha_k, \beta} are structural model paramters.
#' }
#'
#' \subsection{Methodology}{
#'
#' The estimation proceeds in 2 stages (Hu et al. 2025)
#'
#' \emph{Stage 1}: For each endogenous regressor \eqn{P_k}:
#' \enumerate{
#'   \item Estimate the conditional CDF \eqn{\hat{F}(P_k | X)} nonparametrically
#'         using kernel methods (Li and Racine 2008), conditioning on all
#'         exogenous regressors \eqn{X}. A generalised product kernel handles
#'         mixed continuous and discrete \eqn{X} automatically (Hayfield and Racine 2008).
#'   \item Compute the copula correction term
#'         \eqn{C_{i,pk} = \Phi^{-1}(\hat{F}(P_{i,k} | X_i))}
#'         (Equation 21, Hu et al. 2025) where \eqn{\Phi^{-1}} is the
#'         standard normal quantile function.
#'   \item Repeat for all \eqn{K} endogenous regressors.
#' }
#'
#' \emph{Stage 2}: Add \eqn{\{C_{i,pk}\}, k =1, \ldots, K} to the structural model as
#' generated regressors and estimate using OLS (Equation 20):
#' \deqn{Y_i = \mu + \sum_{k=1}^{K} P_{i,k}\alpha_k + X_i'\beta +
#'       \sum_{k=1}^{K} C_{i,pk}\gamma_k + \xi_i}
#' where \eqn{\gamma_k} is the coefficient of the copula correction term \eqn{C_{i, pk}} and
#' \eqn{\xi_i} is the new error term.
#' }
#'
#' \subsection{Parameter \code{bws}}{
#' If given, the bandwidth estimation is skipped for all endogenous regressors, and
#' the bandwidths objects in \code{bws} are passed directly to \code{np::npcdist}
#' instead.
#' They do not serve as starting point for bandwidth estimation but replace it entirely.
#' For each endogenous term in \code{formula}, \code{bws} must contain a
#' \code{condbandwidth} object for the conditional CDF \code{(single endo) ~ (all exo)}.
#' \code{bws} must be a named list, with each each element named after the corresponding
#' endogenous terms.
#' The bandwidths must be estimated using the default \code{np::npcdistbw(xdat=,ydat=)}
#' interface, not the formula interface because this will break downstream usage in
#' \code{np::npcdist()}.
#' }
#'
#' \subsection{Parameter \code{formula}}{
#'
#' The \code{formula} argument follows a two part notation separated by \code{|}.
#' The first part specifies the structural model (e.g \code{y ~ X + P}).
#' The second part lists the endogenous regressors:
#'
#' \preformatted{y ~ X + P | P                       # endogenous P}
#' \preformatted{y ~ X + P1 + log(P2) | P1 + log(P2) # multiple endogenous regressors}
#'
#' At least one exogenous regressor must be present in the model for the
#' nonparametric conditional CDF estimation to be feasible.
#' }
#'
#' If the bootstrap standard errors of the endogenous regressor coefficients are more than
#' 6 times larger than the corresponding OLS standard errors, this may indicate near-multicollinearity
#' between the copula correction term and the original regressors, which may suggest a potential
#' identification issue. (Hu et al. 2025, section 3.5)
#'
#' \subsection{Bootstrap inference}{
#' Standard errors are obtained by resampling the data with replacement and
#' re-running the full two-stage estimation on each resample, including the
#' nonparametric bandwidth selection in Stage 1.
#' Degenerate bootstrap samples are automatically discarded and redrawn.
#' The percentage of discarded samples is reported as a warning if non-zero.
#' }
#'
#' @template template_references_hu2025
#'
#' @references
#' Li, Q. and Racine, J. S. (2008). Nonparametric estimation of conditional CDF and
#' quantile functions with mixed categorical and continuous data.
#' \emph{Journal of Business and Economic Statistics}, 26(4), 423-434
#'
#' Hayfield, T. and Racine, J. S. (2008). Nonparametric econometrics: The np package.
#' \emph{Journal of Statistical Software}, 27(5).
#' \doi{10.18637/jss.v027.i05}
#'
#' @eval doc_rendocopula2scopenp_return()
#'
#' @family copula-based methods
#' @seealso \code{\link[np:npcdistbw]{npcdistbw}} for possible elements of parameter
#'  \code{npcdistbw.args}.
#'
#' @examples
#' # # Set a random number seed because NP procedure is random
#' # set.seed(42)
#'
#' #--------------------------------------------------------------
#' # example 1: Continuous endogenous regressor
#' # (Hu et al. 2025, Section 4.3)
#' #
#' # Demonstrates 2sCOPEnp's unique advantage: it is the only method
#' # that fully eliminates bias when neither the Gaussian copula nor
#' # the mean-dependence assumption holds.
#' # True values: mu = 1, alpha = 1 (P), beta = 2 (X).
#' #--------------------------------------------------------------
#' data("data2sCOPEnpCont")
#' res1 <- copula2sCOPEnp(
#'  y ~ P + X | P,
#'   data = data2sCOPEnpCont)
#'
#' \donttest{
#' #--------------------------------------------------------------
#' # example 2: Binary endogenous regressor
#' # (Hu et al. 2025, Section 4.5)
#' #
#' # To show 2sCOPEnp's ability to handle discrete
#' # endogenous regressors. 2sCOPEnp corrects endogeneity via
#' # the smooth conditional CDF.
#' # True values: mu = 0, alpha = 1 (P), beta = 2 (X).
#' #--------------------------------------------------------------
#' data("data2sCOPEnpBi")
#' res2 <- copula2sCOPEnp(
#'   y ~ P + X | P,
#'   data = data2sCOPEnpBi
#' )
#' summary(res2)
#'
#'
#' #--------------------------------------------------------------
#' # How to tweak the bandwidth selection for the kernel estimate
#' # (parameter npcdistbw.args)
#' #
#' # For a good one-shot fit, we adapt the parameters for the
#' # bandwidth selection to prioritize reliability
#' #   nmulti: Number of restarts at random points. High to
#' #           protect from local minima.
#' #   bwmethod: Bandwidth selection method. Cross-validation ("cv.ls")
#' #           is already default over rule-of-thumb.
#' #   itmax: Max iterations before failing numerical optimization.
#' #          High to prevent silent early stop without convergence.
#' #   ftol,tol: Do not override defaults.
#' #
#' #--------------------------------------------------------------
#' res3 <- copula2sCOPEnp(
#'  y ~ P + X | P,
#'   npcdistbw.args = list(
#'   nmulti = 25,
#'   itmax = 500000
#'   # other common params: bwmethod, bwtype,
#'   ),
#'   data = data2sCOPEnpCont,
#'   num.boots = 100)
#'
#'
#' #--------------------------------------------------------------
#' # How to tweak the bandwidth selection over multiple iterations
#' #
#' # Internally, we are using `np::npcdistbw()` for the bandwidth
#' # estimation which accepts a previously computed `np::condbandwidth`
#' # object as a starting point.
#' #
#' # loosen-then-refine workflow
#' # This allows us to first search broadly using a cheap method
#' # and then further re-fine with a more comprehensive method.
#' #
#' # We can pass a previously estimate
#' # If we passed it in parameter `bws`, it would replace the entire
#' # bandwidth estimation. Instead, we can pass it as part of
#' # `npcdistbw.args` such that it enters `np::npcdistbw(bws=)`.
#'
#' # Exploration (ftol, tol)
#'
#' #--------------------------------------------------------------
#'
#' }
#'
#'
#' @export
#'
#' @importFrom Formula as.Formula
#' @importFrom stats coef terms formula
copula2sCOPEnp <- function(formula, data, bws=NULL, npcdistbw.args=list(), num.boots = 1000, verbose = TRUE) {
  cl <- match.call()

  #Input checks
  # check_err_msg(checkinput_copula2sCOPEnp_formula(formula))
  # check_err_msg(checkinput_copula2sCOPEnp_data(data))
  # check_err_msg(checkinput_copula2sCOPEnp_dataVSformula(data = data, formula = formula))
  check_err_msg(checkinput_copula2scopenp_npcdistbwargs(npcdistbw.args))
  check_err_msg(checkinput_copulashared_numboots(num.boots))
  check_err_msg(checkinput_copulashared_verbose(verbose))

  # checks:
  # - endo are continuous or ordered factors

  F.formula <- as.Formula(formula)

  labels.main <- labels(terms(F.formula, data = data, rhs = 1))
  labels.endo <- labels(terms(F.formula, data = data, rhs = 2))
  labels.exo <- labels.main[!(labels.main %in% labels.endo)]

  if (length(labels.exo) == 0) {
    stop(
      paste0(
        "2sCOPEnp requires at least one exogenous regressor for the ",
        "nonparametric conditional CDF estimation.",
        "Please include exogenous control variables in the formula."
      ),
      call. = FALSE
    )
  }

  if (verbose) {
    message(
      "Fitting 2sCOPEnp model with ",
      length(labels.endo),
      " endogenous regressor(s)."
    )
    message("Note: Nonparametric bandwidth selection could take time.")
  }

  #precomputing the bws once on original data for bootstrap reuse
  #bws estimates are consistent and sampling variability has negligible effect
  #on the bootstrap distribution of the structural coefficients
  if(is.null(bws)){
    bws <- copula2sCOPEnp_bandwidth(
      data = data,
      labels.exo = labels.exo,
      labels.endo = labels.endo,
      npcdistbw.args = npcdistbw.args,
      verbose = verbose
    )
  }

  # TODO: warn if itnmax was hit or bws are at search boundaries (lower bound for
  #   continuous, or lambda at its max). Warn here and not in copula2sCOPEnp_bandwidth
  #   because should also check the user-given bws
  # TODO: warn if the fit diagnostics of the bandwidths are bad?

  fit <- copula2sCOPEnp_fit(
    F.formula = F.formula,
    data = data,
    labels.exo = labels.exo,
    labels.endo = labels.endo,
    verbose = verbose,
    bws = bws
  )

  # Bootstrapping -------------------------------------------------------------------
  fn.fit.boots <- function(data.b) {
    return(copula2sCOPEnp_fit(
      F.formula = F.formula,
      data = data.b,
      labels.exo = labels.exo,
      labels.endo = labels.endo,
      verbose = FALSE,
      bws = bws
    ))
  }

  res.boots <- bootstrap_skip_degenerates(
    fn.fit = fn.fit.boots,
    data = data,
    num.boots = num.boots,
    coef.names = names(coef(fit)),
    verbose = verbose
  )

  # Structural residuals --------------------------------------------------------------

  l.fitted.resid <- copula_compute_structural_fitted_residuals(
    res.lm.aug = fit,
    names.aux.regs = grep("_cop$", names(coef(fit)), value = TRUE)
  )

  # Return object ----------------------------------------------------------------------

  # TODO: return conditional dist object (ie result of npcdist())
  # TODO: summary prints bw fitting: method, kernel type, bwtype (what if diverge for
  # each endo because user-supplied?), scale factors & lambdas (important) of bw estimation
  # TODO: plot: bws

  return(new_rendo_copula2sCOPEnp(
    call = cl,
    F.formula = F.formula,
    res.lm.augmented = fit,
    fitted.values = l.fitted.resid$fitted.values,
    residuals = l.fitted.resid$residuals,
    boots.params = res.boots$boots.params,
    n.boots.attempted = res.boots$n.attempted,
    n.boots.failed = res.boots$n.failed,
    bws = bws,
    names.endo.regs = labels.endo
  ))
}

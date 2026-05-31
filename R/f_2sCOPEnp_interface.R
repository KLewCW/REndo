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
#'
#'
#' @details
#' \strong{Model}
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
#'
#' \strong{Methodology}
#'
#' The estimation proceeds in 2 stages (Hu et al. 2025)
#'
#' \strong{Stage 1}: For each endogenous regressor \eqn{P_k}:
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
#' \strong{Stage 2}: Add \eqn{\{C_{i,pk}\}, k =1, \ldots, K} to the structural model as
#' generated regressors and estimate using OLS (Equation 20):
#' \deqn{Y_i = \mu + \sum_{k=1}^{K} P_{i,k}\alpha_k + X_i'\beta +
#'       \sum_{k=1}^{K} C_{i,pk}\gamma_k + \xi_i}
#' where \eqn{\gamma_k} is the coefficient of the copula correction term \eqn{C_{i, pk}} and
#' \eqn{\xi_i} is the new error term.
#'
#' \strong{Formula interface}
#'
#' The \code{formula} argument follows a two part notation separated by \code{|}.
#' The first part specifies the structural model (e.g \code{y ~ X + P}).
#' The second part identifies the endogenous regressors and their type using \code{continuous()}
#' or \code{discrete()}:
#'
#' \preformatted{y ~ X + P | continuous(P)        # continuous endogenous P}
#' \preformatted{y ~ X + D | discrete(D)          # discrete endogenous D}
#' \preformatted{y ~ X + P + D | continuous(P) + discrete(D)  # mixed}
#'
#' At least one exogenous regressor must be present in the model for the
#' nonparametric conditional CDF estimation to be feasible.
#'
#' If the bootstrap standard errors of the endogenous regressor coefficients are more than
#' 6 times larger than the corresponding OLS standard errors, this may indicate near-multicollinearity
#' between the copula correction term and the original regressors, which may suggest a potential
#' identification issue. (Hu et al. 2025, section 3.5)
#'
#' \strong{Bootstrap inference}
#'
#' Standard errors are obtained by resampling the data with replacement and
#' re-running the full two-stage estimation on each resample, including the
#' nonparametric bandwidth selection in Stage 1.
#' Degenerate bootstrap samples are automatically discarded and redrawn.
#' The percentage of discarded samples is reported as a warning if non-zero.
#'
#' @references
#' Hu, X., Qian, Y., and Xie, H. (2025). Correcting endogeneity via
#' instrument-free two-stage nonparametric copula control functions.
#' NBER Working Paper No. 33607.
#' \url{http://www.nber.org/papers/w33607}
#'
#' Li, Q. and Racine, J. S. (2008). Nonparametric estimation of conditional CDF and
#' quantile functions with mixed categorical and continuous data.
#' \emph{Journal of Business and Economic Statistics}, 26(4), 423-434
#'
#' Hayfield, T. and Racine, J. S. (2008). Nonparametric econometrics: The np package.
#' \emph{Journal of Statistical Software}, 27(5).
#' \doi{10.18637/jss.v027.i05}
#'
#' @family copula-based methods
#'
#' @examples
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
#'  y ~ P + X | continuous(P),
#'   data = data2sCOPEnpCont,
#'   num.boots = 100)
#'
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
#'   y ~ P + X | discrete(P),
#'   data      = data2sCOPEnpBi,
#'   num.boots = 100
#' )
#' summary(res2)
#'
#'
#' @export
#' @importFrom Formula as.Formula
#' @importFrom stats coef model.matrix model.frame formula
#'
copula2sCOPEnp <- function(formula, data, num.boots = 1000, verbose = TRUE) {
  cl <- match.call()

  #Input checks
  # check_err_msg(checkinput_copula2sCOPEnp_formula(formula))
  # check_err_msg(checkinput_copula2sCOPEnp_data(data))
  # check_err_msg(checkinput_copula2sCOPEnp_dataVSformula(data = data, formula = formula))
  # check_err_msg(checkinput_copula2sCOPEnp_numboots(num.boots))
  # check_err_msg(checkinput_copula2sCOPEnp_verbose(verbose))

  F.formula <- Formula::as.Formula(formula)

  # Extracting both continuous and discrete endogenous regressor names

  names.continuous <- formula_readout_special(
    F.formula = F.formula,
    name.special = "continuous",
    from.rhs = 2,
    params.as.chars.only = TRUE
  )

  names.discrete <- formula_readout_special(
    F.formula = F.formula,
    name.special = "discrete",
    from.rhs = 2,
    params.as.chars.only = TRUE
  )

  names.endo.regs <- list(
    continuous = names.continuous,
    discrete = names.discrete,
    all = c(names.continuous, names.discrete)
  )

  if (verbose) {
    message(
      "Fitting 2sCOPEnp model with",
      length(names.endo.regs$all),
      "endogenous regressor(s).",
      if (length(names.discrete) > 0) {
        paste0("(", length(names.discrete), "discrete)")
      } else {
        ""
      }
    )
    message(
      "Note: Nonparametric bandwidth selection could take time."
    )
  }

  fit <- copula2sCOPEnp_fit(
    F.formula = F.formula,
    data = data,
    names.endo.regs = names.endo.regs,
    verbose = verbose,
    bws = NULL
  )

  #precomputing the bws once on original data for bootstrap reuse
  #bws estimates are consistent and sampling variability has negligible effect
  #on the bootstrap distribution of the structural coefficients
  if (verbose) {
    message(
      "Pre-computing bandwidths for ",
      length(names.endo.regs$all),
      " endogenous regressor(s) to speed up bootstrapping..."
    )
  }

  #trying to identify the exogenous columns just like we did in fit.
  #to check if we have correct columns
  F.formula.main <- formula(F.formula, rhs = 1, lhs = 1)
  mf.original <- model.frame(F.formula.main, data = data)
  X.main.original <- model.matrix(F.formula.main, data = mf.original)

  endo.cols.original <- colnames(X.main.original)[
    colnames(X.main.original) %in% names.endo.regs$all
  ]
  exo.cols.original <- colnames(X.main.original)[
    !colnames(X.main.original) %in% c("(Intercept)", endo.cols.original)
  ]

  bws.original <- lapply(endo.cols.original, function(p.var) {
    copula2sCOPEnp_bandwidth(
      y.data = data[, p.var, drop = FALSE],
      x.data = data[, exo.cols.original, drop = FALSE]
    )
  })

  # Bootstrapping -------------------------------------------------------------------
  fn.fit.boots <- function(data.b) {
    return(copula2sCOPEnp_fit(
      F.formula = F.formula,
      data = data.b,
      names.endo.regs = names.endo.regs,
      verbose = FALSE,
      bws = bws.original
    ))
  }

  res.boots <- bootstrap_skip_degenerates(
    fn.fit = fn.fit.boots,
    data = data,
    num.boots = num.boots,
    coef.names = names(coef(fit)),
    verbose = verbose
  )

  # Return value ----------------------------------------------------------------------
  return(new_rendo_copula2sCOPEnp(
    call = cl,
    F.formula = F.formula,
    res.lm = fit,
    boots.params = res.boots$boots.params,
    n.boots.attempted = res.boots$n.attempted,
    n.boots.failed = res.boots$n.failed,
    names.endo.regs = names.endo.regs$all
  ))
}

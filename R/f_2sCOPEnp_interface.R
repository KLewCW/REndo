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
#' @references
#' Hu, X., Qian, Y., and Xie, H. (2025). Correcting endogeneity via
#' instrument-free two-stage nonparametric copula control functions.
#' NBER Working Paper No. 33607.
#' \url{http://www.nber.org/papers/w33607}
#'
#' @export
#' @importFrom Formula as.Formula
#' @importFrom stats coef
copula2sCOPEnp <- function(formula, data, num.boots = 1000, verbose = TRUE){
  cl <- match.call()

  #Input checks
  check_err_msg(checkinput_copula2sCOPEnp_formula(formula))
  check_err_msg(checkinput_copula2sCOPEnp_data(data))
  check_err_msg(checkinput_copula2sCOPEnp_dataVSformula(data = data, formula = formula))
  check_err_msg(checkinput_copula2sCOPEnp_numboots(num.boots))
  check_err_msg(checkinput_copula2sCOPEnp_verbose(verbose))

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

  if (verbose){
    message("Fitting 2sCOPEnp model with", length(names.endo.regs$all), "endogenous regressor(s).",
            if(length(names.discrete) > 0){
              paste0("(", length(names.discrete), "discrete)")
            } else {
              ""
            })
            message(
              "Note: Nonparametric bandwidth selection could take time."
            )
  }

  fit <- copula2sCOPEnp_fit(
    F.formula = F.formula,
    data = data,
    names.endo.regs = names.endo.regs,
    verbose = verbose
  )

  # Bootstrapping -------------------------------------------------------------------
  fn.fit.boots <- function(data.b){
    return(copula2sCOPEnp_fit(
      F.formula = F.formula,
      data = data.b,
      names.endo.regs = names.endo.regs,
      verbose = FALSE
    ))
  }

  res.boots <- bootstrap_skip_degenerated(
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

#' Copula-based Endogeneity Correction with Asymptotic properties (BMW)
#'
#' @description
#' Fits a linear model with endogenous regressors by using the nonparametric
#' control function approach of Breitung, Meyer and Wied (2024). This method
#' corrects endogeneity without any external instrumental variables. It is a
#' copula-based method with asymptotic theory.
#'
#' @details
#' The estimator is done in two steps:
#' \itemsize{
#' \item First: Each endogenous regressor \eqn{P_k} is regressed on the
#' exogenous regressors \eqn{X} using ordinary leas squares in the original
#' variable space to obtain residuals \eqn{\hat{e}_k = P_k - \hat{\delta}' X}.
#' \item Second: The empirical CDF is applied to \eqn{\hat{e}_k} and the result
#' is then transformed through \eqn{\Phi^{-1}} in order to find the correction
#' term \eqn{\hat{\eta}_k}. This is then included as an additional regressor in
#' augmented OLS.
#' }
#'
#' The model is then:
#' \deqn{y = \beta' X + \gamma P + \rho \hat{\eta} + \xi}
#'
#' where \eqn{\hat{\eta} = \Phi^{-1} (\hat{F}_{\hat{e}}(\hat{e}))}
#' \eqn{\hat{F}_{\hat{e}}} is the empirical CDF using.
#'
#' This method requires at least one exogenous regressor for the first-stage
#' regression and supports only continuous regressors.
#'
#' @references
#' Breitung, J., Meyer, M., Wied, D. (2024). Asymptotic properties of endogeneity
#' corrections using nonlinear transformations. \emph{The Econometrics Journal},
#' 27, 362--383/ \doi{10.1093/ectj/utae002}
#'
#' @examples
#'
#'
#'
#' @export
#' @importFrom stats coef
copulaBMW <- function(
    formula,
    data,
    cdf = c("ecdf", "adj.ecdf", "resc.ecdf", "kde"),
    num.boots = 1000,
    verbose = TRUE){

  cl <- match.call()

  check_err_msg(checkinput_copulaBMW_formula(formula))
  check_err_msg(checkinput_copulaBMW_data(data))
  check_err_msg(checkinput_copulaBMW_dataVSformula(data = data, formula = formula))
  check_err_msg(checkinput_copulaBMW_numboots(num.boots))
  check_err_msg(checkinput_copulaBMW_verbose(verbose))
  check_err_msg(checkinput_copulaBMW_cdf(cdf))

  cdf <- match.arg(cdf, choices = c("ecdf", "adj.ecdf", "resc.ecdf", "kde"))

  F.formula <- Formula::as.Formula(formula)
  names.endo.regs <- formula_readout_special(
    F.formula = F.formula,
    name.special = "continuous",
    from.rhs = 2,
    params.as.chars.only = TRUE
  )

  rhs1.vars <- all.vars(formula(F.formula, rhs = 1, lhs = 0))
  exo.vars <- rhs1.vars[!rhs1.vars %in% names.endo.regs]

  if(length(exo.vars) == 0){
    stop("No exogenous regressors were found. BMW method requires at least one",
         "exogeous regressor for the first-stage regression of each endogenous regressor ",
         "P on X.", call. = FALSE) #equation 2.2 & assumption A4
  }

  if(verbose){
    message("Fitting BMW copula model with",
            length(names.endo.regs),
            "endogenous regressors.")
  }

  fit <- copulaBMW_fit(
    F.formula = F.formula,
    data = data,
    names.endo.regs = names.endo.regs,
    cdf = cdf
  )

  fn.fit.boots <- function(data.b){
    return(copulaBMW_fit(F.formula = F.formula, data = data.b, cdf = cdf))
  }

  res.boots <- bootstrap_skip_degenerates(
    fn.fit     = fn.fit.boots,
    data       = data,
    num.boots  = num.boots,
    coef.names = names(coef(fit)),
    verbose    = verbose
  )

  return(new_rendo_copulaBMW(
    call              = cl,
    F.formula         = F.formula,
    res.lm            = fit,
    boots.params      = res.boots$boots.params,
    n.boots.attempted = res.boots$n.attempted,
    n.boots.failed    = res.boots$n.failed,
    cdf               = cdf,
    names.endo.regs   = names.endo.regs
  ))
}

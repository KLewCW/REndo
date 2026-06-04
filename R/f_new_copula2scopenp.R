doc_rendocopula2scopenp_return_list <- function() {
  doc_boots <- doc_rendobootsdegeneratesremoved_return_list()

  doc_boots[['residuals']] <- "\\item{\\code{residuals}}{The structural residuals.}"
  doc_boots[[
    'fitted.values'
  ]] <- "\\item{\\code{fitted.values}}{Fitted values of the structural model.}"

  doc_copula2scopenp <- c(
    bws = "\\item{\\code{bws}}{The bandwidths computed on the input data: a named list with one \\code{np::condbandwidth} object per endogenous regressor (named accordingly), each used to estimate the conditional CDF \\eqn{\\hat{F}(P_k \\mid X)}.}",
    names.endo.regs = "\\item{\\code{names.endo.regs}}{The names of the endogenous regressors.}",
    res.lm.augmented = "\\item{\\code{res.lm.augmented}}{The fitted augmented regression model, including the control function terms.}"
  )

  return(c(doc_boots, doc_copula2scopenp))
}


doc_rendocopula2scopenp_return <- function() {
  doc_intro <- c(
    return = "@return An object of class \\code{rendo.copula.2sCOPE.np} which is a list that contains:"
  )

  return(c(doc_intro, doc_rendocopula2scopenp_return_list()))
}

#' @importFrom stats coef model.frame
new_rendo_copula2sCOPEnp <- function(
    call,
    F.formula,
    fitted.values,
    residuals,
    res.lm.augmented,
    boots.params,
    n.boots.attempted,
    n.boots.failed,
    names.endo.regs,
    bws
) {
  return(.new_rendo_boots_degenerates_removed(
    # Stuff for rendo.boots.degenerates.removed class
    call = call,
    F.formula = F.formula,
    mf = model.frame(res.lm.augmented),
    coefficients = coef(res.lm.augmented),
    names.main.coefs = names(coef(res.lm.augmented)), # OR: row.names(boots.params)
    fitted.values = fitted.values,
    residuals = residuals,
    boots.params = boots.params,
    n.boots.attempted = n.boots.attempted,
    n.boots.failed = n.boots.failed,

    # 2sCOPEnp-specific
    subclass = "rendo.copula.2sCOPE.np",
    res.lm.augmented = res.lm.augmented,
    names.endo.regs = names.endo.regs,
    bws = bws
  ))
}

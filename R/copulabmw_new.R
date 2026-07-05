doc_rendocopulabmw_return_list <- function() {
  doc_boots <- doc_rendobootsdegeneratesremoved_return_list()

  doc_boots[['residuals']] <- "\\item{\\code{residuals}}{The structural residuals.}"
  doc_boots[[
    'fitted.values'
  ]] <- "\\item{\\code{fitted.values}}{Fitted values of the structural model.}"

  doc_copulabmw <- c(
    cdf = "\\item{\\code{cdf}}{The used cdf function.}",
    names.endo.regs = "\\item{\\code{names.endo.regs}}{The names of the continuous endogenous regressors.}",
    res.lm.augmented = "\\item{\\code{res.lm.augmented}}{The fitted augmented regression model, including the control function terms.}"
  )

  return(c(doc_boots, doc_copulabmw))
}


doc_rendocopulabmw_return <- function() {
  doc_intro <- c(
    return = "@return An object of class \\code{rendo.copula.bmw} which is a list that contains:"
  )

  return(c(doc_intro, doc_rendocopulabmw_return_list()))
}

#' @importFrom stats coef model.frame
new_rendo_copulabmw <- function(
    call,
    F.formula,
    res.lm.augmented,
    fitted.values,
    residuals,
    boots.params,
    n.boots.attempted,
    n.boots.failed,
    cdf,
    labels.endo,
    labels.exo,
    labels.pcop
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

    # Stuff specific to copulaBMW
    subclass = "rendo.copula.bmw",
    res.lm.augmented = res.lm.augmented,
    cdf = cdf,
    labels.endo = labels.endo,
    labels.exo = labels.exo,
    labels.pcop = labels.pcop
  ))
}

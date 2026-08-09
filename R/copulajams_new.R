doc_rendocopulajams_return_list <- function() {
  doc_boots <- doc_rendobootsdegeneratesremoved_return_list()

  doc_boots[['residuals']] <- "\\item{\\code{residuals}}{The structural residuals.}"
  doc_boots[[
    'fitted.values'
  ]] <- "\\item{\\code{fitted.values}}{Fitted values of the structural model.}"

  doc_copulajams <- c(
    cdf = "\\item{\\code{cdf}}{The used cdf function.}",
    res.lm.augmented = "\\item{\\code{res.lm.augmented}}{The fitted augmented regression model, including the control function terms.}",
    labels.endo = "\\item{\\code{labels.endo}}{The term labels of the endogenous regressors.}",
    labels.exo = "\\item{\\code{labels.exo}}{The term labels of the exogenous regressors.}",
    labels.pcop = "\\item{\\code{labels.pcop}}{The term labels of the generated, auxiliary regressors.}",
    P  = "\\item{\\code{P}}{The matrix of (continuous) endogenous regressors.}",
    W = "\\item{\\code{W}}{The matrix of continuous exogenous regressors. Only main effects: Discrete (\\code{factor}) regressors and interaction terms are excluded.}",
    df.Z = "\\item{\\code{df.Z}}{\\code{data.frame} of discrete (\\code{factor}) exogenous regressors used to stratify. \\code{NULL} if there are none.}"
  )

  return(c(doc_boots, doc_copulajams))
}


doc_rendocopulajams_return <- function() {
  doc_intro <- c(
    return = "@return An object of class \\code{rendo.copula.jams} which is a list that contains:"
  )

  return(c(doc_intro, doc_rendocopulajams_return_list()))
}

#' @importFrom stats coef model.frame
new_rendo_copulajams <- function(
  call,
  F.formula,
  fitted.values,
  residuals,
  res.lm.augmented,
  boots.params,
  n.boots.attempted,
  n.boots.failed,
  cdf,
  labels.endo,
  labels.exo,
  labels.pcop,
  P,
  W,
  df.Z
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

    # Stuff specific to JAMS
    subclass = "rendo.copula.jams",
    res.lm.augmented = res.lm.augmented,
    cdf = cdf,
    labels.endo = labels.endo,
    labels.exo = labels.exo,
    labels.pcop = labels.pcop,
    P = P,
    W = W,
    df.Z = df.Z
  ))
}

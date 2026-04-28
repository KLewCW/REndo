#' @importFrom Formula as.Formula
#' @importFrom stats lm model.frame model.matrix formula update reformulate
#'

copulaBMW_fit <- function(F.formula, data, names.endo.regs, cdf){

  mf <- model.frame(F.formula, data = data)
  f.main <- formula(mf)

  X.main <- model.matrix(F.formula, data = mf)
  endogenous.cols <- colnames(X.main)[colnames(X.main) %in% names.endo.regs]

  if (length(endogenous.cols)==0)
    stop("No endogenous regressors found in design matrix.")
  if(lenght(endogenous.cols) < length(names.endo.regs))
    stop("Bootstrap sample dropped at least one endogenous regressor. ",
         "This can happen when a regressor becomes constant in a resample."
         )

  # Exogenous columns
  exo.cols <- setdiff(
    colnames(X.main)[colnames(X.main) != "(Intercept)"], endogenous.cols
  )

  # BMW correction
  # step 1: first-stage in the original space
  # step 2: CDF on residuals
  #step 3: then apply qnorm
  cop.terms <- copulaBMW_correction(
    data = data,
    endo.cols = endogenous.cols,
    exo.cols = exo.cols,
    cdf = cdf
  )

  f.pcop <- reformulate(
    termlabels = c(".", colnames(cop.terms)),
    response = NULL,
    intercept = TRUE
  )
  f.final<- update(old = f.main, new = f.pcop)

  return(lm(
    formula = f.final,
    data = cbind(data, cop.terms)
  ))

}

#' @importFrom Formula as.Formula
#' @importFrom stats lm model.frame model.matrix formula update reformulate terms cor
#'
#'
copulaJAMS_fit <- function(f.main, data, names.endo.regs, names.exo.regs, cdf){

  mf <- model.frame(f.main, data= data)
  X.main <- model.matrix(f.main, data= mf)

  endogenous.cols <- colnames(X.main)[colnames(X.main) %in% names.endo.regs]

  if (length(endogenous.cols)==0)
    stop( "No continuous endogenous regressors found in the design matrix",
         call.= FALSE )

  if(length(endogenous.cols) < length(names.endo.regs))
    stop (
      "Bootstrap sample dropped at least one endogenous regressor.",
      "This can happen when a regressor becomes constant in a resample.",
      call. = FALSE
    )

  #All the regressors excluding the intercept.

  P.all <- X.main[, colnames(X.main) != "(Intercept)", drop = FALSE]

  #Checking if there is any factor variables among exogenous regressors
  # if the factor Z is present, then it requires stratified correction per level (see equation 20 and 21)
  # if there is no factor Z, variance-covariance correction (eq. 17 to 19)

  factor.vars <- names(which(sapply(data[names.exo.regs], is.factor)))

  if(length(factor.vars)==0){
    #case no Z
    cop.terms <- copulaJAMS_correction_cont(
      P.all = P.all,
      names.endo.regs = endogenous.cols,
      cdf = cdf
    )
  } else{ #case where Z is present

    cop.terms <- copulaJams_correction_dis(
      data = data,
      names.endo.regs = endogenous.cols,
      names.exo.regs = names.exo.regs,
      factor.vars = factor.vars,
      cdf = cdf
    )
  }

  f.main <- formula(mf)
  has.intercept <- attr(terms(f.main), "intercept") == 1. #equation 19 & 20

  f.pcop <- reformulate(
    termlabels = c(".", colnames(cop.terms)),
    response = NULL,
    intercept = has.intercept
  )

  f.final <- update(old = f.main, new = f.pcop)

  return(lm(formula = f.final, data = cbind(data, cop.terms)))

}

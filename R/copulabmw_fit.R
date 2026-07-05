#' @importFrom stats lm reformulate
copulaBMW_fit <- function(F.formula, data, cdf, labels.endo, labels.exo) {

  # Stage 1 --------------------------------------------------------------------
  # BMW correction
  # step 1: first-stage in the original space
  # step 2: CDF on residuals
  # step 3: apply qnorm

  cop.terms <- copula_create_1ststage_copdata_matrix(
    n = nrow(data),
    labels.endo = labels.endo
  )

  for (k in seq_along(labels.endo)) {
    p.label <- labels.endo[k]

    #case no exo regressors
    #e hat = z - mean(z)
    # first stage with intercept
    # if (length(exo.cols) == 0) {
    # if required, would use formula approach:
    # f.Z <- reformulate(termlabels = p.label, response = NULL, intercept = FALSE)
    # mf.z <- model.frame(f.Z, data=data)
    # Z <- mf.z[, 1, drop=TRUE]
    #   e.hat <- Z - mean(Z)
    # } else {
    #first-stage OLS of Z on X in original space
    #BMW (2024) eq. 2.2, Z = delta'x + e

    # Regress endo ~ (all exo), where (all exo) excludes the intercept
    f.endo.k.on.all.exo <- reformulate(
      response = p.label,
      termlabels = labels.exo,
      intercept = FALSE
    )

    res.lm.first <- lm(formula = f.endo.k.on.all.exo, data = data)
    e.hat <- residuals(res.lm.first)
    # }

    #Apply CDF now, then qnorm to residuals e hat
    P.star <- copulaBMW_pstar(e.hat = e.hat, cdf = cdf)

    #Apply qnorm from eq. 2.3, eta hat =  phi^{-1} (F hat_{e hat} (e hat))
    P.cop <- apply(P.star, 2, qnorm) #eta hat is P_cop

    cop.terms[, k] <- as.vector(P.cop)
  }

  # Stage 2 --------------------------------------------------------------------
  res.2nd.stage <- copula_fit_2ndstage(
    F.formula = F.formula,
    data = data,
    cop.terms = cop.terms
  )

  return(list(
    res.augmented = res.2nd.stage$res.augmented,
    labels.pcop = res.2nd.stage$labels.pcop
  ))
}

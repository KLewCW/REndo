#' @importFrom stats lm reformulate
copulabmw_fit <- function(F.formula, data, cdf, labels.endo, labels.exo) {
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

    # Step 1: Original space --------------------------------------------------------
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

    # Step 2: CDF on residuals ------------------------------------------------------
    # Apply CDF to the first-stage residuals (e hat)
    # Does not apply on the original regressors

    m.ehat <- matrix(e.hat, ncol=1)

    if (cdf == "ecdf") {
      #ecdf: using the theoretical recommendation from BMW (2024) eq. 2.3
      #instead of ecdf() + 10e-7
      #here we use rank/(n+1), so that no arbitary boundary constant is needed.
      #will still keep all the values strictly in (0,1)
      P.star <- copulabmw_ecdf(m.ehat)
    } else {
      # usual cdfs
      P.star <- copula_pstar(P = m.ehat, cdf = cdf)
    }

    # Step 3: Apply qnorm ----------------------------------------------------------
    #Apply qnorm from eq. 2.3, eta hat =  phi^{-1} (F hat_{e hat} (e hat))

    P.cop <- apply(P.star, 2, qnorm) # eta hat is P_cop

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


#According to BMW(2024), they adopt a 'common practice' and rescale by n + 1
# Recommendation from eq. 2.3
#F hat_{e hat} (e hat_i) = rank (e hat_i)/(n+1)
copulabmw_ecdf <- function(P) {
  stopifnot(is.matrix(P))

  n <- nrow(P)
  U <- apply(P, 2, rank, na.last = "keep", ties.method = "average") / (n + 1)

  colnames(U) <- colnames(P)
  return(U)
}

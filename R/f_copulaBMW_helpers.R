#' @importFrom Formula as.Formula
#' @importFrom stats lm residuals qnorm ecdf
#' @importFrom copula pobs
#' @importFrom ks kcde
#'
pobs_adj <- function(x) {
  n <- if (is.matrix(x)) nrow(x) else length(x)

  if (is.matrix(x)) {
    U <- apply(x, 2, rank, na.last = "keep", ties.method = "average") *
      ((n - 1) / n^2) +
      1 / (2 * n)
  } else {
    U <- rank(x, na.last = "keep", ties.method = "average") *
      ((n - 1) / n^2) +
      1 / (2 * n)
  }

  return(U)
}

#According to BMW(2024), they adopt a 'common practice' and rescale by n + 1
# Recommendation from eq. 2.3
#F hat_{e hat} (e hat_i) = rank (e hat_i)/(n+1)
copulaBMW_ecdf <- function(x) {
  n <- if (is.matrix(x)) nrow(x) else length(x)

  if (is.matrix(x)) {
    U <- apply(x, 2, rank, na.last = "keep", ties.method = "average") / (n + 1)
  } else {
    U <- rank(x, na.last = "keep", ties.method = "average") / (n + 1)
  }

  return(U)
}

#Apply CDF to the first-stage residuals (e hat)
# Does not apply on the original regressors

copulaBMW_pstar <- function(e.hat, cdf) {
  if (!is.matrix(e.hat)) {
    e.hat <- matrix(e.hat, ncol = 1)
  }

  if (cdf == "kde") {
    P.star <- apply(e.hat, 2, function(x) {
      Fhat <- ks::kcde(x)
      predict(Fhat, x = x)
    })
  } else if (cdf == "resc.ecdf") {
    P.star <- apply(e.hat, 2, copula::pobs)
  } else if (cdf == "adj.ecdf") {
    P.star <- apply(e.hat, 2, pobs_adj)
  } else {
    #using the theoretical recommendation from BMW (2024) eq. 2.3
    #instead of ecdf() + 10e-7
    #here we use rank/(n+1), so that no arbitary boundary constant is needed.
    #will still keep all the values strictly in (0,1)

    P.star <- copulaBMW_ecdf(e.hat)
    P.star <- as.matrix(P.star)
    colnames(P.star) <- colnames(e.hat)
  }

  return(P.star)
}


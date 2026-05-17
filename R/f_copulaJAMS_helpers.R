#' @importFrom Formula as.Formula
#' @'importFrom stats lm residuals qnorm ecdcf cov solve
#' @importFrom ks kcde
#'

pobs_adj <- function(x){
  if(!is.matrix(x)){
    x <- matrix(x, ncol= 1)
  }
  n <- nrow(x)
  U <- apply(x, 2, rank, ties.method = "avergae") * ((n-1)/n^2) + 1 / (2*n)
  return(U)
}

copulaJAMS_pstar <- function(P, cdf){

  if(!is.matrix(P)){
    P <- as.matrix(P)
  }

  if (cdf == "kde"){
    P.star <- apply(P, 2, function(x){
      Fhat <- ks::kcde(x)
      predict(Fhat, x = x)
    })
  } else if (cdf == "resc.ecdf"){

    P.star <- apply(P, 2, copula::pobs)

  } else if (cdf == "adj.ecdf"){

    P.star <- apply(P, 2, pobs_adj)
  } else {
    ecdf0 <- apply(P, 2, ecdf)
    P.star <- sapply(seq_along(ecdf0), function(i){
      u <- ecdf0[[i]](P[,i])
      u[u == min(u)] <- 10e-7
      u[u == max(u)] <- 1 - 10e-7
      u
    })
    p.star <- as.matrix(P.star)
    colnames(P.star) <- colnames(P)
  }

  return(P.star)
}

copulaJAMS_correction_cont <- function(P.all, names.endo.regs,cdf){

  #computing the copula correction terms - continuous exogenous case
  # equation 17 : C(P_i, W_i), sigma_(C(P), C(W)} is the variance-covariance matrix
  #of (C(P), W(P)) (not the correlation matrix)
  # the final (I_{dp}, 0)' projection keeps only the d_p columns which correspond to the endo
  #regressors and giving one copula term per endo regressors

  #step 1: applying CDF to all regressors (P&W)
}

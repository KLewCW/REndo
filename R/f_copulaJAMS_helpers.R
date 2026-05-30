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
    P.star <- as.matrix(P.star)
    colnames(p.star) <- colnames(P)
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
  P.star <- copulaJAMS_pstar(P= P.all, cdf = cdf)

  #step 2: applying qnorm to obtain normal scores C(P) and C(W)
  C.all <- apply(P.star, 2, qnorm) #using equation C(P_i) = Phi^{-1}(F_hat(P_i)) for each regressor

  #step 3: Estimate variance-covariance matrix \hat{sigma}_{C(P), C(W)}
  Sigma.hat <- cov(C.all) #equation 21: mentioned that it is the variance-covariance so cov should be used
  #instead of cor(). i think ?

  #step 4: inverse of the variance-covariance matrix:
  Sigma.inverse <- solve(Sigma.hat)

  #step 4:finding the d_P copula term to project onto endogenous regressor columns
  # (I_{dp}, 0_{dw x dp})' selects the d_P endo columns of Sigma^{-1}
  endo.index <-  which(colnames(C.all) %in% names.endo.regs) # giving the d_P copula terms (one per endo reg)

  #C(P_i, W_i) = (C(P_i)', C(W_i)') \hat(Sigma^{-1}_{C(P), C(W)}) (I_{dp}, 0_{dw x dp})'
  P.cop <- C.all %*% Sigma.inv[, endo.index, drop = FALSE]
  P.cop <- as.matrix(P.cop)
  colnames(P.cop) <- paste0(names.endo.regs, "_cop")

  return(P.cop)

}

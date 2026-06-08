#' @importFrom Formula as.Formula
#' @importFrom stats lm residuals qnorm ecdf cov
#' @importFrom ks kcde
#'
pobs_adj <- function(x){
  if(!is.matrix(x)){
    x <- matrix(x, ncol= 1)
  }
  n <- nrow(x)
  U <- apply(x, 2, rank, ties.method = "average") * ((n-1)/n^2) + 1 / (2*n)
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
  P.cop <- C.all %*% Sigma.inverse[, endo.index, drop = FALSE]
  P.cop <- as.matrix(P.cop)
  colnames(P.cop) <- paste0(names.endo.regs, "_cop")

  return(P.cop)

}


#Copula correction terms in case of factor or discrete exo case
#' @ importFrom stats cov qnorm

copulaJams_correction_dis <- function(data, names.endo.regs, names.exo.regs, factor.vars, cdf){

  n <- nrow(data)
  result.list <- list()

  for (var in factor.vars){
    levels.var <- levels(as.facto(data[[var]]))

    for (lvl in level.var){

      index <- data[[var]] == lvl
      subdat1 <- data[index, , drop = FALSE]

      #skip if too few observations or no variation in endogenous regressors
      has.variation <- any(sapply(subdat1[, names.endo.regs, drop = FALSE],
                                  function(col) length(unique(col)) >1))
      if (nrow(subdat1) <= 3 || !has.variation) next # 3 is number checked by Haschka's repo line 190. #did not find paper backing this up ?


      ##usually need at least p+1 observations to estimate a pxp cov matrix
      # more generally, instead of '3', we could have tried the most conservative minimum of
      # max(10, 2*p) for reliable estimation
      #min.observation.needed <- max(10L, 2L * length(cols.use))
      #if(nrow(subdat1) <= min.observation.needed || !has.variation){
      # then maybe issue a warning here... saying that the factor level 'lvl' of the variable 'var'
      # does not have enough observation for a reliable covariance estimation and we are skipping this level.


      #now use non-factor columns only for the CDF transformation.
      #factor vairables cannot enter the CDF transformation as they are discrete with no meaningful continuous CDF

      cols.use <- setdiff(c(names.endo.regs, names.exo.regs), factor.vars)
      subdat2 <- as.matrix(subdat1[, cols.use, drop = FALSE])

      subdat2<-subdat2[, apply(subdat2,2, function(x) length(unique(x)) > 1), drop = FALSE] #keeping only columns with variation with subset

      if (ncol(subdat2) ==0) next

      #Using the steps from equation 21 again:

      P.star <- copulaJAMS_pstar(P = subdat2, cdf = cdf)

      C.sub <- apply(P.star, 2, qnorm)

      #checking invertibility before continuing
      is.issue <- tryCatch({
        Sigma.hat <- cov(C.sub)
        solve(Sigma.hat)
        FALSE
        }, error = function(e) TRUE, warning = function(w) TRUE)

      if(is.issue){
        #return zero correction terms
        K.use <- sum(colnames(subdat2) %in% names.endo.regs)
        P.cop <- matrix(0, nrow = nrow(subdat2), ncol = K.use)}
      else{
        Sigma.hat <- cov(C.sub)
        Sigma.inv <- solve(Sigma.hat)

        #project onto endo cols in this subset
        endo.index <- which(colnames(C.sub) %in% names.endo.regs)
        K.use <- length(endo.index)

        P.cop <- C.sub %*% Sigma.inv[, endo.index, drop =FALSE]
        P.cop <- as.matrix(P.cop)
      }

      #naming correction terms with factor level info
      K.actual <- ncol(P.cop)
      endo.for.names <- names.endo.regs[seq_len(K.actual)]
      colnames(P.cop) <- paste0(endo.for.names, "_", var, "_", lvl, "_cop")

      # Expanding back to a full dataset
      # I(Z_i = z) from eq. 20. zero attributed for observations not in this level

      P.cop.full <- amtrix(0, nrow = n, ncol= col(P.cop))
      P.cop.full[idx, ] <- P.cop
      colnames(P.cop.full) <- colnames(P.cop)

      result.list[[paste(var, lvl, sep = "_")]] <- P.cop.full

    }
  }

  if (length(result.list) == 0)
    stop(
      "No valid factor-level subsets found for correction term estimation. ",
      "Check that factor exogenous regressors have sufficient observations ",
      "per level (> 3) and variation in the endogenous regressors.",
       call. = FALSE
      )

  return (do.call(cbind, result.list))
}

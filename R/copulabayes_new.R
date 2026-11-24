new_rendo_copulabayes <- function(
  call,
  F.formula,
  chain, #post burnin thinned chain (all columns)
  chain.struct, #structural parameter draws: alpha, delta, beta, sigma2
  chain.rho, #endo-error copula correlation draws: rho_1,...,rho_K
  post.mean,
  post.sd,
  post.lo,
  post.hi,
  fitted.values,
  residuals,
  names.endo.regs,
  n.obs,
  n.iterations,
  burnin,
  thin,
  n.draws,
  method
) {
  return(structure(
    list(
      call = call,
      F.formula = F.formula,
      chain = chain,
      chain.struct = chain.struct,
      chain.rho = chain.rho,
      post.mean = post.mean,
      post.sd = post.sd,
      post.lo = post.lo,
      post.hi = post.hi,
      fitted.values = fitted.values,
      residuals = residuals,
      names.endo.regs = names.endo.regs,
      n.obs = n.obs,
      n.iterations = n.iterations,
      burnin = burnin,
      thin = thin,
      n.draws = n.draws,
      method = method
    ),
    class = "rendo.copula.bayes"
  ))
}

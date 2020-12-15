#' rmvnorm_prec 
#' @details
#' simualte parameters from the joint precision matrix derived from a TMB objects
#' @param mu vector of MLE both fixed and random effect parameters
#' @param prec precision matrix, derived from sdreport(obj, getJointPrecision = T)
#' @param n.sims integer number of simulations
#' @param random_seed integer seed
#' @return simulated parameters.
rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
  set.seed( random_seed )
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L = Matrix::Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  return(mu + z)
}


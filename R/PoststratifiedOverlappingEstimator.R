#' PoststratifiedOverlappingEstimatorAlt
#' @details a functon that replicates the variacnce estimator from Millar & Olsen 1995
#' Post stratified variance estimator
#' @param y matrix of data which represents data from an areal view of the survey, use NA's for no data/out of bounds
#' @param nx_strat number of sampling units to the right (x-axis) in each strata
#' @param ny_strat number of sampling units in the below (y-axis) in each strata
#' @export
#' @return Post-stratified variance estimator for the poopulation total based on overlapping strata
#' 
PoststratifiedOverlappingEstimator = function(y, nx_strat = 1, ny_strat = 1) {
  ## overlapping from Millar Olsen 1995
  var_over = 0.0
  n_row = nrow(y)
  n_col = ncol(y)
  nh = (nx_strat + 1) * (ny_strat + 1)
  sum_nh = (n_row - nx_strat) * (n_col - nx_strat) * nh;
  strat = vector()
  for(i in 1:(n_row - nx_strat)) {
    for(j in 1:(n_col - nx_strat)) {
      y_h = as.vector(y[c(i:(i + ny_strat)), c(j : (j + nx_strat))])
      var_over = var_over +  var(y_h) / length(y_h)
      strat = c(strat, var(y_h) / length(y_h))
    }
  }
  return(var_over / ((n_row - nx_strat) * (n_col - nx_strat) ))
}

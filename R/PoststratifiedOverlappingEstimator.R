#' PoststratifiedOverlappingEstimator 
#' @detail a functon that replicates the variacnce estimator from Millar & Olsen 1995
#' Post stratified variance estimator
#' @param y matrix of data which represents data from an areal view of the survey, use NA's for no data/out of bounds
#' @param sample_area area of each sampling unit, either a single value of all the same, otherwise a matrix of equal dimensions
#' @param survey_area scalar for survey area
#' @nx_strat number of sampling units to the right (x-axis) in each strata
#' @nt_strat number of sampling units in the below (y-axis) in each strata
#' @export
#' @return Post-stratified variance estimator
#' 
PoststratifiedOverlappingEstimator = function(y, sample_area, survey_area, nx_strat = 1, ny_strat = 1) {
  ## overlapping from Millar Olsen 1995
  var_over = 0.0
  n_row = nrow(y)
  n_col = ncol(y)
  area_mat = NULL;
  if (class(sample_area) == "numeric" & length(area) == 1) {
    area_mat = matrix(sample_area, nrow = n_row, ncol = n_col)
  } else if (class(sample_area) == "matrix" & all(dim(sample_area) %in% dim(y)) ) {
    area_mat = sample_area
  } else {
    stop("PoststratifiedOverlappingEstimator: expects 'area' to be either a scalar or matrix of equal dim to y, this is not the case")
  }
  
  for(i in 1:(n_row - nx_strat)) {
    for(j in 1:(n_col - nx_strat)) {
      y_h = as.vector(y[c(i:(i + ny_strat)), c(j : (j + nx_strat))])
      mean_val = sum(y_h, na.rm = T) / sum(!is.na(y_h))
      var_over = var_over +  sum((y_h - mean_val)^2) / (sum(!is.na(y_h)) - 1)
    }
  }
  result = (1 / (sum(area_mat) / survey_area))^2 * var_over
  return(result)
}

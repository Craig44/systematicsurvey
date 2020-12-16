#' PoststratifiedNonOverlappingEstimator 
#' @details a functon that replicates the variacnce estimator from Lundberg & Strand 2014 eq 7
#' Post stratified variance estimator minus finite population correction.
#' @param y matrix of data which represents data from an areal view of the survey, use NA's for no data/out of bounds
#' @param sample_area area of each sampling unit, either a single value of all the same, otherwise a matrix of equal dimensions
#' @param survey_area scalar for total survey area
#' @param nx_strat number of sampling units in a stratum on the x-axis
#' @param ny_strat number of sampling units in a stratum on the y-axis
#' @export
#' @return Post-stratified variance estimator for population total based on non-overlapping stratum
#' 
PoststratifiedNonOverlappingEstimator = function(y, sample_area, survey_area, nx_strat = 2, ny_strat = 2) {
  if(class(y) != "matrix")
    stop(paste0("expect y, to be a matrix"))
  if(length(survey_area) != 1)
    stop(paste0("expect survey_area to be a scalar, total area of survey"))
  area_mat = NULL;
  if (class(sample_area) == "numeric" & length(sample_area) == 1) {
    area_mat = matrix(sample_area, nrow = n_row, ncol = n_col)
  } else if (class(sample_area) == "matrix" & all(dim(sample_area) %in% dim(y)) ) {
    area_mat = sample_area
  } else {
    stop("PoststratifiedNonOverlappingEstimator: expects 'sample_area' to be either a scalar or matrix of equal dim to y, this is not the case")
  }
  # wh = nh / sum(nh) 
  n_row = nrow(y)
  n_col = ncol(y)
  nh = nx_strat * ny_strat
  ## index for top right sampling unit of each strata
  row_ndx = seq(ny_strat, n_row, by = ny_strat);
  col_ndx = seq(nx_strat, n_col, by = nx_strat);
  sum_nh = length(row_ndx) * length(col_ndx) * nh
  var_mean = 0.0;
  ## first loop over 
  for(i in row_ndx) {
    for(j in col_ndx) {
      this_dat = as.vector(y[c((i - (ny_strat - 1)):i), c(j : (j-1))])
      this_area =  as.vector(area_mat[c((i - (ny_strat - 1)):i), c(j : (j-1))])
      var_mean = var_mean + (nh / sum_nh)^2 * var(this_dat / this_area) / length(this_dat) ## no finite population correction
    }
  }
  ## scale to total
  result = survey_area^2 * var_mean
  return(result)
}

#' CorrelationAdjusted 
#' @details a functon that evaluates the variacnce estimator for population total, with an adjustment based on the spatial autocorrelation, based on Morans I 
#' @param y 
#' @param sample_area area of each sampling unit which has to be a single value i.e. sampling unit area are all the same
#' @param survey_area scalar for survey area
#' @export
#' @import raster
#' @return Spatial autocorrealtion variance estimator for the poopulation total
#' 
CorrelationAdjusted = function(y, sample_area, survey_area) {
  if (class(y) != "RasterLayer")
    stop(paste0("y needs to be of class 'RasterLayer', please reformat it"))

  if (class(sample_area) != "numeric" & length(sample_area) != 1) 
    stop("CorrelationAdjusted: expects 'sample_area' to be a scalar, this is not the case, current form assumes equal area sampling units")
  result = 0.0;
  ## D'Orozia adjusted correlation
  rho_hat = Moran(y)
  result = survey_area^2 * var(y$layer@data@values / sample_area) / length(y$layer@data@values)
  if (rho_hat > 0)
    result = result * (1 + 2/log(rho_hat) + 2 / (1/rho_hat - 1))
  
  return(list(var_est = result, rho = rho_hat))
}

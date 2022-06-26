#' CorrelationAdjusted 
#' @details calculate the correlation adjusted variance estimator for the mean from \insertCite{d2003estimating}{systematicsurvey}. This adjusts the simple random sampling estimator using an adjustment based on the spatial autocorrelation, based on Morans I 
#' @param y sample data (assumed to be densities) of class "RasterLayer" see an example for how this is created.
#' @export
#' @importFrom raster Moran
#' @importFrom Rdpack reprompt
#' @return Spatial autocorrealtion variance estimator for the poopulation total
#' @references
#' \insertAllCited{}
#' 
CorrelationAdjusted = function(y, survey_area) {
  if (class(y) != "RasterLayer")
    stop(paste0("y needs to be of class 'RasterLayer', please reformat it"))
  result = 0.0;
  ## D'Orozia adjusted correlation
  rho_hat = Moran(y)
  var_est = var(y$layer@data@values) / length(y$layer@data@values)
  if (rho_hat > 0)
    var_est = var_est * (1 + 2/log(rho_hat) + 2 / (1/rho_hat - 1))
  
  return(list(var_est = var_est, rho = rho_hat))
}

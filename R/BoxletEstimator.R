#' BoxletEstimator 
#' @details
#' a functon that replicates the Boxlet variance estimator from \insertCite{fewster2011variance}{systematicsurvey}
#' Assumes rectangle/square survey region, assumes width and height of sampling units are constant quadrants or transects
#' @param spatial_df spatial data frame, with coordinates for midpoints. Attributes of spatial_df, should have coordinates with coordinate reference, along with width and height of each sampling unit.
#' @param survey_polygon a spatial polygon that is the same coordinate reference as data. 
#' @param quad_width width of every sampling unit
#' @param quad_height height of every sampling unit
#' @param quad_x_spacing distance between sampling units on the x-axis
#' @param quad_y_spacing distance between sampling units on the spatial_df-axis
#' @param boxlet_sampling_frame a list that is the result of the function BoxletEstimatorSamplingFrame() function, if NULL then we build the frame within this function
#' @param boxlet_frame_settings list, if  build_boxlet_sampling_frame = T. then these are the parameters passed to the function see ?BoxletEstimatorSamplingFrame
#' @param trace print messages as function progresses
#' @export
#' @importFrom rgeos gArea
#' @importFrom stats predict
#' @importFrom sp coordinates
#' @importFrom mgcv gam
#' @return boxlet info, and variance estimator for total
#' @references
#' \insertAllCited{}
#' 
BoxletEstimator = function(spatial_df, survey_polygon, quad_width, quad_height, quad_x_spacing, quad_y_spacing, boxlet_sampling_frame = NULL, boxlet_frame_settings = NULL, trace = F) {
  if(trace)
    cat("Entering: BoxletEstimator()\n")
  if(class(spatial_df) != "SpatialPointsDataFrame")
    stop("spatial_df needs to be of class SpatialPointsDataFrame, please check this.")
  col_names_to_check = c("y_i", "area")
  if (!all(col_names_to_check %in% colnames(spatial_df@data))) 
    stop(paste0("spatial_df@data needs to have colnames ", paste0(col_names_to_check, collapse = ", ")))
  ## midpoints of sampling units
  xy_coords = coordinates(spatial_df)
  
  survey_xlim = survey_polygon@bbox[1,]
  survey_ylim = survey_polygon@bbox[2,]
  survey_area =  rgeos::gArea(survey_polygon)
  h_w = quad_width / 2        # half width
  h_h = quad_height / 2       # half height
  sampling_area = quad_width * quad_height * nrow(spatial_df@data)
  if(trace)
    cat(paste0("Survey spatial coverage ", round(sampling_area / survey_area,2)*100, "%\n"))
  
  if(!all(c("x","y") %in% colnames(spatial_df@data))) {
    spatial_df@data$x = xy_coords[,1]
    spatial_df@data$y = xy_coords[,2]
  }
  n = nrow(spatial_df@data)
  if(is.null(boxlet_sampling_frame)) {
    if(trace)
      cat("Building: boxlet_indicator_matrix, this can be slow\n")
    if(is.null(boxlet_frame_settings))
      stop("If you want to build the Boxlet frame, you need to supply a list with labels that are feed to the function BoxletEstimatorSamplingFrame(), see ?BoxletEstimatorSamplingFrame")
    if(class(boxlet_frame_settings) != "list")
      stop("boxlet_frame_settings param needs to be a list")
    boxlet_params = c("boxlet_per_sample_width", "boxlet_per_sample_height")
    if(!all(boxlet_params %in% names(boxlet_frame_settings)))
      stop(paste("the parameter 'boxlet_frame_settings' needs to have elements", paste(boxlet_params, collapse = ", ")))
    boxlet_sampling_frame = BoxletEstimatorSamplingFrame(survey_polygon, quad_width, quad_height,  quad_x_spacing, quad_y_spacing, 
                                                         boxlet_per_sample_width = boxlet_frame_settings$boxlet_per_sample_width, boxlet_per_sample_height = boxlet_frame_settings$boxlet_per_sample_height, trace = trace)
    
  }
  ##  GAM method used in Boxlet
  gam_fit = mgcv::gam(formula = y_i ~ offset(log(area)) + s(x, y, bs = "tp"), data = spatial_df@data, family = poisson(link = log))
  #gam_fit = gam::gam(formula = count ~ s(x_mid,y_mid, df = 20), data = Data, family = poisson(link = log), offset = log(area))
  ## Predict over all boxlets
  boxlet_sampling_frame$boxlet_df$rate  = stats::predict(gam_fit, newdata = boxlet_sampling_frame$boxlet_df, type = "response")
  musum <- sum(boxlet_sampling_frame$boxlet_df$rate)
  density = sum(spatial_df@data$y_i) / sum(spatial_df@data$area)
  N_hat = density * survey_area
  
  ## predict the striplet distribution X
  if(trace) {
    ## plot gam over survey domain
    #ggplot(data = boxlet_sampling_frame$boxlet_df, aes(x = x, spatial_df = spatial_df, fill = rate)) +
    #  geom_tile()
  }
  pj = boxlet_sampling_frame$boxlet_df$rate / sum(boxlet_sampling_frame$boxlet_df$rate)
  Q_hat = as.vector(pj %*% boxlet_sampling_frame$boxlet_indicator_matrix) ## doesn't need to sum = 1
  ## var(n) is variance in observed count, given constant sampling area of units
  B = ncol(boxlet_sampling_frame$boxlet_indicator_matrix)
  var_n_boxlet = (1/B) * sum(N_hat * Q_hat * (1 - Q_hat) + N_hat^2 * Q_hat^2) - (1/B * sum(N_hat * Q_hat))^2
  # B_alt = ncol(boxlet_sampling_frame$boxlet_indicator_matrix)
  # (1/B_alt) * sum(musum * Q_hat_alt * (1 - Q_hat_alt) + musum^2 * Q_hat_alt^2) - (1/B_alt * sum(musum * Q_hat_alt))^2
  ## covnerte var(n) -> survey_area * var(n) / sample_area
  var_total_boxlet = survey_area^2 * var_n_boxlet / (n * quad_width * quad_height)^2
  return(list(var_total_boxlet = var_total_boxlet, boxlet_df = boxlet_sampling_frame$boxlet_df, N_hat = N_hat, gam_N_hat = musum))
}

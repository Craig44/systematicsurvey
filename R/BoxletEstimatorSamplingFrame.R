#' BoxletEstimatorSamplingFrame 
#' @details 
#' an auxillary functon that replicates the Boxlet variance estimator from Fewster 2011
#' Assumes rectangle/square survey region, assumes width and height of sampling units are constant quadrants or transects
#' @param survey_polygon a spatial polygon that is the same coordinate reference as data. 
#' @param quad_width width of every sampling unit
#' @param quad_height height of every sampling unit
#' @param quad_x_spacing distance between sampling units on the x-axis
#' @param quad_y_spacing distance between sampling units on the y-axis
#' @param boxlet_per_sample_width is building 'boxlet_indicator_matrix' how many boxlets are in a sampling unit for x-axis, ignored if boxlet_indicator_matrix supplied
#' @param boxlet_per_sample_height is building 'boxlet_indicator_matrix' how many boxlets are in a sampling unit for y-axis, ignored if boxlet_indicator_matrix supplied
#' @param trace print messages as function progresses
#' @export
#' @import rgeos
#' @import sp
#' @return returns a list containing a fine lattice which represents all the boxlets.
BoxletEstimatorSamplingFrame = function(survey_polygon, quad_width, quad_height, quad_x_spacing, quad_y_spacing, boxlet_per_sample_width = 2, boxlet_per_sample_height = 2, trace = F) {
  if(trace)
    cat("Entering: BoxletEstimatorSamplingFrame()\n")
  if(class(survey_polygon) != "SpatialPolygons")
    error("survey_polygon needs to be of class SpatialPolygons, please check this.")
  survey_xlim = survey_polygon@bbox[1,]
  survey_ylim = survey_polygon@bbox[2,]
  #survey_area =  rgeos::gArea(survey_polygon)
  h_w = quad_width / 2        # half width
  h_h = quad_height / 2       # half height
  boxlet_width = quad_width / boxlet_per_sample_width
  boxlet_height = quad_height / boxlet_per_sample_height
  ## striplet midpoints
  boxlet_x = seq(from = min(survey_xlim) + boxlet_width/2, to = max(survey_xlim) - boxlet_width/2, by = boxlet_width)
  boxlet_y = seq(from = min(survey_ylim) + boxlet_height/2, to = max(survey_ylim) - boxlet_height/2, by = boxlet_height)
  boxlet_df = expand.grid(boxlet_x, boxlet_y)
  colnames(boxlet_df) = c("x", "y")
  J = nrow(boxlet_df)
  boxlet_df$area = boxlet_width * boxlet_height
  ## Now link each striplet to a b in 1:B
  b_x = seq(from = min(survey_xlim) + h_w, to = quad_x_spacing - h_w, by = boxlet_width)
  b_y = seq(from = min(survey_ylim) + h_h, to = quad_y_spacing - h_h, by = boxlet_height)
  bx1 = sort(rep(b_x, length(b_y)))
  by1 = rep(b_y, length(b_x))
  B = length(b_x) * length(b_y)
  if(trace) 
    cat("n-striplets =", J, "(J), sample frames =", B, "(B)\n")
  ## create an effcient function for this
  boxlet_indicator_matrix = matrix(0, nrow = J, ncol = B)
  for(b in 1:B) {
    if (trace) {
      if(b %% 100 == 0)
        cat("b = ", b, "\n")
    }
    ## midpoints for sampling units with this B
    this_b_x = seq(from = bx1[b], to = survey_xlim[2], by = quad_x_spacing)
    this_b_y = seq(from = by1[b], to = survey_ylim[2], by = quad_y_spacing)
    this_b_x1 = rep(this_b_x, length(this_b_y))
    this_b_y1 = sort(rep(this_b_y, length(this_b_x)))
    
    lb_x = this_b_x1 - h_w
    ub_x = this_b_x1 + h_w
    lb_y = this_b_y1 - h_h
    ub_y = this_b_y1 + h_h
    test_ind = rowSums(sapply(1:length(lb_x), function(i) {findInterval(boxlet_df$x, vec = c(lb_x[i],ub_x[i])) == 1 & findInterval(boxlet_df$y, vec = c(lb_y[i],ub_y[i])) == 1}))
    #plot(this_b_x1, this_b_y1, pch = 16)
    #arrow(x0 = this_b_x1, x1 = this_b_x1, y0 = lb_y, y1 = ub_y, col = "red", lwd = 3)
    boxlet_indicator_matrix[which(test_ind == 1), b] = 1
  }
  return(list(boxlet_indicator_matrix = boxlet_indicator_matrix, boxlet_df = boxlet_df))
}

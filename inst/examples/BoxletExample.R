#' Example
#'
library(rgeos)
library(sp)
library(ggplot2)
set.seed(123)
survey_xlim = c(0,100)      # survey range
survey_ylim = c(0,100)      # survey range
quad_width = 1              # Sampling unit width
quad_height = 1             # Sampling unit height
quad_x_spacing = 10         # distance between samples in x-axis
quad_y_spacing = 10         # distance between samples in ys-axis
N = 10000                   # True population number 
# spatial distribution of population true
beta1 = c(4,1)
beta2 = c(1,4)
## create non-overlapping sampling frame
n_col = max(survey_xlim) /  quad_x_spacing # number of sampling units on x-axis
n_row = max(survey_ylim) /  quad_y_spacing # number of sampling units on y-axis
n = n_row * n_col           # n-sampling units
# Define PSU's, random sample these to define survey
n_l = quad_y_spacing / quad_height # l
n_k = quad_x_spacing / quad_width  # k
total_sampling_units = n_l * n_row * n_col * n_k
within_rows = rep(1:n_l, n_k)
within_cols = sort(rep(1:n_k, n_l))
# Randomly draw PSU
first_quad = sample(1:length(within_cols), size = 1)
PSU_x = within_cols[first_quad] * quad_width
PSU_y = within_rows[first_quad] * quad_height
upper_x1 = seq(from = PSU_x, to = survey_xlim[2], by = quad_x_spacing)
upper_y1 = seq(from = PSU_y, to = survey_ylim[2], by = quad_y_spacing)
upper_x = rep(upper_x1, length(upper_y1))
upper_y = sort(rep(upper_y1, length(upper_x1)))
y_coord = upper_y - quad_height / 2
x_coord = upper_x - quad_width / 2
plot(x = x_coord, y = y_coord, xlim = survey_xlim, ylim = survey_ylim, xlab = "", ylab = "", pch = 16)
## generate True population
pop_x_coord = rbeta(N, beta1[1], beta1[2]) * max(survey_xlim)
pop_y_coord = rbeta(N, beta2[1], beta2[2]) * max(survey_ylim)
points(pop_x_coord, pop_y_coord, pch = 16, cex = 0.7, col = adjustcolor(col = "blue", alpha.f = 0.1))
## sample pop based on systematic survey
Data = data.frame(height = rep(quad_height, n), width = rep(quad_width, n))
Data$count = NA
counter = 1;
#plot(xpp)
for(i in 1:length(upper_x)) {
  #polygon(x = c(upper_x[i] - quad_width, upper_x[i] - quad_width,upper_x[i] ,upper_x[i] ), y = c(upper_y[i]  - quad_height,upper_y[i], upper_y[i],upper_y[i]  - quad_height), col = "blue")
  Data$count[i] = sum(pop_x_coord >= (upper_x[i] - quad_width) & pop_x_coord <= upper_x[i]  & pop_y_coord >= (upper_y[i]  - quad_height) & pop_y_coord <= upper_y[i])
}
Data$area = Data$height * Data$width
head(Data)
## convert to spatial data frame
sp_data = SpatialPointsDataFrame(coords = cbind(x_coord, y_coord), data = Data) ## can add coordinate reference system here as well
## Create survey polgon
s_poly = Polygon(cbind(x = c(survey_xlim[1],survey_xlim[1],survey_xlim[2],survey_xlim[2]), y = c(survey_ylim[1],survey_ylim[2],survey_ylim[2],survey_ylim[1])))
s_polys = Polygons(list(s_poly),1)
sp_survey_poly = SpatialPolygons(list(s_polys)) ## add coordinater reference here
survey_area = rgeos::gArea(sp_survey_poly)
survey_polygon = sp_survey_poly
## debug Boxlet function
y = sp_data
xy_coords = coordinates(y)

survey_polygon = sp_survey_poly
dimensions = 2
## Striplet sampling data frame
quad_width = 1              # Sampling unit width
quad_height = 1             # Sampling unit height
h_w = quad_width / 2        # half width
h_h = quad_height / 2       # half height

quad_x_spacing = 10         # distance between samples in x-axis
quad_y_spacing = 10         # distance between samples in ys-axis
## 
boxlet_per_sample_width = 2
boxlet_per_sample_height = 2
strip_width = quad_width / boxlet_per_sample_width
strip_height = quad_height / boxlet_per_sample_height
## striplet midpoints
strip_x = seq(from = min(survey_xlim) + h_w, to = max(survey_xlim) - h_w, by = strip_width)
strip_y = seq(from = min(survey_ylim) + h_h, to = max(survey_ylim) - h_h, by = strip_height)
strip_df = expand.grid(strip_x, strip_y)
colnames(strip_df) = c("x", "y")
J = nrow(strip_df)
strip_df$area = strip_width * strip_height
## Now link each striplet to a b in 1:B
b_x = seq(from = min(survey_xlim) + h_w, to = quad_x_spacing - h_w, by = strip_width)
b_y = seq(from = min(survey_ylim) + h_h, to = quad_y_spacing - h_h, by = strip_height)
bx1 = sort(rep(b_x, length(b_y)))
by1 = rep(b_y, length(b_x))
B = length(b_x) * length(b_y)
B 
## create an indicator matrix that is J x B
## create an effcient function for this
indicator = matrix(0, nrow = J, ncol = B)
for(b in 1:B) {
  if(b %% 10 == 0)
    cat("b = ", b, "\n")
  ## midpoints for sampling units with this B
  this_b_x = seq(from = bx1[b], to = survey_xlim[2], by = quad_x_spacing)
  this_b_y = seq(from = by1[b], to = survey_ylim[2], by = quad_y_spacing)
  this_b_x1 = rep(this_b_x, length(this_b_y))
  this_b_y1 = sort(rep(this_b_y, length(this_b_x)))
  
  lb_x = this_b_x1 - quad_width / 2
  ub_x = this_b_x1 + quad_width / 2
  lb_y = this_b_y1 - quad_height / 2
  ub_y = this_b_y1 + quad_height / 2
  test_ind = rowSums(sapply(1:length(lb_x), function(i) {findInterval(strip_df$x, vec = c(lb_x[i],ub_x[i])) == 1 & findInterval(strip_df$y, vec = c(lb_y[i],ub_y[i])) == 1}))
  #plot(this_b_x1, this_b_y1, pch = 16)
  #arrow(x0 = this_b_x1, x1 = this_b_x1, y0 = lb_y, y1 = ub_y, col = "red", lwd = 3)
  indicator[which(test_ind == 1), b] = 1
}
table(colSums(indicator))
## check this worked
#plot(this_b_x1, this_b_y1, pch = 16, type = "n")
#for(b in 1:B) 
#  points(strip_df$x[indicator[,b] == 1], strip_df$y[indicator[,b] == 1], col = adjustcolor(col = "blue", alpha.f = 0.2), cex = 0.3, pch = 16)
if(all(c("x","y") %in% colnames(y@data))) {
  y@data$x = xy_coords[,1]
  y@data$y = xy_coords[,2]
}
boxlet_frame = BoxletEstimatorSamplingFrame(survey_polygon, quad_width, quad_height, 
                                            boxlet_per_sample_width = 2, boxlet_per_sample_height = 2, trace = F)
names(boxlet_frame)
nrow(boxlet_frame$boxlet_df)
dim(boxlet_frame$boxlet_indicator_matrix)
table(colSums(boxlet_frame$boxlet_indicator_matrix))
table(colSums(indicator))

##  GAM method used in Boxlet
gam_fit = mgcv::gam(formula = count ~ offset(log(area)) + s(x, y, bs = "tp"), data = y, family = poisson(link = log))
#gam_fit = gam::gam(formula = count ~ s(x_mid,y_mid, df = 20), data = Data, family = poisson(link = log), offset = log(area))
## Predict over all boxlets
strip_df$rate  = predict(gam_fit, newdata = strip_df, type = "response")
musum <- sum(strip_df$rate)
## predict the striplet distribution X
ggplot(data = strip_df, aes(x = x, y = y, fill = rate)) +
  geom_tile()

## Variance estimator
Abvec = as.vector(strip_df$rate %*% indicator)
pj = strip_df$rate / sum(strip_df$rate)
Q_hat = as.vector(pj %*% indicator) ## doesn't need to sum = 1
strip_area = as.vector(strip_df$area %*% indicator)
Lbvec = rep(100, B)

boxplot(cbind(ABvec, Q_hat * musum))
abline(h = sum(Data$count), col = "red", lwd = 2, lty = 2)

## var(n) is variance in observed count, given constant sampling area of units
var_boxlet = (1/B) * sum(musum * Q_hat * (1 - Q_hat) + musum^2 * Q_hat^2) - (1/B * sum(musum * Q_hat))^2
## covnerte var(n) -> survey_area * var(n) / sample_area
sqrt(survey_area^2 * var_boxlet / (n * quad_width * quad_height)^2)
## alternative, if to account for varying sampling areas var(n/a) - here a is a variable
var_boxlet_alt = 1/B * sum((musum * Q_hat * (1 - Q_hat) + musum^2 * Q_hat^2) / strip_area^2) - ((1/B) * sum((musum * Q_hat) / strip_area))^2
sqrt(survey_area^2 * var_boxlet_alt)
## Code from Fewster
var.striplet <- 1/B * sum((Abvec + Abvec^2*(1-1/musum) )/Lbvec^2)- (1/B * sum(Abvec/Lbvec))^2
sqrt(survey_area^2 * var.striplet / quad_width^2)















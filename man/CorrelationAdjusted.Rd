% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CorrelationAdjusted.R
\name{CorrelationAdjusted}
\alias{CorrelationAdjusted}
\title{CorrelationAdjusted}
\usage{
CorrelationAdjusted(y, survey_area)
}
\arguments{
\item{y}{sample data (assumed to be densities) of class "RasterLayer" see an example for how this is created.}
}
\value{
Spatial autocorrealtion variance estimator for the poopulation total
}
\description{
CorrelationAdjusted
}
\details{
calculate the correlation adjusted variance estimator for the mean from \insertCite{d2003estimating}{systematicsurvey}. This adjusts the simple random sampling estimator using an adjustment based on the spatial autocorrelation, based on Morans I
}
\references{
\insertAllCited{}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/app_old2.R
\name{ExtSpecR_app}
\alias{ExtSpecR_app}
\title{ExtSpecR_app()}
\usage{
ExtSpecR_app()
}
\arguments{
\item{...}{plese found the information on our website at https://github.com/Yanjie-Li/ExtSpecR.}
}
\value{
The final data contains essential information such as
the longitude (x), latitude (y), tree ID, point height (Z),
crown area (area), and spectral data in rds data frame format.
}
\description{
This is a shiny app specially used to extract spectral information of a single tree,
and I also provide a sample to show what the data looks like after extraction.
First of all, you need point cloud data with precise positioning information,
which is used to segment each individual plant in a large area of forest land.
In addition, you need to have multi-spectral or hyperspectral information of this forest land.
}
\examples{
library(ExtSpecR)
ExtSpecR::ExtSpecR_app()

}

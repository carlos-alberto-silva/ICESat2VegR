% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICESat2_finder.R
\name{ICESat2_finder}
\alias{ICESat2_finder}
\title{ICESat-2 ATL03 and ATL08 data finder}
\usage{
ICESat2_finder(
  short_name,
  lower_left_lon,
  lower_left_lat,
  upper_right_lon,
  upper_right_lat,
  version = "006",
  daterange = NULL,
  cloud_hosted = TRUE
)
}
\arguments{
\item{short_name}{ICESat-2 ATLAS data level short_name; Options: "ATL03", "ATL08",}

\item{lower_left_lon}{Numeric. Minimum longitude in
(decimal degrees) for the bounding box of the area of interest.}

\item{lower_left_lat}{Numeric. Minimum latitude in
(decimal degrees) for the bounding box of the area of interest.}

\item{upper_right_lon}{Numeric. Maximum longitude in lon
(decimal degrees) for the bounding box of the area of interest.}

\item{upper_right_lat}{Numeric. Maximum latitude in
(decimal degrees) for the bounding box of the area of interest.}

\item{version}{Character. The version of the ICESat-2 ATLAS product files to be
returned (only V005 or V006). Default "006".}

\item{daterange}{Vector. Date range. Specify your start and end dates
using ISO 8601 [YYYY]-[MM]-[DD]T[hh]:[mm]:[ss]Z. Ex.:
c("2019-07-01T00:00:00Z","2020-05-22T23:59:59Z"). If NULL (default),
the date range filter will be not applied.}

\item{cloud_hosted}{Logical. Flag to indicate use of cloud hosted collections.}
}
\value{
Return a vector object pointing out the path saving the downloaded
ICESat-2 ATLAS data within the boundary box coordinates provided
}
\description{
This function finds the exact granule(s) that contain ICESat-2 ATLAS data
for a given region of interest and date range
}
\examples{
\donttest{
# ICESat-2 data finder is a web service provided by NASA
# usually the request takes more than 5 seconds

# Specifying bounding box coordinates
lower_left_lon <- -96.0
lower_left_lat <- 40.0
upper_right_lon <- -100
upper_right_lat <- 42.0

# Specifying the date range
daterange <- c("2019-07-01", "2020-05-22")

# Extracting the path to ICESat-2 ATLAS data for the specified boundary box coordinates
ICESat-2 ATLAS02b_list <- ICESat2_finder(
  short_name = "ATL08",
  ul_lat,
  ul_lon,
  lr_lat,
  lr_lon,
  version = "006",
  daterange = daterange
)
}
}
\seealso{
bbox: Defined by the upper left and lower right corner coordinates,
in lat,lon ordering, for the bounding box of the area of interest
(e.g. lower_left_lon,lower_left_lat,upper_right_lon,upper_right_lat)

This function relies on the existing CMR tool:
\url{https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html}
}

#' @include class.icesat2.R
PAGE_SIZE <- 2000


# #' ICESat-2 ATL03 and ATL08 data finder for cloud computing
# #'
# #' @description This function finds the exact granule(s) that contain ICESat-2 ATLAS data
# #' for a given region of interest and date range for cloud computing. For better processing efficiency, this function should be used
# #' on the cluster hosting ICESat-2 data (e.g. CryoCloud)
# #'
# #' @param short_name ICESat-2 ATLAS data level short_name; Options: "ATL03", "ATL08",
# #' @param lower_left_lon Numeric. Minimum longitude in
# #' (decimal degrees) for the bounding box of the area of interest.
# #' @param lower_left_lat Numeric. Minimum latitude in
# #' (decimal degrees) for the bounding box of the area of interest.
# #' @param upper_right_lon Numeric. Maximum longitude in lon
# #' (decimal degrees) for the bounding box of the area of interest.
# #' @param upper_right_lat Numeric. Maximum latitude in
# #' (decimal degrees) for the bounding box of the area of interest.
# #' @param version Character. The version of the ICESat-2 ATLAS product files to be
# #' returned (only V005 or V006). Default "006".
# #' @param daterange Vector. Date range. Specify your start and end dates
# #' using ISO 8601 \[YYYY\]-\[MM\]-\[DD\]T\[hh\]:\[mm\]:\[ss\]Z. Ex.:
# #' c("2019-07-01T00:00:00Z","2020-05-22T23:59:59Z"). If NULL (default),
# #' the date range filter will be not applied.
# #' @param cloud_hosted Logical. Flag to indicate use of cloud hosted collections.
# #'
# #' @return Return a vector object pointing out the path saving the downloaded
# #' ICESat-2 ATLAS data within the boundary box coordinates provided
# #'
# #' @seealso bbox: Defined by the upper left and lower right corner coordinates,
# #' in lat,lon ordering, for the bounding box of the area of interest
# #' (e.g. lower_left_lon,lower_left_lat,upper_right_lon,upper_right_lat)
# #'
# #' This function relies on the existing CMR tool:
# #' \url{https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html}
# #'
# #' @examples
# #' \donttest{
# #' # ICESat-2 data finder is a web service provided by NASA
# #' # usually the request takes more than 5 seconds
# #'
# #' # Specifying bounding box coordinates
# #' lower_left_lon <- -85
# #' lower_left_lat <- 30.0
# #' upper_right_lon <- -84.0
# #' upper_right_lat <- 31.0
# #'
# #' # Specifying the date range
# #' daterange <- c("2019-07-01", "2020-05-22")
# #'
# #' # Extracting the path to ICESat-2 ATLAS data for the specified boundary box coordinates
# #' \dontrun{
# #' ATLAS02b_list <- ATLAS_dataFinder(
# #'   short_name = "ATL08",
# #'   upper_right_lat,
# #'   lower_left_lon,
# #'   lower_left_lat,
# #'   upper_right_lon,
# #'   version = "006",
# #'   daterange = daterange,
# #'   cloud_computing = TRUE
# #' )
#' 
# #' }
# #' @import jsonlite curl magrittr reticulate
#' @include earthaccess.R
ATLAS_dataFinder_cloud <- function(short_name,
                                 lower_left_lon,
                                 lower_left_lat,
                                 upper_right_lon,
                                 upper_right_lat,
                                 version = "007",
                                 daterange = NULL,
                                 persist = TRUE) {
  if (!reticulate::py_module_available("earthaccess")) {
    tryCatch(expr = {
      reticulate::py_install("earthaccess")
    }, error = function(e) {
      stop("earthaccess module is not available and could not be installed,
      try installing it manually with conda or pip or changing the python environment.")
    })
  }


  if (inherits(earthaccess, "logical")) {
    earthaccess <- reticulate::import("earthaccess", convert = FALSE)
  }

  if (is.null(earthaccess) || is.null(h5py)) {
    tryInitializeCloudCapabilities()
  }

  # Try to refresh login
  earthaccess_login(persist)

  granule_query <- earthaccess$granule_query()
  granule_query$bounding_box(lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)
  if (length(daterange) == 2) {
    granule_query$temporal(daterange[[1]], daterange[[2]])
  }
  granule_query$cloud_hosted(TRUE)
  granule_query$version(version)
  granule_query$short_name(short_name)
  res <- granule_query$get_all()

  res <- new("icesat2.granules_cloud", granules = res)
  return(res)
}

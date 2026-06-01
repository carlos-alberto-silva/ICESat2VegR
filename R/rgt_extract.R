#' Extract reference ground track from ATL03 segments
#'
#' @description
#' Extracts the reference ground track from ICESat-2 ATL03 or ATL08 HDF5
#' data as either a centerline or a bounding polygon, returned as a
#' \code{\link[terra]{SpatVector}} object with orbit metadata attributes.
#'
#' @param h5 An ICESat-2 HDF5 object, output of \code{\link{ATL03_read}}
#'   or \code{\link{ATL08_read}}.
#' @param line logical. If \code{TRUE} (default), extracts the ground track
#'   centerline as a \code{lines} geometry. If \code{FALSE}, extracts the
#'   bounding polygon as a \code{polygons} geometry.
#'
#' @return A \code{\link[terra]{SpatVector}} object in \code{EPSG:4326}
#'   containing the ground track geometry with the following attributes:
#'   \describe{
#'     \item{cycle_number}{ICESat-2 cycle number.}
#'     \item{orbit_number}{Orbit number.}
#'     \item{rgt}{Reference Ground Track number.}
#'   }
#'
#' @details
#' When \code{line = TRUE}, the centerline is computed as the mean of the
#' two sides of the bounding polygon derived from
#' \code{orbit_info/bounding_polygon_lon1} and
#' \code{orbit_info/bounding_polygon_lat1}.
#'
#' When \code{line = FALSE}, the full bounding polygon is returned directly
#' from the \code{orbit_info} group.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#' # Specifying the path to ATL03 H5 file
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Extracting ATL03 RGT track as polygon
#' rgt_polygon <- rgt_extract(h5 = atl03_h5, line = FALSE)
#' head(rgt_polygon)
#' terra::plet(rgt_polygon)
#'
#' # Extracting ATL03 RGT track as centerline
#' rgt_line <- rgt_extract(h5 = atl03_h5, line = TRUE)
#' head(rgt_line)
#' terra::plet(rgt_line)
#'
#' close(atl03_h5)
#'
#' @include class.icesat2.R
#' @importFrom data.table data.table rbindlist
#' @importFrom terra writeVector
#' @export
rgt_extract <- function(h5, line = TRUE) {
  stopifnot(
    "The input file needs to be one ATL03_read or ATL08_read products" = inherits(h5, "icesat2.h5")
  )
  rgt_data <- cbind(
    cycle_number = h5[["orbit_info/cycle_number"]][],
    orbit_number = h5[["orbit_info/orbit_number"]][],
    rgt = h5[["orbit_info/rgt"]][]
  )
  rgt_coords <- cbind(
    h5[["orbit_info/bounding_polygon_lon1"]][],
    h5[["orbit_info/bounding_polygon_lat1"]][]
  )

  if (line == TRUE) {
    n     <- nrow(rgt_coords)
    half  <- floor(n / 2)
    side1 <- rgt_coords[1:half, ]
    side2 <- rgt_coords[(n - 1):(half + 1), ]
    mid   <- (side1 + side2) / 2
    rgt   <- terra::vect(
      mid,
      type = "lines",
      atts = as.data.frame(rgt_data),
      crs  = "EPSG:4326"
    )
  } else {
    rgt <- terra::vect(
      rgt_coords,
      type = "polygons",
      atts = as.data.frame(rgt_data),
      crs  = "EPSG:4326"
    )
  }
  return(rgt)
}

#' Extract Reference Ground Track from ICESat-2 ATL03 or ATL08 Data
#'
#' @description
#' Extracts the reference ground track from ICESat-2 ATL03 or ATL08 data and
#' returns it as a \code{terra::SpatVector}. Optionally, the extracted ground
#' track can also be written to disk, for example as a KML file.
#'
#' @param h5 An ICESat-2 ATL03 or ATL08 object, usually the output of
#' \code{ATL03_read()} or \code{ATL08_read()}.
#' @param line Logical. If \code{TRUE}, returns the reference ground track as a
#' line. If \code{FALSE}, returns the bounding polygon. Default is \code{TRUE}.
#' @param output Character or \code{NULL}. Optional output vector file path.
#' If provided, the extracted ground track is written using
#' \code{terra::writeVector()}. Use a file extension supported by GDAL,
#' such as \code{.kml}, \code{.gpkg}, or \code{.shp}.
#'
#' @return
#' A \code{terra::SpatVector} containing the extracted reference ground track.
#'
#' @details
#' This function uses ICESat-2 orbit information stored in the \code{orbit_info}
#' group to derive the reference ground track. When \code{line = TRUE}, the
#' function estimates the centerline from the bounding polygon coordinates.
#'
#' The returned object can be passed directly to
#' \code{plot_icesat2_orbit_animation()} using the \code{rgt} argument.
#'
#' @seealso
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#' \dontrun{
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' rgt_line <- rgt_extract(h5 = atl03_h5, line = TRUE)
#'
#' plot_icesat2_orbit_animation(
#'   rgt = rgt_line,
#'   launch = TRUE
#' )
#'
#' rgt_extract(
#'   h5 = atl03_h5,
#'   line = TRUE,
#'   output = file.path(tempdir(), "atl03_rgt.kml")
#' )
#'
#' close(atl03_h5)
#' }
#'
#' @include class.icesat2.R
#' @importFrom terra vect writeVector
#' @export
rgt_extract <- function(h5, line = TRUE, output = NULL) {

  stopifnot(
    "The input file needs to be an ATL03_read or ATL08_read product" =
      inherits(h5, "icesat2.h5")
  )

  if (!is.logical(line) || length(line) != 1) {
    stop("line must be TRUE or FALSE.")
  }

  if (!is.null(output) && (!is.character(output) || length(output) != 1)) {
    stop("output must be NULL or a single character file path.")
  }

  rgt_data <- data.frame(
    cycle_number = h5[["orbit_info/cycle_number"]][],
    orbit_number = h5[["orbit_info/orbit_number"]][],
    rgt = h5[["orbit_info/rgt"]][]
  )

  rgt_coords <- cbind(
    h5[["orbit_info/bounding_polygon_lon1"]][],
    h5[["orbit_info/bounding_polygon_lat1"]][]
  )

  rgt_coords <- stats::na.omit(rgt_coords)

  if (nrow(rgt_coords) < 2) {
    stop("Not enough orbit coordinates were found to create a ground track.")
  }

  if (line) {

    n <- nrow(rgt_coords)
    half <- floor(n / 2)

    if (half < 1) {
      stop("Not enough coordinates were found to create a line.")
    }

    side1 <- rgt_coords[seq_len(half), , drop = FALSE]
    side2 <- rgt_coords[(n - 1):(half + 1), , drop = FALSE]

    min_n <- min(nrow(side1), nrow(side2))

    side1 <- side1[seq_len(min_n), , drop = FALSE]
    side2 <- side2[seq_len(min_n), , drop = FALSE]

    mid <- (side1 + side2) / 2

    rgt <- terra::vect(
      mid,
      type = "lines",
      atts = rgt_data,
      crs = "epsg:4326"
    )

  } else {

    rgt <- terra::vect(
      rgt_coords,
      type = "polygons",
      atts = rgt_data,
      crs = "epsg:4326"
    )
  }

  if (!is.null(output)) {
    terra::writeVector(
      rgt,
      output,
      overwrite = TRUE
    )
  }

  return(rgt)
}

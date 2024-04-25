#' Class that represent custom segments created from ATL03 and ATL08 joined data
#'
#' @details
#' This class is actually just a wrap around the [`data.table-class`], but it indicates
#' the output from [`ATL03_ATL08_segment_create()`], which means the dataset will contain
#' the needed structure for computing value for the computing the stats with
#' [`ATL03_ATL08_compute_seg_attributes_dt_segStat()`]

setRefClass("icesat2.atl03_atl08_seg_dt")

#' Compute segments id for a given segment length
#'
#' @description This function reads the ICESat-2 Land and
#' Vegetation Along-Track Products (ATL08) as h5 file.
#'
#' @param atl03_atl08_dt [`ICESat2VegR::icesat2.atl03atl08_dt-class`].
#' The output of the [`ATL03_ATL08_photons_attributes_dt_join()`].
#' @param segment_length [`numeric-class`]. The desired segment length to split the photons.
#' @param centroid character. Method used to calculate the segment centroid, either "mean" or "midpoint",
#' see details. Default 'mean'.
#' @param output Character vector. The output vector file. The GDAL vector format
#' will be inferred by the file extension using [`terra::writeVector()`]
#' @param overwrite logical input to control if the output vector file should be overwritten.
#' Default FALSE.
#'
#' @details
#' The centroid will be computed using either the photons centroid or the approximate segment centroid.
#' - "mean": calculated using the average coordinates from all photons within the segment.
#' This approach will better represent the mean statistics location.
#' - "mid-point": the minimum and maximum coordinates will be averaged to calculate a midpoint within
#' the segment. This will give a better representation of the segment true mid-point
#'
#' @return Returns an S4 object of class [`icesat2.atl03atl08_dt-class`] containing ICESat-2 ATL08 data.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
# Specifying the path to ICESat-2 ATL03 and ATL08 data
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ICESat-2 ATL08 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#'
#' atl03_atl08_dt_seg <- ATL03_ATL08_segment_create(atl03_atl08_dt,
#'   segment_length = 30,
#'   centroid = "mean",
#'   output = NA,
#'   overwrite = FALSE
#' )
#'
#' head(atl03_atl08_dt_seg)
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' @import hdf5r
#' @export
ATL03_ATL08_segment_create <- function(
    atl03_atl08_dt,
    segment_length,
    centroid = "mean",
    output = NA,
    overwrite = FALSE) {
  stopifnot(
    "Object atl03_atl08_dt is not from the icesat2.atl03atl08_dt class" =
      inherits(atl03_atl08_dt, "icesat2.atl03atl08_dt")
  )

  `:=` <- data.table::`:=`

  beam <- lat_ph <- lon_ph <- dist_ph_along <-
    centroid_x <- centroid_y <- segment_id <- NA

  suppressWarnings(atl03_atl08_dt[
    ,
    segment_id := floor(dist_ph_along / segment_length) + 1,
  ])

  dt <- atl03_atl08_dt

  if (centroid == "mean") {
    dt[
      ,
      c(
        "centroid_x",
        "centroid_y"
      ) := list(
        mean(lon_ph),
        mean(lat_ph)
      ),
      by = list(segment_id, beam)
    ]
  } else if (centroid == "midpoint") {
    dt[
      ,
      c(
        "centroid_x",
        "centroid_y"
      ) := list(
        x = mean(range(lon_ph)),
        y = mean(range(lat_ph))
      ),
      by = list(segment_id, beam)
    ]
  } else {
    stop("The centroid extraction method '%s', was not found!\nUse either 'mean' or 'midpoint'.")
  }


  if (is.na(output) == FALSE) {
    centroid <- as.matrix(
      dt[, list(object_id = segment_id, x = centroid_x, y = centroid_y)]
    )

    v <- terra::vect(
      centroid,
      type = "points",
      crs = "epsg:4326"
    )

    terra::writeVector(v, output, overwrite = TRUE)
  }

  prepend_class(dt, "icesat2.atl03_atl08_seg_dt")
  dt
}

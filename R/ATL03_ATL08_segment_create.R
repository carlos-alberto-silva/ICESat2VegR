setRefClass("icesat2.atl03_atl08_seg_dt")

#' Compute segments id for a given segment length
#'
#' @description This function reads the ICESat-2 Land and
#' Vegetation Along-Track Products (ATL08) as h5 file.
#'
#'
#' @usage ATL08_read(atl08_path)
#'
#' @param atl03_atl08_dt [`icesat2.atl03_atl08_dt-class`]. The output of the [`ATL03_ATL08_photons_attributes_dt_join()`].
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
#' @return Returns an S4 object of class ["icesat2.atl03_atl08_dt-class"] containing ICESat-2 ATL08 data.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
#' # Specifying the path to ICESat-2 ATL08 data (zip file)
#' outdir <- tempdir()
#' atl08_fp_zip <- system.file("extdata",
#'   "ATL0802_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ICESat-2 ATL08 data
#' atl08_path <- unzip(atl08_fp_zip, exdir = outdir)
#'
#' # Reading ICESat-2 ATL08 data (h5 file)
#' atl08 <- ATL08_read(atl08_path = atl08_path)
#'
#' close(atl08)
#' @import hdf5r
#' @export
ATL03_ATL08_segment_create <- function(atl03_atl08_dt, segment_length, centroid = "mean", output = NA, overwrite = FALSE) {
  stopifnot(
    "Object atl03_atl08_dt is not from the icesat2.atl03atl08_dt class" =
      inherits(atl03_atl08_dt, "icesat2.atl03atl08_dt")
  )

  .SD <- data.table::.SD
  `:=` <- data.table::`:=`

  lat_ph <- lon_ph <- dist_ph_along <-
    centroid_x <- centroid_y <- segment_id <- NA

  dt <- atl03_atl08_dt[
    ,
    .SD,
    by = list(segment_id = floor(dist_ph_along / segment_length) + 1), .SDcols = names(atl03_atl08_dt)
  ]

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
      by = segment_id
    ]
  } else if (centroid == "midpoint") {
    dt[
      ,
      c(
        "centroid_x",
        "centroid_y"
      ) := list(
        x = mean_red(range(lon_ph)),
        y = mean_red(range(lat_ph))
      ),
      by = segment_id
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

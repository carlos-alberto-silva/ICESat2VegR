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
#' @param centroid character. Method used to calculate the segment centroid, either "photons" or "segment",
#' see details.
#' @param output Character vector. The output vector file. The GDAL vector format
#' will be inferred by the file extension using [`terra::writeVector()`]
#' @param overwrite logical input to control if the output vector file should be overwritten.
#' Default FALSE.
#'
#' @details
#' The centroid will be computed using either the photons centroid or the approximate segment centroid.
#' - Photons centroid: calculated using the average coordinates from all photons within the segment.
#' - Approximate segment centroid: after calculating the centroid from the photons we will use the average
#' photon distance from the beginning of the segment to derive the actual position of the calculated centroid
#' and then adjust it accordingly using the neighbor centroids as nodes for calculating the direction of the
#' offset.
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
#'   package = "rICESat2Veg"
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
ATL03_ATL08_segment_create <- function(atl03_atl08_dt, segment_length, real_centroid = TRUE, output = NA, overwrite = FALSE) {
  stopifnot(
    "Object atl03_atl08_dt is not from the icesat2.atl03atl08_dt class" =
      inherits(atl03_atl08_dt, "icesat2.atl03atl08_dt")
  )

  .SD <- data.table::.SD
  .N <- data.table::.N

  lat_ph <- lon_ph <- dist_ph_along <- NA
  dt2 <- atl03_atl08_dt[, .SD, by = list(segment_id = floor(dist_ph_along / segment_length) + 1), .SDcols = names(atl03_atl08_dt)]

  centroid <- dt2[
    ,
    list(
      y = mean(lat_ph),
      x = mean(lon_ph)
    ),
    by = segment_id
  ]

  if (is.na(output) == FALSE) {
    centroid <- as.matrix(
      centroid[, list(object_id = segment_id, x, y)]
    )

    v <- terra::vect(
      lines,
      type = "points",
      crs = "epsg:4326"
    )

    terra::writeVector(v, output, overwrite = T)
  }

  data.table::setattr(dt, "class", c("icesat2.atl03_atl08_seg_dt", "data.table", "data.frame"))
  dt
}

#' Clips ICESat-2 ATL03 H5 data
#'
#' @param atl03 [`icesat2.atl03_h5-class`] object,
#' obtained through [`ATL03_read()`] for clipping
#' @param output character. Path to the output h5 file, the attribute for polygons
#' will be appended to the file name.
#' @param vect [`terra::SpatVector-class`] for clipping
#' @param polygon_id [`character-class`]. The attribute name used for identifying
#' the different polygons. Default is "id"
#'
#' @return Returns a list of clipped S4 object of class [`icesat2.atl03_h5-class`]
#'
#' @description This function clips ATL03 HDF5 file within beam groups,
#' but keeps metada and ancillary data the same.
#'
#' @examples
##' # Specifying the path to ATL03 file (zip file)
#' outdir <- tempdir()
#' atl03_zip <- system.file("extdata",
#'   "atl03_20220401221822_01501506_005_01.zip",
#'   package = "rICESat2Veg"
#' )
#'
#' # Unzipping ATL03 file
#' atl03_path <- unzip(atl03_zip, exdir = outdir)
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- atl03_read(atl03_path = atl03_path)
#'
#' output <- file.path(outdir, "clipped.h5")
#'
#' vect_path <- system.file("extdata",
#'   "polygons.shp",
#'   package = "rICESat2Veg"
#' )
#'
#' vect <- terra::vect()
#'
#' # Clipping ATL03 photons  by boundary box extent
#' atl03_photons_dt_clip <- ATL03_h5_clipGeometry(
#'   atl03_h5,
#'   output,
#'   vect,
#'   polygon_id = "id"
#' )
#'
#' close(atl03_h5)
#' @import hdf5r
#' @include clipTools.R
#' @export
ATL03_h5_clipGeometry <- function(atl03, output, vect, polygon_id = "id") {
  outputs <- list()

  n_polygons <- nrow(vect)
  n_places <- floor(log10(n_polygons))
  for (geom_idx in seq_along(vect)) {
    current_places <- floor(log10(geom_idx))
    message("=============================", appendLF = FALSE)
    message(paste0(rep("=", current_places + n_places), collapse = ""))
    message(sprintf("== Processing geometry %s/%s ==", geom_idx, n_polygons))
    message("=============================", appendLF = FALSE)
    message(paste0(rep("=", current_places + n_places), collapse = ""))
    geom <- vect[geom_idx]
    sub_output <- gsub(".h5$", sprintf("_%s.h5", geom[[polygon_id]][[1]]), output)

    outputs[[geom[[polygon_id]][[1]]]] <-
      ATL03_h5_clip(atl03, sub_output, geom = geom, ATL03_photons_mask_fn = ATL03_photons_mask_geometry)
  }
}

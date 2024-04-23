#' Generic function for clipping ICESat-2 ATL03 and ATL08 H5 data
#'
#' @description This function clips ATL03 and ATL08 HDF5 file within beam groups,
#' but keeps metada and ancillary data the same.
#'
#' @param x [`icesat2.atl03_h5-class`] or [`icesat2.atl08_h5-class`] object,
#' obtained through [`ATL03_read()`] or [`ATL08_read()`] for clipping
#' @param output character. Path to the output h5 file.
#' @param clip_obj [`numeric-class`], [`terra::SpatExtent`] or [`terra::SpatVector`]
#' object for clipping. The bbox order is the default from NASA's ICESat-2 data searching:
#' [ul_lat, ul_lon, lr_lat, lr_lon].
#' @param ... Forwards to specific method implementation.
#'
#' @return Returns the clipped S4 object of the same class as the input file
#'
#'
#' @seealso
#' [`ATL03_h5_clipBox()`], [`ATL03_h5_clipGeometry()`],
#' [`ATL08_h5_clipBox()`], [`ATL08_h5_clipGeometry()`]
#'
#' @examples
##' # Specifying the path to ATL03 file (zip file)
#' outdir <- tempdir()
#' atl03_zip <- system.file("extdata",
#'   "atl03_20220401221822_01501506_005_01.zip",
#'   package = "ICESat2VegR"
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
#'   package = "ICESat2VegR"
#' )
#'
#' vect <- terra::vect(vect_path)
#' ext <- terra::ext(vect)
#'
#' # Clipping ATL03 photons by boundary box extent
#' atl03_clip <- clip(
#'   atl03_h5,
#'   output,
#'   ext
#' )
#'
#' close(atl03_clip)
#'
#' # Clipping ATL03 photons by geometry
#' atl03_clip_geom <- clip(
#'   atl03_h5,
#'   output,
#'   vect,
#'   polygon_id = "id"
#' )
#'
#' lapply(atl03_clip_geom, close)
#' @import hdf5r
#' @export
setGeneric(
  "clip",
  function(x, output, clip_obj, ...) {
    standardGeneric("clip")
  }
)

#' @include class.icesat2.R
#' @importClassesFrom terra SpatExtent
#' @exportMethod clip
setMethod(
  "clip",
  signature = c("icesat2.atl03_h5", "character", "SpatExtent"),
  function(x, output, clip_obj, ...) {
    ATL03_h5_clipBox(x, output, clip_obj, ...)
  }
)


#' @include class.icesat2.R
#' @importClassesFrom terra SpatVector
#' @exportMethod clip
setMethod(
  "clip",
  signature = c("icesat2.atl03_h5", "character", "SpatVector"),
  function(x, output, clip_obj, polygon_id = "id", ...) {
    ATL03_h5_clipGeometry(x, output, clip_obj, polygon_id, ...)
  }
)


#' @include class.icesat2.R
#' @importClassesFrom terra SpatVector
#' @exportMethod clip
setMethod(
  "clip",
  signature = c("icesat2.atl03_h5", "character", "numeric"),
  function(x, output, clip_obj, ...) {
    print("clipping by bbox")
    bbox_ext <- terra::ext(clip_obj[c(2, 4, 3, 1)])
    ATL03_h5_clipBox(x, output, bbox_ext, ...)
  }
)

#' @include class.icesat2.R ATL08_h5_clipBox.R
#' @importClassesFrom terra SpatExtent
setMethod(
  "clip",
  signature = c("icesat2.atl08_h5", "character", "SpatExtent"),
  function(x, output, clip_obj, ...) {
    ATL08_h5_clipBox(x, output, clip_obj, ...)
  }
)

setMethod(
  "clip",
  signature = c("icesat2.atl08_h5", "character", "numeric"),
  function(x, output, clip_obj, ...) {
    print("clipping by bbox")
    bbox_ext <- terra::ext(clip_obj[c(2, 4, 3, 1)])
    ATL08_h5_clipBox(x, output, bbox_ext, ...)
  }
)

#' @include class.icesat2.R
#' @importClassesFrom terra SpatVector
#' @exportMethod clip
setMethod(
  "clip",
  signature = c("icesat2.atl08_h5", "character", "SpatVector"),
  function(x, output, clip_obj, polygon_id = "id", ...) {
    ATL08_h5_clipGeometry(x, output, clip_obj, polygon_id, ...)
  }
)

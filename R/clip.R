#' Clip ICESat-2 Data
#'
#' Generic function for clipping ICESat-2 ATL03 and ATL08 H5 data.
#'
#' @description This function clips ATL03 and ATL08 HDF5 file within beam groups,
#' but keeps metada and ancillary data the same.
#'
#' @param x [`ICESat2VegR::icesat2.atl03_h5-class`] or
#' [`ICESat2VegR::icesat2.atl08_h5-class`] object, obtained through
#' [`ATL03_read()`] or [`ATL08_read()`] for clipping
#' @param output character. Path to the output h5 file.
#' @param clip_obj [`numeric-class`], [`terra::SpatExtent`] or [`terra::SpatVector`]
#' object for clipping. The bbox order is the default from NASA's ICESat-2 data searching:
#' (ul_lat, ul_lon, lr_lat, lr_lon).
#' @param ... Forwards to specific method implementation.
#'
#' @return Returns the clipped S4 object of the same class as the input file
#'
#' @seealso
#' [`ATL03_h5_clipBox()`], [`ATL03_h5_clipGeometry()`],
#' [`ATL08_h5_clipBox()`], [`ATL08_h5_clipGeometry()`]
#'
#' @examples
#' # ATL08 file path
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' output <- tempfile(fileext = ".h5")
#'
#' vect_path <- system.file("extdata",
#'   "clip_geom.shp",
#'   package = "ICESat2VegR"
#' )
#'
#' vect <- terra::vect(vect_path)
#' ext <- terra::ext(vect)
#'
#' # Clipping ATL08 photons by boundary box extent
#' atl08_clip <- clip(
#'   atl08_h5,
#'   output,
#'   ext
#' )
#'
#' close(atl08_clip)
#'
#' # Clipping ATL08 photons by geometry
#' atl08_clip_geom <- clip(
#'   atl08_h5,
#'   output,
#'   vect,
#'   polygon_id = "id"
#' )
#'
#' close(atl08_clip_geom, close)
#' @import hdf5r
#' @export
setGeneric(
  "clip",
  function(x, output, clip_obj, ...) {
    standardGeneric("clip")
  }
)

#' @describeIn default clip will dispatch to [`graphics::clip()`]
#' @export
setMethod(
  "clip",
  signature = c("ANY", "ANY", "ANY"),
  function(x, output, clip_obj, ...) {
    graphics::clip(x, output, clip_obj, list(...)[[1]])
  }
)

#' @describeIn clip Clips ATL03 data using a `SpatExtent`.
#'
#' @param x An `icesat2.atl03_h5` object.
#' @param output A character string specifying the path to the output H5 file.
#' @param clip_obj A `terra::SpatExtent` object for clipping.
#' @param ... Additional arguments passed to the method implementation.
#' @export
setMethod(
  "clip",
  signature = c("icesat2.atl03_h5", "character", "SpatExtent"),
  function(x, output, clip_obj, ...) {
    ATL03_h5_clipBox(x, output, clip_obj, ...)
  }
)



#' @describeIn clip Clips ATL03 data using a `SpatVector`.
#'
#' This method clips ATL03 HDF5 data within beam groups using a `terra::SpatVector` object,
#' but keeps metadata and ancillary data the same.
#'
#' @param x An `icesat2.atl03_h5` object, obtained through `ATL03_read()`.
#' @param output A character string specifying the path to the output H5 file.
#' @param clip_obj A `terra::SpatVector` object for clipping. The bounding box order is the
#' default from NASA's ICESat-2 data searching: (ul_lat, ul_lon, lr_lat, lr_lon).
#' @param polygon_id A character string specifying the ID of the polygon for clipping.
#' @param ... Additional arguments passed to specific method implementation.
#'
#' @return Returns the clipped `icesat2.atl03_h5` object.
#'
#' @importClassesFrom terra SpatVector
#' @export
setMethod(
  "clip",
  signature = c("icesat2.atl03_h5", "character", "SpatVector"),
  function(x, output, clip_obj, polygon_id = "id", ...) {
    ATL03_h5_clipGeometry(x, output, clip_obj, polygon_id, ...)
  }
)

#' @describeIn clip Clips ATL03 data using numeric coordinates.
#'
#' This method clips ATL03 HDF5 data within beam groups using numeric coordinates,
#' but keeps metadata and ancillary data the same.
#'
#' @param x An `icesat2.atl03_h5` object, obtained through `ATL03_read()`.
#' @param output A character string specifying the path to the output H5 file.
#' @param clip_obj A numeric vector specifying the bounding box for clipping.
#' The bounding box order is the default from NASA's ICESat-2 data searching:
#' (ul_lat, ul_lon, lr_lat, lr_lon).
#' @param ... Additional arguments passed to specific method implementation.
#'
#' @return Returns the clipped `icesat2.atl03_h5` object.
#' @include class.icesat2.R
#' @importClassesFrom terra SpatVector
#' @exportMethod clip
setMethod(
  "clip",
  signature = c("icesat2.atl03_h5", "character", "numeric"),
  function(x, output, clip_obj, ...) {
    bbox_ext <- terra::ext(clip_obj[c(2, 4, 3, 1)])
    ATL03_h5_clipBox(x, output, bbox_ext, ...)
  }
)

#' @describeIn clip Clips ATL08 data using a `SpatExtent`.
#'
#' This method clips ATL08 HDF5 data within beam groups using a `terra::SpatExtent` object,
#' but keeps metadata and ancillary data the same.
#'
#' @param x An `icesat2.atl08_h5` object, obtained through `ATL08_read()`.
#' @param output A character string specifying the path to the output H5 file.
#' @param clip_obj A `terra::SpatExtent` object for clipping.
#' @param ... Additional arguments passed to specific method implementation.
#' @return Returns the clipped `icesat2.atl08_h5` object.
#' @export
setMethod(
  "clip",
  signature = c("icesat2.atl08_h5", "character", "SpatExtent"),
  function(x, output, clip_obj, ...) {
    ATL08_h5_clipBox(x, output, clip_obj, ...)
  }
)

#' @describeIn clip Clips ATL08 data using numeric coordinates.
#'
#' This method clips ATL08 HDF5 data within beam groups using numeric coordinates,
#' but keeps metadata and ancillary data the same.
#'
#' @param x An `icesat2.atl08_h5` object, obtained through `ATL08_read()`.
#' @param output A character string specifying the path to the output H5 file.
#' @param clip_obj A numeric vector specifying the bounding box for clipping.
#' The bounding box order is the default from NASA's ICESat-2 data searching:
#' (ul_lat, ul_lon, lr_lat, lr_lon).
#' @param ... Additional arguments passed to specific method implementation.
#'
#' @return Returns the clipped `icesat2.atl08_h5` object.
setMethod(
  "clip",
  signature = c("icesat2.atl08_h5", "character", "numeric"),
  function(x, output, clip_obj, ...) {
    print("clipping by bbox")
    bbox_ext <- terra::ext(clip_obj[c(2, 4, 3, 1)])
    ATL08_h5_clipBox(x, output, bbox_ext, ...)
  }
)

#' @describeIn clip Clips ATL08 data using a `SpatVector`.
#'
#' This method clips ATL08 HDF5 data within beam groups using a `terra::SpatVector` object,
#' but keeps metadata and ancillary data the same.
#'
#' @param x An `icesat2.atl08_h5` object, obtained through `ATL08_read()`.
#' @param output A character string specifying the path to the output H5 file.
#' @param clip_obj A `terra::SpatVector` object for clipping. The bounding box order is the
#' default from NASA's ICESat-2 data searching: (ul_lat, ul_lon, lr_lat, lr_lon).
#' @param polygon_id A character string specifying the ID of the polygon for clipping.
#' @param ... Additional arguments passed to specific method implementation.
#'
#' @return Returns the clipped `icesat2.atl08_h5` object.
#'
#' @include class.icesat2.R
#' @importClassesFrom terra SpatVector
#' @export
setMethod(
  "clip",
  signature = c("icesat2.atl08_h5", "character", "SpatVector"),
  function(x, output, clip_obj, polygon_id, ...) {
    ATL08_h5_clipGeometry(x, output, clip_obj, polygon_id, ...)
  }
)

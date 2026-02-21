# =============================================================================
# Generic: clip
# =============================================================================
#' Clip ICESat-2 data (h5, attributes, or ATL03-ATL08 join) by box or geometry
#'
#' @description
#' Unified clipping interface for ICESat-2 ATL03/ATL08 data. The generic
#' `clip()` dispatches on the class of `x` and internally calls appropriate
#' `_clipBox()` or `_clipGeometry()` helpers.
#'
#' @param x ICESat-2 object (ATL03/ATL08 h5, attributes, or joined table)
#' Can be either:
#'
#' - [`icesat2.atl03_h5-class`]
#' - [`icesat2.atl03_dt-class`]
#' - [`icesat2.atl03_seg_dt-class`]
#' - [`icesat2.atl08_h5-class`]
#' - [`icesat2.atl08_dt-class`]
#' - [`icesat2.atl03atl08_dt-class`]
#'
#' @param clip_obj Clipping object. Supported for box-based clipping:
#'   * numeric bounding box: c(xmin, ymin, xmax, ymax)
#'   * `terra::SpatExtent`
#'
#'   For geometry-based clipping:
#'   * `sf`, `sfc`, `terra::SpatVector`, etc.
#'   * Will also accept additional argument `split_by`:
#'       - `split_by` Character. The SpatVector attribute
#'   whose values are used to identify and name each clipping
#'   geometry. Defaults to `id`.
#'
#'   To clip by a spatial object's extent, extract its bbox first.
#'
#' @param ... Additional arguments passed to `_clipBox()` or `_clipGeometry()`.
#'
#' For hdf5 based clip, the additional arguments are:
#'   - `output` Character. Path to the output HDF5 file that will store the
#'   clipped ATL03 dataset.
#'   - `beam` Character vector specifying which ATL03 beams to include.
#'   Defaults to:
#'   `c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")`.
#'   - `additional_groups` Character vector of additional non-beam HDF5 groups
#'   to copy unchanged into the output file. Defaults to:
#'   `c("orbit_info")`.
#' 
#' @export
setGeneric(
  "clip",
  function(x, clip_obj, ...) standardGeneric("clip")
)

# =============================================================================
# Methods: clip() for ATL03/ATL08 HDF5 handles
# =============================================================================

#' Clip ICESat-2 ATL03 HDF5 Data Using a Bounding Extent
#'
#' @description
#' Clips an ICESat-2 ATL03 HDF5 file to a specified spatial extent. Only
#' geolocated photon and segment datasets within the ATL03 beam groups are
#' clipped; all metadata, orbit information, and ancillary groups are preserved
#' unchanged. This function provides bounding-box-based clipping, complementing
#' geometry-based clipping workflows.
#'
#' @param x An [`ICESat2VegR::icesat2.atl03_h5-class`] object created by
#'   [`ATL03_read()`], representing the ATL03 HDF5 file to be clipped.
#'
#' @param clip_obj Bounding extent used for clipping. Supported inputs:
#'   \itemize{
#'     \item A numeric vector of length 4:
#'       `c(ul_lat, ul_lon, lr_lat, lr_lon)`, following the ICESat-2 CMS
#'       bounding-box convention (upper-left latitude/longitude and
#'       lower-right latitude/longitude).
#'     \item A [`terra::SpatExtent`] object.
#'   }
#' 
#' @param ... Additional arguments:
#' 
#'   - `output` Character. Path to the output HDF5 file that will store the
#'   clipped ATL03 dataset.
#'   - `beam` Character vector specifying which ATL03 beams to include.
#'   Defaults to:
#'   `c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")`.
#'   - `additional_groups` Character vector of additional non-beam HDF5 groups
#'   to copy unchanged into the output file. Defaults to:
#'   `c("orbit_info")`.
#'
#' @return
#' An S4 object of class [`ICESat2VegR::icesat2.atl03_h5-class`] representing
#' the clipped ATL03 HDF5 file saved at the `output` location.
#'
#' @details
#' The function performs spatial subsetting at the photon and segment level.
#' Internally, it:
#'
#' \enumerate{
#'   \item Copies file- and group-level attributes.
#'   \item Extracts beam-level datasets and evaluates segment-level and
#'         photon-level masks using the supplied bounding box.
#'   \item Applies these masks to all datasets whose dimensions correspond to
#'         segments or photons.
#'   \item Recomputes dependent indexing datasets (e.g., `ph_index_beg`)
#'         to maintain ATL03 structural integrity.
#'   \item Writes the clipped data into a new, valid ATL03 HDF5 file.
#' }
#'
#' This function preserves the full ATL03 metadata, ancillary information, and
#' file structure while trimming the spatial domain to the user-specified
#' bounding extent.
#'
#' @examples
#' # ATL03 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ATL03 data (HDF5)
#' atl03_h5 <- ATL03_read(atl03_path)
#'
#' # Bounding rectangle coordinates (ul_lat, ul_lon, lr_lat, lr_lon)
#' xmin <- -106.5723
#' xmax <- -106.5693
#' ymin <- 41.533
#' ymax <- 41.537
#'
#' bbox <- c(ymax, xmin, ymin, xmax)
#'
#' # Clip ATL03 using the bounding extent
#' output <- tempfile(fileext = ".h5")
#' atl03_clip <- ATL03_h5_clipBox(
#'   atl03_h5,
#'   output = output,
#'   clip_obj = bbox
#' )
#'
#' close(atl03_h5)
#'
#' @import hdf5r
#' @include clipTools.R
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_h5", clip_obj = "numeric"),
  function(x, clip_obj, ...) {
    ATL03_h5_clipBox(x, clip_obj = clip_obj, ...)
  }
)

#' Clip ICESat-2 ATL03 HDF5 Data Using a Bounding Extent
#'
#' @description
#' Clips an ICESat-2 ATL03 HDF5 file to a specified spatial extent. Only
#' geolocated photon and segment datasets within the ATL03 beam groups are
#' clipped; all metadata, orbit information, and ancillary groups are preserved
#' unchanged. This function provides bounding-box-based clipping, complementing
#' geometry-based clipping workflows.
#'
#' @param x An [`ICESat2VegR::icesat2.atl03_h5-class`] object created by
#'   [`ATL03_read()`], representing the ATL03 HDF5 file to be clipped.
#'
#' @param clip_obj Bounding extent used for clipping. Supported inputs:
#'   \itemize{
#'     \item A numeric vector of length 4:
#'       `c(ul_lat, ul_lon, lr_lat, lr_lon)`, following the ICESat-2 CMS
#'       bounding-box convention (upper-left latitude/longitude and
#'       lower-right latitude/longitude).
#'     \item A [`terra::SpatExtent`] object.
#'   }
#'
#' @param ... Additional arguments:
#' 
#'   - `output` Character. Path to the output HDF5 file that will store the
#'   clipped ATL03 dataset.
#'   - `beam` Character vector specifying which ATL03 beams to include.
#'   Defaults to:
#'   `c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")`.
#'   - `additional_groups` Character vector of additional non-beam HDF5 groups
#'   to copy unchanged into the output file. Defaults to:
#'   `c("orbit_info")`.
#'
#' @return
#' An S4 object of class [`ICESat2VegR::icesat2.atl03_h5-class`] representing
#' the clipped ATL03 HDF5 file saved at the `output` location.
#'
#' @details
#' The function performs spatial subsetting at the photon and segment level.
#' Internally, it:
#'
#' \enumerate{
#'   \item Copies file- and group-level attributes.
#'   \item Extracts beam-level datasets and evaluates segment-level and
#'         photon-level masks using the supplied bounding box.
#'   \item Applies these masks to all datasets whose dimensions correspond to
#'         segments or photons.
#'   \item Recomputes dependent indexing datasets (e.g., `ph_index_beg`)
#'         to maintain ATL03 structural integrity.
#'   \item Writes the clipped data into a new, valid ATL03 HDF5 file.
#' }
#'
#' This function preserves the full ATL03 metadata, ancillary information, and
#' file structure while trimming the spatial domain to the user-specified
#' bounding extent.
#'
#' @examples
#' # ATL03 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ATL03 data (HDF5)
#' atl03_h5 <- ATL03_read(atl03_path)
#'
#' # Bounding rectangle coordinates (ul_lat, ul_lon, lr_lat, lr_lon)
#' xmin <- -106.5723
#' xmax <- -106.5693
#' ymin <- 41.533
#' ymax <- 41.537
#'
#' bbox <- c(ymax, xmin, ymin, xmax)
#'
#' # Clip ATL03 using the bounding extent
#' output <- tempfile(fileext = ".h5")
#' atl03_clip <- ATL03_h5_clipBox(
#'   atl03_h5,
#'   output = output,
#'   clip_obj = bbox
#' )
#'
#' close(atl03_h5)
#'
#' @import hdf5r
#' @include clipTools.R
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_h5", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) {
    ATL03_h5_clipBox(x, clip_obj = clip_obj, ...)
  }
)

#' Clip ICESat-2 ATL03 HDF5 Data Using Geometry-Based Boundaries
#'
#' @description
#' Clips an ICESat-2 ATL03 HDF5 file using one or more geometry-based clipping
#' objects provided as a [`terra::SpatVector`]. Each geometry in `vect`
#' defines an individual clipping region. The function iteratively clips the
#' ATL03 data for each region and writes a separate output HDF5 file for every
#' geometry. The name of each output file is automatically appended with the
#' value of the attribute specified in `split_by`.
#'
#' Only datasets within the ATL03 beam groups are spatially clipped; metadata
#' and ancillary non-beam groups remain unchanged in the resulting files.
#'
#' @param x An [`ICESat2VegR::icesat2.atl03_h5-class`] object obtained using
#'   [`ATL03_read()`], representing the ATL03 HDF5 file to be clipped.
#'
#'
#' @param clip_obj A [`terra::SpatVector`] object containing the geometric clip
#'   boundaries. Each feature (row) in the SpatVector defines an independent
#'   clipping region.
#' 
#' @param ... Additional arguments:
#'
#'   - `output` Character. Path to the output filename. The final written files
#'   will append the unique value of `split_by` for each clipping geometry,
#'   for example:
#'   `output_id1.h5`, `output_id2.h5`, etc.
#'   - `split_by` Character. The SpatVector attribute whose values are used to
#'   identify and name each clipping geometry. Defaults to `id`.
#'   - `beam` Character vector specifying which ATL03 beams to include. The
#'   default includes all six beams:
#'   `c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")`.
#'   - `additional_groups` Character vector of non-beam HDF5 groups that should
#'   be copied unchanged into each clipped file. Defaults to:
#'   `c("orbit_info")`.
#'
#' @return
#' Returns a list of clipped S4 objects of class
#' [`ICESat2VegR::icesat2.atl03_h5-class`], one for each clipping geometry
#' contained in `vect`.
#'
#' @examples
#' # ATL03 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ATL03 file
#' atl03_h5 <- ATL03_read(atl03_path)
#'
#' # Output base filename
#' output <- tempfile(fileext = ".h5")
#'
#' # Load clipping geometries
#' vect_path <- system.file("extdata", "clip_geom.shp",
#'   package = "ICESat2VegR"
#' )
#' vect <- terra::vect(vect_path)
#'
#' # Clip ATL03 data using polygon geometries
#' atl03_clipped_list <- ATL03_h5_clipGeometry(
#'   atl03_h5,
#'   output,
#'   vect,
#'   split_by = "id"
#' )
#'
#' close(atl03_h5)
#'
#' @import hdf5r
#' @include clipTools.R
#' @export
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_h5", clip_obj = "ANY"),
  function(x, clip_obj, ...) {
    ATL03_h5_clipGeometry(x, clip_obj = clip_obj, ...)
  }
)

#' Clip ICESat-2 ATL08 HDF5 Data Using a Bounding Extent
#'
#' @description
#' Clips an ICESat-2 ATL08 HDF5 file to a specified spatial extent. Only
#' segment-level datasets within the ATL08 beam groups are clipped; all metadata,
#' orbit information, and ancillary groups are preserved unchanged. This
#' function is the bounding-box-based counterpart to geometry-based clipping
#' functions.
#'
#' @param x An [`ICESat2VegR::icesat2.atl08_h5-class`] object obtained via
#'   [`ATL08_read()`], representing the ATL08 HDF5 file to be clipped.
#'
#' @param clip_obj Bounding extent used for clipping. Supported inputs:
#'   \itemize{
#'     \item A numeric vector of length 4:
#'       `c(ul_lat, ul_lon, lr_lat, lr_lon)`, following the ICESat-2 CMS
#'       bounding-box convention (upper-left latitude/longitude and
#'       lower-right latitude/longitude).
#'     \item A [`terra::SpatExtent`] object.
#'   }
#'
#' @param ... Additional arguments: 
#'
#'   - `output` Character path specifying where the clipped HDF5 file should be
#'   written.
#'   - `beam` Character vector specifying which ATL08 beams to include.
#'   Defaults to:
#'   `c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")`.
#'   - `additional_groups` Character vector specifying additional non-beam
#'   HDF5 groups to copy unchanged into the output file. Defaults to:
#'   `c("orbit_info")`.
#'
#' @return
#' An S4 object of class [`ICESat2VegR::icesat2.atl08_h5-class`] representing
#' the clipped ATL08 HDF5 file.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Copies file-level attributes and all requested non-beam groups.
#'   \item Clips beam-level datasets based on latitude/longitude coordinates and
#'         the supplied spatial extent.
#'   \item Reconstructs dependent indexing datasets (e.g., photon index ranges)
#'         when necessary to maintain valid ATL08 structure.
#' }
#'
#' The resulting HDF5 file remains fully compliant with the ATL08 product
#' structure, but contains only data within the specified bounding extent.
#'
#' @examples
#' # ATL08 file path
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ATL08 data
#' atl08_h5 <- ATL08_read(atl08_path)
#'
#' # Bounding rectangle coordinates (ul_lat, ul_lon, lr_lat, lr_lon)
#' ul_lon <- -106.5723
#' lr_lon <- -106.5693
#' lr_lat <- 41.533
#' ul_lat <- 41.537
#'
#' # Clip ATL08 data using the bounding extent
#' atl08_clip <- ATL08_h5_clipBox(
#'   atl08_h5,
#'   output = tempfile(fileext = ".h5"),
#'   clip_obj = c(ul_lat, ul_lon, lr_lat, lr_lon)
#' )
#'
#' close(atl08_h5)
#' close(atl08_clip)
#'
#' @import hdf5r
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_h5", clip_obj = "numeric"),
  function(x, clip_obj, ...) {
    ATL08_h5_clipBox(x, clip_obj = clip_obj, ...)
  }
)

#' Clip ICESat-2 ATL08 HDF5 Data Using a Bounding Extent
#'
#' @description
#' Clips an ICESat-2 ATL08 HDF5 file to a specified spatial extent. Only
#' segment-level datasets within the ATL08 beam groups are clipped; all metadata,
#' orbit information, and ancillary groups are preserved unchanged. This
#' function is the bounding-box-based counterpart to geometry-based clipping
#' functions.
#'
#' @param x An [`ICESat2VegR::icesat2.atl08_h5-class`] object obtained via
#'   [`ATL08_read()`], representing the ATL08 HDF5 file to be clipped.
#' @param clip_obj Bounding extent used for clipping. Supported inputs:
#'   \itemize{
#'     \item A numeric vector of length 4:
#'       `c(ul_lat, ul_lon, lr_lat, lr_lon)`, following the ICESat-2 CMS
#'       bounding-box convention (upper-left latitude/longitude and
#'       lower-right latitude/longitude).
#'     \item A [`terra::SpatExtent`] object.
#'   }
#'
#' @param ... Additional arguments:
#'   - `output` Character path specifying where the clipped HDF5 file should be
#'   written.
#'   - `beam` Character vector specifying which ATL08 beams to include.
#'   Defaults to:
#'   `c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")`.
#'
#'   - `additional_groups` Character vector specifying additional non-beam
#'   HDF5 groups to copy unchanged into the output file. Defaults to:
#'   `c("orbit_info")`.
#'
#' @return
#' An S4 object of class [`ICESat2VegR::icesat2.atl08_h5-class`] representing
#' the clipped ATL08 HDF5 file.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Copies file-level attributes and all requested non-beam groups.
#'   \item Clips beam-level datasets based on latitude/longitude coordinates and
#'         the supplied spatial extent.
#'   \item Reconstructs dependent indexing datasets (e.g., photon index ranges)
#'         when necessary to maintain valid ATL08 structure.
#' }
#'
#' The resulting HDF5 file remains fully compliant with the ATL08 product
#' structure, but contains only data within the specified bounding extent.
#'
#' @examples
#' # ATL08 file path
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ATL08 data
#' atl08_h5 <- ATL08_read(atl08_path)
#'
#' # Bounding rectangle coordinates (ul_lat, ul_lon, lr_lat, lr_lon)
#' ul_lon <- -106.5723
#' lr_lon <- -106.5693
#' lr_lat <- 41.533
#' ul_lat <- 41.537
#'
#' # Clip ATL08 data using the bounding extent
#' atl08_clip <- ATL08_h5_clipBox(
#'   atl08_h5,
#'   output = tempfile(fileext = ".h5"),
#'   clip_obj = c(ul_lat, ul_lon, lr_lat, lr_lon)
#' )
#'
#' close(atl08_h5)
#' close(atl08_clip)
#'
#' @import hdf5r
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_h5", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) {
    ATL08_h5_clipBox(x, clip_obj = clip_obj, ...)
  }
)

#' Clips ICESat-2 ATL08 data
#'
#' @param x [`ICESat2VegR::icesat2.atl08_h5-class`] object, obtained through [`ATL08_read()`]
#' for clipping
#' @param clip_obj [`terra::SpatVector-class`] for clipping
#' @param ... Additional arguments:
#'   - `output` character. Path to the output h5 file.
#'   - `split_by` [`character-class`]. The attribute name used for identifying
#' the different clip_objs. Default is "id"
#'   - `beam` [`character-class`]. The vector of beams to include, default
#' all c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")
#'   - `additional_groups` [`character-class`]. Other addional groups that should be included, default
#' c("orbit_info")
#'
#' @return Returns a list of clipped S4 object of class [`ICESat2VegR::icesat2.atl08_h5-class`]
#'
#' @description This function clips ATL08 HDF5 file within beam groups,
#' but keeps metada and ancillary data the same.
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
#'
#' # Clipping ATL08 photons by boundary box extent
#' atl08_photons_dt_clip <- ATL08_h5_clipGeometry(
#'   atl08_h5,
#'   output,
#'   vect,
#'   split_by = "id"
#' )
#'
#' close(atl08_h5)
#' @import hdf5r
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_h5", clip_obj = "ANY"),
  function(x, clip_obj, ...) {
    ATL08_h5_clipGeometry(x, clip_obj = clip_obj, ...)
  }
)

# =============================================================================
# Methods: clip() for ATL03/ATL08 attribute tables
# =============================================================================

#' Clip ATL03 photons by bounding extent
#'
#' @description
#' Clips ATL03 photon attributes within a given bounding extent, defined by a
#' single `clip_obj` argument. The clipping extent can be provided as:
#'
#' \itemize{
#'   \item (i) a numeric bounding box, or
#'   \item (ii) a spatial extent object.
#' }
#'
#' @param x An `icesat2.atl03_dt` object (output of
#'   [ATL03_photons_attributes_dt()]).
#'
#' @param clip_obj Bounding extent used to perform the clipping. Supported
#'   inputs:
#'   \itemize{
#'     \item Numeric vector of length 4: `c(xmin, ymin, xmax, ymax)` in
#'           decimal degrees.
#'     \item `SpatExtent` (package **terra**): the extent is used directly.
#'   }
#' @param ... not used
#'
#' @return
#' Returns a subset of the original `icesat2.atl03_dt` object containing
#' only photons within the bounding box.
#'
#' @details
#' When `clip_obj` is a `SpatExtent`, the package \pkg{terra} must be
#' installed. If not found, the function stops with an informative message.
#'
#' @seealso
#'  \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#' \dontrun{
#' atl03_path <- system.file("extdata", "atl03_clip.h5",
#'                           package = "ICESat2VegR")
#'
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5)
#'
#' # 1) Using a numeric bbox: xmin ymin xmax ymax
#' bbox <- c(-106.57, 41.53, -106.5698, 41.54)
#' atl03_clip_bbox <- ATL03_photons_attributes_dt_clipBox(
#'   atl03_photons_dt = atl03_photons_dt,
#'   clip_obj = bbox
#' )
#'
#' # 2) Using a SpatExtent
#' # library(terra)
#' # ext <- terra::ext(-106.57, -106.5698, 41.53, 41.54)
#' # atl03_clip_ext <- ATL03_photons_attributes_dt_clipBox(
#' #   atl03_photons_dt = atl03_photons_dt,
#' #   clip_obj = ext
#' # )
#'
#' close(atl03_h5)
#' }
#'
#' @import hdf5r stats
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_dt", clip_obj = "numeric"),
  function(x, clip_obj, ...) {
    ATL03_photons_attributes_dt_clipBox(x, clip_obj = clip_obj, ...)
  }
)

#' Clip ATL03 photons by bounding extent
#'
#' @description
#' Clips ATL03 photon attributes within a given bounding extent, defined by a
#' single `clip_obj` argument. The clipping extent can be provided as:
#'
#' \itemize{
#'   \item (i) a numeric bounding box, or
#'   \item (ii) a spatial extent object.
#' }
#'
#' @param x An `icesat2.atl03_dt` object (output of
#'   [ATL03_photons_attributes_dt()]).
#'
#' @param clip_obj Bounding extent used to perform the clipping. Supported
#'   inputs:
#'   \itemize{
#'     \item Numeric vector of length 4: `c(xmin, ymin, xmax, ymax)` in
#'           decimal degrees.
#'     \item `SpatExtent` (package **terra**): the extent is used directly.
#'   }
#' @param ... not used
#'
#' @return
#' Returns a subset of the original `icesat2.atl03_dt` object containing
#' only photons within the bounding box.
#'
#' @details
#' When `clip_obj` is a `SpatExtent`, the package \pkg{terra} must be
#' installed. If not found, the function stops with an informative message.
#'
#' @seealso
#'  \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#' \dontrun{
#' atl03_path <- system.file("extdata", "atl03_clip.h5",
#'                           package = "ICESat2VegR")
#'
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5)
#'
#' # 1) Using a numeric bbox: xmin ymin xmax ymax
#' bbox <- c(-106.57, 41.53, -106.5698, 41.54)
#' atl03_clip_bbox <- ATL03_photons_attributes_dt_clipBox(
#'   atl03_photons_dt = atl03_photons_dt,
#'   clip_obj = bbox
#' )
#'
#' # 2) Using a SpatExtent
#' # library(terra)
#' # ext <- terra::ext(-106.57, -106.5698, 41.53, 41.54)
#' # atl03_clip_ext <- ATL03_photons_attributes_dt_clipBox(
#' #   atl03_photons_dt = atl03_photons_dt,
#' #   clip_obj = ext
#' # )
#'
#' close(atl03_h5)
#' }
#'
#' @import hdf5r stats
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_dt", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) {
    ATL03_photons_attributes_dt_clipBox(x, clip_obj = clip_obj, ...)
  }
)

#' Clip ATL03 photons by Coordinates
#'
#' @description This function clips ATL03 photon attributes within given bounding coordinates.
#'
#' @param x An ATL03 photon data table. An S4 object of class [`ICESat2VegR::icesat2.atl03_dt-class`].
#' @param clip_obj Spatial clip_obj. An object of class [`terra::SpatVector`],
#'   which can be loaded as an ESRI shapefile using the [terra::vect] function in the
#'   \emph{sf} package.
#' @param ... . The `split_by` parameter can be specified to clip_obj
#' by attribute. If defined, GEDI data will be clipped by each clip_obj using
#' the clip_obj id from the attribute table defined by the user.
#'
#' @return Returns an S4 object of class [`ICESat2VegR::icesat2.atl03_dt-class`]
#'   containing the ATL03 photon attributes.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#' # ATL03 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Extracting ATL03 photon attributes
#' atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5)
#'
#' # Specifying the path to shapefile
#' clip_obj_filepath <-
#'   system.file(
#'     "extdata",
#'     "clip_geom.shp",
#'     package = "ICESat2VegR"
#'   )
#'
#' # Reading shapefile as sf object
#' clip_obj <- terra::vect(clip_obj_filepath)
#'
#' # Clipping ATL03 photon attributes by Geometry
#' atl03_photons_dt_clip <-
#'   ATL03_photons_attributes_dt_clipGeometry(atl03_photons_dt, clip_obj, split_by = "id")
#'
#' head(atl03_photons_dt_clip)
#'
#' close(atl03_h5)
#' @import hdf5r stats
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_dt", clip_obj = "ANY"),
  function(x, clip_obj, ...) {
    ATL03_photons_attributes_dt_clipGeometry(x, clip_obj = clip_obj, ...)
  }
)

# =============================================================================
# Clip ATL08 terrain and canopy attributes by bounding extent
# =============================================================================

#' Clip ATL08 terrain and canopy attributes by bounding extent
#'
#' @description
#' Clips ATL08 terrain and canopy segment attributes within a given bounding
#' extent, defined by a single `clip_obj` argument. The clipping extent can
#' be provided as:
#'
#' \itemize{
#'   \item (i) a numeric bounding box, or
#'   \item (ii) a spatial extent object.
#' }
#'
#' @param x An `icesat2.atl08_dt` object (output of
#'   [ICESat2VegR::ATL08_seg_attributes_dt()] function).
#'
#' @param clip_obj Bounding extent used to perform the clipping. Supported
#'   inputs:
#'   \itemize{
#'     \item Numeric vector of length 4: `c(xmin, ymin, xmax, ymax)` in
#'           decimal degrees.
#'     \item `SpatExtent` (package **terra**): the extent is used directly.
#'   }
#' @param ... not used.
#'
#' @return
#' Returns an S4 object of class
#'   [`ICESat2VegR::icesat2.atl08_dt-class`]
#' containing the clipped ATL08 terrain and canopy attributes.
#'
#' @details
#' When `clip_obj` is a `SpatExtent`, the package \pkg{terra} must be
#' installed. If not found, the function stops with an informative message.
#'
#' @examples
#' \dontrun{
#' # Specifying the path to the ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path)
#'
#' # Extracting ATL08-derived segment attributes
#' atl08_seg_att_dt <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)
#'
#' # 1) Using a numeric bounding box: xmin ymin xmax ymax
#' bbox <- c(-103.7604, 59.4672, -103.7600, 59.4680)
#'
#' atl08_seg_att_dt_clip <- ATL08_seg_attributes_dt_clipBox(
#'   atl08_seg_att_dt = atl08_seg_att_dt,
#'   clip_obj = bbox
#' )
#'
#' # 2) Using a SpatExtent
#' # library(terra)
#' # ext <- terra::ext(-103.7604, -103.7600, 59.4672, 59.4680)
#' # atl08_seg_att_dt_clip_ext <- ATL08_seg_attributes_dt_clipBox(
#' #   atl08_seg_att_dt = atl08_seg_att_dt,
#' #   clip_obj = ext
#' # )
#'
#' close(atl08_h5)
#' }
#'
#' @import hdf5r stats
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_dt", clip_obj = "numeric"),
  function(x, clip_obj, ...) {
    ATL08_seg_attributes_dt_clipBox(x, clip_obj = clip_obj, ...)
  }
)

# =============================================================================
# Clip ATL08 terrain and canopy attributes by bounding extent
# =============================================================================

#' Clip ATL08 terrain and canopy attributes by bounding extent
#'
#' @description
#' Clips ATL08 terrain and canopy segment attributes within a given bounding
#' extent, defined by a single `clip_obj` argument. The clipping extent can
#' be provided as:
#'
#' \itemize{
#'   \item (i) a numeric bounding box, or
#'   \item (ii) a spatial extent object.
#' }
#'
#' @param x An `icesat2.atl08_dt` object (output of
#'   [ICESat2VegR::ATL08_seg_attributes_dt()] function).
#'
#' @param clip_obj Bounding extent used to perform the clipping. Supported
#'   inputs:
#'   \itemize{
#'     \item Numeric vector of length 4: `c(xmin, ymin, xmax, ymax)` in
#'           decimal degrees.
#'     \item `SpatExtent` (package **terra**): the extent is used directly.
#'   }
#' @param ... not used.
#'
#' @return
#' Returns an S4 object of class
#'   [`ICESat2VegR::icesat2.atl08_dt-class`]
#' containing the clipped ATL08 terrain and canopy attributes.
#'
#' @details
#' When `clip_obj` is a `SpatExtent`, the package \pkg{terra} must be
#' installed. If not found, the function stops with an informative message.
#'
#' @import hdf5r stats
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_dt", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) {
    ATL08_seg_attributes_dt_clipBox(x, clip_obj = clip_obj, ...)
  }
)

#' Clip ATL08 Terrain and Canopy Attributes by Geometry
#'
#' @description This function clips ATL08 Terrain and Canopy Attributes within a given geometry
#'
#' @param x An [icesat2.atl08_dt-class] object (output of
#' [`ATL08_seg_attributes_dt()`] function).
#' @param clip_obj clip_obj. An object of class [`terra::SpatVector`],
#' which can be loaded as an ESRI shapefile using [terra::vect] function in the
#' \emph{sf} package.
#' @param ... The `split_by` parameter can be specified to clip_obj id. 
#' If defined, ATL08 data will be clipped by each clip_obj using
#' the clip_obj id from table of attribute defined by the user
#'
#' @return Returns an S4 object of class [`ICESat2VegR::icesat2.atl08_dt-class`]
#' containing the clipped ATL08 Terrain and Canopy Attributes.
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL08-derived terrain and canopy attributes
#' atl08_seg_att_dt <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)
#'
#' clip_obj_path <- system.file("extdata",
#'   "clip_geom.shp",
#'   package = "ICESat2VegR"
#' )
#'
#' if (require(terra)) {
#'   polygon <- terra::vect(clip_obj_path)
#'
#'   head(atl08_seg_att_dt)
#'   # Clipping ATL08 Terrain and Canopy Attributes by Geometry
#'   atl08_seg_att_dt_clip <- ATL08_seg_attributes_dt_clipGeometry(atl08_seg_att_dt,
#'       polygon, split_by = "id")
#'
#'   hasLeaflet <- require(leaflet)
#'
#'   if (hasLeaflet) {
#'     leaflet() %>%
#'       addCircleMarkers(atl08_seg_att_dt_clip$longitude,
#'         atl08_seg_att_dt_clip$latitude,
#'         radius = 1,
#'         opacity = 1,
#'         color = "red"
#'       ) %>%
#'       addScaleBar(options = list(imperial = FALSE)) %>%
#'       addPolygons(
#'         data = polygon, weight = 1, col = "white",
#'         opacity = 1, fillOpacity = 0
#'       ) %>%
#'       addProviderTiles(providers$Esri.WorldImagery,
#'         options = providerTileOptions(minZoom = 3, maxZoom = 17)
#'       )
#'   }
#' }
#' close(atl08_h5)
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_dt", clip_obj = "ANY"),
  function(x, clip_obj, ...) {
    ATL08_seg_attributes_dt_clipGeometry(x, clip_obj = clip_obj, ...)
  }
)

# =============================================================================
# Methods: clip() for joined ATL03-ATL08 photons
# =============================================================================

#' Clip joined ATL03 and ATL08 photons by bounding extent
#'
#' @description
#' Clips joined ATL03 and ATL08 photon attributes within a given bounding
#' extent, defined by a single `clip_obj` argument. The clipping extent can
#' be provided as:
#'
#' \itemize{
#'   \item (i) a numeric bounding box, or
#'   \item (ii) a spatial extent object.
#' }
#'
#' @param x An S4 object of class
#'   [`ICESat2VegR::icesat2.atl03atl08_dt-class`] containing joined ATL03 and
#'   ATL08 data (output of
#'   [ICESat2VegR::ATL03_ATL08_photons_attributes_dt_join()] function).
#'
#' @param clip_obj Bounding extent used to perform the clipping. Supported
#'   inputs:
#'   \itemize{
#'     \item Numeric vector of length 4: `c(xmin, ymin, xmax, ymax)` in
#'           decimal degrees.
#'     \item `SpatExtent` (package **terra**): the extent is used directly.
#'   }
#' @param ... not used.
#'
#' @return
#' Returns an S4 object of class
#' [`ICESat2VegR::icesat2.atl03atl08_dt-class`] containing a subset of the
#' joined ATL03 and ATL08 photon attributes restricted to the clipping extent.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @details
#' When `clip_obj` is a `SpatExtent`, the package \pkg{terra} must be
#' installed. If not found, the function stops with an informative message.
#'
#' @examples
#' \dontrun{
#' # Specifying the path to ATL03 and ATL08 files
#' atl03_path <- system.file("extdata", "atl03_clip.h5",
#'                           package = "ICESat2VegR")
#'
#' atl08_path <- system.file("extdata", "atl08_clip.h5",
#'                           package = "ICESat2VegR")
#'
#' # Reading ATL03 and ATL08 data (h5 files)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Joining ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#' head(atl03_atl08_dt)
#'
#' # 1) Using a numeric bounding box: xmin ymin xmax ymax
#' bbox <- c(-103.7604, 59.4672, -103.7600, 59.4680)
#'
#' atl03_atl08_dt_clip <- ATL03_ATL08_photons_attributes_dt_clipBox(
#'   atl03_atl08_dt = atl03_atl08_dt,
#'   clip_obj = bbox
#' )
#' head(atl03_atl08_dt_clip)
#'
#' # 2) Using a SpatExtent (example)
#' # library(terra)
#' # ext <- terra::ext(-103.7604, -103.7600, 59.4672, 59.4680)
#' # atl03_atl08_dt_clip_ext <- ATL03_ATL08_photons_attributes_dt_clipBox(
#' #   atl03_atl08_dt = atl03_atl08_dt,
#' #   clip_obj = ext
#' # )
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' }
#'
#' @import hdf5r stats
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03atl08_dt", clip_obj = "numeric"),
  function(x, clip_obj, ...) {
    ATL03_ATL08_photons_attributes_dt_clipBox(
      atl03_atl08_dt = x,
      clip_obj       = clip_obj,
      ...
    )
  }
)

#' Clip joined ATL03 and ATL08 photons by bounding extent
#'
#' @description
#' Clips joined ATL03 and ATL08 photon attributes within a given bounding
#' extent, defined by a single `clip_obj` argument. The clipping extent can
#' be provided as:
#'
#' \itemize{
#'   \item (i) a numeric bounding box, or
#'   \item (ii) a spatial extent object.
#' }
#'
#' @param x An S4 object of class
#'   [`ICESat2VegR::icesat2.atl03atl08_dt-class`] containing joined ATL03 and
#'   ATL08 data (output of
#'   [ICESat2VegR::ATL03_ATL08_photons_attributes_dt_join()] function).
#'
#' @param clip_obj Bounding extent used to perform the clipping. Supported
#'   inputs:
#'   \itemize{
#'     \item Numeric vector of length 4: `c(xmin, ymin, xmax, ymax)` in
#'           decimal degrees.
#'     \item `SpatExtent` (package **terra**): the extent is used directly.
#'   }
#' @param ... not used.
#'
#' @return
#' Returns an S4 object of class
#' [`ICESat2VegR::icesat2.atl03atl08_dt-class`] containing a subset of the
#' joined ATL03 and ATL08 photon attributes restricted to the clipping extent.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @details
#' When `clip_obj` is a `SpatExtent`, the package \pkg{terra} must be
#' installed. If not found, the function stops with an informative message.
#'
#' @examples
#' \dontrun{
#' # Specifying the path to ATL03 and ATL08 files
#' atl03_path <- system.file("extdata", "atl03_clip.h5",
#'                           package = "ICESat2VegR")
#'
#' atl08_path <- system.file("extdata", "atl08_clip.h5",
#'                           package = "ICESat2VegR")
#'
#' # Reading ATL03 and ATL08 data (h5 files)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Joining ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#' head(atl03_atl08_dt)
#'
#' # 1) Using a numeric bounding box: xmin ymin xmax ymax
#' bbox <- c(-103.7604, 59.4672, -103.7600, 59.4680)
#'
#' atl03_atl08_dt_clip <- ATL03_ATL08_photons_attributes_dt_clipBox(
#'   atl03_atl08_dt = atl03_atl08_dt,
#'   clip_obj = bbox
#' )
#' head(atl03_atl08_dt_clip)
#'
#' # 2) Using a SpatExtent (example)
#' # library(terra)
#' # ext <- terra::ext(-103.7604, -103.7600, 59.4672, 59.4680)
#' # atl03_atl08_dt_clip_ext <- ATL03_ATL08_photons_attributes_dt_clipBox(
#' #   atl03_atl08_dt = atl03_atl08_dt,
#' #   clip_obj = ext
#' # )
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' }
#'
#' @import hdf5r stats
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03atl08_dt", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) {
    ATL03_ATL08_photons_attributes_dt_clipBox(
      atl03_atl08_dt = x,
      clip_obj       = clip_obj,
      ...
    )
  }
)

#' Clip Joined ATL03 and ATL08 by Geometry
#'
#' @description This function clips joined ATL03 and ATL08 photon attributes within a given geometry.
#'
#' @param x An S4 object of class [`ICESat2VegR::icesat2.atl03atl08_dt-class`]
#' containing ATL03 and ATL08 data (output of
#' [`ICESat2VegR::ATL03_ATL08_photons_attributes_dt_join()`] function).
#' @param clip_obj An object of class [`terra::SpatVector`].
#' which can be loaded as an ESRI shapefile using [terra::vect] function.
#' @param ... `split_by` Optional. clip_obj ID.
#' If defined, the data will be clipped by each clip_obj using
#' the clip_obj ID from the attribute table.
#'
#' @return Returns an S4 object of class [`ICESat2VegR::icesat2.atl03atl08_dt-class`]
#' containing the clipped ATL08 attributes.
#'
#' @examples
#' # Specifying the path to ATL03 and ATL08 files
#' atl03_path <- system.file("extdata", "atl03_clip.h5", package = "ICESat2VegR")
#' atl08_path <- system.file("extdata", "atl08_clip.h5", package = "ICESat2VegR")
#'
#' # Reading ATL03 and ATL08 data (h5 files)
#' atl03_h5 <- ATL03_read(atl03_path)
#' atl08_h5 <- ATL08_read(atl08_path)
#'
#' # Joining ATL03 and ATL08 photon attributes
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#' head(atl03_atl08_dt)
#'
#' # Specifying the path to the shapefile
#' clip_obj_filepath <- system.file("extdata", "clip_geom.shp", package = "ICESat2VegR")
#'
#' # Reading shapefile as a SpatVector object
#' clip_obj <- terra::vect(clip_obj_filepath)
#'
#' # Clipping ATL08 terrain attributes by geometry
#' atl03_atl08_dt_clip <- ATL03_ATL08_photons_attributes_dt_clipGeometry(
#'   atl03_atl08_dt,
#'   clip_obj,
#'   split_by = "id"
#' )
#' head(atl03_atl08_dt_clip)
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' @export
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03atl08_dt", clip_obj = "ANY"),
  function(x, clip_obj, ...) {
    ATL03_ATL08_photons_attributes_dt_clipGeometry(
      atl03_atl08_dt = x,
      clip_obj       = clip_obj,
      ...
    )
  }
)

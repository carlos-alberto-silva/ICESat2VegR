#' Clips ICESat-2 ATL03 H5 data
#'
#' @param x [`icesat2.atl03_h5-class`] object,
#' obtained through [`ATL03_read()`] for clipping
#' @param output character. Path to the output h5 file.
#' @param bbox [`numeric-class`] or [`terra::SpatExtent`] for clipping, the
#' order of the bbox is the default from NASA's ICESat-2 CMS searching:
#' [ul_lat, ul_lon, lr_lat, lr_lon].
#'
#' @return Returns the clipped S4 object of class [`icesat2.atl03_h5-class`]
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
#'
#' # Bounding rectangle coordinates
#' xmin <- -107.7
#' xmax <- -106.5
#' ymin <- 32.75
#' ymax <- 42.75
#'
#' # Clipping ATL03 photons  by boundary box extent
#' atl03_photons_dt_clip <- ATL03_h5_clipBox(atl03_h5, outdir, xmin, xmax, ymin, ymax)
#'
#' close(atl03_h5)
#' @import hdf5r
#' @export
ATL03_h5_clip <- function(atl03, output, clipObj, mask_fn) {
  dataset.rank <- dataset.dims <- obj_type <- name <- NA

  # Create a new HDF5 file
  newFile <- hdf5r::H5File$new(output, mode = "w")


  # Create all groups
  groups <- atl03$ls_groups(recursive = TRUE)


  for (group in groups) {
    grp <- newFile$create_group(group)

    # Create all atributes within group
    attributes <- atl03[[group]]$ls_attrs()
    for (attribute in attributes) {
      grp$create_attr(attribute, atl03[[group]]$attr(attribute))
    }
  }

  # Create root attributes
  attributes <- atl03$ls_attrs()
  for (attribute in attributes) {
    hdf5r::h5attr(newFile, attribute) <- atl03$attr(attribute)
  }

  # Get all beams
  beams <- atl03$beams

  nBeam <- 0
  nBeams <- length(beams)

  # Loop the beams
  # beamName = beams[nBeam + 1]
  for (beamName in beams) {
    nBeam <- nBeam + 1
    message(sprintf("Clipping %s (%d/%d)", beamName, nBeam, nBeams))

    # Get the reference beam
    beam <- atl03[[beamName]]

    # Get the beam to update
    updateBeam <- newFile[[beamName]]

    # Get the masks
    segmentsMask <- mask_fn(beam, clipObj)
    photonsMask <- ATL03_photons_segment_mask(beam, segmentsMask)

    # Get sizes of clipping datasets
    photonsSize <- beam[["heights/h_ph"]]$dims
    segmentsSize <- beam[["geolocation/ph_index_beg"]]$dims

    # Get all datasets
    datasets_dt <- beam$dt_datasets(recursive = TRUE)

    # Get all types of clipping photons/segment/no cut
    photonsCut <- datasets_dt[dataset.dims == photonsSize]$name
    photonsCut2D <- datasets_dt[grepl(photonsSize, dataset.dims) & dataset.rank == 2]$name
    specialCuts <- c(
      "geolocation/ph_index_beg"
    )
    segmentsCut <- datasets_dt[
      dataset.dims == segmentsSize &
        !name %in% specialCuts
    ]$name
    segmentsCut2D <- datasets_dt[
      grepl(segmentsSize, dataset.dims) & dataset.rank == 2
    ]$name
    allCuts <- c(photonsCut, photonsCut2D, segmentsCut, segmentsCut2D, specialCuts)
    nonCuts <- datasets_dt[
      !name %in% allCuts
    ]$name

    qtyList <- lapply(datasets_dt$dataset.dims, function(x) eval(parse(text = gsub("x", "*", x))))
    qty <- sum(unlist(qtyList))

    pb <- utils::txtProgressBar(min = 0, max = qty, style = 3)

    # Do clipping and copying

    clipByMask(beam, updateBeam, segmentsCut, segmentsMask, pb)
    clipByMask2D(beam, updateBeam, segmentsCut2D, segmentsMask, pb)
    clipByMask(beam, updateBeam, photonsCut, photonsMask, pb)
    clipByMask2D(beam, updateBeam, photonsCut2D, photonsMask, pb)

    photons_per_segment <- beam[["geolocation/segment_ph_cnt"]][segmentsMask]
    copyDataset(beam, updateBeam, "geolocation/ph_index_beg", cumsum(photons_per_segment), pb)

    for (dataset in nonCuts) {
      copyDataset(beam, updateBeam, dataset, beam[[dataset]][], pb)
    }
    close(pb)
  }
  newFile$close_all()

  ATL03_read(output)
}

#' Generic function for clipping ICESat-2 ATL03 and ATL08 H5 data
#'
#' @param x [`icesat2.atl03_h5-class`] or [`icesat2.atl08_h5-class`] object,
#' obtained through [`ATL03_read()`] or [`ATL08_read()`] for clipping
#' @param output character. Path to the output h5 file.
#' @param clip_obj [`numeric-class`], [`terra::SpatExtent`] or [`terra::SpatVector`]
#' object for clipping. The bbox order is the default from NASA's ICESat-2 CMS searching:
#' [ul_lat, ul_lon, lr_lat, lr_lon].
#' @param ... For making it compatible with passing ul_lat, ul_lon, lr_lat, lr_lon.
#'
#' @return Returns the clipped S4 object of the same class as input
#'
#' @details 
#' This function clips ATL03 and ATL08 HDF5 file within beam groups,
#' but keeps metada and ancillary data the same.
#' 
#' The clipping process will keep the entire segments if the reference
#' photon is within the clip_obj.
#' 
#' This function will dispatch to one of the specifics functions for clipping.
#'
#' @seealso
#' [`ATL03_h5_clipBox()`], [`ATL03_h5_clipGeometry()`],
#' [`ATL08_h5_clipBox()`], [`ATL08_h5_clip()`]
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
setMethod(
  "clip",
  signature = c("icesat2.atl03_h5", "character", "SpatExtent"),
  function(x, output, clip_obj) {
    ATL03_h5_clipBox(x, output, clip_obj)
  }
)


#' @include class.icesat2.R
#' @importClassesFrom terra SpatVector
setMethod(
  "clip",
  signature = c("icesat2.atl03_h5", "character", "SpatVector"),
  function(x, output, clip_obj, polygon_id = "id") {
    ATL03_h5_clipGeometry(x, output, clip_obj, polygon_id)
  }
)

#' @include class.icesat2.R
#' @importClassesFrom terra SpatVector
setMethod(
  "clip",
  signature = c("icesat2.atl03_h5", "character", "numeric"),
  function(x, output, clip_obj) {
    print("clipping by bbox")
    bbox_ext <- terra::ext(clip_obj[c(2, 4, 3, 1)])
    ATL03_h5_clipBox(x, output, bbox_ext)
  }
)



#' Clips ICESat-2 ATL03 H5 data
#'
#' @param x [`icesat2.atl03_h5-class`] object,
#' obtained through [`ATL03_read()`] for clipping
#' @param output character. Path to the output h5 file.
#' @param bbox [`numeric-class`] or [`terra::SpatExtent`] for clipping, the
#' order of the bbox is the default from NASA's ICESat-2 CMS searching:
#' [ul_lat, ul_lon, lr_lat, lr_lon].
#'
#' @return Returns the clipped S4 object of class [`icesat2.atl03_h5-class`]
#'
#' @description This function clips ATL03 HDF5 file within beam groups,
#' but keeps metada and ancillary data the same.
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
#'
#' # Bounding rectangle coordinates
#' xmin <- -107.7
#' xmax <- -106.5
#' ymin <- 32.75
#' ymax <- 42.75
#'
#' # Clipping ATL03 photons  by boundary box extent
#' atl03_photons_dt_clip <- ATL03_h5_clipBox(atl03_h5, outdir, xmin, xmax, ymin, ymax)
#'
#' close(atl03_h5)
#' @import hdf5r
#' @include clipTools.R
#' @export
ATL03_h5_clipBox <- function(atl03, output, bbox) {
  ATL03_h5_clip(atl03, output, bbox, ATL03_segments_mask)
}



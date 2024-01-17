#' @include class.icesat2.R
#' @importClassesFrom terra SpatExtent
setMethod(
  "clip",
  signature = c("icesat2.atl08_h5", "character", "SpatExtent"),
  function(x, output, clip_obj) {
    ATL08_h5_clipBox(x, output, clip_obj)
  }
)

setMethod(
  "clip",
  signature = c("icesat2.atl08_h5", "character", "numeric"),
  function(x, output, clip_obj) {
    print("clipping by bbox")
    bbox_ext <- terra::ext(clip_obj[c(2, 4, 3, 1)])
    ATL08_h5_clipBox(x, output, bbox_ext)
  }
)

#' @include class.icesat2.R
#' @importClassesFrom terra SpatVector
#' @exportMethod clip
setMethod(
  "clip",
  signature = c("icesat2.atl08_h5", "character", "SpatVector"),
  function(x, output, clip_obj, polygon_id = "id") {
    ATL08_h5_clipGeometry(x, output, clip_obj, polygon_id)
  }
)

#' Clips ICESat-2 ATL08 data
#'
#' @param atl08 [`icesat2.atl08_h5-class`] object, obtained through [`ATL08_read()`]
#' for clipping
#' @param output character. Path to the output h5 file.
#' @param bbox [`numeric-class`] or [`terra::SpatExtent`] for clipping, the
#' order of the bbox is the default from NASA's ICESat-2 CMS searching:
#' [ul_lat, ul_lon, lr_lat, lr_lon].
#'
#' @return Returns the clipped S4 object of class [`icesat2.atl08_h5-class`]
#'
#' @description This function clips the ATl08 HDF5 file. This function
#' will only clip the beam groups within hdf5, it won't change metadata
#' or ancillary data.
#'
#' @examples
#' # Specifying the path to ATL08 file (zip file)
#' outdir <- tempdir()
#' atl08_zip <- system.file("extdata",
#'   "ATL08_20220401221822_01501506_005_01.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ATL08 file
#' atl08_path <- unzip(atl08_zip, exdir = outdir)
#'
#' # Reading ATL08 data (h5 file)
# atl08_h5<-ATL08_read(atl08_path=atl08_path)
#'
#' # Bounding rectangle coordinates
#' xmin <- -107.7
#' xmax <- -106.5
#' ymin <- 32.75
#' ymax <- 42.75
#'
#' # Clipping ATL08 terrain and canopy attributes by boundary box
#' atl08_seg_att_dt_clip <- ATL08_h5_clipBox(atl08_h5, outdir, xmin, xmax, ymin, ymax)
#'
#' close(atl08_h5)
#' @import hdf5r
#' @export
ATL08_h5_clipBox <- function(atl08, output, bbox) {
  ATL08_h5_clip(atl08, output, bbox, landsegmentsMask_bbox)
}

#' Clips ICESat-2 ATL08 data
#'
#' @param atl08 [`icesat2.atl08_h5-class`] object, obtained through [`ATL08_read()`]
#' for clipping
#' @param output character. Path to the output h5 file.
#' @param vect [`terra::SpatVector-class`] for clipping
#' @param polygon_id [`character-class`]. The attribute name used for identifying
#' the different polygons. Default is "id"
#'
#' @return Returns a list of clipped S4 object of class [`icesat2.atl08_h5-class`]
#'
#' @description This function clips ATL08 HDF5 file within beam groups,
#' but keeps metada and ancillary data the same.
#'
#' @examples
##' # Specifying the path to ATL08 file (zip file)
#' outdir <- tempdir()
#' atl08_zip <- system.file("extdata",
#'   "atl08_20220401221822_01501506_005_01.zip",
#'   package = "rICESat2Veg"
#' )
#'
#' # Unzipping ATL08 file
#' atl08_path <- unzip(atl08_zip, exdir = outdir)
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- atl08_read(atl08_path = atl08_path)
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
#' # Clipping ATL08 photons  by boundary box extent
#' atl08_photons_dt_clip <- ATL08_h5_clipGeometry(
#'   atl08_h5,
#'   output,
#'   vect,
#'   polygon_id = "id"
#' )
#'
#' close(atl08_h5)
#' @import hdf5r
#' @export
ATL08_h5_clipGeometry <- function(atl08, output, vect, polygon_id = "id") {
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
      ATL08_h5_clip(atl08, sub_output, clip_obj = geom, landSegmentsMask_fn = landsegmentsMask_geom)
  }
}

landSegments_bbox <- function(beam, bbox) {
  latitude <- longitude <- NA
  .I <- data.table::.I

  land_segments_dt <- data.table::data.table(
    latitude = beam[["land_segments/latitude"]][],
    longitude = beam[["land_segments/longitude"]][]
  )

  landSegmentsDt <- land_segments_dt[, list(latitude, longitude, .I)][
    longitude >= bbox$xmin &
      longitude <= bbox$xmax &
      latitude >= bbox$ymin &
      latitude <= bbox$ymax
  ]

  landSegmentsDt
}

landsegmentsMask_bbox <- function(beam, bbox) {
  landSegmentsBbox <- landSegments_bbox(beam, bbox)

  landSegmentsBbox$I
}

landsegmentsMask_geom <- function(beam, geom) {
  bbox <- terra::ext(geom)
  landSegmentsBbox <- landSegments_bbox(beam, bbox)

  pts <- terra::vect(
    landSegmentsBbox,
    geom = c("longitude", "latitude"),
    crs = "epsg:4326"
  )

  res <- terra::intersect(pts, geom)
  res$I
}

ATL08_h5_clip <- function(atl08, output, clip_obj, landSegmentsMask_fn) {
  obj_type <- dataset.rank <- dataset.dims <- NA

  # Create a new HDF5 file
  newFile <- hdf5r::H5File$new(output, mode = "w")


  # Create all groups
  structure_dt <- data.table::as.data.table(atl08@h5$ls(recursive = T))
  groups <- structure_dt[obj_type == "H5I_GROUP"]$name


  for (group in groups) {
    grp <- newFile$create_group(group)

    # Create all atributes within group
    attributes <- hdf5r::list.attributes(atl08[[group]])
    for (attribute in attributes) {
      grp$create_attr(attribute, hdf5r::h5attr(atl08[[group]], attribute))
    }
  }

  # Create root attributes
  attributes <- hdf5r::list.attributes(atl08@h5)
  for (attribute in attributes) {
    hdf5r::h5attr(newFile, attribute) <- hdf5r::h5attr(atl08@h5, attribute)
  }

  # Get all beams
  beams <- getBeams(atl08)


  nBeam <- 0
  nBeams <- length(beams)

  # Loop the beams
  # beamName <- beams[2]
  for (beamName in beams) {
    nBeam <- nBeam + 1
    message(sprintf("Clipping %s (%d/%d)", beamName, nBeam, nBeams))

    # Get the reference beam
    beam <- atl08[[beamName]]

    # Get the beam to update
    updateBeam <- newFile[[beamName]]


    # Get land segments mask
    landSegmentsMask <- landSegmentsMask_fn(beam, clip_obj)

    if (length(landSegmentsMask) == 0) {
      next
    }

    # Get photons for the masked land segments
    ph_ndx_beg <- beam[["land_segments/ph_ndx_beg"]][landSegmentsMask]
    n_seg_ph <- beam[["land_segments/n_seg_ph"]][landSegmentsMask]

    # Create mask for photons
    photonsMask <- unlist(
      Vectorize(seq.default, vectorize.args = c("from", "to"))(from = ph_ndx_beg, to = ph_ndx_beg + n_seg_ph)
    )

    # Get sizes of clipping datasets
    photonsSize <- beam[["signal_photons/ph_h"]]$dims
    segmentsSize <- beam[["land_segments/segment_watermask"]]$dims

    # Get all datasets
    datasets_dt <- data.table::as.data.table(beam$ls(recursive = TRUE))[obj_type == 5]

    # Get all types of clipping photons/segment/no cut
    photonsCut <- datasets_dt[dataset.dims == photonsSize]$name

    segmentsCut <- datasets_dt[
      dataset.dims == segmentsSize
    ]$name
    segmentsCut2D <- datasets_dt[
      grepl(segmentsSize, dataset.dims) & dataset.rank == 2
    ]$name

    qtyList <- lapply(datasets_dt$dataset.dims, function(x) eval(parse(text = gsub("x", "*", x))))
    qty <- sum(unlist(qtyList))

    pb <- utils::txtProgressBar(min = 0, max = qty, style = 3)

    # Do clipping and copying


    clipByMask(beam, updateBeam, segmentsCut, landSegmentsMask, pb)
    clipByMask2D(beam, updateBeam, segmentsCut2D, landSegmentsMask, pb)
    clipByMask(beam, updateBeam, photonsCut, photonsMask, pb)
    updateBeam[["land_segments/ph_ndx_beg"]][] <- c(1, cumsum(updateBeam[["land_segments/n_seg_ph"]][])[-1])

    close(pb)
  }

  newFile$close_all()
  ATL08_read(output)
}

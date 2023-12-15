#' Clips ICESat-2 ATL03 data
#'
#' @param atl03 [`icesat2.atl03_h5-class`] object, obtained through [`ATL03_read()`]
#' for clipping
#' @param output character. Path to the output h5 file.
#' @param bbox [`numeric-class`] or [`terra::SpatExtent`] for clipping, the
#' order of the bbox is the default from NASA's ICESat-2 CMS searching:
#' [ul_lat, ul_lon, lr_lat, lr_lon].
#'
#' @return Returns the clipped S4 object of class [`icesat2.atl03_h5-class`]
#'
#' @description This function clips the ATl03 HDF5 file. This function
#' will only clip the beam groups within hdf5, it won't change metadata
#' or ancillary data.
#'
#' @examples
#' # Specifying the path to ICESat-2 ATL03 data (zip file)
#' outdir <- tempdir()
#' atl03_fp_zip <- system.file("extdata",
#'   "ATL0302_A_2019103030338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rICESat2Veg"
#' )
#'
#' # Unzipping ICESat-2 ATL03 data
#' atl03_path <- unzip(atl03_fp_zip, exdir = outdir)
#'
#' # Reading ICESat-2 ATL03 data (h5 file)
#' atl03 <- ATL03_read(atl03_path = atl03_path)
#'
#' close(atl03)
#' @import hdf5r
#' @export
setGeneric(
  "ATL03_h5_clipBox",
  function(atl03, output, bbox) {
    standardGeneric("ATL03_h5_clipBox")
  }
)

#' @include class.icesat2.R
#' @importClassesFrom terra SpatExtent
setMethod(
  "ATL03_h5_clipBox",
  signature = c("icesat2.atl03_h5", "character", "SpatExtent"),
  function(atl03, output, bbox) {
    clipBoxATL03(atl03, output, bbox)
  }
)

setMethod(
  "ATL03_h5_clipBox",
  signature = c("icesat2.atl03_h5", "character", "numeric"),
  function(atl03, output, bbox) {
    print("clipping by bbox")
    bbox_ext <- terra::ext(bbox[c(2, 4, 3, 1)])
    ATL03_h5_clipBox(atl03, bbox_ext)
  }
)



ATL03_photons_mask <- function(beam, bbox) {
  x <- y <- 0

  xy <- data.table::data.table(
    x = beam[["heights/lon_ph"]][],
    y = beam[["heights/lat_ph"]][]
  )

  mask <- xy[
    , x >= bbox$xmin &
      x <= bbox$xmax &
      y >= bbox$ymin &
      y <= bbox$ymax
  ]

  mask <- seq_along(mask)[mask]
  return(mask)
}

# Count number of photons per segment
ATL03_photons_per_segment <- function(beam, photonsMask) {
  seg_indices <- beam[["geolocation/ph_index_beg"]][]
  seg_indices <- seg_indices[seg_indices != 0]
  photons_segment <- findInterval(photonsMask, seg_indices)
  table(photons_segment)
}



clipBoxATL03 <- function(atl03, output, bbox) {
  dataset.rank <- dataset.dims <- obj_type <- name <- NA

  # Create a new HDF5 file
  newFile <- hdf5r::H5File$new(output, mode = "w")


  # Create all groups
  structure_dt <- data.table::as.data.table(atl03@h5$ls(recursive = T))
  groups <- structure_dt[obj_type == "H5I_GROUP"]$name


  for (group in groups) {
    grp <- newFile$create_group(group)

    # Create all atributes within group
    attributes <- hdf5r::list.attributes(atl03[[group]])
    for (attribute in attributes) {
      grp$create_attr(attribute, hdf5r::h5attr(atl03[[group]], attribute))
    }
  }

  # Create root attributes
  attributes <- hdf5r::list.attributes(atl03@h5)
  for (attribute in attributes) {
    hdf5r::h5attr(newFile, attribute) <- hdf5r::h5attr(atl03@h5, attribute)
  }

  # Get all beams
  beams <- getBeams(atl03)

  nBeam <- 0
  nBeams <- length(beams)

  # Loop the beams
  for (beamName in beams) {
    nBeam <- nBeam + 1
    message(sprintf("Clipping %s (%d/%d)", beamName, nBeam, nBeams))

    # Get the reference beam
    beam <- atl03[[beamName]]

    # Get the beam to update
    updateBeam <- newFile[[beamName]]

    # Get the masks
    photonsMask <- ATL03_photons_mask(beam, bbox)
    photons_per_segment <- ATL03_photons_per_segment(beam, photonsMask)
    segmentsMask <- as.integer(names(photons_per_segment))

    # Get sizes of clipping datasets
    photonsSize <- beam[["heights/h_ph"]]$dims
    segmentsSize <- beam[["geolocation/ph_index_beg"]]$dims

    # Get all datasets
    datasets_dt <- data.table::as.data.table(beam$ls(recursive = TRUE))[obj_type == 5]

    # Get all types of clipping photons/segment/no cut
    photonsCut <- datasets_dt[dataset.dims == photonsSize]$name
    photonsCut2D <- datasets_dt[grepl(photonsSize, dataset.dims) & dataset.rank == 2]$name
    specialCuts <- c(
      "geolocation/segment_ph_cnt",
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

    copyDataset(beam, updateBeam, "geolocation/segment_ph_cnt", photons_per_segment, pb)
    copyDataset(beam, updateBeam, "geolocation/ph_index_beg", cumsum(photons_per_segment), pb)

    for (dataset in nonCuts) {
      copyDataset(beam, updateBeam, dataset, beam[[dataset]][], pb)
    }
    close(pb)
  }
  newFile$close_all()

  ATL03_read(output)
}

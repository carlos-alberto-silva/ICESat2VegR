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
#' # Specifying the path to ICESat-2 ATL08 data (zip file)
#' outdir <- tempdir()
#' atl08_fp_zip <- system.file("extdata",
#'   "ATL0802_A_2019103030338_O01964_T05337_02_001_01_sub.zip",
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
setGeneric(
  "ATL08_h5_clipBox",
  function(atl08, atl03, output, bbox) {
    standardGeneric("ATL08_h5_clipBox")
  }
)

#' @include class.icesat2.R
#' @importClassesFrom terra SpatExtent
setMethod(
  "ATL08_h5_clipBox",
  signature = c("icesat2.atl08_h5", "icesat2.atl03_h5", "character", "SpatExtent"),
  function(atl08, atl03, output, bbox) {
    clipBoxATL08(atl08, atl03, output, bbox)
  }
)

setMethod(
  "ATL08_h5_clipBox",
  signature = c("icesat2.atl08_h5", "icesat2.atl03_h5", "character", "numeric"),
  function(atl08, atl03, output, bbox) {
    print("clipping by bbox")
    bbox_ext <- terra::ext(bbox[c(2, 4, 3, 1)])
    ATL08_h5_clipBox(atl08, bbox_ext)
  }
)

ATL08_photons_mask <- function(joined_dt, bbox) {
  lon_ph <- lat_ph <- NA

  mask <- joined_dt[
    , lon_ph >= bbox$xmin &
      lon_ph <= bbox$xmax &
      lat_ph >= bbox$ymin &
      lat_ph <= bbox$ymax
  ]

  mask <- seq_along(mask)[mask]
  return(mask)
}

clipBoxATL08 <- function(atl08, atl03, output, bbox) {
  latitude <- longitude <- obj_type <- dataset.rank <- dataset.dims <- ph_segment_id <- NA
  .N <- data.table::.N
  .I <- data.table::.I
  . <- list

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
  for (beamName in beams) {
    nBeam <- nBeam + 1
    message(sprintf("Clipping %s (%d/%d)", beamName, nBeam, nBeams))

    # Get the reference beam
    beam <- atl08[[beamName]]

    # Get the beam to update
    updateBeam <- newFile[[beamName]]

    # Get the masks
    joined_data <- ATL03_ATL08_photons_attributes_dt_join(atl03, atl08, beam = beamName)
    joined_dt <- joined_data@dt[!is.na(ph_segment_id)]
    photonsMask <- ATL08_photons_mask(joined_dt, bbox)

    photons_per_segment <- joined_dt[, .N, by = .(ph_segment_id)]
    segmentsMask <- photons_per_segment$ph_segment_id

    # Get sizes of clipping datasets

    photonsSize <- nrow(joined_dt)
    segmentsSize <- beam[["land_segments/segment_watermask"]]$dims

    # Get all datasets
    datasets_dt <- data.table::as.data.table(beam$ls(recursive = TRUE))[obj_type == 5]

    seg_dt <- ATL08_seg_attributes_dt(atl08, beam = beamName, attribute = c("n_toc_photons"))
    segmentsMask <- seg_dt[, .(latitude, longitude, .I)][
      longitude >= bbox$xmin &
        longitude <= bbox$xmax &
        latitude >= bbox$ymin &
        latitude <= bbox$ymax, I
    ]

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

    clipByMask(beam, updateBeam, segmentsCut, segmentsMask, pb)
    clipByMask2D(beam, updateBeam, segmentsCut2D, segmentsMask, pb)
    clipByMask(beam, updateBeam, photonsCut, photonsMask, pb)

    close(pb)
  }

  newFile$close_all()
  ATL08_read(output)
}

#' @include class.icesat2.R ATL08_read.R
#' @import data.table hdf5r
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

#' @include class.icesat2.R ATL08_read.R
#' @import data.table hdf5r
ATL08_h5_clip <- function(atl08, output, clip_obj, landSegmentsMask_fn) {
  dataset.rank <- dataset.dims <- name <- NA

  # Create a new HDF5 file
  newFile <- hdf5r::H5File$new(output, mode = "w")

  # Create all groups
  groups <- atl08$ls_groups(recursive = TRUE)


  for (group in groups) {
    grp <- newFile$create_group(group)

    # Create all atributes within group
    attributes <- atl08[[group]]$ls_attrs()
    for (attribute in attributes) {
      grp$create_attr(attribute, atl08[[group]]$attr(attribute))
    }
  }

  # Create root attributes
  attributes <- atl08$ls_attrs()
  for (attribute in attributes) {
    hdf5r::h5attr(newFile, attribute) <- atl08$attr(attribute)
  }

  # Get all beams
  beams <- atl08$beams


  nBeam <- 0
  nBeams <- length(beams)

  # Loop the beams
  # beamName <- beams[nBeam + 1]
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
      Vectorize(seq.default, vectorize.args = c("from", "to"))(from = ph_ndx_beg, to = ph_ndx_beg + n_seg_ph - 1)
    )

    # Get sizes of clipping datasets
    photonsSize <- beam[["signal_photons/ph_h"]]$dims
    segmentsSize <- beam[["land_segments/segment_watermask"]]$dims

    # Get all datasets
    datasets_dt <- beam$dt_datasets(recursive = TRUE)

    # Get all types of clipping photons/segment/no cut
    photonsCut <- datasets_dt[dataset.dims == photonsSize]$name

    segmentsCut <- datasets_dt[
      dataset.dims == segmentsSize
    ]$name
    segmentsCut2D <- datasets_dt[
      grepl(segmentsSize, dataset.dims) & dataset.rank == 2
    ]$name

    allCuts <- c(photonsCut, segmentsCut, segmentsCut2D)

    nonCuts <- datasets_dt[
      !name %in% allCuts
    ]$name


    qtyList <- lapply(datasets_dt$dataset.dims, function(x) eval(parse(text = gsub("x", "*", x))))
    qty <- sum(unlist(qtyList))

    pb <- utils::txtProgressBar(min = 0, max = qty, style = 3)

    # Do clipping and copying
    if (length(landSegmentsMask) == 0) {
      utils::setTxtProgressBar(pb, qty)
      close(pb)
      next
    }

    clipByMask(beam, updateBeam, segmentsCut, landSegmentsMask, pb)
    clipByMask2D(beam, updateBeam, segmentsCut2D, landSegmentsMask, pb)
    clipByMask(beam, updateBeam, photonsCut, photonsMask, pb)
    updateBeam[["land_segments/ph_ndx_beg"]][] <- c(1, cumsum(updateBeam[["land_segments/n_seg_ph"]][])[-1])

    for (dataset in nonCuts) {
      copyDataset(beam, updateBeam, dataset, beam[[dataset]][], pb)
    }

    close(pb)
  }

  newFile$close_all()
  ATL08_read(output)
}

# #' Clips ICESat-2 ATL03 H5 data
# #'
# #' @param atl03 [`ICESat2VegR::icesat2.atl03_h5-class`] object,
# #' obtained through [`ATL03_read()`] for clipping
# #' @param output character. Path to the output h5 file.
# #' @param clip_obj [`numeric-class`], [`terra::SpatExtent`] or [`terra::SpatVector`]
# #' object for clipping. The bbox order is the default from NASA's ICESat-2 CMS searching:
# #' (ul_lat, ul_lon, lr_lat, lr_lon).
# #' @param beam [`character-class`]. The vector of beams to include, default
# #' all c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")
# #' @param additional_groups [`character-class`]. Other addional groups that should be included, default
# #' c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")
# #'
# #' @return Returns the clipped S4 object of class [`ICESat2VegR::icesat2.atl03_h5-class`]
# #'
# #' @description This function clips ATL03 HDF5 file within beam groups,
# #' but keeps metada and ancillary data the same.
# #'
# #' @export
#' @import hdf5r
ATL03_h5_clip <- function(
  atl03,
  output,
  clipObj,
  mask_fn,
  beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
  additional_groups = c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")
) {
  # Create a new HDF5 file
  newFile <- hdf5r::H5File$new(output, mode = "w")

  all_groups <- c(beam, additional_groups)
  all_groups <- intersect(all_groups, atl03$ls_groups())
  starts_with_regex <- paste0("^(", paste(all_groups, collapse = "|"), ")")
  groups <- atl03$ls_groups(recursive = TRUE)
  groups <- grep(starts_with_regex, groups, value = TRUE)

  # Copy groups and attributes
  copyGroupsAndAttributes(atl03, newFile, beam, groups)

  # Copy non-beam groups and datasets
  copyNonBeamGroupsAndDatasets(atl03, newFile, beam, groups)

  # Loop through beams
  clipBeams(atl03, newFile, beam, mask_fn, clipObj)

  newFile$close_all()
  ATL03_read(output)
}

copyGroupsAndAttributes <- function(atl03, newFile, beam, groups) {
  for (group in groups) {
    grp <- newFile$create_group(group)
    attributes <- atl03[[group]]$ls_attrs()
    for (attribute in attributes) {
      grp$create_attr(attribute, atl03[[group]]$attr(attribute))
    }
  }

  attributes <- atl03$ls_attrs()
  for (attribute in attributes) {
    hdf5r::h5attr(newFile, attribute) <- atl03$attr(attribute)
  }
}

copyNonBeamGroupsAndDatasets <- function(atl03, newFile, beam, groups) {
  non_beams_groups <- grep("^gt[1-3][rl]", groups, value = TRUE, invert = TRUE)

  # non_beam_group <- "orbit_info"
  for (non_beam_group in non_beams_groups) {
    datasets_dt <- atl03[[non_beam_group]]$dt_datasets()
    # dataset <- "bounding_polygon_lat1"
    for (dataset in datasets_dt$name) {
      ds_dims <- length(atl03[[non_beam_group]][[dataset]]$dims)
      args <- list(atl03[[non_beam_group]][[dataset]])
      alist_text <- gettextf("alist(%s)", paste0(rep(",", max(0, ds_dims - 1)), collapse = ""))
      args <- c(args, eval(parse(text = alist_text)))
      ds_data <- do.call("[", args)
      copyDataset(
        beam = atl03[[non_beam_group]],
        updateBeam = newFile[[non_beam_group]],
        dataset = dataset,
        data = ds_data,
        NULL
      )
    }
  }
}

clipBeams <- function(atl03, newFile, beam, mask_fn, clipObj) {
  beams <- intersect(beam, atl03$beams)
  nBeams <- length(beams)

  for (nBeam in seq_along(beams)) {
    beamName <- beams[nBeam]
    message(sprintf("Clipping %s (%d/%d)", beamName, nBeam, nBeams))

    beam <- atl03[[beamName]]
    updateBeam <- newFile[[beamName]]

    segmentsMask <- mask_fn(beam, clipObj)
    photonsMask <- ATL03_photons_segment_mask(beam, segmentsMask)

    photonsSize <- beam[["heights/h_ph"]]$dims
    segmentsSize <- beam[["geolocation/ph_index_beg"]]$dims

    datasets_dt <- beam$dt_datasets(recursive = TRUE)

    photonsCut <- datasets_dt[dataset.dims == photonsSize]$name
    photonsCut2D <- datasets_dt[grepl(photonsSize, dataset.dims) & dataset.rank == 2]$name
    specialCuts <- c("geolocation/ph_index_beg")
    segmentsCut <- datasets_dt[
      dataset.dims == segmentsSize &
        !(name %in% specialCuts)
    ]$name
    segmentsCut2D <- datasets_dt[
      grepl(segmentsSize, dataset.dims) & dataset.rank == 2
    ]$name
    allCuts <- c(photonsCut, photonsCut2D, segmentsCut, segmentsCut2D, specialCuts)
    nonCuts <- datasets_dt[
      !(name %in% allCuts)
    ]$name

    qtyList <- lapply(datasets_dt$dataset.dims, function(x) eval(parse(text = gsub("x", "*", x))))
    qty <- sum(unlist(qtyList))

    pb <- utils::txtProgressBar(min = 0, max = qty, style = 3)

    if (length(segmentsMask) == 0) {
      utils::setTxtProgressBar(pb, qty)
      close(pb)
      next
    }
    clipByMask(beam, updateBeam, photonsCut, photonsMask, pb)
    clipByMask2D(beam, updateBeam, photonsCut2D, photonsMask, pb)
    clipByMask(beam, updateBeam, segmentsCut, segmentsMask, pb)
    clipByMask2D(beam, updateBeam, segmentsCut2D, segmentsMask, pb)

    photons_per_segment <- beam[["geolocation/segment_ph_cnt"]][segmentsMask]
    createDatasetClip(
      beam,
      updateBeam,
      "geolocation/ph_index_beg",
      c(1, cumsum(photons_per_segment)[-length(photons_per_segment)]),
      pb
    )

    for (dataset in nonCuts) {
      copyDataset(beam, updateBeam, dataset, beam[[dataset]][], pb)
    }
    close(pb)
  }
}

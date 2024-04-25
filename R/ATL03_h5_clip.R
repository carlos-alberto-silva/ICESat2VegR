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
# #' @examples
# ##' # Specifying the path to ATL03 file (zip file)
# #' outdir <- tempdir()
# #' atl03_zip <- system.file("extdata",
# #'   "atl03_20220401221822_01501506_005_01.zip",
# #'   package = "rICESat2Veg"
# #' )
# #'
# #' # Unzipping ATL03 file
# #' atl03_path <- unzip(atl03_zip, exdir = outdir)
# #'
# #' # Reading ATL03 data (h5 file)
# #' atl03_h5 <- atl03_read(atl03_path = atl03_path)
# #'
# #' # Bounding rectangle coordinates
# #' xmin <- -107.7
# #' xmax <- -106.5
# #' ymin <- 32.75
# #' ymax <- 42.75
# #'
# #' # Clipping ATL03 photons  by boundary box extent
# #' atl03_photons_dt_clip <- ATL03_h5_clipBox(atl03_h5, outdir, xmin, xmax, ymin, ymax)
# #'
# #' close(atl03_h5)
# #' @export
#' @import hdf5r
ATL03_h5_clip <- function(
    atl03,
    output,
    clipObj,
    mask_fn,
    beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")) {
  dataset.rank <- dataset.dims <- name <- NA

  # Create a new HDF5 file
  newFile <- hdf5r::H5File$new(output, mode = "w")

  all_groups <- c(beam, additional_groups)
  all_groups <- intersect(all_groups, atl03$ls_groups())
  starts_with_regex <- paste0("^(", paste(all_groups, collapse = "|"), ")")

  # Create all groups
  groups <- atl03$ls_groups(recursive = TRUE)
  groups <- grep(starts_with_regex, groups, value = TRUE)

  # Remove unselected beams groups
  not_beam <- setdiff(atl03$beams, beam)
  not_beam_regex <- paste(not_beam, collapse = "|")

  groups <- grep(not_beam_regex, groups, value = TRUE, invert = TRUE)

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

  non_beams_groups <- grep("^gt[1-3][rl]", groups, value = TRUE, invert = TRUE)

  for (non_beam_group in non_beams_groups) {
    datasets_dt <- atl03[[non_beam_group]]$dt_datasets()
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

  # Get all beams
  beams <- intersect(beam, atl03$beams)


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

    # Do clipping and copying
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
  newFile$close_all()

  ATL03_read(output)
}
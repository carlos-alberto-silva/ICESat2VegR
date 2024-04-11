#' Clips ICESat-2 ATL08 H5 data
#'
#' @param x [`icesat2.atl08_h5-class`] object,
#' obtained through [`ATL08_read()`] for clipping
#' @param output character. Path to the output h5 file.
#' @param bbox [`numeric-class`] or [`terra::SpatExtent`] for clipping, the
#' order of the bbox is the default from NASA's ICESat-2 CMS searching:
#' [ul_lat, ul_lon, lr_lat, lr_lon].
#' @param beam [`character-class`]. The vector of beams to include, default
#' all c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")
#' @param additional_groups [`character-class`]. Other addional groups that should be included, default
#' c("METADATA", "orbit_info", "quality_assessment", "ancillary_data")
#'
#' @return Returns the clipped S4 object of class [`icesat2.atl08_h5-class`]
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
#'
#' # Bounding rectangle coordinates
#' xmin <- -107.7
#' xmax <- -106.5
#' ymin <- 32.75
#' ymax <- 42.75
#'
#' # Clipping ATL08 photons  by boundary box extent
#' atl08_photons_dt_clip <- ATL08_h5_clipBox(atl08_h5, outdir, xmin, xmax, ymin, ymax)
#'
#' close(atl08_h5)
#' @include class.icesat2.R ATL08_read.R
#' @import data.table hdf5r
#' @export
ATL08_h5_clip <- function(
    atl08, output, clip_obj, landSegmentsMask_fn,
    beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")) {
  dataset.rank <- dataset.dims <- name <- NA

  # Create a new HDF5 file
  newFile <- hdf5r::H5File$new(output, mode = "w")

  all_groups <- c(beam, additional_groups)
  all_groups <- intersect(all_groups, atl08$ls_groups())
  starts_with_regex <- paste0("^(", paste(all_groups, collapse="|"), ")")

  # Create all groups
  groups <-  atl08$ls_groups(recursive = TRUE)
  groups <- grep(starts_with_regex, groups, value = TRUE)

  # Remove unselected beams groups
  not_beam <- setdiff(atl08$beams, beam)
  not_beam_regex <- paste(not_beam, collapse = "|")
  
  groups <- grep(not_beam_regex, groups, value = TRUE, invert = TRUE)

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

  non_beams_groups <- grep("^gt[1-3][rl]", groups, value = TRUE, invert = TRUE)

  for (non_beam_group in non_beams_groups) {
    datasets_dt <- atl08[[non_beam_group]]$dt_datasets()
    for (dataset in datasets_dt$name) {
      ds_dims <- length(atl08[[non_beam_group]][[dataset]]$dims)
      args <- list(atl08[[non_beam_group]][[dataset]])
      alist_text <- gettextf("alist(%s)", paste0(rep(",", max(0, ds_dims - 1)), collapse = ""))
      args <- c(args, eval(parse(text = alist_text)))
      ds_data <- do.call("[", args)
      copyDataset(
        beam = atl08[[non_beam_group]],
        updateBeam = newFile[[non_beam_group]],
        dataset = dataset,
        data = ds_data,
        NULL
      )
    }
  }
  # Get all beams
  beams <- intersect(beam, atl08$beams)


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
      createDatasetClip(beam, updateBeam, dataset, beam[[dataset]][], pb)
    }

    close(pb)
  }

  newFile$close_all()
  ATL08_read(output)
}

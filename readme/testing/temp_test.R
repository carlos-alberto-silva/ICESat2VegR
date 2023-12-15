library(rICESat2Veg)

atl03_path <- "../inst/exdata/ATL03_20220401221822_01501506_005_01.h5"
atl08_path <- "../inst/exdata/ATL08_20220401221822_01501506_005_01.h5"

atl03 <- ATL03_read(atl03_path)
atl08 <- ATL08_read(atl08_path)


ul_lat <- 59.50
ul_lon <- -106.3
lr_lat <- 26.99
lr_lon <- -105.8

output <- file.path("C:/Users/caiohamamura/Desktop/saida/", "output2.h5")

bbox <- terra::ext(c(ul_lon, lr_lon, lr_lat, ul_lat))

ATL03_h5_clipBox(atl03, output, bbox)



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

# Count number of photons per segment
ATL08_photons_per_segment <- function(beam, photonsMask) {
  seg_indices <- beam[["geolocation/ph_index_beg"]][]
  seg_indices <- seg_indices[seg_indices != 0]
  photons_segment <- findInterval(photonsMask, seg_indices)
  table(photons_segment)
}


clipByMask <- function(beam, updateBeam, datasets, mask, pb) {
  maskSize <- length(mask)

  for (dataset in datasets) {
    h5_dtype <- beam[[dataset]]$get_type()
    h5_pl <- beam[[dataset]]$get_create_plist()
    plist <- hdf5r::H5P_DATASET_CREATE$new()
    plist$set_filter(h5_pl$get_filter(0)$filter)
    plist$set_fill_value(h5_dtype, h5_pl$get_fill_value(dtype = h5_dtype))
    updateBeam$create_dataset(
      name = dataset,
      robj = beam[[dataset]][mask],
      dims = c(maskSize),
      dtype = beam[[dataset]]$get_type(),
      dataset_create_pl = plist,
      chunk_dim = beam[[dataset]]$chunk_dims
    )
    utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + exp(sum(log(beam[[dataset]]$dims))))
  }
}

clipByMask2D <- function(beam, updateBeam, datasets, mask, pb) {
  maskSize <- length(mask)
  for (dataset in datasets) {
    h5_dtype <- beam[[dataset]]$get_type()
    h5_pl <- beam[[dataset]]$get_create_plist()
    plist <- hdf5r::H5P_DATASET_CREATE$new()
    plist$set_filter(h5_pl$get_filter(0)$filter)
    plist$set_fill_value(h5_dtype, h5_pl$get_fill_value(dtype = h5_dtype))
    updateBeam$create_dataset(
      name = dataset,
      robj = beam[[dataset]][, mask],
      dims = c(beam[[dataset]]$dims[1], length(mask)),
      dtype = beam[[dataset]]$get_type(),
      dataset_create_pl = plist,
      chunk_dim = beam[[dataset]]$chunk_dims
    )

    utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + exp(sum(log(beam[[dataset]]$dims))))
  }
}

copyDataset <- function(beam, updateBeam, dataset, data, pb) {
  h5_dtype <- beam[[dataset]]$get_type()
  h5_pl <- beam[[dataset]]$get_create_plist()
  plist <- hdf5r::H5P_DATASET_CREATE$new()
  plist$set_filter(h5_pl$get_filter(0)$filter)
  plist$set_fill_value(h5_dtype, h5_pl$get_fill_value(dtype = h5_dtype))
  updateBeam$create_dataset(
    name = dataset,
    robj = data,
    dims = c(length(data)),
    dtype = beam[[dataset]]$get_type(),
    dataset_create_pl = plist,
    chunk_dim = beam[[dataset]]$chunk_dims
  )

  utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + exp(sum(log(beam[[dataset]]$dims))))
}

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


# Loop the beams

# Comment this
# beamName <- beams[6]
nBeam <- 0
nBeams <- length(beams)

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
    latitude <= bbox$ymax, I]

  # Get all types of clipping photons/segment/no cut
  photonsCut <- datasets_dt[dataset.dims == photonsSize]$name

  segmentsCut <- datasets_dt[
    dataset.dims == segmentsSize
  ]$name
  segmentsCut2D <- datasets_dt[
    grepl(segmentsSize, dataset.dims) & dataset.rank == 2
  ]$name

  allCuts <- c(photonsCut, segmentsCut, segmentsCut2D)

  qtyList <- lapply(datasets_dt$dataset.dims, function(x) eval(parse(text = gsub("x", "*", x))))
  qty <- sum(unlist(qtyList))

  pb <- utils::txtProgressBar(min = 0, max = qty, style = 3)

  # Do clipping and copying

  datasets <- photonsCut
  mask <- photonsMask

  clipByMask(beam, updateBeam, segmentsCut, segmentsMask, pb)
  clipByMask2D(beam, updateBeam, segmentsCut2D, segmentsMask, pb)
  clipByMask(beam, updateBeam, photonsCut, photonsMask, pb)

  close(pb)
}

newFile$close_all()

# # Copy the rest of the datasets
# newFile <- hdf5r::H5File$new(output, mode = "r")
# head(newFile[[beamName]][[dataset]][])

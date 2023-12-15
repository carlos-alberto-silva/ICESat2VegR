library(rICESat2Veg)

atl03_path <- "../inst/exdata/ATL03_20220401221822_01501506_005_01.h5"

atl03 <- ATL03_read(atl03_path)

ul_lat <- 59.50
ul_lon <- -106.3
lr_lat <- 26.99
lr_lon <- -104.8


bbox_numeric <- c(ul_lat, ul_lon, lr_lat, lr_lon)
bbox_spatextent <- terra::ext(c(ul_lon, lr_lon, lr_lat, ul_lat))
bbox <- bbox_spatextent

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

  return(mask)
}

# Count number of photons per segment
ATL03_photons_per_segment <- function(beam, photonsMask) {
  photons_within <- seq_along(photonsMask)[photonsMask]
  seg_indices <- beam[["geolocation/ph_index_beg"]][]
  seg_indices <- seg_indices[seg_indices != 0]
  photons_segment <- findInterval(photons_within, seg_indices)
  table(photons_segment)
}

clipByMask <- function(beam, updateBeam, datasets, mask) {
  for (dataset in datasets) {
    beam[[dataset]][] <- beam[[dataset]][]
  }
}

clipByMask2D <- function() {
  
}

# Create a new HDF5 file

# Create all groups

# Create all atributes

# Get all beams
beams <- getBeams(atl03)


# Loop the beams
{
  # beamName = beams[[1]]
  # Get the reference beam
  beam <- atl03[[beamName]]

  # Get the beam to update
  # updateBeam <- newATL03[[beamName]]

  # Get the masks
  photonsMask <- ATL03_photons_mask(beam, bbox)
  photons_per_segment <- ATL03_photons_per_segment(beam, photonsMask)
  segmentsMask <- as.integer(names(photons_per_segment))

  # Get sizes of clipping datasets
  photonsSize <- beam[["heights/h_ph"]]$dims
  segmentsSize <- beam[["geolocation/ph_index_beg"]]$dims

  # Get all datasets
  datasets_dt <- data.table::as.data.table(beam$ls(recursive = T))[obj_type == 5]

  # Get all types of clipping photons/segment/no cut
  photonsCut <- datasets_dt[dataset.dims == photonsSize]$name
  photonsCut2D <- datasets_dt[grepl(photonsSize, dataset.dims) & dataset.rank == 2]$name
  specialCuts <- c(
    "geolocation/segment_ph_cnt",
    "geolocation/ph_index_beg",
    "geolocation/reference_photon_index"
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


  # Do clipping and copying
  clipByMask(beam, updateBeam, segmentsCut, segmentsMask)
  clipByMask2D(beam, updateBeam, segmentsCut2D, segmentsMask)
  clipByMask(beam, updateBeam, photonsCut, photonsMask)
  clipByMask2D(beam, updateBeam, photonsCut2D, photonsMask)

  # Copy the rest of the datasets
}

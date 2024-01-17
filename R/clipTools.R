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

ATL03_photons_mask_geometry <- function(beam, geom) {
  ext <- terra::ext(geom)
  mask_ext <- ATL03_photons_mask(beam, ext)

  xy <- terra::vect(
    data.frame(
      longitude = beam[["heights/lon_ph"]][mask_ext],
      latitude = beam[["heights/lat_ph"]][mask_ext],
      row = mask_ext
    ),
    geom = c("longitude", "latitude"),
    crs = "epsg:4326"
  )

  res <- terra::intersect(xy, geom)

  return(res$row)
}

ATL03_photons_per_segment <- function(beam, photonsMask) {
  seg_indices <- beam[["geolocation/ph_index_beg"]][]
  seg_indices <- seg_indices[seg_indices != 0]
  photons_segment <- findInterval(photonsMask, seg_indices)
  table(photons_segment)
}

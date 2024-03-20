clipByMask <- function(beam, updateBeam, datasets, mask, pb) {
  maskSize <- length(mask)

  for (dataset in datasets) {
    h5_pl <- beam[[dataset]]$get_creation_property_list()
    updateBeam$create_dataset(
      name = dataset,
      robj = beam[[dataset]][mask],
      dims = c(maskSize),
      dtype = beam[[dataset]]$get_type(),
      dataset_create_pl = h5_pl,
      chunk_dim = beam[[dataset]]$chunk_dims
    )
    utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + exp(sum(log(beam[[dataset]]$dims))))
  }
}

clipByMask2D <- function(beam, updateBeam, datasets, mask, pb) {
  maskSize <- length(mask)
  for (dataset in datasets) {
    h5_pl <- beam[[dataset]]$get_creation_property_list()
    updateBeam$create_dataset(
      name = dataset,
      robj = beam[[dataset]][, mask],
      dims = c(beam[[dataset]]$dims[1], maskSize),
      dtype = beam[[dataset]]$get_type(),
      dataset_create_pl = h5_pl,
      chunk_dim = beam[[dataset]]$chunk_dims
    )

    utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + exp(sum(log(beam[[dataset]]$dims))))
  }
}

copyDataset <- function(beam, updateBeam, dataset, data, pb) {
  h5_pl <- beam[[dataset]]$get_creation_property_list()
  updateBeam$create_dataset(
    name = dataset,
    robj = data,
    dims = c(length(data)),
    dtype = beam[[dataset]]$get_type(),
    dataset_create_pl = h5_pl,
    chunk_dim = beam[[dataset]]$chunk_dims
  )

  utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + exp(sum(log(beam[[dataset]]$dims))))
}

clipByMask2D <- function(beam, updateBeam, datasets, mask, pb) {
  maskSize <- length(mask)
  for (dataset in datasets) {
    h5_pl <- beam[[dataset]]$get_creation_property_list()
    updateBeam$create_dataset(
      name = dataset,
      robj = beam[[dataset]][, mask],
      dims = c(beam[[dataset]]$dims[1], length(mask)),
      dtype = beam[[dataset]]$get_type(),
      dataset_create_pl = h5_pl,
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

ATL03_segments_mask <- function(beam, bbox) {
  x <- y <- 0

  xy <- data.table::data.table(
    x = beam[["geolocation/reference_photon_lon"]][],
    y = beam[["geolocation/reference_photon_lat"]][]
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

ATL03_segments_mask_geometry <- function(beam, geom) {
  ext <- terra::ext(geom)
  mask_ext <- ATL03_segments_mask(beam, ext)

  xy <- terra::vect(
    data.frame(
      longitude = beam[["geolocation/reference_photon_lon"]][mask_ext],
      latitude = beam[["geolocation/reference_photon_lat"]][mask_ext],
      row = mask_ext
    ),
    geom = c("longitude", "latitude"),
    crs = "epsg:4326"
  )

  res <- terra::intersect(xy, geom)

  return(res$row)
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

seq_lens <- Vectorize(function(from, length.out) {seq.int(from, length.out = length.out)}, vectorize.args = c("from", "length.out"), SIMPLIFY = TRUE)
seq_lens_simplify <- function(from, length.out) {
  unlist(seq_lens(from, length.out))
}

ATL03_photons_segment_mask <- function(beam, segmentsMask) {
  seg_indices <- beam[["geolocation/ph_index_beg"]][segmentsMask]
  ph_cnt <- beam[["geolocation/segment_ph_cnt"]][segmentsMask]
  photons_mask <- seq_lens_simplify(seg_indices, ph_cnt)
  photons_mask
}


ATL03_photons_per_segment <- function(beam, photonsMask) {
  seg_indices <- beam[["geolocation/ph_index_beg"]][]
  seg_indices <- seg_indices[seg_indices != 0]
  photons_segment <- findInterval(photonsMask, seg_indices)
  table(photons_segment)
}

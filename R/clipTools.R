clipByMask <- function(beam, updateBeam, datasets, mask, pb) {
  maskSize <- length(mask)

  for (dataset in datasets) {
    h5_pl <- beam[[dataset]]$get_creation_property_list()
    robj <- beam[[dataset]][mask]

    updateBeam$create_dataset(
      name = dataset,
      robj = robj,
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
    robj <- beam[[dataset]][, mask]
    updateBeam$create_dataset(
      name = dataset,
      robj = robj,
      dims = c(beam[[dataset]]$dims[1], maskSize),
      dtype = beam[[dataset]]$get_type(),
      dataset_create_pl = h5_pl,
      chunk_dim = beam[[dataset]]$chunk_dims
    )

    utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + exp(sum(log(beam[[dataset]]$dims))))
  }
}

createDatasetClip <- function(beam, updateBeam, dataset, data, pb) {
  h5_pl <- beam[[dataset]]$get_creation_property_list()
  res <- try({
    updateBeam$create_dataset(
      name = dataset,
      robj = data,
      dims = if (is.null(dim(data))) {
        length(data)
      } else {
        dim(data)
      },
      dtype = beam[[dataset]]$get_type(),
      dataset_create_pl = h5_pl,
      chunk_dim = if (any(is.na(beam[[dataset]]$chunk_dims))) {
        beam[[dataset]]$dims
      } else {
        beam[[dataset]]$chunk_dims
      }
    )
  })
}

copyDataset <- function(beam, updateBeam, dataset, data, pb) {
  h5_pl <- beam[[dataset]]$get_creation_property_list()
  res <- try(
    {
      updateBeam$create_dataset(
        name = dataset,
        robj = data,
        dims = beam[[dataset]]$dims,
        dtype = beam[[dataset]]$get_type(),
        dataset_create_pl = h5_pl,
        chunk_dim = if (any(is.na(beam[[dataset]]$chunk_dims))) {
          beam[[dataset]]$dims
        } else {
          beam[[dataset]]$chunk_dims
        }
      )
    },
    silent = TRUE
  )
  if (inherits(res, "try-error")) {
    res2 <- try({
      updateBeam$create_dataset(
        name = dataset,
        robj = data,
        dims = beam[[dataset]]$dims,
        dtype = beam[[dataset]]$get_type(),
        chunk_dim = if (any(is.na(beam[[dataset]]$chunk_dims))) {
          beam[[dataset]]$dims
        } else {
          beam[[dataset]]$chunk_dims
        }
      )
    })

    if (inherits(res2, "try-error")) {
      message("Still can't create ", dataset)
    } else {
      message("Success after failure in ", dataset)
    }
  }

  if (!is.null(pb)) {
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



ATL03_photons_segment_mask <- function(beam, segmentsMask) {
  if (is.null(segmentsMask)) {
    return(integer(0))
  }

  seg_indices <- beam[["geolocation/ph_index_beg"]][segmentsMask]
  ph_cnt <- beam[["geolocation/segment_ph_cnt"]][segmentsMask]
  mask2 <- seg_indices > 0
  photons_mask <- pkg_module$seq_lens_simplify(seg_indices[mask2], ph_cnt[mask2])
  photons_mask
}


ATL03_photons_per_segment <- function(beam, photonsMask) {
  seg_indices <- beam[["geolocation/ph_index_beg"]][]
  seg_indices <- seg_indices[seg_indices != 0]
  photons_segment <- findInterval(photonsMask, seg_indices)
  table(photons_segment)
}

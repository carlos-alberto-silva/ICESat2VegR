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

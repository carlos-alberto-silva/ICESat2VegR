def_co <- c(
  "COMPRESS=DEFLATE",
  "BIGTIFF=IF_SAFER",
  "TILED=YES",
  "BLOCKXSIZE=512",
  "BLOCKYSIZE=512"
)

finalizer_default <- list(
  sd = ~ ifelse(n > 1, sqrt(variance / (n - 1)), NA)
  # skew = ~ sqrt((n * (n - 1))) * ((sqrt(n) * M3) / (variance^1.5)) / (n - 2),
  # kur = ~ ((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * ((n * M4) / (variance^2) - 3.0) + 6)
)

agg_function_default <- function(x) {
  list(
    n = length(x),
    mean = mean(x, na.rm = TRUE),
    variance = var(x) * (length(x) - 1),
    # M3 = e1071::moment(x, order = 3, center = TRUE, na.rm = TRUE) * length(x),
    # M4 = e1071::moment(x, order = 4, center = TRUE, na.rm = TRUE) * length(x),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

agg_join_default <- function(x1, x2) {
  combined <- data.table::data.table()
  x1$n[is.na(x1$n)] <- 0
  x1$mean[is.na(x1$mean)] <- 0
  x1$variance[is.na(x1$variance)] <- 0
  # x1$M3[is.na(x1$M3)] <- 0
  # x1$M4[is.na(x1$M4)] <- 0
  x1$max[is.na(x1$max)] <- -Inf
  x1$min[is.na(x1$min)] <- Inf

  combined$n <- x1$n + x2$n

  delta <- x2$mean - x1$mean
  delta2 <- delta * delta

  combined$mean <- (x1$n * x1$mean + x2$n * x2$mean) / combined$n

  combined$variance <- x1$variance + x2$variance +
    delta2 * x1$n * x2$n / combined$n

  # combined$M3 <- x1$M3 + x2$M3 +
  #   delta3 * x1$n * x2$n * (x1$n - x2$n) / (combined$n * combined$n)
  # combined$M3 <- combined$M3 + 3.0 * delta * (x1$n * x2$variance - x2$n * x1$variance) / combined$n

  # combined$M4 <- x1$M4 + x2$M4 + delta4 * x1$n * x2$n * (x1$n * x1$n - x1$n * x2$n + x2$n * x2$n) /
  #   (combined$n * combined$n * combined$n)
  # combined$M4 <- combined$M4 + 6.0 * delta2 * (x1$n * x1$n * x2$variance + x2$n * x2$n * x1$variance) / (combined$n * combined$n) +
  #   4.0 * delta * (x1$n * x2$M3 - x2$n * x1$M3) / combined$n

  combined$min <- pmin(x1$min, x2$min, na.rm = F)
  combined$max <- pmax(x1$max, x2$max, na.rm = F)
  return(combined)
}


# Function to convert from x, y integers to single index
xy_to_index <- function(x, y, ysize) {
  index <- (x - 1) * ysize + (y - 1) + 1
  return(index)
}

# Function convert back from single index to x, y integers
index_to_xy <- function(index, ysize) {
  index <- index - 1
  x <- floor(index / ysize) + 1
  y <- index %% ysize + 1
  return(list(x = x, y = y))
}

#' Rasterizes the model prediction saved in the HDF5 file
#'
#' This is used after running the prediction using [`predict_h5()`]
#' function to rasterize and aggregate the prediction within raster cells.
#' By default it will calculate (n, mean, variance * (n - 1), min, max, sd)
#' in a single raster file with those 6 bands in that order.
#' You can modify this behavior by changing the `agg_function`, `agg_join`
#' and `finalizer` parameters, see details section.
#'
#' @param h5_input The input HDF5 file path
#' @param output The output raster file path
#' @param bbox The bounding box of the raster
#' @param res The resolution of the raster
#' @param ... Additional parameters (see details section)
#'
#' @details
#' This function will create five different aggregate statistics
#' (n, mean, variance, min, max).
#'
#' Within `...` additional parameters we can use:
#'
#' `agg_function`: is a formula which return a data.table with the
#' aggregate function to perform over the data.
#' The default is:
#'
#' ```
#' ~data.table(
#'     n = length(x),
#'     mean = mean(x,na.rm = TRUE),
#'     variance = var(x) * (length(x) - 1),
#'     min = min(x, na.rm=T),
#'     max = max(x, na.rm=T)
#'   )
#' ```
#'
#' `agg_join`:  is a function to merge two data.table aggregates
#' from the `agg_function`. Since the h5 files will be aggregated
#' in chunks to avoid memory overflow, the statistics from the
#' different chunks should have a function to merge them.
#'
#' The default function is:
#'
#' ```
#' function(x1, x2) {
#'     combined = data.table()
#'     x1$n[is.na(x1$n)] = 0
#'     x1$mean[is.na(x1$mean)] = 0
#'     x1$variance[is.na(x1$variance)] = 0
#'     x1$max[is.na(x1$max)] = -Inf
#'     x1$min[is.na(x1$min)] = Inf
#'
#'     combined$n = x1$n + x2$n
#'
#'     delta = x2$mean - x1$mean
#'     delta2 = delta * delta
#'
#'     combined$mean = (x1$n * x1$mean + x2$n * x2$mean) / combined$n
#'     combined$variance = x1$variance + x2$variance +
#'       delta2 * x1$n * x2$n / combined$n
#'
#'     combined$min = pmin(x1$min, x2$min, na.rm=F)
#'     combined$max = pmax(x1$max, x2$max, na.rm=F)
#'     return(combined)
#' }
#' ```
#'
#' `finalizer`: is a list of formulas to generate the final
#' rasters based on the intermediate statistics from the previous
#' functions. The default `finalizer` will calculate the `sd`,
#' based on the `variance` and `n` values. It is defined as:
#'
#' ```
#' list(
#'   sd = ~sqrt(variance/(n - 1)),
#' )
#' ```
#'
#' @examples
#' atl08_path <- system.file(
#'   "extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#' atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)
#'
#' xmin <- min(atl08_dt$longitude)
#' xmax <- max(atl08_dt$longitude)
#' ymin <- min(atl08_dt$latitude)
#' ymax <- max(atl08_dt$latitude)
#'
#' linear_model <- stats::lm(h_canopy ~ canopy_openness, data = atl08_dt)
#' output_h5 <- tempfile(fileext = ".h5")
#' predicted_h5 <- predict_h5(linear_model, atl08_dt, output_h5)
#' output_raster <- tempfile(fileext = ".tif")
#'
#' rasterize_h5(
#'   predicted_h5,
#'   output = output_raster,
#'   bbox = terra::ext(xmin, xmax, ymin, ymax),
#'   res = 0.003
#' )
#'
#' @export
setGeneric("rasterize_h5", function(
    h5_input,
    output,
    bbox,
    res,
    ...) {
  standardGeneric("rasterize_h5")
})


#' Rasterizes the model prediction saved in the HDF5 file
#'
#' @param h5_input The input HDF5 file path
#' @param output The output raster file path
#' @param bbox The bounding box of the raster
#' @param res The resolution of the raster
#' @param chunk_size The chunk size to read the HDF5 file
#' @param agg_function The function to aggregate the data
#' @param agg_join The function to join the aggregated data
#' @param finalizer The function to finalize the raster
#'
#' @export
setMethod("rasterize_h5",
  signature = c("icesat2.predict_h5", "character", "SpatExtent", "numeric"),
  function(h5_input,
           output,
           bbox,
           res,
           chunk_size = 512 * 512,
           agg_function = agg_function_default,
           agg_join = agg_join_default,
           finalizer = finalizer_default) {
    if (length(res) == 1) {
      res <- c(res, -res)
    }

    `.` <- NA

    max_size <- length(h5_input)
    xsize <- ceiling((bbox$xmax - bbox$xmin) / abs(res[1]))
    ysize <- ceiling((bbox$ymax - bbox$ymin) / abs(res[2]))

    temp_h5_path <- tempfile(fileext = ".h5")
    temp_h5 <- hdf5r::H5File$new(temp_h5_path, "w")

    nodata <- list("numeric" = NA_real_, "integer" = NA_integer_)[[
      class(h5_input[["prediction"]]$get_fill_value())
    ]]

    stats_names <- names(agg_function(0))
    for (name in stats_names) {
      temp_h5[[name]] <- rep(nodata, xsize * ysize)
    }

    temp_h5[["x"]] <- rep(seq(ysize), xsize)
    temp_h5[["y"]] <- rep(seq(xsize), each = ysize)

    for (chunk in seq(1, max_size, by = chunk_size)) {
      chunk_lenght <- as.integer(min(chunk_size, max_size - chunk + 1))
      chunk_range <- seq(from = chunk, length.out = chunk_lenght)

      x <- as.integer(floor((h5_input[["longitude"]][chunk_range] - bbox$xmin) / res[1]) + 1)
      y <- ysize + as.integer(ceiling((h5_input[["latitude"]][chunk_range] - bbox$ymin) / res[2]))

      index <- xy_to_index(x, y, ysize)
      values <- h5_input[["prediction"]][chunk_range]
      dt <- data.table::data.table(index, values)
      agg_dt <- dt[, agg_function(values), by = .(index)]

      stats_names <- setdiff(names(agg_dt), "index")

      prev_values <- lapply(stats_names, function(x) {
        temp_h5[[x]][agg_dt$index]
      })
      names(prev_values) <- stats_names
      dt_agg_join <- agg_join(prev_values, agg_dt)

      for (name in stats_names) {
        temp_h5[[name]][agg_dt$index] <- dt_agg_join[[name]]
      }

      temp_h5$close_all()
      temp_h5 <- hdf5r::H5File$new(temp_h5_path)
    }

    finalizer_bands <- names(finalizer)
    output_raster <- createDataset(
      output,
      nbands = length(stats_names) + length(finalizer_bands),
      projstring = "epsg:4326",
      datatype = GDALDataType$GDT_Float32,
      lr_lat = bbox$ymin,
      ul_lat = bbox$ymax,
      lr_lon = bbox$xmax,
      ul_lon = bbox$xmin,
      res = res,
      nodata = -3.402823466E38,
      co = def_co
    )

    band1 <- output_raster[[1]]

    block_xsize <- band1$GetBlockXSize()
    block_ysize <- band1$GetBlockYSize()
    xsize <- band1$GetXSize()
    ysize <- band1$GetYSize()
    n_x_blocks <- ceiling(xsize / block_xsize)
    n_y_blocks <- ceiling(ysize / block_ysize)

    band_index <- 0
    for (stat_name in stats_names) {
      band_index <- band_index + 1
      for (x_block_index in seq_len(n_x_blocks)) {
        for (y_block_index in seq_len(n_y_blocks)) {
          x_off <- (x_block_index - 1) * block_xsize
          y_off <- (y_block_index - 1) * block_ysize

          y_range <- seq(y_off + 1, min(y_off + block_ysize, ysize))
          x_range <- seq(x_off + 1, min(x_off + block_xsize, xsize))

          y_range_rep <- rep(y_range, length(x_range))
          x_range_rep <- rep(x_range, each = length(y_range))

          indexes_h5 <- xy_to_index(x_range_rep, y_range_rep, ysize)
          xy_indexes <- index_to_xy(indexes_h5, ysize)

          indexes_write_raster <- (xy_indexes$y + y_off - 1) * block_xsize + ((xy_indexes$x - 1) - x_off)
          values_write <- rep(NA, block_xsize * block_ysize)
          values_write[indexes_write_raster + 1] <- temp_h5[[stat_name]][indexes_h5]

          band <- output_raster[[band_index]]
          band[[x_block_index - 1, y_block_index - 1]] <- values_write
        }
      }
    }

    bands <- lapply(seq_along(stats_names), function(ii) output_raster[[ii]])
    names(bands) <- stats_names

    for (finalizer_band in finalizer_bands) {
      band_index <- band_index + 1
      formula <- finalizer[[finalizer_band]]
      formulaCalculate(
        formula,
        data = bands,
        updateBand = output_raster[[band_index]]
      )
    }

    close(output_raster)
  }
)
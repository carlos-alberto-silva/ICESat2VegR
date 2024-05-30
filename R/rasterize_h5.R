def_co <- c(
  "COMPRESS=DEFLATE",
  "BIGTIFF=IF_SAFER",
  "TILED=YES",
  "BLOCKXSIZE=512",
  "BLOCKYSIZE=512"
)

finalizer <- list(
  sd = ~ sqrt(M2 / (n - 1))
  # skew = ~ sqrt((n * (n - 1))) * ((sqrt(n) * M3) / (M2^1.5)) / (n - 2),
  # kur = ~ ((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * ((n * M4) / (M2^2) - 3.0) + 6)
)

agg_function <- function(x) {
  list(
    n = length(x),
    M1 = mean(x, na.rm = TRUE),
    M2 = var(x) * (length(x) - 1),
    # M3 = e1071::moment(x, order = 3, center = TRUE, na.rm = TRUE) * length(x),
    # M4 = e1071::moment(x, order = 4, center = TRUE, na.rm = TRUE) * length(x),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

agg_join <- function(x1, x2) {
  combined <- data.table()
  x1$n[is.na(x1$n)] <- 0
  x1$M1[is.na(x1$M1)] <- 0
  x1$M2[is.na(x1$M2)] <- 0
  # x1$M3[is.na(x1$M3)] <- 0
  # x1$M4[is.na(x1$M4)] <- 0
  x1$max[is.na(x1$max)] <- -Inf
  x1$min[is.na(x1$min)] <- Inf

  combined$n <- x1$n + x2$n

  delta <- x2$M1 - x1$M1
  delta2 <- delta * delta
  delta3 <- delta * delta2
  delta4 <- delta2 * delta2

  combined$M1 <- (x1$n * x1$M1 + x2$n * x2$M1) / combined$n

  combined$M2 <- x1$M2 + x2$M2 +
    delta2 * x1$n * x2$n / combined$n

  # combined$M3 <- x1$M3 + x2$M3 +
  #   delta3 * x1$n * x2$n * (x1$n - x2$n) / (combined$n * combined$n)
  # combined$M3 <- combined$M3 + 3.0 * delta * (x1$n * x2$M2 - x2$n * x1$M2) / combined$n

  # combined$M4 <- x1$M4 + x2$M4 + delta4 * x1$n * x2$n * (x1$n * x1$n - x1$n * x2$n + x2$n * x2$n) /
  #   (combined$n * combined$n * combined$n)
  # combined$M4 <- combined$M4 + 6.0 * delta2 * (x1$n * x1$n * x2$M2 + x2$n * x2$n * x1$M2) / (combined$n * combined$n) +
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


#' @export
setGeneric("rasterize_h5", function(
    h5_input,
    output,
    bbox,
    res,
    chunk_size = 512 * 512,
    agg_function = agg_function,
    agg_join = agg_join,
    finalizer = finalizer) {
  standardGeneric("rasterize_h5")
})



#' @export
setMethod("rasterize_h5",
  signature = c("icesat2.predict_h5", "character", "SpatExtent", "numeric"),
  function(h5_input, output, bbox, res, chunk_size = 512 * 512) {
    if (length(res) == 1) {
      res <- c(res, -res)
    }

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
      # temp_h5$link_delete(name)
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

    output_raster <- gdalBindings::createDataset(
      output,
      nbands = length(stats_names),
      projstring = "epsg:4326",
      datatype = gdalBindings::GDALDataType$GDT_Float32,
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
    # stat_name <- "M1"
    # x_block_index <- 1
    # y_block_index <- 1
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
    close(output_raster)
  }
)


#' @export
"close.icesat2.predict_h5" <- function(con, ...) {
  con$close_all()
}

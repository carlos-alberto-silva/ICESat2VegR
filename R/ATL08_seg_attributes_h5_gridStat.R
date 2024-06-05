def_co <- c(
  "COMPRESS=DEFLATE",
  "BIGTIFF=IF_SAFER",
  "TILED=YES",
  "BLOCKXSIZE=512",
  "BLOCKYSIZE=512"
)

default_finalizer <- list(
  sd = ~ sqrt(variance / (n - 1))
  # skew = ~ sqrt((n * (n - 1))) * ((sqrt(n) * M3) / (variance^1.5)) / (n - 2),
  # kur = ~ ((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * ((n * M4) / (variance^2) - 3.0) + 6)
)

default_agg_function <- ~ data.table::data.table(
  n = length(x),
  mean = mean(x, na.rm = TRUE),
  variance = var(x) * (length(x) - 1),
  # M3 = e1071::moment(x, order = 3, center = TRUE, na.rm = TRUE) * length(x),
  # M4 = e1071::moment(x, order = 4, center = TRUE, na.rm = TRUE) * length(x),
  min = min(x, na.rm = TRUE),
  max = max(x, na.rm = TRUE)
)

default_agg_join <- function(x1, x2) {
  combined <- data.table()
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
  # delta3 <- delta * delta2
  # delta4 <- delta2 * delta2

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

#' Rasterize ATL08 canopy attributes from h5 files at large scale
#'
#' @description
#' This function will read multiple ATL08 H5 files and create a stack of raster layers: count, and 1st, 2nd, 3rd and 4th moments (count, m1, m2, m3 and m4) for each metric selected, from which we can calculate statistics such as Mean, SD, Skewness and Kurtosis.
#'
#' @param atl08_dir CharacterVector. The directory paths where the ATL08 H5 files are stored;
#' @param metrics CharacterVector. A vector of canopy attributes available from ATL08 product (e.g. "h_canopy")
#' @param out_root Character. The root name for the raster output files, the pattern is
#' \{out_root\}_\{metric\}_\{count/m1/m2\}.tif. This should include the full path for the file.
#' @param beam Character vector indicating beams to process (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#' @param ul_lat Numeric. Upper left latitude for the bounding box
#' @param ul_lon Numeric. Upper left longitude for the bounding box
#' @param lr_lat Numeric. Lower right latitude for the bounding box
#' @param lr_lon Numeric. Lower right longitude for the bounding box
#' @param res NumericVector. Resolution lon lat for the output raster in coordinates decimal degrees
#' @param creation_options CharacterVector. The GDAL creation options for the tif file. Default c("COMPRESS=PACKBITS", "BIGTIFF=IF_SAFER", "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512") will create BIGTIFF if needed, with DEFLATE compression and tiled by 512x512 pixels.
#' @param agg_function Formula function-like. An aggregate function which should return a data.table with the aggregate statistics
#' @param agg_join Function. A function to merge two different agg objects.
#' @param finalizer List<name, formula>. A list with the final raster names and the formula which uses the base statistics.
#'
#' @details
#' This function will create five different aggregate statistics
#' (n, mean, variance, min, max).
#' One can calculate mean and standard deviation with the following
#' formulas according to Terriberry (2007) and
#' \insertCite{Joanes1998;textual}{ICESat2VegR}:
#'
#' The `agg_function` is a formula which return a data.table with the
#' aggregate function to perform over the data.
#' The default is:
#'
#' ```
#' ~data.table(
#'     n = length(x),
#'     mean = mean(x,na.rm = TRUE),
#'     var = var(x) * (length(x) - 1),
#'     min = min(x, na.rm=T),
#'     max = max(x, na.rm=T)
#'   )
#' ```
#'
#' The `agg_join` is a function to merge two data.table aggregates
#' from the `agg_function`. Since the h5 files will be aggregated
#' one by one, the statistics from the different h5 files should
#' have a function to merge them. The default function is:
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
#' The `finalizer` is a list of formulas to generate the final
#' rasters based on the intermediate statistics from the previous
#' functions. The default `finalizer` will calculate the `sd`,
#' `skewness` and `kurtosis` based on the `variance`, `M3`, `M4` and `n`
#' values. It is defined as:
#'
#' ```
#' list(
#'   sd = ~sqrt(variance/(n - 1)),
#' )
#' ```
#'
#' @references
#' \insertAllCited{}
#'
#' Terriberry, Timothy B. (2007), Computing Higher-Order Moments Online, archived from the original on 23 April 2014, retrieved 5 May 2008
#'
#' @return Nothing. It outputs multiple raster tif files to the out_root specified path.
#'
#' @examples
#' # Specifying the path to GEDI leveatl08_canopy_dt data (zip file)
#' library(ICESat2VegR)
#' library(data.table)
#'
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Bounding rectangle coordinates
#' ul_lat <- 41.5386848449707031
#' ul_lon <- -106.5708541870117188
#' lr_lat <- 41.5314979553222656
#' lr_lon <- -106.5699081420898438
#'
#' res <- 100 # meters
#' lat_to_met_factor <- 1 / 110540
#' lon_to_met_factor <- 1 / 111320
#' xres <- lon_to_met_factor * res
#' yres <- lat_to_met_factor * res
#'
#' agg_function <- ~ data.table(
#'   min = min(x),
#'   max = max(x),
#'   sum = sum(x),
#'   n = length(x)
#' )
#'
#' agg_join <- function(agg1, agg2) {
#'   agg1[is.na(agg1)] <- 0
#'   data.table(
#'     min = pmin(agg1$min, agg2$min),
#'     max = pmax(agg1$max, agg2$max),
#'     sum = agg1$sum + agg2$sum,
#'     n = agg1$n + agg2$n
#'   )
#' }
#'
#' finalizer <- list(
#'   mean = "sum/n",
#'   range = "max-min"
#' )
#'
#' ATL08_seg_attributes_h5_gridStat(
#'   atl08_dir = dirname(atl08_path),
#'   metrics = c("h_canopy"),
#'   out_root = tempdir(),
#'   ul_lat = ul_lat,
#'   ul_lon = ul_lon,
#'   lr_lat = lr_lat,
#'   lr_lon = lr_lon,
#'   res = c(xres, -yres),
#'   creation_options = c(
#'     "COMPRESS=DEFLATE",
#'     "BIGTIFF=IF_SAFER",
#'     "TILED=YES",
#'     "BLOCKXSIZE=512",
#'     "BLOCKYSIZE=512"
#'   ),
#'   agg_function = agg_function,
#'   agg_join = agg_join,
#'   finalizer = finalizer
#' )
#'
#' close(atl08_h5)
#'
#' @import data.table
#' @export
ATL08_seg_attributes_h5_gridStat <- function(
    atl08_dir,
    metrics = c(
      "h_canopy",
      "canopy_rh_conf",
      "h_median_canopy_abs",
      "h_min_canopy",
      "h_mean_canopy_abs",
      "h_median_canopy",
      "h_canopy_abs",
      "toc_roughness",
      "h_min_canopy_abs",
      "h_dif_canopy",
      "h_canopy_quad",
      "h_canopy_20m",
      "n_ca_photons",
      "photon_rate_can",
      "centroid_height",
      "canopy_h_metrics_abs",
      "h_mean_canopy",
      "subset_can_flag",
      "canopy_h_metrics",
      "n_toc_photons",
      "h_max_canopy_abs",
      "h_canopy_uncertainty",
      "canopy_openness",
      "h_max_canopy",
      "segment_cover"
    ),
    beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
    out_root,
    ul_lat,
    ul_lon,
    lr_lat,
    lr_lon,
    res,
    creation_options = def_co,
    agg_function = default_agg_function,
    agg_join = default_agg_join,
    finalizer = default_finalizer) {

      stopifnot("gdalBindings package is not installed. Please install it to use this function." = requireNamespace("gdalBindings", quietly = TRUE))

  block_inds <-
    block_xind <-
    block_yind <-
    inds <-
    latitude <-
    longitude <-
    x_block <-
    y_block <-
    x_ind <-
    y_ind <-
    NULL
  projstring <- 'GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.01745329251994328,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]'
  atl08_list <- sapply(atl08_dir, function(search_path) {
    list.files(search_path,
      pattern = "08.*h5",
      recursive = TRUE, full.names = TRUE
    )
  })
  total_files <- length(atl08_list)

  xres <- res[1]
  yres <- res[2]

  cols.coord <- c("latitude", "longitude")

  metricCounter <- 0
  nMetrics <- length(metrics)

  call <- parse(text = as.character(eval.parent(substitute(agg_function)))[2])
  tempenv <- new.env()
  tempenv$x <- c(1,1)
  stats <- eval(call, envir = tempenv)
  classes <- lapply(stats, class)
  stats_names <- names(stats)
  # metric = metrics[1]
  for (metric in metrics) {
    metricCounter <- metricCounter + 1
    message(sprintf("Metric %s (%d/%d)", metric, metricCounter, nMetrics), appendLF = TRUE)
    cols <- c(cols.coord, metric)

    rast_paths <- sprintf("%s_%s_%s.tif", out_root, metric, stats_names)
    rasts <- list()

    # stat_ind = 1
    for (stat_ind in seq_along(stats_names)) {
      datatype <- gdalBindings::GDALDataType$GDT_Float64
      nodata <- -9999.0
      if (classes[[stat_ind]] == "integer") {
        datatype <- gdalBindings::GDALDataType$GDT_Int32
        nodata <- 0
      }
      rasts[[stats_names[[stat_ind]]]] <- gdalBindings::createDataset(
        raster_path = rast_paths[[stat_ind]],
        nbands = 1,
        datatype = datatype,
        projstring = projstring,
        lr_lat = lr_lat,
        ul_lat = ul_lat,
        ul_lon = ul_lon,
        lr_lon = lr_lon,
        res = c(xres, yres),
        nodata = nodata,
        co = creation_options
      )
    }


    bands <- lapply(rasts, function(x) x[[1]])

    block_x_size <- bands[[1]]$GetBlockXSize()
    block_y_size <- bands[[1]]$GetBlockYSize()

    file_index <- 0
    # atl08_path = atl08_list[1]
    for (atl08_path in atl08_list) {
      file_index <- file_index + 1
      message(sprintf("Reading file %s (%d/%d)", basename(atl08_path), file_index, total_files), appendLF = TRUE)

      # read H5
      atl08_h5 <- ATL08_read(atl08_path = atl08_path)

      vals <- ATL08_seg_attributes_dt(atl08_h5, beam = beam, attributes = cols[-c(1:2)])
      vals <- ATL08_seg_attributes_dt_clipBox(vals, lower_left_lon = ul_lon, upper_right_lon = lr_lon, lower_left_lat = lr_lat, upper_right_lat = ul_lat)
      # Use only night photons
      # cols_night_flag = c(setdiff(cols, "night_flag"))
      # vals = vals[night_flag == 1, cols_night_flag, with = FALSE]

      # vals = vals[vals@data$agbd >= 0]


      # Goto next file if no data is available after clipping
      if (nrow(vals) == 0) next

      # Compute x/y indices (0-based) from lon/lat
      vals[, x_ind := as.integer(vals[, floor((longitude - ul_lon) / xres)])]
      vals[, y_ind := as.integer(vals[, floor((latitude - ul_lat) / yres)])]
      # Compute single index (1-based) aggregating x/y indices based on block size
      vals[, inds := 1 + x_ind + y_ind * block_x_size]
      names(vals) <- gsub(metric, "x", names(vals))

      # Compute statistics from the function provided by the user
      aggs <- vals[, eval(call), by = list(inds, x_ind, y_ind)]

      # Calculate x_block/y_block indices
      aggs[,
        c("x_block", "y_block") := lapply(.SD, function(x) as.integer(floor(x / block_x_size))),
        .SDcols = c("x_ind", "y_ind")
      ]
      aggs[, `:=`(
        block_xind = x_ind - x_block * block_x_size,
        block_yind = y_ind - y_block * block_y_size
      )]
      aggs[, block_inds := 1 + block_xind + block_yind * block_x_size]

      # Get calculated stats aggregated by each block
      blocks <- aggs[, list(vals = list(.SD)), by = list(x_block, y_block), .SDcols = c(stats_names, "block_inds")]

      thisEnv <- new.env()
      assign("ii", 0, thisEnv)
      total_rows <- nrow(blocks)

      # Join aggregates for each block of raster
      # ii = 1
      invisible(apply(blocks, 1, function(row) {
        ii <- get("ii", thisEnv) + 1
        assign("ii", ii, thisEnv)

        message(sprintf("\rProcessing blocks...%.2f%%", (100.0 * ii) / total_rows), appendLF = F)
        agg1 <- data.table::as.data.table(
          lapply(bands, function(x) x[[row$x_block, row$y_block]])
        )

        agg1[row$vals$block_inds] <- agg_join(agg1[row$vals$block_inds], row$vals[, 1:(ncol(row$vals) - 1)])
        lapply(stats_names, function(x) bands[[x]][[row$x_block, row$y_block]] <- agg1[[x]])
      }))
      message()
      rm(list = ls(envir = thisEnv), envir = thisEnv)
      rm(thisEnv)

      close(atl08_h5)
    }
    # Update statistics for bands
    lapply(bands, function(x) x$CalculateStatistics())

    invisible(lapply(names(finalizer), function(x) {
      rast_name <- sprintf("%s_%s_%s.tif", out_root, metric, x)
      message(sprintf("Writing raster: %s", rast_name))
      rast <- gdalBindings::createDataset(
        raster_path = rast_name,
        nbands = 1,
        datatype = gdalBindings::GDALDataType$GDT_Float64,
        projstring = projstring,
        lr_lat = lr_lat,
        ul_lat = ul_lat,
        ul_lon = ul_lon,
        lr_lon = lr_lon,
        res = c(xres, yres),
        nodata = -9999.0,
        co = creation_options
      )

      band <- rast[[1]]
      formula <- finalizer[[x]]
      gdalBindings::formulaCalculate(formula, bands, band)
      rast$Close()
    }))

    lapply(rasts, function(x) x$Close())
  }
}

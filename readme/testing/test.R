# Specifying the path to GEDI leveatl08_canopy_dt data (zip file)
library(rICESat2Veg)
library(data.table)

# Specifying the path to ATL08 file (zip file)
outdir <- tempdir()

# atl08_zip <- system.file("extdata",
#   "ATL08_20220401221822_01501506_005_01.h5",
#   package = "rICESat2Veg"
# )

# # Unzipping ATL08 file
# atl08_path <- unzip(atl08_zip, exdir = outdir)

# Bounding rectangle coordinates
ul_lat <- -13.72016
ul_lon <- -44.14000
lr_lat <- -13.74998
lr_lon <- -44.11009

res <- 100 # meters
lat_to_met_factor <- 1 / 110540
lon_to_met_factor <- 1 / 111320
xres <- lon_to_met_factor * res
yres <- lat_to_met_factor * res

agg_function <- ~ data.table(
  min = min(x),
  max = max(x),
  sum = sum(x),
  n = length(x)
)

agg_join <- function(agg1, agg2) {
  agg1[is.na(agg1)] <- 0
  data.table(
    min = pmin(agg1$min, agg2$min),
    max = pmax(agg1$max, agg2$max),
    sum = agg1$sum + agg2$sum,
    n = agg1$n + agg2$n
  )
}

finalizer <- list(
  mean = "sum/n",
  range = "max-min"
)

atl08_path <- "inst/extdat/ATL08_20220401221822_01501506_005_01.h5"
metrics <- c("h_canopy")
out_root <- file.path(outdir, "output")
ul_lat <- ul_lat
ul_lon <- ul_lon
lr_lat <- lr_lat
lr_lon <- lr_lon
res <- c(xres, -yres)
creation_options <- c(
  "COMPRESS=DEFLATE",
  "BIGTIFF=IF_SAFER",
  "TILED=YES",
  "BLOCKXSIZE=512",
  "BLOCKYSIZE=512"
)
agg_function <- agg_function
agg_join <- agg_join
finalizer <- finalizer



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
atl08_list = atl08_path
atl08_list <- sapply(atl08_path, function(search_path) {
  list.files(search_path,
    pattern = "ATL08_.*h5",
    recursive = TRUE, full.name = TRUE
  )
})
total_files <- length(atl08_list)

xres <- res[1]
yres <- res[2]


cols.coord <- c("latitude", "longitude")

metricCounter <- 0
nMetrics <- length(metrics)


func <- lazyeval::f_interp(agg_function)
call <- lazyeval::as_call(func)
x <- 1
stats <- eval(call)
classes <- lapply(stats, class)
stats <- names(stats)
# metric = metrics[1]
for (metric in metrics) {
  metricCounter <- metricCounter + 1
  message(sprintf("Metric %s (%d/%d)", metric, metricCounter, nMetrics), appendLF = T)
  cols <- c(cols.coord, metric)

  rast_paths <- sprintf("%s_%s_%s.tif", out_root, metric, stats)
  rasts <- list()

  # stat_ind = 1
  for (stat_ind in seq_along(stats)) {
    datatype <- GDALDataType$GDT_Float64
    nodata <- -9999.0
    if (classes[[stat_ind]] == "integer") {
      datatype <- GDALDataType$GDT_Int32
      nodata <- 0
    }
    rasts[[stats[[stat_ind]]]] <- createDataset(
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


  xsize <- rasts[[1]]$GetRasterXSize()
  ysize <- rasts[[1]]$GetRasterYSize()

  bands <- lapply(rasts, function(x) x[[1]])

  block_x_size <- bands[[1]]$GetBlockXSize()
  block_y_size <- bands[[1]]$GetBlockYSize()

  file_index <- 0
  # atl08_path = atl08_list[1]
  for (atl08_path in atl08_list) {
    file_index <- file_index + 1
    message(sprintf("Reading file %s (%d/%d)", basename(atl08_path), file_index, total_files), appendLF = T)
    atl08 <- ATL08_read(atl08_path)
    
    atl03 <- ATL03_read(gsub('ATL08','ATL03',atl08_path))

    atl08_canopy_dt <- ATL03_ATL08_join_dt(atl03, atl08)
    vals <- ATL08_canopy_attributes_dt(atl08_canopy_dt, beam = beam, canopy_attribute = cols[-c(1:2)])

    ## Clip metrics by extent
    vals <- ATL08_canopy_dt_clipBox(vals, ul_lon, lr_lon, lr_lat, ul_lat)

    # vals = ATL08_canopy_dt_clipGeometry(vals, polygon=polygon, split_by=NULL)

    vals <- vals@dt
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
    blocks <- aggs[, list(vals = list(.SD)), by = list(x_block, y_block), .SDcols = c(stats, "block_inds")]

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
      lapply(stats, function(x) bands[[x]][[row$x_block, row$y_block]] <- agg1[[x]])
    }))
    message()
    rm(list = ls(envir = thisEnv), envir = thisEnv)
    rm(thisEnv)

    close(atl08_canopy_dt)
  }
  # Update statistics for bands
  lapply(bands, function(x) x$CalculateStatistics())

  finalize_rasts <- lapply(names(finalizer), function(x) {
    rast_name <- sprintf("%s_%s_%s.tif", out_root, metric, x)
    message(sprintf("Writing raster: %s", rast_name))
    rast <- createDataset(
      raster_path = rast_name,
      nbands = 1,
      datatype = GDALDataType$GDT_Float64,
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
    formulaCalculate(formula, bands, band)
    rast$Close()
  })

  lapply(rasts, function(x) x$Close())
}

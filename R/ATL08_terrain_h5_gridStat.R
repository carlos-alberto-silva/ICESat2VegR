def_co = c("COMPRESS=DEFLATE",
        "BIGTIFF=IF_SAFER",
        "TILED=YES",
        "BLOCKXSIZE=512",
        "BLOCKYSIZE=512"
        )

default_finalizer = list(
  sd = ~sqrt(M2/(n - 1)),
  skew = ~sqrt((n * (n - 1))) * ((sqrt(n) * M3) / (M2^1.5)) / (n - 2),
  kur = ~((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * ((n * M4) / (M2^2) - 3.0) + 6)
)

default_agg_function = ~data.table::data.table(
    n = length(x),
    M1 = mean(x,na.rm = T),
    M2 = e1071::moment(x, order = 2, center = TRUE, na.rm = T) * length(x),
    M3 = e1071::moment(x, order = 3, center = TRUE, na.rm = T) * length(x),
    M4 = e1071::moment(x, order = 4, center = TRUE, na.rm = T) * length(x),
    min = min(x, na.rm=T),
    max = max(x, na.rm=T)
  )

default_agg_join = function(x1, x2) {
  combined = data.table()
  x1$n[is.na(x1$n)] = 0
  x1$M1[is.na(x1$M1)] = 0
  x1$M2[is.na(x1$M2)] = 0
  x1$M3[is.na(x1$M3)] = 0
  x1$M4[is.na(x1$M4)] = 0
  x1$max[is.na(x1$max)] = -Inf
  x1$min[is.na(x1$min)] = Inf

  combined$n = x1$n + x2$n

  delta = x2$M1 - x1$M1
  delta2 = delta * delta
  delta3 = delta * delta2
  delta4 = delta2 * delta2

  combined$M1 = (x1$n * x1$M1 + x2$n * x2$M1) / combined$n

  combined$M2 = x1$M2 + x2$M2 +
    delta2 * x1$n * x2$n / combined$n

  combined$M3 = x1$M3 + x2$M3 +
    delta3 * x1$n * x2$n * (x1$n - x2$n) / (combined$n * combined$n)
  combined$M3 = combined$M3 + 3.0 * delta * (x1$n * x2$M2 - x2$n * x1$M2) / combined$n

  combined$M4 = x1$M4 + x2$M4 + delta4 * x1$n * x2$n * (x1$n * x1$n - x1$n * x2$n + x2$n * x2$n) /
    (combined$n * combined$n * combined$n)
  combined$M4 = combined$M4 + 6.0 * delta2 * (x1$n * x1$n * x2$M2 + x2$n * x2$n * x1$M2) / (combined$n * combined$n) +
    4.0 * delta * (x1$n * x2$M3 - x2$n * x1$M3) / combined$n

  combined$min = pmin(x1$min, x2$min, na.rm=F)
  combined$max = pmax(x1$max, x2$max, na.rm=F)
  return(combined)
}

#' Rasterize ATL08 terrain attributes from h5 files at large scale
#'
#' @description
#' This function will read multiple ATL08 H5 files and create a stack of raster layers: count, and 1st, 2nd, 3rd and 4th moments (count, m1, m2, m3 and m4) for each metric selected, from which we can calculate statistics such as Mean, SD, Skewness and Kurtosis.
#'
#' @param atl08_path CharacterVector. The directory paths where the ATL08 H5 files are stored;
#' @param metrics CharacterVector. A vector of terrain attributes available from ATL08 product (e.g. "h_te_best_fit")
#' @param out_root Character. The root name for the raster output files, the pattern is {out_root}_{metric}_{count/m1/m2/m3/m4}.tif. This should include the full path for the file.
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
#' This function will create seven different aggregate statistics
#' (n, m1, m2, m3, m4, min, max). m1 to m4 are the central moments.
#' One can calculate mean, standard deviation, skewness and kurtosis
#' with the following formulas according to Terriberry (2007) and
#' \insertCite{Joanes1998;textual}{rGEDI}:
#'
#' \deqn{ \bar{x} = m_1 }{mean = m1}
#'
#' \deqn{ s = \sqrt{\frac{m_2}{n - 1}} }{sd = sqrt(m2/(n - 1))}
#'
#' \deqn{ g_1 = \frac{\sqrt{n} \cdot m_3}{m_2^{1.5}} }{ g1 = (sqrt(n) * m3) / (m2^1.5)}
#'
#' \deqn{ g_2 = \frac{n \cdot m_4}{m_2^2} - 3 }{g2 = (n * m4) / (m2 * m2) - 3.0}
#'
#' \deqn{ skewness = \frac{\sqrt{n(n - 1)}}{n-2} g_1 }{skewness = sqrt((n * (n - 1))) * g1 / (n - 2)}
#'
#' \deqn{ kurtosis = \frac{n - 1}{(n - 2)(n - 3)}[(n + 1)g_2 + 6] }{kurtosis = ((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * g2 + 6)}
#'
#' The `agg_function` is a formula which return a data.table with the
#' aggregate function to perform over the data.
#' The default is:
#'
#' ```{r, eval=FALSE}
#' ~data.table(
#'     n = length(x),
#'     M1 = mean(x,na.rm = T),
#'     M2 = e1071::moment(x, order = 2, center = TRUE, na.rm = T) * length(x),
#'     M3 = e1071::moment(x, order = 3, center = TRUE, na.rm = T) * length(x),
#'     M4 = e1071::moment(x, order = 4, center = TRUE, na.rm = T) * length(x),
#'     min = min(x, na.rm=T),
#'     max = max(x, na.rm=T)
#'   )
#' ```
#'
#' The `agg_join` is a function to merge two data.table aggregates
#' from the `agg_function`. Since the h5 files will be aggregated
#' one by one, the statistics from the different h5 files should
#' have a function to merge. The default function is:
#'
#' ```{r, eval=FALSE}
#' function(x1, x2) {
#'     combined = data.table()
#'     x1$n[is.na(x1$n)] = 0
#'     x1$M1[is.na(x1$M1)] = 0
#'     x1$M2[is.na(x1$M2)] = 0
#'     x1$M3[is.na(x1$M3)] = 0
#'     x1$M4[is.na(x1$M4)] = 0
#'     x1$max[is.na(x1$max)] = -Inf
#'     x1$min[is.na(x1$min)] = Inf
#'
#'     combined$n = x1$n + x2$n
#'
#'     delta = x2$M1 - x1$M1
#'     delta2 = delta * delta
#'     delta3 = delta * delta2
#'     delta4 = delta2 * delta2
#'
#'     combined$M1 = (x1$n * x1$M1 + x2$n * x2$M1) / combined$n
#'
#'     combined$M2 = x1$M2 + x2$M2 +
#'       delta2 * x1$n * x2$n / combined$n
#'
#'     combined$M3 = x1$M3 + x2$M3 +
#'       delta3 * x1$n * x2$n * (x1$n - x2$n) / (combined$n * combined$n)
#'     combined$M3 = combined$M3 + 3.0 * delta * (x1$n * x2$M2 - x2$n * x1$M2) / combined$n
#'
#'     combined$M4 = x1$M4 + x2$M4 + delta4 * x1$n * x2$n * (x1$n * x1$n - x1$n * x2$n + x2$n * x2$n) /
#'       (combined$n * combined$n * combined$n)
#'     combined$M4 = combined$M4 + 6.0 * delta2 * (x1$n * x1$n * x2$M2 + x2$n * x2$n * x1$M2) / (combined$n * combined$n) +
#'       4.0 * delta * (x1$n * x2$M3 - x2$n * x1$M3) / combined$n
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
#' `skewness` and `kurtosis` based on the `M2`, `M3`, `M4` and `n`
#' values. It is defined as:
#'
#' ```{r, eval=FALSE}
#' list(
#'   sd = ~sqrt(M2/(n - 1)),
#'   skew = ~sqrt((n * (n - 1))) * ((sqrt(n) * M3) / (M2^1.5)) / (n - 2),
#'   kur = ~((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * ((n * M4) / (M2^2) - 3.0) + 6)
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
#' # Specifying the path to GEDI leveatl08_terrain_dt data (zip file)
#' library(rICESat2Veg)
#' library(data.table)
#'
#'# Specifying the path to ATL08 file (zip file)
#'outdir = tempdir()
#'
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#' # Bounding rectangle coordinates
#' ul_lat = -13.72016
#' ul_lon = -44.14000
#' lr_lat = -13.74998
#' lr_lon = -44.11009
#'
#' res = 100 # meters
#' lat_to_met_factor = 1 / 110540
#' lon_to_met_factor = 1 / 111320
#' xres = lon_to_met_factor * res
#' yres = lat_to_met_factor * res
#'
#' agg_function = ~data.table(
#'     min = min(x),
#'     max = max(x),
#'     sum = sum(x),
#'     n = length(x))
#'
#' agg_join = function(agg1, agg2) {
#' agg1[is.na(agg1)] = 0
#' data.table(
#'     min = pmin(agg1$min, agg2$min),
#'     max = pmax(agg1$max, agg2$max),
#'     sum = agg1$sum + agg2$sum,
#'     n = agg1$n + agg2$n
#' )
#' }
#'
#' finalizer = list(
#'     mean = "sum/n",
#'     range = "max-min"
#' )
#'
#' ATL08_terrain_h5_gridStat(
#'   atl08_path = outdir,
#'   metrics = c("h_te_best_fit"),
#'   out_root = file.path(outdir, "output"),
#'   ul_lat = ul_lat,
#'   ul_lon = ul_lon,
#'   lr_lat = lr_lat,
#'   lr_lon = lr_lon,
#'   res = c(xres, -yres),
#'   creation_options = c("COMPRESS=DEFLATE" ,
#'     "BIGTIFF=IF_SAFER",
#'     "TILED=YES",
#'     "BLOCKXSIZE=512",
#'     "BLOCKYSIZE=512"),
#'   agg_function = agg_function,
#'   agg_join = agg_join,
#'   finalizer = finalizer
#'   )
#'
#' close(leveatl08_terrain_dt)
#'
#' @import e1071
#' @import data.table
#' @export
ATL08_terrain_h5_gridStat = function(
    atl08_path,
    metrics=c(
      "h_te_best_fit",
      "h_te_best_fit_20m",
      "h_te_interp",
      "h_te_max",
      "h_te_mean",
      "h_te_median",
      "h_te_min",
      "h_te_mode",
      "h_te_rh25",
      "h_te_skew",
      "h_te_std",
      "h_te_uncertainty",
      "n_te_photons",
      "photon_rate_te",
      "subset_te_flag",
      "terrain_slope"),
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
    finalizer=default_finalizer
    ) {
  block_inds =
  block_xind =
  block_yind =
  inds =
  latitude =
  longitude =
  x_block =
  y_block =
  x_ind =
  y_ind =
    NULL
  projstring = 'GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.01745329251994328,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]'
  atl08_list = sapply(atl08_path, function(search_path) list.files(search_path,
                                                                 pattern = "ATL08_.*h5",
                                                                 recursive = TRUE, full.name = TRUE))
  total_files = length(atl08_list)

  xres = res[1]
  yres = res[2]


  cols.coord = c("latitude", "longitude")

  metricCounter = 0
  nMetrics = length(metrics)


  func = lazyeval::f_interp(agg_function)
  call = lazyeval::as_call(func)
  x = 1
  stats = eval(call)
  classes = lapply(stats, class)
  stats = names(stats)
  # metric = metrics[1]
  for (metric in metrics) {
    metricCounter = metricCounter + 1
    message(sprintf("Metric %s (%d/%d)", metric, metricCounter, nMetrics), appendLF = T)
    cols = c(cols.coord, metric)

    rast_paths = sprintf("%s_%s_%s.tif", out_root, metric, stats)
    rasts = list()
    for (stat_ind in 1:length(stats)) {
      datatype = GDALDataType$GDT_Float64
      nodata = -9999.0
      if (classes[[stat_ind]] == "integer") {
        datatype = GDALDataType$GDT_Int32
        nodata = 0
      }
      rasts[[stats[[stat_ind]]]] = createDataset(
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
        co = creation_options)
    }


    xsize = rasts[[1]]$GetRasterXSize()
    ysize = rasts[[1]]$GetRasterYSize()

    bands = lapply(rasts, function(x) x[[1]])

    block_x_size = bands[[1]]$GetBlockXSize()
    block_y_size = bands[[1]]$GetBlockYSize()

    file_index = 0
    # atl08_path = atl08_list[1]
    for (atl08_path in atl08_list) {
      file_index = file_index + 1
      message(sprintf("Reading file %s (%d/%d)", basename(atl08_path), file_index, total_files), appendLF = T)
      atl08_terrain_dt = ATL03_ATL08_joined_dt_LAS(atl08_path)

      vals = ATL08_terrain_attributes(atl08_terrain_dt, beam=beam, cols = cols[-c(1:2)])

      ## Clip metrics by extent
      vals = ATL08_terrain_dt_clipBox(vals, ul_lon, lr_lon, lr_lat, ul_lat)

      #vals = ATL08_terrain_dt_clipGeometry(vals, polygon=polygon, split_by=NULL)

      vals<-vals@dt
      # Use only night photons
      #cols_night_flag = c(setdiff(cols, "night_flag"))
      #vals = vals[night_flag == 1, cols_night_flag, with = FALSE]

      #vals = vals[vals@data$agbd >= 0]


      # Goto next file if no data is available after clipping
      if (nrow(vals) == 0) next

      # Compute x/y indices (0-based) from lon/lat
      vals[, x_ind := as.integer(vals[, floor((longitude - ul_lon) / xres)])]
      vals[, y_ind := as.integer(vals[, floor((latitude - ul_lat) / yres)])]
      # Compute single index (1-based) aggregating x/y indices based on block size
      vals[, inds := 1 + x_ind + y_ind * block_x_size]
      names(vals) = gsub(metric, "x", names(vals))

      # Compute statistics from the function provided by the user
      aggs = vals[,eval(call), by = list(inds, x_ind, y_ind)]

      # Calculate x_block/y_block indices
      aggs[,
        c("x_block", "y_block") := lapply(.SD, function(x) as.integer(floor(x / block_x_size))),
        .SDcols = c("x_ind", "y_ind")
      ]
      aggs[, `:=`(block_xind = x_ind - x_block * block_x_size,
                  block_yind = y_ind - y_block * block_y_size)]
      aggs[, block_inds := 1 + block_xind + block_yind * block_x_size]

      # Get calculated stats aggregated by each block
      blocks = aggs[,list(vals=list(.SD)) , by=list(x_block, y_block), .SDcols=c(stats, "block_inds")]

      thisEnv = new.env()
      assign("ii", 0, thisEnv)
      total_rows = nrow(blocks)

      # Join aggregates for each block of raster
      # ii = 1
      invisible(apply(blocks, 1, function(row) {
        ii = get("ii", thisEnv) + 1
        assign("ii", ii, thisEnv)

        message(sprintf("\rProcessing blocks...%.2f%%", (100.0 * ii) / total_rows), appendLF = F)
        agg1 = data.table::as.data.table(
          lapply(bands, function(x)x[[row$x_block, row$y_block]])
        )

        agg1[row$vals$block_inds] = agg_join(agg1[row$vals$block_inds], row$vals[,1:(ncol(row$vals)-1)])
        lapply(stats, function(x) bands[[x]][[row$x_block, row$y_block]] = agg1[[x]])
      }))
      message()
      rm(list = ls(envir=thisEnv), envir= thisEnv)
      rm(thisEnv)

      close(atl08_terrain_dt)
    }
    # Update statistics for bands
    lapply(bands, function(x) x$CalculateStatistics())

    finalize_rasts = lapply(names(finalizer), function(x) {
        rast_name = sprintf("%s_%s_%s.tif", out_root, metric, x)
        message(sprintf("Writing raster: %s", rast_name))
        rast = createDataset(
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
          co = creation_options)

        band = rast[[1]]
        formula = finalizer[[x]]
        formulaCalculate(formula, bands, band)
        rast$Close()
      })

    lapply(rasts, function(x) x$Close())
  }
}

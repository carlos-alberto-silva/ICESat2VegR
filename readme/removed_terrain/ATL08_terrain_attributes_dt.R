# ATL08_terrain.var.map
ATL08_terrain.var.map <- list()
ATL08_terrain.var.map[["h_te_best_fit"]] <- "h_te_best_fit"
ATL08_terrain.var.map[["h_te_best_fit_20m"]] <- "h_te_best_fit_20m"
ATL08_terrain.var.map[["h_te_interp"]] <- "h_te_interp"
ATL08_terrain.var.map[["h_te_max"]] <- "h_te_max"
ATL08_terrain.var.map[["h_te_mean"]] <- "h_te_mean"
ATL08_terrain.var.map[["h_te_median"]] <- "h_te_median"
ATL08_terrain.var.map[["h_te_mode"]] <- "h_te_mode"
ATL08_terrain.var.map[["h_te_rh25"]] <- "h_te_rh25"
ATL08_terrain.var.map[["h_te_skew"]] <- "h_te_skew"
ATL08_terrain.var.map[["h_te_std"]] <- "h_te_std"
ATL08_terrain.var.map[["h_canopy"]] <- "h_canopy"
ATL08_terrain.var.map[["h_te_uncertainty"]] <- "h_te_uncertainty"
ATL08_terrain.var.map[["n_te_photons"]] <- "n_te_photons"
ATL08_terrain.var.map[["photon_rate_te"]] <- "photon_rate_te"
ATL08_terrain.var.map[["subset_te_flag"]] <- "subset_te_flag"
ATL08_terrain.var.map[["terrain_slope"]] <- "terrain_slope"

#' ATL08 Terrain Metrics
#'
#' @description This function extracts terrain metrics from ICESat-2 ATL08 data
#'
#' @usage ATL08_terrain_attributes_dt(atl08_h5, beam)
#'
#' @param atl08_h5 A ICESat-2 ATL08 object (output of [ATL08_read()] function). An S4 object of class [ICESat2VegR::icesat2.atl08_dt].
#' An S4 object of class [ICESat2VegR::icesat2.atl08_dt].
#' @param beam Character vector indicating beams to process (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#' @param terrain_attribute A character vector containing the list of metrics to be extracted. See the default columns in the description.
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the ATL08-derived terrain metrics.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#'
#' @examples
#'
#' # Specifying the path to ATL08 file (zip file)
#' outdir <- tempdir()
#' atl08_zip <- system.file("extdata",
#'   "ATL08_20220401221822_01501506_005_01.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ATL08 file
#' atl08_path <- unzip(atl08_zip, exdir = outdir)
#'
#' # Reading ATL08 data (h5 file)
# atl08_h5<-ATL08_read(atl08_path=atl08_path)
#'
#' # Extracting ATL08-derived Terrain Metrics
#' terrain_metrics <- ATL08_terrain_attributes_dt(atl08_h5 = atl08_h5)
#' head(terrain_metrics)
#'
#' close(atl08_h5)
#' @export
ATL08_terrain_attributes_dt <- function(atl08_h5,
                                        beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                        terrain_attribute = c(
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
                                          "terrain_slope"
                                        )) {
  # Check file input
  if (!class(atl08_h5) == "icesat2.atl08_h5") {
    stop("atl08_h5 must be an object of class 'icesat2.atl08_h5' - output of [ATL08_read()] function ")
  }

  # Check beams to select
  groups_id <- getBeams(atl08_h5)


  beam <- intersect(groups_id, beam)
  terrain.dt <- data.table::data.table()

  pb <- utils::txtProgressBar(min = 0, max = length(beam), style = 3)

  i_s <- 0

  if (length(terrain_attribute) > 1) {
    for (i in beam) {
      i_s <- i_s + 1
      utils::setTxtProgressBar(pb, i_s)

      atl08_h5v2_i <- atl08_h5[[paste0(i, "/land_segments/terrain")]]

      lat_i <- atl08_h5[[paste0(i, "/land_segments/latitude")]][]
      lon_i <- atl08_h5[[paste0(i, "/land_segments/longitude")]][]

      m <- data.table::data.table(latitude = lat_i, longitude = lon_i, beam = i)

      for (col in terrain_attribute) {
        # print(col)
        metric_address <- ATL08_terrain.var.map[[col]]

        if (is.null(metric_address)) {
          if (h5exists(atl08_h5v2_i, col)) {
            metric_address <- col
          } else {
            if (!col %in% names(atl08_h5v2_i)) {
              warning(
                sprintf(
                  "The column '%s' is not available in the ATL08 product!",
                  col
                )
              )
            }
            m[, eval(col) := NA]
            next
          }
        }
        base_addr <- gsub("^(.*)/.*", "\\1", metric_address)
        if (h5exists(atl08_h5v2_i, base_addr) && h5exists(atl08_h5v2_i, metric_address)) {
          if (metric_address %in% c("h_te_best_fit_20m", "subset_te_flag")) {
            m <- cbind(m, t(atl08_h5v2_i[[metric_address]][, ]))

            if (col == "subset_te_flag") {
              colnames(m)[(ncol(m) - 4):ncol(m)] <- paste0("subset_te_flag_geo_", 1:5)
            }
            if (col == "h_te_best_fit_20m") {
              colnames(m)[(ncol(m) - 4):ncol(m)] <- paste0("h_te_best_fit_20m_", 1:5)
            }
          } else {
            m[, eval(col) := atl08_h5v2_i[[metric_address]][]]
          }
        }
      }
      terrain_dt <- data.table::rbindlist(list(terrain.dt, m), fill = TRUE)
    }
  }
  setattr(terrain_dt, "class", c("icesat2.atl08_dt", "data.table", "data.frame"))
  close(pb)

  return(terrain_dt)
}

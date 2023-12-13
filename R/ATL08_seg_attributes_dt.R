# ATL08_canopy.var.map
ATL08_canopy.var.map <- list()
ATL08_canopy.var.map[["h_canopy"]] <- "h_canopy"
ATL08_canopy.var.map[["canopy_rh_conf"]] <- "canopy_rh_conf"
ATL08_canopy.var.map[["h_median_canopy_abs"]] <- "h_median_canopy_abs"
ATL08_canopy.var.map[["h_min_canopy"]] <- "h_min_canopy"
ATL08_canopy.var.map[["h_mean_canopy_abs"]] <- "h_mean_canopy_abs"
ATL08_canopy.var.map[["h_median_canopy"]] <- "h_median_canopy"
ATL08_canopy.var.map[["h_canopy_abs"]] <- "h_canopy_abs"
ATL08_canopy.var.map[["toc_roughness"]] <- "toc_roughness"
ATL08_canopy.var.map[["h_min_canopy_abs"]] <- "h_min_canopy_abs"
ATL08_canopy.var.map[["h_dif_canopy"]] <- "h_dif_canopy"
ATL08_canopy.var.map[["h_canopy_quad"]] <- "h_canopy_quad"
ATL08_canopy.var.map[["h_canopy_20m"]] <- "h_canopy_20m"
ATL08_canopy.var.map[["n_ca_photons"]] <- "n_ca_photons"
ATL08_canopy.var.map[["photon_rate_can"]] <- "photon_rate_can"
ATL08_canopy.var.map[["centroid_height"]] <- "centroid_height"
ATL08_canopy.var.map[["canopy_h_metrics_abs"]] <- "canopy_h_metrics_abs"
ATL08_canopy.var.map[["h_mean_canopy"]] <- "h_mean_canopy"
ATL08_canopy.var.map[["subset_can_flag"]] <- "subset_can_flag"
ATL08_canopy.var.map[["canopy_h_metrics"]] <- "canopy_h_metrics"
ATL08_canopy.var.map[["n_toc_photons"]] <- "n_toc_photons"
ATL08_canopy.var.map[["h_max_canopy_abs"]] <- "h_max_canopy_abs"
ATL08_canopy.var.map[["h_canopy_uncertainty"]] <- "h_canopy_uncertainty"
ATL08_canopy.var.map[["canopy_openness"]] <- "canopy_openness"
ATL08_canopy.var.map[["h_max_canopy"]] <- "h_max_canopy"
ATL08_canopy.var.map[["segment_cover"]] <- "segment_cover"


#' ATL08 Canopy Height Metrics
#'
#' @description This function extracts canopy height metrics from ICESat-2 ATL08 data
#'
#' @usage ATL08_seg_attributes_dt(atl08_h5, beam)
#'
#' @param atl08_h5 A ICESat-2 ATL08 object (output of [ATL08_read()] function). An S4 object of class [rICESat2Veg::icesat2.atl08_dt].
#' An S4 object of class [rICESat2Veg::icesat2.atl08_dt].
#' @param beam Character vector indicating beams to process (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#' @param canopy_attribute A character vector containing the list of metrics to be extracted. See the default columns in the description.
#'
#' @return Returns an S4 object of class [rICESat2Veg::icesat2.atl08_dt]
#' containing the ATL08-derived vegetation relative heights.
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
#'   package = "rICESat2Veg"
#' )
#'
#' # Unzipping ATL08 file
#' atl08_path <- unzip(atl08_zip, exdir = outdir)
#'
#' # Reading ATL08 data (h5 file)
# atl08_h5<-ATL08_read(atl08_path=atl08_path)
#'
#' # Extracting ATL08-derived Canopy Metrics
#' canopy_metrics <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)
#' head(canopy_metrics)
#'
#' close(atl08_h5)
#' @export
ATL08_seg_attributes_dt <- function(atl08_h5,
                                       beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                       canopy_attribute = c(
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
                                       )) {
  # Check file input
  if (!class(atl08_h5) == "icesat2.atl08_h5") {
    stop("atl08_h5 must be an object of class 'icesat2.atl08_h5' - output of [ATL08_read()] function ")
  }

  # Check beams to select
  groups_id <- getBeams(atl08_h5, recursive = F)

  check_beams <- groups_id %in% beam
  beam <- groups_id[check_beams]

  canopy.dt <- data.table::data.table()

  pb <- utils::txtProgressBar(min = 0, max = length(beam), style = 3)

  i_s <- 0

  if (length(canopy_attribute) > 0) {
    for (i in beam) {
      i_s <- i_s + 1
      utils::setTxtProgressBar(pb, i_s)
      atl08_canopy_h5 <- atl08_h5[[paste0(i, "/land_segments/canopy")]]
      lat_i <- atl08_h5[[paste0(i, "/land_segments/latitude")]][]
      lon_i <- atl08_h5[[paste0(i, "/land_segments/longitude")]][]

      m <- data.table::data.table(latitude = lat_i, longitude = lon_i, beam = i)

      for (col in canopy_attribute) {
        # print(col)
        metric_address <- ATL08_canopy.var.map[[col]]

        if (is.null(metric_address)) {
          if (atl08_canopy_h5$exists(col)) {
            metric_address <- col
          } else {
            if (i.s == 1) {
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

        if (atl08_canopy_h5$exists(base_addr) && atl08_canopy_h5$exists(metric_address)) {
          if (metric_address %in% c("h_canopy_20m", "canopy_h_metrics_abs", "subset_can_flag", "canopy_h_metrics")) {
            m <- cbind(m, t(atl08_canopy_h5[[metric_address]][, ]))

            if (col == "h_canopy_20m") {
              colnames(m)[(ncol(m) - 4):ncol(m)] <- paste0("h_canopy_20m_geo_", 1:5)
            }
            if (col == "canopy_h_metrics_abs") {
              colnames(m)[(ncol(m) - 17):ncol(m)] <- paste0("AH", seq(10, 95, 5))
            }
            if (col == "subset_can_flag") {
              colnames(m)[(ncol(m) - 4):ncol(m)] <- paste0("subset_can_flag_geo_", 1:5)
            }
            if (col == "canopy_h_metrics") {
              colnames(m)[(ncol(m) - 17):ncol(m)] <- paste0("RH", seq(10, 95, 5))
            }
          } else {
            m[, eval(col) := atl08_canopy_h5[[metric_address]][]]
          }
        }
      }
      canopy.dt <- data.table::rbindlist(list(canopy.dt, m), fill = TRUE)
    }
  }


  canopy_dt <- new("icesat2.atl08_dt", dt = canopy.dt)

  close(pb)

  return(canopy_dt)
}

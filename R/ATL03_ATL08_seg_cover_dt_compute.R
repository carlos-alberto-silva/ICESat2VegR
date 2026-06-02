#' Compute canopy cover from ATL03/ATL08 classified photons
#'
#' @description
#' Computes canopy cover metrics from ATL03 and ATL08 classified photons
#' within each segment. Cover is calculated using the ratio of ground
#' to vegetation photons, optionally weighted by a reflectance ratio.
#' The function adds three new columns (\code{cover}, \code{n_g},
#' \code{n_v}) to the input object and returns it.
#'
#' @param atl03_atl08_dt An object of class
#'   [`ICESat2VegR::icesat2.atl03_atl08_seg_dt-class`] containing ATL03
#'   and ATL08 data (output of [ATL03_ATL08_segment_create()] function).
#' @param reflectance_ratio numeric. The reflectance ratio
#'   \eqn{\rho_v / \rho_g}, where \eqn{\rho_v} is the vegetation
#'   reflectance and \eqn{\rho_g} is the ground reflectance.
#'   Default is \code{1} (same reflectance for ground and vegetation).
#'
#' @return Returns the input object with three additional columns:
#'   \describe{
#'     \item{cover}{Numeric. Canopy cover fraction per segment (0-1).}
#'     \item{n_g}{Integer. Number of ground photons per segment
#'       (\code{classed_pc_flag == 1}).}
#'     \item{n_v}{Integer. Number of vegetation photons per segment
#'       (\code{classed_pc_flag > 1}).}
#'   }
#'
#' @details
#' Cover is calculated with the formula:
#' \deqn{cover = \frac{1}{1 + \frac{\rho_v \cdot N_g}{\rho_g \cdot N_v}}}
#'
#' where \eqn{N_g} is the number of ground photons
#' (\code{classed_pc_flag == 1}) and \eqn{N_v} is the number of
#' vegetation photons (\code{classed_pc_flag > 1}) within each segment.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
#' # Specifying the path to ATL03 and ATL08 files
#' atl03_path <- system.file("extdata", "atl03_clip.h5", package = "ICESat2VegR")
#' atl08_path <- system.file("extdata", "atl08_clip.h5", package = "ICESat2VegR")
#'
#' # Reading ATL03 and ATL08 data (h5 files)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#'
#' # Converting to seg_dt class (required input)
#' atl03_atl08_seg_dt <- ATL03_ATL08_segment_create(
#'   atl03_atl08_dt,
#'   segment_length = 30
#' )
#'
#' # Computing canopy cover metrics per segment
#' cover <- ATL03_ATL08_seg_cover_dt_compute(atl03_atl08_seg_dt)
#' head(cover[, c("segment_id", "cover", "n_g", "n_v")])
#'
#' # Computing cover with a custom reflectance ratio (rho_v/rho_g = 1.5)
#' cover2 <- ATL03_ATL08_seg_cover_dt_compute(
#'   atl03_atl08_seg_dt,
#'   reflectance_ratio = 1.5
#' )
#' head(cover2[, c("segment_id", "cover", "n_g", "n_v")])
#'
#' close(atl03_h5)
#' close(atl08_h5)
#'
#' @import data.table
#' @export
ATL03_ATL08_seg_cover_dt_compute <- function(
    atl03_atl08_dt,
    reflectance_ratio = 1) {
  segment_id <-
    classed_pc_flag <- NA

  # ensure proper data.table reference to avoid shallow copy warning
  data.table::setDT(atl03_atl08_dt)

  atl03_atl08_seg_dt <- na.omit(atl03_atl08_dt)  # fix: use correct param name
  atl03_atl08_seg_dt[
    ,
    c(
      "cover",
      "n_g",
      "n_v"
    ) := list(
      1 / (
        1 +
          reflectance_ratio *
          sum(classed_pc_flag == 1) /
          sum(classed_pc_flag > 1)
      ),
      sum(classed_pc_flag == 1),
      sum(classed_pc_flag > 1)
    ),
    by = list(segment_id)
  ]
  return(atl03_atl08_seg_dt)
}

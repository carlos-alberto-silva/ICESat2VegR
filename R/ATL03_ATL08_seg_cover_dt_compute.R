#' Compute segments id for a given segment length
#'
#' @description This function reads the ICESat-2 Land and
#' Vegetation Along-Track Products (ATL08) as h5 file.
#'
#' @param atl03_atl08_dt [`ICESat2VegR::icesat2.atl03atl08_dt-class`].
#' The output of the [`ATL03_ATL08_photons_attributes_dt_join()`].
#' @param reflectance_ratio Numeric. The reflectance ratio to use to compute the coverage metric.
#' \eqn{\rho_v / \rho_g}, where \eqn{\rho_v} is the vegetation reflectance
#' and \eqn{\rho_v} is the ground reflectance. Default is 1 (same reflectance).
#'
#' @details
#' Coverage is calculated with the formula
#' \deqn{\frac{1}{1 + \frac{\rho_v N_g}{\rho_g N_v}}}
#'
#' @return Returns an S4 object of class [`data.table::data.table-class`] containing ICESat-2 ATL08 data.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
#' # ATL08 file path
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ICESat-2 ATL08 data (h5 file)
#' atl08 <- ATL08_read(atl08_path = atl08_path)
#'
#' close(atl08)
#' @import hdf5r
#' @export
ATL03_ATL08_seg_cover_dt_compute <- function(
    atl03_atl08_dt,
    reflectance_ratio = 1) {
  segment_id <-
    classed_pc_flag <- NA

  atl03_atl08_seg_dt <- na.omit(atl03_atl08_seg_dt)

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

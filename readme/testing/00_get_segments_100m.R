# devtools::load_all()
all_data = list()
all_data <- parallel::mclapply(seq_along(granules03), mc.cores = 4, function(ii) {
  atl08_h5 <- ATL08_read(granules08[ii])
  atl08_dt <- ATL08_seg_attributes_dt(atl08_h5, power_beam_filter = TRUE, attribute = c(
    "h_canopy",
    "night_flag",
    "n_ca_photons",
    "n_te_photons",
    "n_toc_photons",
    "h_canopy_20m",
    "longitude_20m",
    "latitude_20m"
  ))

  atl08_dt

  atl08_dt <- atl08_dt[
    h_canopy < 50
  ]
  if (nrow(atl08_dt) > 0) {
    atl08_seg_dt <- ATL08_seg_attributes_dt_clipGeometry(atl08_dt, geom)
    atl08_seg_dt$nid <- NULL

    if (nrow(atl08_seg_dt) > 0) {
      close(atl08_h5)
      return(atl08_seg_dt)
    }
  }

  close(atl08_h5)
  return(list())
})

all_data <- parallel::mclapply(seq_along(granules03), mc.cores = 4, function(ii) {
  atl03_h5 <- ATL03_read(granules03[ii])
  atl08_h5 <- ATL08_read(granules08[ii])
  atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5, strong_beam_filter = TRUE)
  out_name <- gsub("ATL03", "ATL03_ATL08", basename(granules03[ii]))

  atl03_atl08_dt <- atl03_atl08_dt[
    ph_h < 50 & 
    night_flag == 1 & 
    ph_h >= 0 
    ]
  if (nrow(atl03_atl08_dt) > 0) {
    atl03_atl08_dt_seg <- ATL03_ATL08_segment_create(
      atl03_atl08_dt,
      segment_length = 30
    )

    atl03_atl08_seg_dt <- ATL03_ATL08_compute_seg_attributes_dt_segStat(
      atl03_atl08_dt_seg,
      list(
        h_canopy_ge0 = quantile(ph_h, 0.98),
        h_canopy_gt0 = quantile(ph_h[ph_h > 0], 0.98),
        n_ground = sum(classed_pc_flag == 1),
        n_mid_canopy = sum(classed_pc_flag == 2),
        n_top_canopy = sum(classed_pc_flag == 3),
        n_canopy_total = sum(classed_pc_flag >= 2)
      ),
      ph_class = c(1, 2, 3)
    )

    prepend_class(atl03_atl08_seg_dt, "icesat2.atl08_dt")
    atl03_atl08_seg_dt <- ATL08_seg_attributes_dt_clipGeometry(atl03_atl08_seg_dt, geom)
    atl03_atl08_seg_dt$nid <- NULL

    if (nrow(atl03_atl08_seg_dt) > 0) {
      close(atl03_h5)
      close(atl08_h5)
      return(atl03_atl08_seg_dt)
    }
  }

  close(atl03_h5)
  close(atl08_h5)
  return(list())
})
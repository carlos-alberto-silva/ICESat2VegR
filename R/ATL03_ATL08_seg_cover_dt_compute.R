ATL03_ATL08_seg_cover_dt_compute <- function(atl03_atl08_dt, atl08_dt, ground_veg_reflectance_ratio = 1) {
  atl03_atl08_dt <- na.omit(atl03_atl08_dt)

  atl03_atl08_dt[
    ,
    cover := list(
      cover = 1 / (
        1 +
          ground_veg_reflectance_ratio *
            sum(classed_pc_flag == 1) /
            sum(classed_pc_flag > 1)
      )
    ),
    by = list(ph_segment_id)
  ]
}
